#include "Submergence.hpp"

namespace moordyn {

SidewaysCylinderBuoyancy::SidewaysCylinderBuoyancy(real phi1,
                                                   real phi2,
                                                   real diameter,
                                                   real length)
  : phi1(phi1)
  , phi2(phi2)
  , d(diameter)
  , l(length)
  , phi1Values(phi1)
  , phi2Values(phi2)
{
}

BuoyancyResult2D
SidewaysCylinderBuoyancy::calculate()
{
	real cosDiff = phi1Values.cosPhi - phi2Values.cosPhi;
	if (abs(cosDiff) < 1e-6) {
		real wettedVolume = 0.5 * d * d * l * phi2Values.A;
		real yV = (-(d * d * d * l)) / 4 * phi2Values.B;
		real zV = wettedVolume * l / 2.0;
		return BuoyancyResult2D{ wettedVolume, yV, zV };
	}
	// Calculate the wetted volume
	real denom = cosDiff;
	real numerator = 0.5 * d * d * l;
	real wettedVolume = (numerator / denom) * (phi1Values.E - phi2Values.E);

	// Calculate the center of buoyancy
	real yV =
	    (-(d * d * d * l)) / (4 * cosDiff) * (phi1Values.F - phi2Values.F);

	real zV = ((d * d * l * l) / (4 * cosDiff * cosDiff)) *
	          ((phi1Values.E - 2 * phi2Values.E) * phi1Values.cosPhi +
	           phi2Values.E * phi2Values.cosPhi - phi1Values.F + phi2Values.F);

	return BuoyancyResult2D{ wettedVolume, yV, zV };
}

AnalyticalBuoyancyCalculator::AnalyticalBuoyancyCalculator(real length,
                                                           real diameter,
                                                           vec3 cylinderBase,
                                                           vec3 cylinderDir)
  : cylinderLength(length)
  , cylinderDiameter(diameter)
  , bottom_pos(cylinderBase)
  , q(cylinderDir)
  , r(diameter / 2.0)
{
	if (cylinderDir.z() < 0.0) {
		q = -cylinderDir;
		bottom_pos = cylinderBase + length * cylinderDir;
	}
	// the horizontal component of the cylinderDir is the magnitude of the
	// x and y components
	h = sqrt(q.x() * q.x() + q.y() * q.y());

	// This computes the vector that points from bottom_pos to the point on
	// the edge of the bottom face with the highest z value.
	// This vector should be orthogonal to q.
	if (h > 1e-6) {
		leftDir = vec3(-q.x() * q.z() / h, -q.y() * q.z() / h, h);
	} else {
		// just choose an arbitrary direction since we're basically vertical
		leftDir = vec3(1.0, 0.0, 0.0);
	}
	// std::cout << "leftDir = " << leftDir.transpose() << std::endl;

	// Compute the "corners" where top and bottom refer to the top and
	// bottom face and left and right refer to the highest point and lowest
	// point on those faces
	bottomLeft = bottom_pos + r * leftDir;
	bottomRight = bottom_pos - r * leftDir;
	topLeft = bottomLeft + length * q;
	topRight = bottomRight + length * q;

	area = pi * r * r;
}

AnalyticalBuoyancyCalculator::BuoyancyResult
AnalyticalBuoyancyCalculator::calculateBuoyancy()
{
	bool bottom_partial_submerged =
	    bottomRight.z() < 0.0 && bottomLeft.z() > 0.0;

	bool top_partial_submerged = topRight.z() < 0.0 && topLeft.z() > 0.0;
	// std::cout << "bottom partial sub? " << bottom_partial_submerged
	//           << ", top partial sub? " << top_partial_submerged
	//           << std::endl;
	if (bottom_partial_submerged && top_partial_submerged) {
		// this is the case that the external function wetted_volume can
		// handle. we just need to calculate psi1 and psi2

		vec3 bottomIntersect =
		    bottomRight + leftDir * (-bottomRight.z() / leftDir.z());
		real bottom_b = (bottomIntersect - bottomLeft).norm();
		real phi1 = acos((r - bottom_b) / r);

		vec3 topIntersect = topRight + leftDir * (-topRight.z() / leftDir.z());
		real top_b = (topIntersect - topLeft).norm();
		real phi2 = acos((r - top_b) / r);

		SidewaysCylinderBuoyancy buoyancyCalculation(
		    phi1, phi2, cylinderDiameter, cylinderLength);
		auto result = buoyancyCalculation.calculate();
		real V = result.wettedVolume;
		vec3 centerOfBuoyancy =
		    bottom_pos + q * result.zV / V - leftDir * result.yV / V;

		return BuoyancyResult{ V, centerOfBuoyancy };
	}
	if (!bottom_partial_submerged && !top_partial_submerged) {
		// in this case we are some sort of cylindrical segment
		// or we are fully in or out of the water
		if (bottomRight.z() > 0) {
			return { 0, bottom_pos };
		}
		if (topLeft.z() < 0) {
			return { area * cylinderLength,
				     bottom_pos + (cylinderLength / 2.0) * q };
		}

		// since neither base intersects the water and the cylinder is
		// neither fully in or out of the water, these calculations should
		// work properly
		vec3 leftIntersect = bottomLeft + q * (-bottomLeft.z() / q.z());
		vec3 rightIntersect = bottomRight + q * (-bottomRight.z() / q.z());

		real hl = (leftIntersect - bottomLeft).norm();
		real hr = (rightIntersect - bottomRight).norm();
		real cylinderVolume = hl * area;
		real cylinderZ = hl / 2.0;

		// the length of the cylinder the intersects the water surface
		real newL = hr - hl;
		if (newL > 1e-6) {
			// only consider the wedge if it has a meaningful length
			SidewaysCylinderBuoyancy buoyancyCalculation(
			    0, pi, cylinderDiameter, newL);
			auto wedgeResult = buoyancyCalculation.calculate();

			real totalVolume = cylinderVolume + wedgeResult.wettedVolume;
			real wedgeCobY = wedgeResult.yV / wedgeResult.wettedVolume;
			real wedgeCobZ = hl + wedgeResult.zV / wedgeResult.wettedVolume;

			real totalCobY =
			    (wedgeResult.wettedVolume * wedgeCobY) / totalVolume;
			real totalCobZ = (cylinderVolume * cylinderZ +
			                  wedgeResult.wettedVolume * wedgeCobZ) /
			                 totalVolume;
			vec3 centerOfBuoyancy =
			    bottom_pos + q * totalCobZ - leftDir * totalCobY;
			return { totalVolume, centerOfBuoyancy };
		} else {
			// just the result for the bottom cylinder
			vec3 centerOfBuoyancy = bottom_pos + q * cylinderZ;
			return { cylinderVolume, centerOfBuoyancy };
		}
	}
	if (bottom_partial_submerged && !top_partial_submerged) {
		// we just have a part of our bottom corner submerged in the
		// water

		vec3 rightIntersect = bottomRight + q * (-bottomRight.z() / q.z());
		vec3 bottomIntersect =
		    bottomRight + leftDir * (-bottomRight.z() / leftDir.z());
		real bottom_b = (bottomIntersect - bottomLeft).norm();
		real phi1 = acos((r - bottom_b) / r);

		real newL = (rightIntersect - bottomRight).norm();

		SidewaysCylinderBuoyancy buoyancyCalculation(
		    phi1, pi, cylinderDiameter, newL);
		auto result = buoyancyCalculation.calculate();

		real V = result.wettedVolume;
		vec3 centerOfBuoyancy =
		    bottom_pos + q * result.zV / V - leftDir * result.yV / V;

		return BuoyancyResult{ V, centerOfBuoyancy };
	}
	if (!bottom_partial_submerged && top_partial_submerged) {
		// The bottom is fully submerged and the top is partially
		// submerged. To find the volume we break it down into a fully
		// cylindrical part and a truncated cylinder part

		vec3 leftIntersect = bottomLeft + q * (-bottomLeft.z() / q.z());
		real hl = (leftIntersect - bottomLeft).norm();
		real cylinderVolume = hl * area;
		real cylinderZ = hl / 2.0;

		real newL = cylinderLength - hl;

		vec3 topIntersect = topRight + leftDir * (-topRight.z() / leftDir.z());
		real top_b = (topIntersect - topLeft).norm();
		real phi2 = acos((r - top_b) / r);

		SidewaysCylinderBuoyancy buoyancyCalculation(
		    0, phi2, cylinderDiameter, newL);
		auto wedgeResult = buoyancyCalculation.calculate();

		real totalVolume = cylinderVolume + wedgeResult.wettedVolume;
		real wedgeCobY = wedgeResult.yV / wedgeResult.wettedVolume;
		real wedgeCobZ = hl + wedgeResult.zV / wedgeResult.wettedVolume;
		real totalCobY = (wedgeResult.wettedVolume * wedgeCobY) / totalVolume;
		real totalCobZ = (cylinderVolume * cylinderZ +
		                  wedgeResult.wettedVolume * wedgeCobZ) /
		                 totalVolume;

		vec3 centerOfBuoyancy =
		    bottom_pos + q * totalCobZ - leftDir * totalCobY;
		return { totalVolume, centerOfBuoyancy };
	}
	throw moordyn::unhandled_error(
	    "AnalyticalBuoyancy Calculator reached an impossible state");
}

} // namespace moordyn