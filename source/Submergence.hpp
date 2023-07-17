#pragma once

#include "Misc.hpp"

namespace moordyn {

struct BuoyancyResult2D
{
	real wettedVolume;
	// Dimensions for center of buoyancy scaled by wetted volume
	real yV;
	real zV;
};

/**
 * @brief Does the computation of the wetted volume and 2D center of mass for a
 * cylinder where both bases intersect the water plane
 *
 * This probably shouldn't be used publicly, but is here so it can be tested
 * directly.
 */
class SidewaysCylinderBuoyancy
{

	// Input angles the water makes at the cylinder bases
	real phi1;
	real phi2;
	// cylinder diameter
	real d;
	// cylinder length
	real l;

	// Used to avoid redundant trig calculations
	struct TrigMemo
	{
		real cosPhi;
		real sinPhi;
		real A;
		real B;
		real C;
		real D;
		real E;
		real F;
		TrigMemo(real phi)
		{
			sinPhi = sin(phi);
			cosPhi = cos(phi);
			// this is the original equation, but can be expanded to avoid more
			// trig computation
			// A = 0.5 * (pi - phi + 0.5 * sin(2 * phi));
			A = 0.5 * (pi - phi + sinPhi * cosPhi);
			B = (-1.0 / 3.0) * sinPhi * sinPhi * sinPhi;

			// this is the original equation, but can be expanded to avoid trig
			// computation
			// C = (1.0 / 8.0) * (pi - phi + (1.0 / 4.0) * sin(4 *phi));
			C = (1.0 / 8.0) * (pi - phi +
			                   (cosPhi * cosPhi * cosPhi * sinPhi -
			                    sinPhi * sinPhi * sinPhi * cosPhi));
			E = A * cosPhi - B;
			F = B * cosPhi - C;
		}
	};

	TrigMemo phi1Values;
	TrigMemo phi2Values;

  public:
	SidewaysCylinderBuoyancy(real phi1, real phi2, real diameter, real length);
	BuoyancyResult2D calculate();
};

class AnalyticalBuoyancyCalculator
{

	/// Input Cylinder Data
	real cylinderLength;
	real cylinderDiameter;
	vec3 bottom_pos;
	vec3 q;

	/// Derived values

	/// Corner positions, (right means lower, left means upper)
	vec3 bottomRight, bottomLeft, topRight, topLeft;

	/// horizontal component of cylinderDir
	real h;

	/// radial vector that is the most up pointing
	vec3 leftDir;

	/// Cylinder radius
	real r;
	real area;

  public:
	AnalyticalBuoyancyCalculator(real length,
	                             real diameter,
	                             vec3 cylinderBase,
	                             vec3 cylinderDir);

	struct BuoyancyResult
	{
		real wettedVolume;
		vec3 centerOfBuoyancy;
	};

	BuoyancyResult calculateBuoyancy();
};

} // namespace moordyn