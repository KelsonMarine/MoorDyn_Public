#include <iostream>
#include <sstream>

#include "Misc.hpp"
#include "Submergence.hpp"
#include <catch2/catch_test_macros.hpp>
#include <tuple>
#include "catch2/catch_tostring.hpp"
#include "catch2/matchers/catch_matchers_templated.hpp"
#include "catch2/generators/catch_generators.hpp"
#include "util.h"

namespace Catch {
template<typename T, int N>
struct StringMaker<Eigen::Vector<T, N>>
{
	static std::string convert(const Eigen::Vector<T, N>& value)
	{
		Eigen::IOFormat testFmt(4, Eigen::DontAlignCols, ", ", "\n", "[", "]");
		std::stringstream ss;
		ss << (value.transpose()).format(testFmt);
		return ss.str();
	}
};
}

template<typename DerivedA>
struct IsCloseMatcher : Catch::Matchers::MatcherGenericBase
{
	IsCloseMatcher(
	    const Eigen::Ref<const DerivedA> a,
	    const typename DerivedA::RealScalar rtol =
	        Eigen::NumTraits<typename DerivedA::RealScalar>::dummy_precision(),
	    const typename DerivedA::RealScalar atol =
	        Eigen::NumTraits<typename DerivedA::RealScalar>::epsilon())
	  : a(a)
	  , rtol(rtol)
	  , atol(atol)
	{
	}

	template<typename DerivedB>
	bool match(const Eigen::DenseBase<DerivedB>& b) const
	{
		return ((a.derived() - b.derived()).array().abs() <=
		        (atol + rtol * b.derived().array().abs()))
		    .all();
	}

	std::string describe() const override
	{
		std::stringstream ss;
		ss << "Is close to: " << Catch::StringMaker<DerivedA>::convert(a)
		   << "\nrtol = " << rtol << ", atol = " << atol;
		return ss.str();
	}

  private:
	const Eigen::Ref<const DerivedA> a;
	const typename DerivedA::RealScalar rtol;
	const typename DerivedA::RealScalar atol;
};

template<typename T>
IsCloseMatcher<T>
IsClose(T value, double rtol = 1e-10, double atol = 1e-12)
{
	return IsCloseMatcher<T>(value, rtol, atol);
}

// md::real is faster to type
namespace md = moordyn;
using namespace md;

TEST_CASE("getH gives the cross product matrix")
{

	vec testVec;
	testVec << 1, 2, 3;

	vec v{ 0.3, 0.2, 0.1 };
	// getH() should create a matrix that replicates the behavior of the cross
	// product such that getH(v) * a == v.cross(a)
	REQUIRE_THAT(getH(v) * testVec, IsClose(v.cross(testVec)));
}

TEST_CASE("translateMass linear acceleration")
{
	/**
	 * This test imagines that we have some point whose center of mass
	 * is 1 meter in the x direction away from our reference point.
	 *
	 * A force applied in the y direction through this center of mass should
	 * result in no rotation and a acceleration in the y direction according to
	 * F = ma
	 *
	 * In our local coordinate system, this force will result in a torque.
	 *
	 * The mass matrix produced by translateMass should be such that we can
	 * correctly predict the acceleration in this situation
	 *
	 */

	md::real m = 1.0;
	mat mass = m * mat::Identity();
	mat sphereI = ((2.0 / 5.0) * m) * mat::Identity();

	vec offset{ 1.0, 1.0, 1.0 };

	mat6 M6 = translateMass(offset, mass);
	mat sphereIRef = sphereI - mass * getH(offset) * getH(offset);
	M6.bottomRightCorner<3, 3>() += sphereIRef;

	vec6 F = vec6::Zero();

	vec f3{ 0, 10, 0 };
	F.head<3>() = f3;
	F.tail<3>() = offset.cross(f3);

	// std::cout << "F = " << F.transpose() << std::endl;
	// std::cout << "M6 = \n" << M6 << std::endl;

	vec6 acc = solveMat6(M6, F);

	// linear acceleration by F = ma
	// no angular acceleration
	vec6 expectedAcc = vec6::Zero();
	expectedAcc.head<3>() = f3 / m;
	REQUIRE_THAT(acc, IsClose(expectedAcc));
}

TEST_CASE("translateMass6 linear acceleration")
{
	/**
	 * Like the testTranslateMass test except we do a series of two offsets
	 * We verify both that the acceleration is computed correctly,
	 * and that the mass matrix we get matched what we get by translating
	 * the mass by both offsets simultaneously.
	 */
	md::real m = 1.0;
	mat mass = m * mat::Identity();

	mat6 mass6 = mat6::Zero();
	mass6.topLeftCorner<3, 3>() = mass;

	vec offset{ 1.0, 0.0, 0.0 };
	vec offset2{ 2, 1, 0.2 };

	mat6 M6 = translateMass(offset, mass);

	M6 = translateMass6(offset2, M6);

	REQUIRE_THAT(translateMass((offset + offset2), mass), IsClose(M6));

	// we add some moment of inertia to prevent a singular matrix.
	// presume it's a sphere with radius 1
	mat sphereI = ((2.0 / 5.0) * m) * mat::Identity();
	mat sphereIRef =
	    sphereI - mass * getH(offset + offset2) * getH(offset + offset2);
	M6.bottomRightCorner<3, 3>() += sphereIRef;
	vec6 F = vec6::Zero();

	vec f3{ 0, 10, 0 };
	F.head<3>() = f3;
	F.tail<3>() = (offset + offset2).cross(f3);
	// std::cout << "F = " << F.transpose() << std::endl;
	// std::cout << "M6 = \n" << M6 << std::endl;
	// std::cout << "det(M6) = " << M6.determinant() << std::endl;
	vec6 acc = solveMat6(M6, F);

	// linear acceleration by F = ma
	// no angular acceleration
	vec6 expectedAcc = vec6::Zero();
	expectedAcc.head<3>() = f3 / m;

	REQUIRE_THAT(acc, IsClose(expectedAcc));
}

TEST_CASE("rotateMass simple")
{
	md::real m = 1.0;
	mat mass = m * mat::Identity();
	mat sphereI = ((2.0 / 5.0) * m) * mat::Identity();

	vec offset{ 1.0, 1.0, 0.0 };

	mat6 M6 = translateMass(offset, mass);

	vec3 axis{ 0, 0, 1 };
	axis.normalize();
	// rotate -90 degrees around the z axis
	Eigen::AngleAxisd rot(-pi / 2, axis);
	mat6 rotatedMass = rotateMass6(rot.toRotationMatrix(), M6);

	// this offset represents the offset after rotation
	vec newOffset{ 1, -1, 0.0 };

	REQUIRE_THAT(rot.toRotationMatrix() * offset, IsClose(newOffset));
	REQUIRE_THAT(translateMass(newOffset, mass), IsClose(rotatedMass));

	// add some moment of inertia to prevent singular mass matrix
	mat sphereIRef = sphereI - mass * getH(newOffset) * getH(newOffset);
	rotatedMass.bottomRightCorner<3, 3>() += sphereIRef;
	vec6 F = vec6::Zero();

	vec f3{ 0, 10, 0 };
	F.head<3>() = f3;
	F.tail<3>() = newOffset.cross(f3);

	vec6 acc = solveMat6(rotatedMass, F);

	vec6 expectedAcc = vec6::Zero();
	expectedAcc.head<3>() = f3 / m;

	REQUIRE_THAT(acc, IsClose(expectedAcc));
}

TEST_CASE("SidewaysCylinder Volume")
{
	md::real D = 0.5;
	md::real L = 1.0;
	using tup = std::tuple<md::real, md::real, md::real>;
	// precalculated results using existing, tested, implementation
	auto [phi1, phi2, expected_result] = GENERATE(
	    tup{ 0.0, 0.0, 0.19634954084936207 },
	    tup{ 0.0, 1.0471975511965976, 0.1806095061943583 },
	    tup{ 0.0, 2.0943951023931953, 0.12565301568124013 },
	    tup{ 0.0, 3.141592653589793, 0.09817477042468103 },
	    tup{ 1.0471975511965976, 0.0, 0.1806095061943583 },
	    tup{ 1.0471975511965976, 1.0471975511965976, 0.15796298776783843 },
	    tup{ 1.0471975511965976, 2.0943951023931953, 0.09817477042468106 },
	    tup{ 1.0471975511965976, 3.141592653589793, 0.07069652516812194 },
	    tup{ 2.0943951023931953, 0.0, 0.12565301568124013 },
	    tup{ 2.0943951023931953, 1.0471975511965976, 0.09817477042468106 },
	    tup{ 2.0943951023931953, 2.0943951023931953, 0.03838655308152367 },
	    tup{ 2.0943951023931953, 3.141592653589793, 0.015740034655003753 },
	    tup{ 3.141592653589793, 0.0, 0.09817477042468103 },
	    tup{ 3.141592653589793, 1.0471975511965976, 0.07069652516812194 },
	    tup{ 3.141592653589793, 2.0943951023931953, 0.015740034655003753 },
	    tup{ 3.141592653589793, 3.141592653589793, -7.654042494670958e-18 });
	SidewaysCylinderBuoyancy buoyancyCalculation(phi1, phi2, D, L);
	auto V = buoyancyCalculation.calculate().wettedVolume;
	// std::cout << "phi1 = " << phi1 << ", phi2 = " << phi2 << ", V = " << V
	//           << std::endl;
	// std::cout << "expected_result = " << expected_result << std::endl;
	if (!isclose(V, expected_result, 0.0, 1e-10)) {

		std::cout << "phi1 = " << phi1 << ", phi2 = " << phi2 << std::endl;
		std::cout << "calculated V = " << V
		          << ", expected V = " << expected_result << std::endl;
	}
	REQUIRE(isclose(V, expected_result, 0.0, 1e-10));
}

TEST_CASE("SidewaysCylinder COB")
{
	md::real D = 0.5;
	md::real L = 1.0;
	using t = std::tuple<md::real, md::real, md::real, md::real>;

	// Precalculated expected values based on existing, tested, implementation
	auto [phi1, phi2, expected_yV, expected_zV] = GENERATE(
	    t{ 0.0, 0.0, 0.0, 0.09817477042 },
	    t{ 0.0, 1.047197551, 0.003106863268, 0.08699217152 },
	    t{ 0.0, 2.094395102, 0.007145609779, 0.04878847916 },
	    t{ 0.0, 3.141592654, 0.006135923152, 0.03067961576 },
	    t{ 1.047197551, 0.0, 0.003106863268, 0.09361733468 },
	    t{ 1.047197551, 1.047197551, 0.006765823467, 0.07898149388 },
	    t{ 1.047197551, 2.094395102, 0.009164983035, 0.03893865001 },
	    t{ 1.047197551, 3.141592654, 0.007145609779, 0.0213102339 },
	    t{ 2.094395102, 0.0, 0.007145609779, 0.07686453652 },
	    t{ 2.094395102, 1.047197551, 0.009164983035, 0.05923612041 },
	    t{ 2.094395102, 2.094395102, 0.006765823467, 0.01919327654 },
	    t{ 2.094395102, 3.141592654, 0.003106863268, 0.004557435746 },
	    t{ 3.141592654, 0.0, 0.006135923152, 0.06749515467 },
	    t{ 3.141592654, 1.047197551, 0.007145609779, 0.04938629127 },
	    t{ 3.141592654, 2.094395102, 0.003106863268, 0.01118259891 },
	    t{ 3.141592654, 3.141592654, 1.913204185e-50, -3.827021247e-18 });

	SidewaysCylinderBuoyancy buoyancyCalculation(phi1, phi2, D, L);
	auto result = buoyancyCalculation.calculate();
	auto yV = result.yV;
	auto zV = result.zV;
	// std::cout << "phi1 = " << phi1 << ", phi2 = " << phi2 << ", yV = " << yV
	//           << ", zV = " << zV << std::endl;
	// std::cout << "expected_result = " << expected_result << std::endl;

	// if (!isclose(zV, expected_zV, 1e-5, 1e-8)) {

	// 	std::cout << "phi1 = " << phi1 << ", phi2 = " << phi2 << std::endl;
	// 	std::cout << "calculated zV = " << zV
	// 	          << ", expected zV = " << expected_zV << std::endl;
	// }
	REQUIRE(isclose(yV, expected_yV, 0.0, 1e-10));
	REQUIRE(isclose(zV, expected_zV, 0.0, 1e-10));
}

TEST_CASE("Rod buoyancy")
{

	struct P
	{
		md::real length;
		md::real diameter;
		vec3 bottom_pos;
		vec3 q;
		md::real expectedV;
		vec3 expectedCOB;
	};

	// Results from an existing, tested implementation
	// clang-format off
	P p = GENERATE(
P{1, 1, {0, 0., -2}, {-1.2246467991473532e-16, 0., -1.0}, 0.7853981633974483,  {-6.123233995736766e-17, 0., -2.5}},
P{1, 1, {0, 0., -2}, {-0.7818314824680299, 0., -0.6234898018587335}, 0.7853981633974483,  {-0.39091574123401496, 0., -2.3117449009293667}},
P{1, 1, {0, 0., -2}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.7853981633974483,  {-0.4874639560909118, 0., -1.8887395330218428}},
P{1, 1, {0, 0., -2}, {-0.4338837391175581, 0., 0.9009688679024191}, 0.7853981633974483,  {-0.21694186955877906, 0., -1.5495155660487905}},
P{1, 1, {0, 0., -2}, {0.4338837391175581, 0., 0.9009688679024191}, 0.7853981633974483,  {0.21694186955877906, 0., -1.5495155660487905}},
P{1, 1, {0, 0., -2}, {0.9749279121818236, 0., 0.22252093395631445}, 0.7853981633974483,  {0.4874639560909118, 0., -1.8887395330218428}},
P{1, 1, {0, 0., -2}, {0.7818314824680299, 0., -0.6234898018587335}, 0.7853981633974483,  {0.39091574123401496, 0., -2.3117449009293667}},
P{1, 1, {0, 0., -1}, {-1.2246467991473532e-16, 0., -1.0}, 0.7853981633974483,  {-6.123233995736766e-17, 0., -1.5}},
P{1, 1, {0, 0., -1}, {-0.7818314824680299, 0., -0.6234898018587335}, 0.7853981633974483,  {-0.39091574123401496, 0., -1.3117449009293667}},
P{1, 1, {0, 0., -1}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.7853981633974483,  {-0.4874639560909118, 0., -0.8887395330218428}},
P{1, 1, {0, 0., -1}, {-0.4338837391175581, 0., 0.9009688679024191}, 0.7761097654125233,  {-0.21870275033704578, 0., -0.556501233942026}},
P{1, 1, {0, 0., -1}, {0.4338837391175581, 0., 0.9009688679024191}, 0.7761097654125233,  {0.21870275033704578, 0., -0.556501233942026}},
P{1, 1, {0, 0., -1}, {0.9749279121818236, 0., 0.22252093395631445}, 0.7853981633974483,  {0.4874639560909118, 0., -0.8887395330218428}},
P{1, 1, {0, 0., -1}, {0.7818314824680299, 0., -0.6234898018587335}, 0.7853981633974483,  {0.39091574123401496, 0., -1.3117449009293667}},
P{1, 1, {0.0, 0., -0.2}, {-1.2246467991473532e-16, 0., -1.0}, 0.7853981633974483,  {-6.123233995736766e-17, 0., -0.7}},
P{1, 1, {0.0, 0., -0.2}, {-0.7818314824680299, 0., -0.6234898018587335}, 0.7667610911937136,  {-0.3927142092874629, 0., -0.5255266393634712}},
P{1, 1, {0.0, 0., -0.2}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.4824048601167302,  {-0.48548125428935995, 0., -0.25330251656977065}},
P{1, 1, {0.0, 0., -0.2}, {-0.4338837391175581, 0., 0.9009688679024191}, 0.1744219686871725,  {-0.18421034140937215, 0., -0.12935579518985824}},
P{1, 1, {0.0, 0., -0.2}, {0.4338837391175581, 0., 0.9009688679024191}, 0.1744219686871726,  {0.18421034140937217, 0., -0.12935579518985818}},
P{1, 1, {0.0, 0., -0.2}, {0.9749279121818236, 0., 0.22252093395631445}, 0.48240486011673017,  {0.4854812542893598, 0., -0.25330251656977065}},
P{1, 1, {0.0, 0., -0.2}, {0.7818314824680299, 0., -0.6234898018587335}, 0.7667610911937135,  {0.39271420928746287, 0., -0.5255266393634712}},
P{1, 1, {0.0, 0., 0.2}, {-1.2246467991473532e-16, 0., -1.0}, 0.6283185307179586,  {-7.960204194457795e-17, 0., -0.4}},
P{1, 1, {0.0, 0., 0.2}, {-0.7818314824680299, 0., -0.6234898018587335}, 0.5148249212496222,  {-0.4133211426933896, 0., -0.26414857114020945}},
P{1, 1, {0.0, 0., 0.2}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.10038946060967581,  {-0.43071079296204434, 0., -0.0851653619652063}},
P{1, 1, {0.0, 0., 0.2}, {-0.4338837391175581, 0., 0.9009688679024191}, 7.672960561261578e-05,  {-0.4377807498816102, 0., -0.00484969277273728}},
P{1, 1, {0.0, 0., 0.2}, {0.4338837391175581, 0., 0.9009688679024191}, 7.672960561261598e-05,  {0.4377807498816101, 0., -0.00484969277273728}},
P{1, 1, {0.0, 0., 0.2}, {0.9749279121818236, 0., 0.22252093395631445}, 0.10038946060967581,  {0.4307107929620443, 0., -0.0851653619652063}},
P{1, 1, {0.0, 0., 0.2}, {0.7818314824680299, 0., -0.6234898018587335}, 0.5148249212496222,  {0.4133211426933896, 0., -0.26414857114020945}},
P{1, 1, {0, 0., 1}, {-1.2246467991473532e-16, 0., -1.0}, 1.5308084989341912e-17,  {0.2945243112740431, 0., -2.7051619130462695e-17}},
P{1, 1, {0, 0., 1}, {-0.7818314824680299, 0., -0.6234898018587335}, 3.069737100296182e-05,  {-0.469838516420322, 0., -0.004119510363701007}},
P{1, 1, {0, 0., 1}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.0,  {-0.4874639560909118, 0., 1.1112604669781572}},
P{1, 1, {0, 0., 1}, {-0.4338837391175581, 0., 0.9009688679024191}, 0.0,  {-0.21694186955877906, 0., 1.4504844339512095}},
P{1, 1, {0, 0., 1}, {0.4338837391175581, 0., 0.9009688679024191}, 0.0,  {0.21694186955877906, 0., 1.4504844339512095}},
P{1, 1, {0, 0., 1}, {0.9749279121818236, 0., 0.22252093395631445}, 0.0,  {0.4874639560909118, 0., 1.1112604669781572}},
P{1, 1, {0, 0., 1}, {0.7818314824680299, 0., -0.6234898018587335}, 3.069737100296182e-05,  {0.4698385164203219, 0., -0.004119510363701007}},
P{1, 1, {0, 0., 2}, {-1.2246467991473532e-16, 0., -1.0}, 0.0,  {-6.123233995736766e-17, 0., 1.5}},
P{1, 1, {0, 0., 2}, {-0.7818314824680299, 0., -0.6234898018587335}, 0.0,  {-0.39091574123401496, 0., 1.6882550990706333}},
P{1, 1, {0, 0., 2}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.0,  {-0.4874639560909118, 0., 2.111260466978157}},
P{1, 1, {0, 0., 2}, {-0.4338837391175581, 0., 0.9009688679024191}, 0.0,  {-0.21694186955877906, 0., 2.4504844339512095}},
P{1, 1, {0, 0., 2}, {0.4338837391175581, 0., 0.9009688679024191}, 0.0,  {0.21694186955877906, 0., 2.4504844339512095}},
P{1, 1, {0, 0., 2}, {0.9749279121818236, 0., 0.22252093395631445}, 0.0,  {0.4874639560909118, 0., 2.111260466978157}},
P{1, 1, {0, 0., 2}, {0.7818314824680299, 0., -0.6234898018587335}, 0.0,  {0.39091574123401496, 0., 1.6882550990706333}},
P{2, 0.2, {0.0, 0., -0.1}, {-1.0, 0., 6.123233995736766e-17}, 0.06283185307179587,  {-1.0, 0., -0.09999999999999995}},
P{2, 0.2, {0.0, 0., -0.1}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.014118189231609806,  {-0.27654216509041113, 0., -0.06188105542439016}},
P{2, 0.2, {0.0, 0., -0.1}, {-0.9009688679024191, 0., 0.4338837391175582}, 0.007240632386867575,  {-0.1346690340511963, 0., -0.0601468112616171}},
P{2, 0.2, {0.0, 0., -0.1}, {-0.7818314824680298, 0., 0.6234898018587336}, 0.0050387233988818236,  {-0.08446582055505547, 0., -0.05764075583722696}},
P{2, 0.2, {0.0, 0., -0.1}, {-0.6234898018587335, 0., 0.7818314824680298}, 0.004018247824547358,  {-0.055935386256286655, 0., -0.05485924416277308}},
P{2, 0.2, {0.0, 0., -0.1}, {-0.4338837391175581, 0., 0.9009688679024191}, 0.0034869047816311984,  {-0.03498486044089578, 0., -0.052353188738382905}},
P{2, 0.2, {0.0, 0., -0.1}, {-0.2225209339563144, 0., 0.9749279121818236}, 0.003222384562320222,  {-0.01697699051886916, 0., -0.05061894457560985}},
P{2, 0.2, {0.0, 0., -0.2}, {-1.0, 0., 6.123233995736766e-17}, 0.06283185307179587,  {-1.0, 0., -0.19999999999999996}},
P{2, 0.2, {0.0, 0., -0.2}, {-0.9749279121818236, 0., 0.22252093395631445}, 0.028236378463219594,  {-0.4668675526103172, 0., -0.10594052771219511}},
P{2, 0.2, {0.0, 0., -0.2}, {-0.9009688679024191, 0., 0.4338837391175582}, 0.01448126477373515,  {-0.2230736217685235, 0., -0.10507340563080855}},
P{2, 0.2, {0.0, 0., -0.2}, {-0.7818314824680298, 0., 0.6234898018587336}, 0.01007744679776365,  {-0.13627993560223056, 0., -0.10382037791861351}},
P{2, 0.2, {0.0, 0., -0.2}, {-0.6234898018587335, 0., 0.7818314824680298}, 0.008036495649094717,  {-0.08777819729432366, 0., -0.10242962208138659}},
P{2, 0.2, {0.0, 0., -0.2}, {-0.4338837391175581, 0., 0.9009688679024191}, 0.006973809563262401,  {-0.05361052663101258, 0., -0.10117659436919148}},
P{2, 0.2, {0.0, 0., -0.2}, {-0.2225209339563144, 0., 0.9749279121818236}, 0.006444769124640448,  {-0.025606755838695854, 0., -0.10030947228780497}}
	);
	// clang-format on
	md::AnalyticalBuoyancyCalculator buoyancyCalculator(
	    p.length, p.diameter, p.bottom_pos, p.q);

	auto result = buoyancyCalculator.calculateBuoyancy();

	REQUIRE(isclose(result.wettedVolume, p.expectedV, 0, 1e-10));
	// In cases where the cylinder isn't submerged, the two versions may
	// disagree about COB location, this is fine
	if (p.expectedV > 1e-8) {
		REQUIRE_THAT(result.centerOfBuoyancy, IsClose(p.expectedCOB, 0, 1e-10));
	}
}