/*
 * Copyright (c) 2022 Jose Luis Cercos-Pita <jlc@core-marine.com>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/** @file rods.cpp
 * Minimal tests for rods in different situations
 */

#define _USE_MATH_DEFINES

#include "MoorDyn2.h"
#include "Log.hpp"
#include "Misc.hpp"
#include "Rod.hpp"
#include "Submergence.hpp"
#include "Waves.hpp"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>
#include "util.h"

#define TOL 1e-2

using namespace std;

bool
added_mass()
{
	// Do the math for our expected acceleration given the mass and that
	// Ca = 1.0
	double length = 2.0;
	double mass = 100 * length;
	double r = 0.25;
	double V = length * moordyn::pi * r * r;
	double rho_w = 1025;
	double g = 9.8;
	double buoyancy = rho_w * V * g;
	double weight = mass * g;
	double Fnet = buoyancy - weight;
	double added_mass = 1.0 * V * rho_w;
	double acceleration = Fnet / (mass + added_mass);

	int err;
	cout << endl << " => " << __PRETTY_FUNC_NAME__ << "..." << endl;

	MoorDyn system = MoorDyn_Create("Mooring/rod_tests/AddedMass.txt");
	if (!system) {
		cerr << "Failure Creating the Mooring system" << endl;
		return false;
	}

	err = MoorDyn_Init(system, NULL, NULL);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure during the mooring initialization: " << err << endl;
		return false;
	}

	double t = 0, Tmax = 2.5, dt = 0.1;

	while (t < Tmax) {
		err = MoorDyn_Step(system, NULL, NULL, NULL, &t, &dt);
		if (err != MOORDYN_SUCCESS) {
			cerr << "Failure during the mooring initialization: " << err
			     << endl;
			return false;
		}

		const auto rod = MoorDyn_GetRod(system, 1);

		unsigned int n_nodes;
		err = MoorDyn_GetRodN(rod, &n_nodes);
		if (err != MOORDYN_SUCCESS) {
			cerr << "Failure getting the number of nodes for rod " << 1 << ": "
			     << err << endl;
			return false;
		}

		double expected_z = -10 + 0.5 * acceleration * t * t;
		double pos[3];
		for (unsigned int i = 0; i < n_nodes; i++) {

			err = MoorDyn_GetRodNodePos(rod, i, pos);
			if (err != MOORDYN_SUCCESS) {
				cerr << "Failure getting the position of nodes " << i
				     << " for rod 1"
				     << ": " << err << endl;
				return false;
			}
			if (!isclose(pos[2], expected_z, 1e-8, 1e-10)) {
				cerr << "Node " << i << " of Rod 1 should have a z position of "
				     << expected_z << " but has a z pos of " << pos[2]
				     << " at t=" << t << endl;
				return false;
			}
		}
		// when the rod get higher than z = -1, stop simulating
		if (expected_z > -1.0) {
			break;
		}
	}

	err = MoorDyn_Close(system);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure closing Moordyn: " << err << endl;
		return false;
	}

	return true;
}

/** @brief Pinned rod horizontally lying, which should float towards the
 * vertical direction.
 *
 * Since there is no damping forces, the rods will be oscillating from one side
 * to the other
 */
bool
pinned_floating()
{
	int err;
	cout << endl << " => " << __PRETTY_FUNC_NAME__ << "..." << endl;

	MoorDyn system = MoorDyn_Create("Mooring/RodPinnedFloating.txt");
	if (!system) {
		cerr << "Failure Creating the Mooring system" << endl;
		return false;
	}

	err = MoorDyn_Init(system, NULL, NULL);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure during the mooring initialization: " << err << endl;
		return false;
	}

	double t = 0, T = 10.0;
	err = MoorDyn_Step(system, NULL, NULL, NULL, &t, &T);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure during the mooring initialization: " << err << endl;
		return false;
	}

	unsigned int n_rods;
	err = MoorDyn_GetNumberRods(system, &n_rods);
	for (unsigned int i_rod = 1; i_rod <= n_rods; i_rod++) {
		const auto rod = MoorDyn_GetRod(system, i_rod);
		if (!rod) {
			cerr << "Failure getting the rod " << i_rod << endl;
			return false;
		}

		unsigned int n_nodes;
		err = MoorDyn_GetRodN(rod, &n_nodes);
		if (err != MOORDYN_SUCCESS) {
			cerr << "Failure getting the number of nodes for rod " << i_rod
			     << ": " << err << endl;
			return false;
		}
		// Check that the rod is floating upwards
		double pos[3];
		err = MoorDyn_GetRodNodePos(rod, 0, pos);
		if (err != MOORDYN_SUCCESS) {
			cerr << "Failure getting first node position for rod " << i_rod
			     << ": " << err << endl;
			return false;
		}
		cout << pos[0] << ", " << pos[1] << ", " << pos[2] << endl;
		const double z0 = pos[2];
		err = MoorDyn_GetRodNodePos(rod, n_nodes, pos);
		if (err != MOORDYN_SUCCESS) {
			cerr << "Failure getting last node position for rod " << i_rod
			     << ": " << err << endl;
			return false;
		}
		const double z1 = pos[2];
		cout << pos[0] << ", " << pos[1] << ", " << pos[2] << endl;

		if (z1 < z0) {
			cerr << "The last node is below the first one for rod " << i_rod
			     << ": " << z1 << " vs. " << z0 << endl;
			return false;
		}
	}

	err = MoorDyn_Close(system);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure closing Moordyn: " << err << endl;
		return false;
	}

	return true;
}

/** @brief Rod hanging from 2 identical ropes, which should move horizontally
 */
bool
hanging()
{
	int err;
	cout << endl << " => " << __PRETTY_FUNC_NAME__ << "..." << endl;

	MoorDyn system = MoorDyn_Create("Mooring/RodHanging.txt");
	if (!system) {
		cerr << "Failure Creating the Mooring system" << endl;
		return false;
	}

	err = MoorDyn_Init(system, NULL, NULL);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure during the mooring initialization: " << err << endl;
		return false;
	}

	double t = 0, T = 10.0;
	err = MoorDyn_Step(system, NULL, NULL, NULL, &t, &T);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure during the mooring initialization: " << err << endl;
		return false;
	}

	const auto rod = MoorDyn_GetRod(system, 1);
	if (!rod) {
		cerr << "Failure getting the rod" << endl;
		return false;
	}

	unsigned int n_nodes;
	err = MoorDyn_GetRodN(rod, &n_nodes);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure getting the number of nodes: " << err << endl;
		return false;
	}
	// Check that the rod is is not rotating
	const double l = 20.0;
	double posa[3], posb[3];
	err = MoorDyn_GetRodNodePos(rod, 0, posa);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure getting first node position: " << err << endl;
		return false;
	}
	cout << posa[0] << ", " << posa[1] << ", " << posa[2] << endl;
	err = MoorDyn_GetRodNodePos(rod, n_nodes, posb);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure getting last node position: " << err << endl;
		return false;
	}
	cout << posb[0] << ", " << posb[1] << ", " << posb[2] << endl;

	if (((posa[0] + 0.5 * l) / l > TOL) || ((posb[0] - 0.5 * l) / l > TOL) ||
	    (fabs(posa[1]) / l > TOL) || (fabs(posb[1]) / l > TOL) ||
	    (fabs(posa[2] - posb[2]) / l > TOL)) {
		cerr << "The rod is rotating" << endl;
		return false;
	}

	err = MoorDyn_Close(system);
	if (err != MOORDYN_SUCCESS) {
		cerr << "Failure closing Moordyn: " << err << endl;
		return false;
	}

	return true;
}
template<typename T>
bool
valuesAreClose(const T& actual,
               const T& expected,
               const double rtol = 1e-8,
               const double atol = 1e-10)
{
	auto a_begin = std::cbegin(actual);
	const auto& a_end = std::cend(actual);
	auto e_begin = std::cbegin(expected);
	const auto& e_end = std::cend(expected);
	while (a_begin != a_end && e_begin != e_end) {
		if (!isclose(*a_begin, *e_begin, rtol, atol)) {
			return false;
		}
		++a_begin;
		++e_begin;
	}
	return true;
}
bool
buoyancyTest()
{
	const moordyn::real length = 2;
	const moordyn::real diameter = 0.5;
	const unsigned int NumSegs = 5;
	const double w = 10;
	const double totalMass = length * w;
	const double g = 9.8;
	const double rho_w = 1025;
	moordyn::Log* log = new moordyn::Log(MOORDYN_NO_OUTPUT, MOORDYN_NO_OUTPUT);
	EnvCondRef env = std::make_shared<EnvCond>();
	env->g = g;
	env->WtrDpth = 20.;
	env->rho_w = rho_w;
	env->kb = 3.0e6;
	env->cb = 3.0e5;
	env->waterKinOptions = moordyn::waves::WaterKinOptions();
	env->FrictionCoefficient = 0.0;
	env->FricDamp = 200.0;
	env->StatDynFricScale = 1.0;
	const moordyn::WavesRef waves = std::make_shared<moordyn::Waves>(log);
	waves->setup(env, nullptr, nullptr);

	RodProps rodProps{
		"",  diameter,
		w, // linear weight in air
		0.0, 0.0,      0.0, 0.0, 0.0, 0.0,
	};
	moordyn::Rod testRod(log, 0);

	testRod.setup(1,
	              moordyn::Rod::types::FREE,
	              &rodProps,
	              moordyn::vec6::Zero(),
	              NumSegs,
	              nullptr,
	              "");
	testRod.setEnv(env, waves, nullptr);
	waves->addRod(&testRod);

	moordyn::real phi = 0.0;
	const moordyn::real MAX_PHI = moordyn::pi;

	const moordyn::vec3 center_pos{ 0, 0, -0.1 };
	std::vector<moordyn::real> angles;
	std::vector<moordyn::vec6> actualFnets;
	std::vector<moordyn::vec6> actualAcc;
	std::vector<moordyn::vec6> expectedFnets;
	std::vector<moordyn::vec6> expectedAcc;

	for (; phi < MAX_PHI; phi += 0.02) {
		// std::cout << "phi = " << rad2deg * phi << std::endl;
		angles.push_back(phi);

		moordyn::vec3 q{ sin(phi), 0., cos(phi) };
		moordyn::vec6 endCoords{};
		endCoords.head<3>() = center_pos - ((length / 2.0) * q);
		endCoords.tail<3>() = center_pos + ((length / 2.0) * q);
		// std::cout << "endCoords = " << endCoords.transpose() << std::endl;
		testRod.setup(1,
		              moordyn::Rod::types::FREE,
		              &rodProps,
		              endCoords,
		              NumSegs,
		              nullptr,
		              "");
		testRod.initialize();
		auto [dpos, dvel] = testRod.getStateDeriv();

		moordyn::vec6 fnet = testRod.getFnet();

		actualFnets.push_back(fnet);
		actualAcc.push_back(dvel);
		Eigen::IOFormat testFmt(4, Eigen::DontAlignCols, ",", "\n", "[", "]");

		moordyn::vec3 rod_base = endCoords.head<3>();
		moordyn::AnalyticalBuoyancyCalculator buoyancyCalc(
		    length, diameter, rod_base, q);
		auto result = buoyancyCalc.calculateBuoyancy();
		// std::cout << "wetted volume = " << result.wettedVolume << std::endl;
		// std::cout << "COB = " << result.centerOfBuoyancy.transpose()
		//           << std::endl;

		moordyn::vec3 weight = moordyn::vec3::UnitZ() * totalMass * -g;
		moordyn::vec3 buoyancyForce =
		    moordyn::vec3::UnitZ() * rho_w * g * result.wettedVolume;
		moordyn::vec3 buoyancyMoment =
		    (result.centerOfBuoyancy - rod_base).cross(buoyancyForce);
		// the center of gravity also creates a moment because the torques are
		// about the end of the rod
		buoyancyMoment += (center_pos - rod_base).cross(weight);
		moordyn::vec6 expectedFnet{};
		expectedFnet.head<3>() = weight + buoyancyForce;
		expectedFnet.tail<3>() = buoyancyMoment;
		if (!valuesAreClose(fnet, expectedFnet)) {
			std::cout << "phi = " << moordyn::rad2deg * phi << std::endl;
			std::cout << "rod calculated fnet is "
			          << fnet.transpose().format(testFmt) << std::endl;
			std::cout << "rod expected fnet is "
			          << expectedFnet.transpose().format(testFmt) << std::endl;
			return false;
		}
		expectedFnets.push_back(expectedFnet);
		moordyn::vec6 centerAcc{};
		centerAcc.head<3>() = expectedFnet.head<3>() / totalMass;
		moordyn::real moment =
		    (1.0 / 12.0) * totalMass *
		    ((3.0 / 4.0) * diameter * diameter + length * length);
		// we can ignore the q axis moment because no torque there
		centerAcc.tail<3>() =
		    (result.centerOfBuoyancy - center_pos).cross(buoyancyForce) /
		    moment;
		moordyn::vec6 trueAcc = centerAcc;
		// because we are looking at the acceleration of the end of the rod, we
		// also must account for the acceleration due to rotation about the rod
		// center
		trueAcc.head<3>() += centerAcc.tail<3>().cross(-((length / 2.0) * q));

		if (!valuesAreClose(dvel, trueAcc)) {
			std::cout << "phi = " << moordyn::rad2deg * phi << std::endl;
			std::cout << "rod calculated acceleration is "
			          << dvel.transpose().format(testFmt) << std::endl;
			std::cout << "rod expected acc is "
			          << trueAcc.transpose().format(testFmt) << std::endl;
			return false;
		}
		expectedAcc.push_back(trueAcc);

		// std::cout << "expected Fnet is " << expectedFnet.transpose()
		//           << std::endl;
	}
	// Possible CSV output
	// std::cout <<
	// "phi,fx,fy,fz,tx,ty,tz,ax,ay,az,aax,aay,aaz,exp_fx,exp_fy,exp_"
	//              "fz,exp_tx,"
	//              "exp_ty,exp_tz,exp_ax,exp_ay,exp_az,exp_aax,exp_aay,exp_aaz"
	//           << std::endl;
	// for (int i = 0; i < angles.size(); i++) {
	// 	std::cout << angles[i] << ",";
	// 	auto&& actualF = actualFnets[i];
	// 	auto&& acc = actualAcc[i];
	// 	auto&& expectedF = expectedFnets[i];
	// 	auto&& trueAcc = expectedAcc[i];

	// 	Eigen::IOFormat testFmt(4, Eigen::DontAlignCols, ",", "", "", "");
	// 	std::cout << (actualF.transpose()).format(testFmt) << ",";
	// 	std::cout << (acc.transpose()).format(testFmt) << ",";
	// 	std::cout << (expectedF.transpose()).format(testFmt) << ",";
	// 	std::cout << (trueAcc.transpose()).format(testFmt) << std::endl;
	// }

	return true;
}
/** @brief Runs all the test
 * @return 0 if the tests have ran just fine. The index of the failing test
 * otherwise
 */
int
main(int, char**)
{
	if (!pinned_floating())
		return 1;
	if (!hanging())
		return 2;
	if (!added_mass())
		return 3;
	if (!buoyancyTest()) {
		return 4;
	}

	cout << "rods.cpp passed successfully" << endl;
	return 0;
}
