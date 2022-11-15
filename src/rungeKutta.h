#pragma once
#include <fstream>
#include <vector>
#include <array>
#include <cfenv>
#include <iostream>

namespace rk
{
	template<typename State>
	struct sensor {
		virtual bool measure(double t, State const& y, State const& deriv, State const& erry,
		                     int step, double dt, double& forceDt) = 0; // return true to stop integration
	};

	struct butcher7TableauFSAL {
		const double c2, a21;
		const double c3, a31, a32;
		const double c4, a41, a42, a43;
		const double c5, a51, a52, a53, a54;
		const double c6, a61, a62, a63, a64, a65;
		const double c7, a71, a72, a73, a74, a75, a76;
		const double     e1, e2, e3, e4, e5, e6, e7;
		const int stepper_order, error_order;
	};
	constexpr butcher7TableauFSAL rkdp54 = { // Dormand–Prince 5(4)
		1. / 5,       1. / 5,
		3. / 10,      3. / 40,         9. / 40,
		4. / 5,      44. / 45,       -56. / 15,      32. / 9,
		8. / 9,   19372. / 6561,  -25360. / 2187, 64448. / 6561,   -212. / 729,
		1,         9017. / 3168,    -355. / 33,   46732. / 5247,     49. / 176,   -5103. / 18656,
		1,           35. / 384,          0,         500. / 1113,    125. / 192,   -2187. / 6784,    11. / 84,
				     71. / 57600,        0,         -71. / 16695,    71. / 1920, -17253. / 339200,  22. / 525,   -1. / 40,
		5,4
	};

	struct butcher8TableauFSAL {
		const double c2, a21;
		const double c3, a31, a32;
		const double c4, a41, a42, a43;
		const double c5, a51, a52, a53, a54;
		const double c6, a61, a62, a63, a64, a65;
		const double c7, a71, a72, a73, a74, a75, a76;
		const double c8, a81, a82, a83, a84, a85, a86, a87;
		const double     e1, e2, e3, e4, e5, e6, e7, e8;
		const int stepper_order, error_order;
	};
	constexpr butcher8TableauFSAL rkbs54 = { // Bogacki-Shampine 5(4)
		1. / 6,       1. / 6,
		2. / 9,       2. / 27,          4. / 27,
		3. / 7,     183. / 1372,     -162. / 343,       1053. / 1372,
		2. / 3,      68. / 297,        -4. / 11,          42. / 143,           1960. / 3861,
		3. / 4,     597. / 22528,      81. / 352,      63099. / 585728,       58652. / 366080,     4617. / 20480,
		1,       174197. / 959244, -30942. / 79937,  8152137. / 19744439,    666106. / 1039181,  -29421. / 29068, 482048. / 414219,
		1,          587. / 8064,          0,         4440339. / 15491840,     24353. / 124800,      387. / 44800,   2152. / 5985,        7267. / 94080,
				   3817. / 1959552,       0,         -140181. / 15491840,   4224731. / 272937600, -8557. / 403200, 57928. / 4363065, 23930231. / 4366535040, -3293. / 556956,
		5,4
	};

	template<typename System, typename State, butcher7TableauFSAL tbl>
	void rkStep(System const& sys, double t, State const& y, double dt, State& out_y, State& out_erry, State const& k1, State& out_k7) {
		State k2, k3, k4, k5, k6;
		k2 = sys(t + tbl.c2 * dt, y + (k1 * tbl.a21) * dt);
		k3 = sys(t + tbl.c3 * dt, y + (k1 * tbl.a31 + k2 * tbl.a32) * dt);
		k4 = sys(t + tbl.c4 * dt, y + (k1 * tbl.a41 + k2 * tbl.a42 + k3 * tbl.a43) * dt);
		k5 = sys(t + tbl.c5 * dt, y + (k1 * tbl.a51 + k2 * tbl.a52 + k3 * tbl.a53 + k4 * tbl.a54) * dt);
		k6 = sys(t + tbl.c6 * dt, y + (k1 * tbl.a61 + k2 * tbl.a62 + k3 * tbl.a63 + k4 * tbl.a64 + k5 * tbl.a65) * dt);
		out_y = y + (k1 * tbl.a71 + k2 * tbl.a72 + k3 * tbl.a73 + k4 * tbl.a74 + k5 * tbl.a75 + k6 * tbl.a76) * dt;
		out_k7 = sys(t + tbl.c7 * dt, out_y); // FSAL property allows a reuse in the next step
		out_erry = (k1 * tbl.e1 + k2 * tbl.e2 + k3 * tbl.e3 + k4 * tbl.e4 + k5 * tbl.e5 + k6 * tbl.e6 + out_k7 * tbl.e7) * dt;
	}

	template<typename System, typename State, butcher8TableauFSAL tbl>
	void rkStep(System const& sys, double t, State const& y, double dt, State& out_y, State& out_erry, State const& k1, State& out_k8) {
		State k2, k3, k4, k5, k6, k7;
		k2 = sys(t + tbl.c2 * dt, y + (k1 * tbl.a21) * dt);
		k3 = sys(t + tbl.c3 * dt, y + (k1 * tbl.a31 + k2 * tbl.a32) * dt);
		k4 = sys(t + tbl.c4 * dt, y + (k1 * tbl.a41 + k2 * tbl.a42 + k3 * tbl.a43) * dt);
		k5 = sys(t + tbl.c5 * dt, y + (k1 * tbl.a51 + k2 * tbl.a52 + k3 * tbl.a53 + k4 * tbl.a54) * dt);
		k6 = sys(t + tbl.c6 * dt, y + (k1 * tbl.a61 + k2 * tbl.a62 + k3 * tbl.a63 + k4 * tbl.a64 + k5 * tbl.a65) * dt);
		k7 = sys(t + tbl.c7 * dt, y + (k1 * tbl.a71 + k2 * tbl.a72 + k3 * tbl.a73 + k4 * tbl.a74 + k5 * tbl.a75 + k6 * tbl.a76) * dt);
		out_y = y + (k1 * tbl.a81 + k2 * tbl.a82 + k3 * tbl.a83 + k4 * tbl.a84 + k5 * tbl.a85 + k6 * tbl.a86 + k7 * tbl.a87) * dt;
		out_k8 = sys(t + tbl.c8 * dt, out_y); // FSAL property allows a reuse in the next step
		out_erry = (k1 * tbl.e1 + k2 * tbl.e2 + k3 * tbl.e3 + k4 * tbl.e4 + k5 * tbl.e5 + k6 * tbl.e6 + k7 * tbl.e7 + out_k8 * tbl.e8) * dt;
	}

	struct integrationInfo {
		int stopReason;
		double tFrom, tTo;
		int nsteps;
	};

	template<typename System, typename State, auto tbl>
	integrationInfo rkIntegrate(System const* sys, double t0, State y0, double t1,
		std::vector<sensor<State>*> const& sensors, const double maxdt) {

		int maxSteps = 20'000;
		const double startStepSize = 1e-8;
		const double absTol = 1e-8;
		const double relTol = 1e-8;
		const double stepSizeSafetyFactor = 0.9;
		const double stepSizeRatioBound = 5;

		double t = t0, dt = startStepSize, forcedt = std::numeric_limits<double>::max();
		State y = y0, newY, erry = State{ 0 };
		State deriv = (*sys)(t, y), newDeriv;
		double relErrNorm;
		int nSteps = 0;

		for (int i = 0; i < sensors.size(); ++i) {
			sensors[i]->measure(t, y, deriv, erry, nSteps, dt, forcedt);
		}
		while (t < t1) {
			if (dt > maxdt)   dt = maxdt;
			if (dt > forcedt) dt = forcedt;
			if (t + dt > t1)  dt = t1 - t;

			while (true) {
				rkStep<System, State, tbl>(*sys, t, y, dt, newY, erry, deriv, newDeriv);

				relErrNorm = reduceMax(abs(erry) / ((abs(y) + deriv * dt) * relTol + absTol));

				if (relErrNorm <= 1)
					break;

				// decrease stepsize
				dt *= std::max(stepSizeSafetyFactor * pow(relErrNorm, -1. / (tbl.error_order - 1.)), 1. / tbl.stepper_order);
			};

			if (newY.isbad())
				return { -2,t0,t,nSteps };

			y = newY;
			deriv = newDeriv;
			t += dt;
			nSteps++;

			// increase stepsize
			if (relErrNorm < 0.5) {
				dt *= stepSizeSafetyFactor * std::max(pow(relErrNorm, -1. / tbl.stepper_order), stepSizeRatioBound);
			}

			forcedt = std::numeric_limits<double>::max();
			int sensorbreak = -1;
			for (int i = 0; i < sensors.size(); ++i) {
				if (sensors[i]->measure(t, y, deriv, erry, nSteps, dt, forcedt)) sensorbreak = i;
			}
			if (sensorbreak != -1) return { sensorbreak + 1,t0,t,nSteps };

			if (nSteps == maxSteps)
				return { -1,t0,t,nSteps };
		}
		return { 0,t0,t,nSteps };
	}

}
