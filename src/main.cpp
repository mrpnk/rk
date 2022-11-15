#include "rungeKutta.h"
#include "odesystem.h"
#include "sensor.h"
#include <iostream>
#include <fstream>
#include <string>

// m*y'(t) + k*y(t) == f
struct TestSystem : public rk::FirstOrderSystem {
	using FirstOrderSystem::FirstOrderSystem;
	vec<1> operator()(double z, vec<1> const& y) const override {
		vec<1> y_t;
		y_t[0] = (params.f - params.k * y[0]) / params.m;
		return y_t;
	}
};

int main() {
	using namespace rk;

	logger<1> lg;


	double t0 = 0, t1 = 6;
	vec<1> initialState = {-1.};
	parameters params{ .m = 3, .k = 7, .f = 1 };


	TestSystem sys(params);
	rkIntegrate<FirstOrderSystem,vec<1>,rkdp54>(&sys, t0, initialState, t1, { &lg }, 1e0);


	std::ofstream file("rksol.txt");
	lg.print(file);
	file.close();

	return 0;
}
