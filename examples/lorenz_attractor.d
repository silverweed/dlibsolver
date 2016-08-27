import solver;
import std.stdio;

enum TIME = 1000;

void main() {
	immutable a = 10., b = 28., c = 8./3.;
	PosVec X0 = [-0.1, -0.01, .2];
	ForceVec lorenz_attractor = [
		(real t, PosVec X) => -a * X[0] + a * X[1],
		(real t, PosVec X) => X[0] * X[2]-X[1] + b * X[0],
		(real t, PosVec X) => -X[0] * X[1] - c * X[2]
	];

	auto rksolver = new RK4Solver(0., X0, lorenz_attractor, 0.001);
	
	auto output = File("lorenz.txt", "w");

	for (PosVec X = X0; rksolver.time < TIME; X = rksolver.step()) {
		output.writefln("%f %f %f", X[0], X[1], X[2]);
	}
}
