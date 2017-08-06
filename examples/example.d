import solver;
import std.stdio;

enum STEPS = 100;

void main() {
	// define initial conditions vector { x1, x2, ..., v1, v2, ...}
	PosVec X0 = [ 0., 1. ];
	// define force vector. this is a vector of double function(double, PosVec); 
	// the differential equation that will be solved is dX/dt = F(t,X(t))
	// In our example F = { v, -x } (simple harmonic oscillator with frequency=1)
	ForceVec double_osc_f = [
		(real t, PosVec x) => x[1],
		(real t, PosVec x) => -x[0]
	];

	// Solvers can be declared as instance of a derived class...
	// Solver(t0, X0, Force, step)
	auto eulersolver = new EulerSolver(0., X0, double_osc_f, 0.1);
	// ...or as pointers to the base class
	Solver rksolver = new RK2Solver(0., X0, double_osc_f, 0.1);
	
	// we want all output printed on file
	auto eulerOutput = File("euler_example.txt", "w");
	auto rkOutput = File("rk2_example.txt", "w");

	// do STEPS steps
	for (int i=0; i < STEPS; ++i) {
		PosVec X = eulersolver.step();
		double t = eulersolver.time;
		// print current time
		eulerOutput.writef("%f ", t);
		// print coordinates
		for (int j = 0; j < 2; ++j)
			eulerOutput.writef("%f ", X[j]);
		eulerOutput.writeln();
	}

	// same for the Runge-Kutta solver
	/* we use a compact (although not very readable) way to do enough steps
		to arrive to the desired time */
	// we use a different timestap than the one specified in the constructor
	for (PosVec X = X0; rksolver.time < STEPS; X = rksolver.step(0.2))
		rkOutput.writefln("%f %f %f", rksolver.time, X[0], X[1]);
}
