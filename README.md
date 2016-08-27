# Libsolver
A partial rewrite of [bluehood/diff_eq_solver](https://github.com/bluehood/diff_eq_solver),
a tiny C++ library for solving differential equations, in D.

Currently only includes the core solver itself, not the ROOT plotters.

## Building / Installing

Currently the library is a single file `solver.d`. To use, it:

- download the source  
- use `import solver` in your source  
- compile your files along with `solver.d`

## Usage

Sample usage:

```d
// example.d
import solver;
import std.stdio;

enum STEPS = 100;

void main() {
	// Define initial conditions vector [ x1, x2, ..., v1, v2, ... ]
	PosVec X0 = [0., 1.];

	// Define force vector. This is a vector of `real delegate(real, PosVec)`; 
	// the differential equation that will be solved is dX/dt = F(t, x(t))
	// In our example F = [ v, -x ] (simple harmonic oscillator with frequency=1)
	ForceVec double_osc_f = [
		(double t, PosVec x) => x[1],
		(double t, PosVec x) => -x[0]
	];

	// Solver(t0, X0, Force, step)
	// Use Euler's algorithm
	auto eulersolver = new EulerSolver(0., X0, double_osc_f, 0.1);
	// Use Runge-Kutta's algorithm (2nd order - 4th order is available as RK4Solver)
	auto rksolver = new RK2Solver(0., X0, double_osc_f, 0.1);
	
	// We want all output printed on file
	auto eulerOutput = File("euler_example.txt", "w");
	auto rkOutput = File("rk2_example.txt", "w");

	// Do `STEPS` steps
	for (int i = 0; i < STEPS; ++i) {
		PosVec X = eulersolver.step();
		real t = eulersolver.time;
		// Print current time
		eulerOutput.writef("%f ", t);
		// Print coordinates
		for (int j = 0; j < 2; ++j)
			eulerOutput.writef("%f ", X[j]);
		eulerOutput.writeln();
	}

	// Same for the Runge-Kutta solver.
	// We use a compact (although not very readable) way to do enough steps
	// to arrive to the desired time .
	// We use a different time step than the one specified in the constructor
	for (PosVec X = X0; rksolver.time < STEPS; X = rksolver.step(0.2))
		rkOutput.writefln("%f %f %f", rksolver.time, X[0], X[1]);
}
```

Then, compile it:  
```
dmd -release -O example.d /path/to/solver.d
```

## License
DLibsolver, like the C++ library, is licensed under the GNU GPL v2.
