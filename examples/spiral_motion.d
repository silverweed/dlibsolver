import solver;
import std.math;
import std.stdio;

enum TIME = 100000;

void main() {
	PosVec X0 = [ 1./sqrt(2.), 1./sqrt(2.), 0., 0., 0., 1. ];
	ForceVec double_osc_f = [
		(real t, PosVec X) => X[3],
		(real t, PosVec X) => X[4],
		(real t, PosVec X) => X[5],
		(real t, PosVec X) => -0.7*X[0],
		(real t, PosVec X) => -X[1],
		(real t, PosVec X) => 0
	];

	auto eulersolver = new EulerSolver(0.,X0,double_osc_f,0.1);
	auto rksolver = new RK2Solver(0.,X0,double_osc_f,0.1);
	
	auto eulerOutput = File("euler_spiral.txt", "w");
	auto rkOutput = File("rk2_spiral.txt", "w");
	for(PosVec X = X0; eulersolver.time < TIME; X = eulersolver.step())
	{
		eulerOutput.writef("%f ", eulersolver.time);
		for(int i=0; i<3; ++i)
			eulerOutput.writef("%f ", X[i]);
		eulerOutput.writeln();
	}
	for(PosVec X = X0; rksolver.time < TIME; X = rksolver.step())
	{
		rkOutput.writef("%f ", rksolver.time);
		for(int i=0; i<3; ++i)
			rkOutput.writef("%f ", X[i]);
		rkOutput.writeln();
	}
}
