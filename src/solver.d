module solver;

alias PosVec = real[];
alias Force = real delegate(real, PosVec);
alias ForceVec = Force[];

class Solver {
	this(real t0, in PosVec x0, ForceVec f, real step) {
		reset(t0, x0, f, step);
	}

	abstract PosVec step(real step);
	PosVec step() {
		return step(mstep);	
	}
	
	void reset(real t0, in PosVec x0, ForceVec f = mF, real step = mstep) 
	in {
		assert(x0.length > 0);
		assert(x0.length == f.length);
	} body {
		mt = t0;
		mX = x0.dup;
		mF = f;
		mstep = step;
	}

	@property real time() const { return mt; }
	@property const(PosVec) pos() const { return mX; }

	protected real mt, mstep;
	protected PosVec mX;
	protected ForceVec mF;
}

class EulerSolver : Solver {
	this(real t0, PosVec x0, ForceVec f, real step = 0) {
		super(t0, x0, f, step);
	}

	alias step = Solver.step;
	override PosVec step(real step) {
		PosVec x = mX;
		for (int i = 0; i < x.length; ++i)
			x[i] += step * mF[i](mt, mX);
		mX = x;
		mt += step;
		return mX;
	}
}

class RK2Solver : Solver {
	this(real t0, in PosVec x0, ForceVec f, real step = 0) {
		super(t0, x0, f, step);
	}

	alias step = Solver.step;
	override PosVec step(real step) {
		PosVec x = mX;
		auto dim = mX.length;

		auto k = new real[][2];
		for (int i = 0; i < k.length; ++i)
			k[i] = new real[](dim);

		auto t_coef = [0, 0.5];
		real[2][2] x_coefs;
		x_coefs[1][0] = 0.5;
		auto weights = [0, 1];
		
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < dim; ++j) {
				PosVec midPoint = mX.dup;
				for (int l = i - 1; l >= 0; --l)
					for (int m = 0; m < dim; ++m)
						midPoint[m] += x_coefs[i][l] * k[l][m];
				k[i][j] = step * mF[j](mt + t_coef[i] * step, midPoint);
			}
		}

		// evaluate new position
		for (int i = 0; i < dim; ++i)
			x[i] += weights[0] * k[0][i] + weights[1] * k[1][i];

		mX = x;
		mt += step;
		return mX;
	}
}

class RK4Solver : Solver {
	this(real t0, in PosVec x0, ForceVec f, real step = 0) {
		super(t0, x0, f, step);
	}

	alias step = Solver.step;
	override PosVec step(real step) {
		PosVec x = mX;
		auto dim = mX.length;


		auto k = new real[][4];
		for (int i = 0; i < k.length; ++i)
			k[i] = new real[](dim);
		auto t_coef = [0, 1/3., 2/3., 1];
		real[4][4] x_coefs;
		x_coefs[1][0] = 1./3.;
		x_coefs[2][1] = 1.;
		x_coefs[2][0] = -1./3.;
		x_coefs[3][2] = 1.;
		x_coefs[3][1] = -1.;
		x_coefs[3][0] = 1.;
		auto weights = [1/8., 3/8., 3/8., 1/8.];
		
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < dim; ++j) {
				PosVec midPoint = mX.dup;
				for (int l = i - 1; l >= 0; --l)
					for (int m = 0; m < dim; ++m)
						midPoint[m] += x_coefs[i][l] * k[l][m];
				k[i][j] = step * mF[j](mt + t_coef[i] * step, midPoint);
			}
		}

		// evaluate new position
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < dim; ++j)
				x[j] += weights[i] * k[i][j];

		mX = x;
		mt += step;
		return mX;
	}
}
