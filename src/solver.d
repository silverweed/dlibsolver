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

		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < dim; ++j) {
				PosVec midPoint = mX.dup;
				for (int l = i - 1; l >= 0; --l)
					for (int m = 0; m < dim; ++m)
						midPoint[m] += X_COEFS[i][l] * k[l][m];
				k[i][j] = step * mF[j](mt + T_COEF[i] * step, midPoint);
			}
		}

		// evaluate new position
		for (int i = 0; i < dim; ++i)
			x[i] += WEIGHTS[0] * k[0][i] + WEIGHTS[1] * k[1][i];

		mX = x;
		mt += step;
		return mX;
	}

	private static immutable T_COEF = [0, 0.5];
	private static immutable WEIGHTS = [0, 1];
	private static immutable X_COEFS = [
		[0.,  0],
		[0.5, 0]
	];
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

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < dim; ++j) {
				PosVec midPoint = mX.dup;
				for (int l = i - 1; l >= 0; --l)
					for (int m = 0; m < dim; ++m)
						midPoint[m] += X_COEFS[i][l] * k[l][m];
				k[i][j] = step * mF[j](mt + T_COEF[i] * step, midPoint);
			}
		}

		// evaluate new position
		for(int i = 0; i < 4; ++i)
			for(int j = 0; j < dim; ++j)
				x[j] += WEIGHTS[i] * k[i][j];

		mX = x;
		mt += step;
		return mX;
	}

	private static immutable T_COEF = [0, 1/3., 2/3., 1];
	private static immutable WEIGHTS = [1/8., 3/8., 3/8., 1/8.];
	private static immutable X_COEFS = [
		[0.,    0, 0, 0],
		[1/3.,  0, 0, 0],
		[-1/3., 1, 0, 0],
		[1.,   -1, 1, 0]
	];
}
