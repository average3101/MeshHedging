package numerics;

public class BisectionSolver extends OneDimSolver {
	// Assumes convexity of fct to optimize
	// Output of solver : [0] optimal objective, [1] minimizer, [2] nb fct evals

	protected double C = (3 - Math.sqrt(5)) / 2;

	protected double R = 1 - this.C;

	public BisectionSolver(OneDimFunction f) {
		this.f = f;
	}

	@Override
	public OneDimSolverResults minimize(double ax, double bx, double cx) {
		// MATLAB Algorithm for the Golden Section Search:
		// function [xmin, fmin] = golden(f, ax, bx, cx, tol)
		// %GOLDEN Minimize function of one variable using golden section search
		// %
		// % [xmin, fmin] = golden(f, ax, bx, cx, tol) computes a local minimum
		// % of f. xmin is the computed local minimizer of f and fmin is
		// % f(xmin). xmin is computed to an relative accuracy of TOL.
		// %
		// % The parameters ax, bx and cx must satisfy the following conditions:
		// % ax < bx < cx, f(bx) < f(ax) and f(bx) < f(cx).
		// %
		// % xmin satisfies ax < xmin < cx. golden is guaranteed to succeed if f
		// % is continuous between ax and cx
		// %
		// % Roman Geus, ETH Zuerich, 9.12.97

		double x1, x2;
		double x0 = ax;
		double x3 = cx;
		if (Math.abs(cx - bx) > Math.abs(bx - ax)) {
			x1 = bx;
			x2 = bx + this.C * (cx - bx);
		} else {
			x2 = bx;
			x1 = bx - this.C * (bx - ax);
		}

		double f1 = this.f.evaluate(x1);
		double f2 = this.f.evaluate(x2);

		int k = 1;
		while (Math.abs(x3 - x0) > this.tol * (Math.abs(x1) + Math.abs(x2))) {
			// fprintf(1,'k=%4d, |a-b|=%e\n', k, abs(x3-x0));
			if (f2 < f1) {
				x0 = x1;
				x1 = x2;
				x2 = this.R * x1 + this.C * x3; // x2 = x1+c*(x3-x1)
				f1 = f2;
				f2 = this.f.evaluate(x2);
			} else {
				x3 = x2;
				x2 = x1;
				x1 = this.R * x2 + this.C * x0; // % x1 = x2+c*(x0-x2)
				f2 = f1;
				f1 = this.f.evaluate(x1);
			}

			k++;
		}

		if (f1 < f2) {
			this.results.optim_f = f1;
			this.results.optim_v = x1;
		} else {
			this.results.optim_f = f2;
			this.results.optim_v = x2;
		}
		this.results.nbSteps = k;

		return this.results;
	}

	public void setFunctionToMinimize(OneDimFunction f) {
		this.f = f;
	}

}
