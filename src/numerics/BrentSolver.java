package numerics;

public class BrentSolver extends OneDimSolver {
	// Assumes convexity of fct to optimize
	// Output of solver : [0] optimal objective, [1] minimizer, [2] nb fct evals

	int maxIter = 15;

	public BrentSolver(OneDimFunction f) {
		this.f = f;
	}

	@Override
	public OneDimSolverResults minimize(double ax, double bx, double cx) {
		// PA Tremblay : The original method was modified in the following way :
		// 1) Extra condition to exit sooner when tolerance is reached (to avoid
		// further function evaluations)
		// 2) Extra condition to exit if max number of iterations is reached

		/**
		 *
		 * <p>
		 * This method performs a 1-dimensional minimization. It implements
		 * Brent's method which combines a golden-section search and parabolic
		 * interpolation. The introductory comments from the FORTRAN version are
		 * provided below.
		 *
		 * This method is a translation from FORTRAN to Java of the Netlib
		 * function fmin. In the Netlib listing no author is given.
		 *
		 * Translated by Steve Verrill, March 24, 1998.
		 *
		 * @param a
		 *            Left endpoint of initial interval
		 * @param b
		 *            Right endpoint of initial interval
		 * @param minclass
		 *            A class that defines a method, f_to_minimize, to minimize.
		 *            The class must implement the Fmin_methods interface (see
		 *            the definition in Fmin_methods.java). See FminTest.java
		 *            for an example of such a class. f_to_minimize must have
		 *            one double valued argument.
		 * @param tol
		 *            Desired length of the interval in which the minimum will
		 *            be determined to lie (This should be greater than,
		 *            roughly, 3.0e-8.)
		 *
		 */

		/*
		 *
		 * Here is a copy of the Netlib documentation:
		 *
		 * c c An approximation x to the point where f attains a minimum on c
		 * the interval (ax,bx) is determined. c c input.. c c ax left endpoint
		 * of initial interval c bx right endpoint of initial interval c f
		 * function subprogram which evaluates f(x) for any x c in the interval
		 * (ax,bx) c tol desired length of the interval of uncertainty of the
		 * final c result (.ge.0.) c c output.. c c fmin abcissa approximating
		 * the point where f attains a c minimum c c The method used is a
		 * combination of golden section search and c successive parabolic
		 * interpolation. Convergence is never much slower c than that for a
		 * Fibonacci search. If f has a continuous second c derivative which is
		 * positive at the minimum (which is not at ax or c bx), then
		 * convergence is superlinear, and usually of the order of c about
		 * 1.324. c The function f is never evaluated at two points closer
		 * together c than eps*abs(fmin)+(tol/3), where eps is approximately the
		 * square c root of the relative machine precision. If f is a unimodal c
		 * function and the computed values of f are always unimodal when c
		 * separated by at least eps*abs(x)+(tol/3), then fmin approximates c
		 * the abcissa of the global minimum of f on the interval (ax,bx) with c
		 * an error less than 3*eps*abs(fmin)+tol. If f is not unimodal, c then
		 * fmin may approximate a local, but perhaps non-global, minimum to c
		 * the same accuracy. c This function subprogram is a slightly modified
		 * version of the c Algol 60 procedure localmin given in Richard Brent,
		 * Algorithms For c Minimization Without Derivatives, Prentice-Hall,
		 * Inc. (1973). c
		 *
		 */

		double c, d, e, eps, xm, p, q, r, tol1, t2, u, v, w, fu, fv, fw, fx, x, tol3;

		double a = ax;
		double b = cx;

		c = .5 * (3.0 - Math.sqrt(5.0));
		d = 0.0;

		// 1.1102e-16 is machine precision

		eps = 1.2e-16;
		tol1 = eps + 1.0;
		eps = Math.sqrt(eps);

		v = a + c * (b - a);
		w = v;
		x = v;
		e = 0.0;
		fx = this.f.evaluate(x);
		fv = fx;
		fw = fx;
		tol3 = this.tol / 3.0;

		xm = .5 * (a + b);
		tol1 = eps * Math.abs(x) + tol3;
		t2 = 2.0 * tol1;

		// main loop
		int k = 1;
		while (Math.abs(x - xm) > (t2 - .5 * (b - a))) {

			p = q = r = 0.0;

			if (Math.abs(e) > tol1) {

				// fit the parabola

				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);

				if (q > 0.0) {
					p = -p;
				} else {
					q = -q;
				}

				r = e;
				e = d;

				// brace below corresponds to statement 50
			}

			if ((Math.abs(p) < Math.abs(.5 * q * r)) && (p > q * (a - x)) && (p < q * (b - x))) {

				// a parabolic interpolation step

				d = p / q;
				u = x + d;

				// f must not be evaluated too close to a or b

				if (((u - a) < t2) || ((b - u) < t2)) {

					d = tol1;
					if (x >= xm) {
						d = -d;
					}

				}

				// brace below corresponds to statement 60
			} else {

				// a golden-section step

				if (x < xm) {
					e = b - x;
				} else {
					e = a - x;
				}
				d = c * e;
			}

			// f must not be evaluated too close to x

			if (Math.abs(d) >= tol1) {
				u = x + d;
			} else {

				if (d > 0.0) {
					u = x + tol1;
				} else {
					u = x - tol1;
				}

			}

			fu = this.f.evaluate(u);

			// Update a, b, v, w, and x

			if (fx <= fu) {

				if (u < x) {
					a = u;
				} else {
					b = u;
				}

				// brace below corresponds to statement 140
			}

			if (fu <= fx) {

				// PA Tremblay : Condition added to exit ASAP.
				if (Math.abs(x - u) < this.tol) {
					x = u;
					fx = fu;
					break;
				}

				if (u < x) {
					b = x;
				} else {
					a = x;
				}

				v = w;
				fv = fw;
				w = x;
				fw = fx;
				x = u;
				fx = fu;

				xm = .5 * (a + b);
				tol1 = eps * Math.abs(x) + tol3;
				t2 = 2.0 * tol1;

				// brace below corresponds to statement 170
			} else {

				if ((fu <= fw) || (w == x)) {

					v = w;
					fv = fw;
					w = u;
					fw = fu;

					xm = .5 * (a + b);
					tol1 = eps * Math.abs(x) + tol3;
					t2 = 2.0 * tol1;

				} else if ((fu > fv) && (v != x) && (v != w)) {

					xm = .5 * (a + b);
					tol1 = eps * Math.abs(x) + tol3;
					t2 = 2.0 * tol1;

				} else {

					v = u;
					fv = fu;

					xm = .5 * (a + b);
					tol1 = eps * Math.abs(x) + tol3;
					t2 = 2.0 * tol1;

				}

			}
			k++;

			// PA Tremblay : Condition added to exit if max number of iterations
			// reached
			if (k >= this.maxIter) {
				break;
			}

		}

		this.results.optim_v = x;
		this.results.optim_f = fx;
		this.results.nbSteps = k;

		return this.results;
	}
}
