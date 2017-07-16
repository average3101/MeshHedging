package derivatives;

import org.apache.commons.math3.util.FastMath;

import umontreal.ssj.probdist.NormalDist;

public class OptionPricing {
	// Créé par Pierre-Alexande Tremblay, 09-11-2004
	// Formule de Black-Scholes tirée de Hull 2e Edition, p.286-287

	static double sqrtTwoPi = FastMath.sqrt(2 * Math.PI);

	public static double binaryCallBSFast(double S, double X, double T, double sigma, double r, double q,
			double sigmaSqrtT, double exp_rT, double exp_qT) {
		if (T <= 0.0) {
			return FastMath.max(0.0, S - X);
		} else {
			double d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			double d2 = d1 - sigmaSqrtT;

			return exp_rT * NormalDist.cdf01(d2);
		}
	}

	public static double callBS(double S, double X, double T, double sigma, double r, double q) {
		if (T <= 0.0) {
			return FastMath.max(0.0, S - X);
		} else {
			double sigmaSqrtT, d1, d2;
			sigmaSqrtT = sigma * FastMath.sqrt(T);
			d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			d2 = d1 - sigmaSqrtT;

			return S * FastMath.exp(-q * T) * NormalDist.cdf01(d1) - X * FastMath.exp(-r * T) * NormalDist.cdf01(d2);
		}
	}

	public static double callBS(double S, double X, double T, double sigma, double r, double q, double sqrT,
			double exp_rT, double exp_qT) {
		if (T <= 0.0) {
			return FastMath.max(0.0, S - X);
		} else {
			double sigmaSqrtT = sigma * sqrT;
			double d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			double d2 = d1 - sigmaSqrtT;

			return S * exp_qT * NormalDist.cdf01(d1) - X * exp_rT * NormalDist.cdf01(d2);
		}
	}

	public static double callBSFast(double S, double X, double T, double sigma, double r, double q, double sigmaSqrtT,
			double exp_rT, double exp_qT) {
		if (T <= 0.0) {
			return FastMath.max(0.0, S - X);
		} else {
			double d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			double d2 = d1 - sigmaSqrtT;

			return S * exp_qT * NormalDist.cdf01(d1) - X * exp_rT * NormalDist.cdf01(d2);
		}
	}

	public static double callBSmx(double F, double X, double dtExp, double dtDel, double sigma, double r) {
		double sigmaSqrtTexp, d1, d2;
		sigmaSqrtTexp = sigma * FastMath.sqrt(dtExp);
		d1 = (Math.log(F / X) + (sigma * sigma / 2.0) * dtExp) / sigmaSqrtTexp;
		d2 = d1 - sigmaSqrtTexp;

		return FastMath.exp(-r * dtDel) * (F * normalCdfMx(d1) - X * normalCdfMx(d2));
	}

	public static double deltaBinaryCallBS(double S, double X, double T, double sigma, double r, double q) {
		return gammaBS(S, X, T, sigma, r, q);
	}

	public static double deltaCallBS(double S, double X, double T, double sigma, double r, double q) {
		double sigmaSqrtT, d1;
		sigmaSqrtT = sigma * FastMath.sqrt(T);
		d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;

		return FastMath.exp(-q * T) * NormalDist.cdf01(d1);
	}

	public static double deltaCallBS(double S, double X, double T, double sigma, double r, double q, double sqrT,
			double exp_qT) {
		double sigmaSqrtT = sigma * sqrT;
		double d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2) * T) / sigmaSqrtT;

		return exp_qT * NormalDist.cdf01(d1); // exp_qT = FastMath.exp(-q * T);
	}

	public static double deltaCallBSFast(double S, double X, double T, double sigma, double r, double q,
			double sigmaSqrtT, double exp_rT, double exp_qT) {
		double d1;
		sigmaSqrtT = sigma * FastMath.sqrt(T);
		d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;

		return FastMath.exp(-q * T) * NormalDist.cdf01(d1);
	}

	public static double deltaPutBS(double S, double X, double T, double sigma, double r, double q) {
		if (T <= 0.0) {
			if (S > X) {
				return 0.0;
			} else {
				return -1.0;
			}
		} else {
			double sigmaSqrtT, d1; // , d2;
			sigmaSqrtT = sigma * FastMath.sqrt(T);
			d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			// d2 = d1 - sigmaSqrtT;

			return -Math.exp(-q * T) * NormalDist.cdf01(-d1);
		}
	}

	public static double deltaPutBSFast(double S, double X, double T, double sigma, double r, double q,
			double sigmaSqrtT, double exp_rT, double exp_qT) {
		if (T <= 0.0) {
			if (S > X) {
				return 0.0;
			} else {
				return -1.0;
			}
		} else {
			double d1; // , d2;
			d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			// d2 = d1 - sigmaSqrtT;

			return -exp_qT * NormalDist.cdf01(-d1);
		}
	}

	public static double gammaBS(double S, double X, double T, double sigma, double r, double q) {
		double sigmaSqrtT, d1;
		sigmaSqrtT = sigma * FastMath.sqrt(T);
		d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;

		return FastMath.exp(-q * T) * NormalDist.density(0.0, 1.0, d1) / (S * sigmaSqrtT);
	}

	public static double gammaFDCallBS(double S, double X, double T, double sigma, double r, double q) {

		double precision = 0.001;
		return (callBS(S * (1.0 + precision), X, T, sigma, r, q) + callBS(S * (1.0 - precision), X, T, sigma, r, q)
				- 2.0 * callBS(S, X, T, sigma, r, q)) / Math.pow(precision * S, 2.0);
	}

	public static double impVolDeltaCallBS(double S, double X, double T, double delta, double r, double q,
			double tolerance) {
		double diff = 2 * tolerance;
		double sigmaL = 0.0;
		// First find upper bound.
		double sigmaU = 0.20;
		while (deltaCallBS(S, X, T, sigmaU, r, q) < delta) {
			sigmaU += 0.10;
		}
		// Bisection method
		double sigmaM = .10, deltaM;
		while (Math.abs(diff) > tolerance) {
			sigmaM = (sigmaU + sigmaL) / 2;
			deltaM = deltaCallBS(S, X, T, sigmaM, r, q);
			diff = deltaM - delta;
			if (diff > 0) {
				sigmaU = sigmaM;
			} else {
				sigmaL = sigmaM;
			}
		}
		return sigmaM;

	}

	public static void main(String[] args) {
	}

	public static double normalCdfMx(double x) {
		double y, z, nu;

		y = 1.0 / (1.0 + 0.2316419 * Math.abs(x));
		z = 1.330274429 * Math.pow(y, 5) - 1.821255978 * Math.pow(y, 4) + 1.781477937 * Math.pow(y, 3)
				- 0.356563782 * Math.pow(y, 2) + 0.31938153 * y;
		nu = 1.0 - z * FastMath.exp(-x * x / 2.0) / sqrtTwoPi;

		if (x > 0.0) {
			return nu;
		} else {
			return (1.0 - nu);
		}
	}

	public static double putBS(double S, double X, double T, double sigma, double r, double q) {
		if (T <= 0.0) {
			return FastMath.max(0.0, X - S);
		} else {
			double sigmaSqrtT, d1, d2;

			sigmaSqrtT = sigma * FastMath.sqrt(T);
			d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			d2 = d1 - sigmaSqrtT;

			return X * FastMath.exp(-r * T) * NormalDist.cdf01(-d2) - S * FastMath.exp(-q * T) * NormalDist.cdf01(-d1);
		}
	}

	public static double putBSFast(double S, double X, double T, double sigma, double r, double q, double sigmaSqrtT,
			double exp_rT, double exp_qT) {
		if (T <= 0.0) {
			return FastMath.max(0.0, X - S);
		} else {
			double d1, d2;

			d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			d2 = d1 - sigmaSqrtT;

			return X * exp_rT * NormalDist.cdf01(-d2) - S * exp_qT * NormalDist.cdf01(-d1);
		}
	}

	public static double putBSmx(double F, double X, double dtExp, double dtDel, double sigma, double r) {

		double sigmaSqrtTexp, d1, d2;

		sigmaSqrtTexp = sigma * FastMath.sqrt(dtExp);
		d1 = (Math.log(F / X) + (sigma * sigma / 2.0) * dtExp) / sigmaSqrtTexp;
		d2 = d1 - sigmaSqrtTexp;

		return FastMath.exp(-r * dtDel) * (F * (normalCdfMx(d1) - 1.0) - X * (normalCdfMx(d2) - 1.0));
	}

	public static double straddleBS(double S, double X, double T, double sigma, double r, double q) {
		if (T <= 0.0) {
			return FastMath.max(0.0, X - S) + FastMath.max(0.0, S - X);
		} else {

			double exp_rT = FastMath.exp(-r * T);
			double exp_qT = FastMath.exp(-q * T);

			double sigmaSqrtT = sigma * FastMath.sqrt(T);
			double d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;
			double d2 = d1 - sigmaSqrtT;
			double cdfd1 = NormalDist.cdf01(d1);
			double cdfd2 = NormalDist.cdf01(d2);

			double put = X * exp_rT * (1 - cdfd2) - S * exp_qT * (1 - cdfd1);
			double call = S * exp_qT * cdfd1 - X * exp_rT * cdfd2;

			return put + call;

		}
	}

	public static double vegaBS(double S, double X, double T, double sigma, double r, double q) {
		double sqrtT, sigmaSqrtT, d1;

		sqrtT = FastMath.sqrt(T);
		sigmaSqrtT = sigma * sqrtT;
		d1 = (Math.log(S / X) + (r - q + sigma * sigma / 2.0) * T) / sigmaSqrtT;

		return S * sqrtT * FastMath.exp(-d1 * d1 / 2.0 - q * T) / FastMath.sqrt(2 * Math.PI);
	}

	public static double vegaBSmx(double F, double X, double dtExp, double dtDel, double sigma, double r) {
		double sqrtTexp, sigmaSqrtTexp, d1;

		sqrtTexp = FastMath.sqrt(dtExp);
		sigmaSqrtTexp = sigma * sqrtTexp;
		d1 = (Math.log(F / X) + (sigma * sigma / 2.0) * dtExp) / sigmaSqrtTexp;

		return FastMath.exp(-r * dtDel) * F * sqrtTexp * FastMath.exp(-d1 * d1 / 2.0) / FastMath.sqrt(2 * Math.PI);
	}

}
