package policies;

import org.apache.commons.math3.util.FastMath;

import derivatives.Derivative;
import dynprog.State;
import stochprocess.MarketVector;

public class ZQF2006Hedging extends DeltaBoundHedging {

	final private double c0 = 1.12; // 1.08;
	final private double c1 = 0.31; // 0.31;
	final private double c2 = 0.25;
	final private double c3 = 0.5;
	final private double c4 = 4.76; // 6.85;
	final private double c5 = 0.78;
	final private double c6 = 0.15;

	private double dt;
	private double e_dt_sigma14;
	private double gamma;
	private double sigma;

	public ZQF2006Hedging() {
		this.setPolicyHedgingMethod(HedgingMethod.ZQF);
	}

	@Override
	public double computeDelta(State x) {
		int tIndex = x.timeStep;
		MarketVector marketInfo = x.marketVector;
		double s = marketInfo.values[MarketVector.PRICE];

		double h_s = this.c4 * (FastMath.pow(this.getbTC(), this.c5) * this.e_dt_sigma14)
				* FastMath.pow(this.getRiskAversionCoeff() * s * s * this.gamma, this.c6);
		double sigma_m = this.sigma * FastMath.sqrt(1 + h_s);

		double delta = -this.getDerivative().delta(s, sigma_m, tIndex);

		return delta;
	}

	@Override
	public double computeDeltaBound(State x) {
		MarketVector marketInfo = x.marketVector;
		double s = marketInfo.values[MarketVector.PRICE];
		double ra = this.getRiskAversionCoeff();

		double h_0 = this.getbTC() / (ra * s * this.sigma * this.sigma * this.dt);
		double h_w = this.c0 * (FastMath.pow(this.getbTC(), this.c1) * this.e_dt_sigma14)
				* FastMath.pow(this.gamma / ra, this.c3);

		double bound = h_w + h_0;

		return bound;
	}

	@Override
	public void precomputeCommonQuantities(State x) {
		int tIndex = x.timeStep;
		MarketVector marketVector = x.marketVector;

		this.dt = this.getMaturity() - this.getObservationTimes()[tIndex];
		this.sigma = marketVector.values[MarketVector.VOL];
		Derivative derivative = this.getDerivative();
		this.gamma = FastMath.abs(derivative.gamma(tIndex, marketVector.values));
		this.e_dt_sigma14 = FastMath.pow(FastMath.exp(-this.getRfRate() * this.dt) / this.sigma, this.c2);
	}

}