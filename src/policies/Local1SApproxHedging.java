package policies;

import dynprog.State;
import stochprocess.MarketVector;

public class Local1SApproxHedging extends DeltaBoundHedging {

	private double phi;
	private double variance;
	private double z;

	public Local1SApproxHedging() {
		this.setPolicyHedgingMethod(HedgingMethod.LOC_1SA);
	}

	@Override
	public double computeDelta(State x) {
		int tIndex = x.timeStep;
		MarketVector marketInfo = x.marketVector;
		double u0 = x.stockQuantity;

		double s = marketInfo.values[MarketVector.PRICE];
		double sigma = marketInfo.values[MarketVector.VOL];
		double delta = -this.getDerivative().delta(s, sigma, tIndex);
		double d = (delta + this.z * u0) / (1.0 + this.z);

		return d;
	}

	@Override
	public double computeDeltaBound(State x) {
		double bound = 1.0 / ((1.0 + this.z) * this.phi * this.getRiskAversionCoeff());

		return bound;
	}

	@Override
	public void precomputeCommonQuantities(State x) {
		MarketVector marketInfo = x.marketVector;
		int k = x.timeStep;
		double s = marketInfo.values[MarketVector.PRICE];

		this.phi = this.getaTC() + this.getbTC() * s;
		this.variance = marketInfo.getCurrentExpectedSecondMoment(k);
		this.z = this.phi * this.phi / this.variance;
	}
}
