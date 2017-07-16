package policies;

import org.apache.commons.math3.util.FastMath;

import dynprog.State;
import stochprocess.MarketVector;

public class WWHedging extends DeltaBoundHedging {

	public WWHedging() {
		this.setPolicyHedgingMethod(HedgingMethod.WW);
	}

	@Override
	public double computeDelta(State x) {
		int tIndex = x.timeStep;
		MarketVector marketInfo = x.marketVector;
		double s = marketInfo.values[MarketVector.PRICE];
		double sigma = marketInfo.values[MarketVector.VOL];
		double delta = -this.getDerivative().delta(s, sigma, tIndex);

		return delta;
	}

	@Override
	public double computeDeltaBound(State x) {
		int tIndex = x.timeStep;
		MarketVector marketVector = x.marketVector;
		double s = marketVector.values[MarketVector.PRICE];
		double horizon = this.getMaturity() - this.getObservationTimes()[tIndex];
		double df = FastMath.exp(-this.getRfRate() * horizon);
		double gamma = this.getDerivative().gamma(tIndex, marketVector.values);
		double unitCosts = this.getaTC() + this.getbTC() * s;

		double bound = FastMath.pow(1.5 * (df * unitCosts * gamma * gamma) / this.getRiskAversionCoeff(), 1.0 / 3.0);

		return bound;
	}

	@Override
	public void precomputeCommonQuantities(State x) {
		// TODO Auto-generated method stub

	}

}