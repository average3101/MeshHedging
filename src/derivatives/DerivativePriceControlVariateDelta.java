package derivatives;

import stochprocess.MarketVector;

public class DerivativePriceControlVariateDelta extends DerivativePriceControlVariate {
	private Derivative derivative;

	public DerivativePriceControlVariateDelta(Derivative derivative) {
		this.derivative = derivative;
	}

	@Override
	public double evaluate(int k0, float[] marketVector0, int k1, float[] marketVector1) {
		// double expectation =
		// this.marketProcess.expectedPriceConditionalOnMarketVector(k0,
		// marketVector0, k1,
		// marketVector1);
		return marketVector1[MarketVector.PRICE] - marketVector0[MarketVector.PRICE]; // -
																						// expectation;
																						// //
																						// marketVector0[MarketVector.PRICE];
	}

	@Override
	public double getExpectedValue(int k, float[] marketVector) {
		return 0;
	}

	@Override
	public double getMultiplicativeFactor(double s, double sigma, int k) {
		return this.derivative.delta(s, sigma, k);
	}

}
