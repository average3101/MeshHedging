package derivatives;

import stochprocess.MarketVector;

public class DerivativePriceControlVariateCondExpect extends DerivativePriceControlVariate {
	private Derivative derivative;

	public DerivativePriceControlVariateCondExpect(Derivative derivative) {
		this.derivative = derivative;
	}

	@Override
	public double evaluate(int k0, float[] marketVector0, int k1, float[] marketVector1) {
		float sigma0 = marketVector0[MarketVector.VOL];
		// float s0 =
		// this.marketProcess.expectedPriceConditionalOnMarketVector(k0,
		// marketVector0, k1, marketVector1);
		// float sigma0 =
		// this.marketProcess.expectedVolatilityConditionalOnMarketVector(k0,
		// marketVector0, k1,
		// marketVector1);
		float sigma1Backup = marketVector1[MarketVector.VOL];
		// float s1Backup = marketVector1[MarketVector.PRICE];

		marketVector1[MarketVector.VOL] = sigma0;
		// marketVector1[MarketVector.PRICE] = s0;
		double p = this.derivative.evaluate(k1, marketVector1);

		marketVector1[MarketVector.VOL] = sigma1Backup;
		// marketVector1[MarketVector.PRICE] = s1Backup;

		return p;
	}

	@Override
	public double getExpectedValue(int k, float[] marketVector) {
		//
		// float s0 =
		// this.marketProcess.expectedPriceConditionalOnMarketVector(k,
		// marketVector, -1, null);
		// float sigma0 =
		// this.marketProcess.expectedVolatilityConditionalOnMarketVector(k,
		// marketVector, -1, null);
		// float sigmaBackup = marketVector[MarketVector.VOL];
		// float sBackup = marketVector[MarketVector.PRICE];
		//
		// marketVector[MarketVector.VOL] = sigma0;
		// marketVector[MarketVector.PRICE] = s0;
		double p = this.derivative.evaluate(k, marketVector);

		// marketVector[MarketVector.VOL] = sigmaBackup;
		// marketVector[MarketVector.PRICE] = sBackup;

		return p;
	}

}
