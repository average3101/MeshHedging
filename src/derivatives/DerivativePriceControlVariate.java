package derivatives;

import stochprocess.MarketVectorProcess;

public abstract class DerivativePriceControlVariate {

	public enum BetaType {
		ANALYTIC, ESTIMATED
	}

	public enum Label {
		COND_EXPECT, DELTA
	};

	public static DerivativePriceControlVariate createControl(Label label, Derivative derivative) {
		DerivativePriceControlVariate control = null;
		switch (label) {
		case COND_EXPECT:
			control = new DerivativePriceControlVariateCondExpect(derivative);
			break;
		case DELTA:
			control = new DerivativePriceControlVariateDelta(derivative);
			break;
		}
		return control;
	}

	protected MarketVectorProcess marketProcess;

	public abstract double evaluate(int k0, float[] marketVector0, int k1, float[] marketVector1);

	public abstract double getExpectedValue(int k, float[] marketVector);

	public double getMultiplicativeFactor(double s, double sigma, int k) {
		return 1.0;
	}

}
