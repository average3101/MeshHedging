package stochprocess;

public abstract class MarketVector {
	public final static int PRICE = 0;
	public final static int VOL = 1;
	public final static int INDIC0 = 2;
	public final static int INDIC1 = 3;

	public int simulationIndex;
	public float[] values;

	public static void copy(MarketVector originMarketInfo, MarketVector destinationMarketInfo) {
		int length = originMarketInfo.values.length;
		for (int i = 0; i < length; i++) {
			destinationMarketInfo.values[i] = originMarketInfo.values[i];
		}
		destinationMarketInfo.simulationIndex = originMarketInfo.simulationIndex;
	}

	@Override
	public abstract MarketVector clone();

	public float[] cloneValues() {
		float[] values = this.getValues();
		int l = values.length;
		float[] newY = new float[l];
		for (int i = 0; i < l; i++) {
			newY[i] = values[i];
		}
		return newY;
	}

	public void copyValues(int index, float[] valuesToCopy) {
		this.simulationIndex = index;
		for (int i = 0; i < valuesToCopy.length; i++) {
			this.values[i] = valuesToCopy[i];
		}
	}

	public double getCurrentExpectedSecondMoment(int k) {
		return Double.NaN;
	}

	public abstract MarketVector getNewMarketVector(float[] y);

	public float[] getValues() {
		return this.values;
	}

	public void setSimulationIndex(int simulationIndex) {
		this.simulationIndex = simulationIndex;
	}

	public void setValue(int j, float v) {
		this.values[j] = v;
	}

	public void setValues(float[] y) {
		this.values = y;
	}
}