package stochprocess;

public class MarketVectorGBM extends MarketVector {

	public MarketVectorGBM(float[] y) {
		this.setValues(y);
	}

	public static void copy(MarketVector originMarketInfo, MarketVector destinationMarketInfo) {
		float[] originvalues = originMarketInfo.getValues();
		int length = originvalues.length;
		for (int i = 0; i < length; i++) {
			destinationMarketInfo.setValue(i, originvalues[i]);
		}
		destinationMarketInfo.setSimulationIndex(originMarketInfo.simulationIndex);
	}

	@Override
	public MarketVectorGBM clone() {
		int simulationIndex = this.simulationIndex;
		float[] newY = this.cloneValues();

		MarketVectorGBM newMarketInfo = new MarketVectorGBM(newY);
		newMarketInfo.setSimulationIndex(simulationIndex);

		return newMarketInfo;
	}

	@Override
	public MarketVector getNewMarketVector(float[] y) {
		return new MarketVectorGBM(y);
	}

}