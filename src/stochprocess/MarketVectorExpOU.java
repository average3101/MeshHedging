package stochprocess;

public class MarketVectorExpOU extends MarketVector {

	public MarketVectorExpOU(float[] y) {
		this.setValues(y);
	}

	@Override
	public MarketVectorExpOU clone() {

		int simulationIndex = this.simulationIndex;
		float[] newY = this.cloneValues();

		MarketVectorExpOU newMarketInfo = new MarketVectorExpOU(newY);
		newMarketInfo.setSimulationIndex(simulationIndex);

		return newMarketInfo;
	}

	@Override
	public MarketVector getNewMarketVector(float[] y) {
		return new MarketVectorExpOU(y);
	}

}