package stochprocess;

public class MarketVectorExpOUWithIndicators extends MarketVectorExpOU {

	public MarketVectorExpOUWithIndicators(float[] y) {
		super(y);
	}

	@Override
	public MarketVectorExpOUWithIndicators clone() {
		int simulationIndex = this.simulationIndex;
		float[] newY = this.cloneValues();

		MarketVectorExpOUWithIndicators newMarketInfo = new MarketVectorExpOUWithIndicators(newY);
		newMarketInfo.setSimulationIndex(simulationIndex);

		return newMarketInfo;
	}

	@Override
	public MarketVector getNewMarketVector(float[] y) {
		return new MarketVectorExpOUWithIndicators(y);
	}

}
