package meshmethods;

import stochprocess.MarketVectorProcess;

public class MeshWeightsAverageDensity extends MeshWeights {

	public MeshWeightsAverageDensity() {
		super();
		this.setUsingNormalizedWeights(false);
	}

	@Override
	public void computeWeightsForIntermediateStage(int timeStep) {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		int nbPaths = stochasticMesh.getNbPaths();

		for (int i = 0; i < nbPaths; i++) {
			this.initializeWeightsListForStepAndOriginIndex(timeStep, i);
		}

		for (int destinationIndex = 0; destinationIndex < nbPaths; destinationIndex++) {
			this.copyStatePointToTempMarketVectorDestination(timeStep + 1, destinationIndex);
			this.computeWeightsForAllOrigins(timeStep, destinationIndex);
		}
	}

	@Override
	public String getDescription() {
		return "AD";
	}

	@Override
	public void precomputeInverseUnconditionalDensity(int k, MarketVectorProcess process) {

	}

	private void computeWeightsForAllOrigins(int timeStep, int destinationIndex) {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		int nbPaths = stochasticMesh.getNbPaths();
		float[] tempWeights = this.getTempWeights();
		float[] tempMarketVectorOrigin = stochasticMesh.getTempMarketVectorOrigin();
		double weightSum = 0.0;
		for (int originIndex = 0; originIndex < nbPaths; originIndex++) {
			this.copyStatePointToTempMarketVectorOrigin(timeStep, originIndex);

			float w = (float) this.computeConditionalDensity(timeStep, destinationIndex, tempMarketVectorOrigin);
			tempWeights[originIndex] = w;
			weightSum += w;
		}

		double g = 1 / weightSum;
		this.inverseUncondDensity[timeStep][destinationIndex] = g;

		for (int i = 0; i < nbPaths; i++) {
			float w = (float) (tempWeights[i] * g);
			this.updateWeightInformation(w, destinationIndex, timeStep, i);
		}
	}

}
