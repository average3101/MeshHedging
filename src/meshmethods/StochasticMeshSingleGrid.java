package meshmethods;

import java.util.List;

import cern.colt.list.FloatArrayList;
import cern.colt.list.IntArrayList;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.rng.RandomStream;

public class StochasticMeshSingleGrid extends StochasticMesh {
	private MarketVectorProcess processForGeneratingOnlyOneStep;
	private MarketVectorProcess processForGeneratingOnlyOneStepForPricing;

	public StochasticMeshSingleGrid(MarketVectorProcess multiSp) {
		super(multiSp);
		this.processForGeneratingOnlyOneStep = multiSp.getProcessForGeneratingOnlyOneStep();
		this.setBasicLabel("SG");
	}

	@Override
	public void computeAllWeights(MeshWeights weights, MarketVectorProcess process) {
		weights.init();

		int nbSteps = this.getNbSteps();
		weights.precomputeInverseUnconditionalDensity(0, process);
		this.assignIdenticalUnconditionalDensityToAllPeriods(weights);
		if (nbSteps > 1) {
			int timeStep = nbSteps - 1;
			weights.computeWeightsForIntermediateStage(timeStep);
			this.assignIdenticalWeightsToAllPeriods(weights);
		}
		weights.computeWeightsForIntermediateStage(0);

		this.applyAdditionalTransformations(weights);
	}

	@Override
	public void defineDefaultMeshWeights() {
		super.defineDefaultMeshWeights();
		MeshWeights weights = this.getMeshWeights();
		weights.setRemoveTransitionsToIdenticalIndex(true);
		weights.setUsingIdenticalWeightsAtFirstStep(false);
	}

	public MarketVectorProcess getProcessForGeneratingOnlyOneStep() {
		return this.processForGeneratingOnlyOneStep;
	}

	@Override
	public void resetNextRandomSubstream() {
		RandomStream stream = this.processForGeneratingOnlyOneStep.getStream();
		this.setRandomStreamForMeshStateGeneration(stream);
		stream.resetNextSubstream();
	}

	@Override
	public void setMarketVectorProcessForPricing(MarketVectorProcess processForPricing) {
		if (!processForPricing.getClass().equals(this.marketProcess.getClass())) { // processForGeneratingOnlyOneStep.getClass()))
																					// {
			this.setUsingDifferentProcessForPricing(true);

			// this.setProcessForGeneratingOnlyOneStepForPricing(oneStepProcess);
			this.marketProcessForPricing = processForPricing;
			this.processForGeneratingOnlyOneStepForPricing = processForPricing.getProcessForGeneratingOnlyOneStep();

			MeshWeights weights = this.newDefaultMeshWeights();
			this.setMeshWeightsForDerivative(weights);
		}
	}

	@Override
	public void setObservationTimes(double[] obsTimes, int nbSteps) {
		this.marketProcess.setObservationTimes(obsTimes, nbSteps);

		this.setNbSteps(nbSteps);
		double dt = obsTimes[nbSteps] - obsTimes[0];
		this.processForGeneratingOnlyOneStep.setObservationTimes(dt, 1);

		if (this.isUsingDifferentProcessForPricing()) {
			this.marketProcessForPricing.setObservationTimes(obsTimes, nbSteps);
			this.processForGeneratingOnlyOneStepForPricing.setObservationTimes(dt, 1);
		}
	}

	@Override
	public void simulatePaths() {
		this.resetNextRandomSubstream();
		int nbPaths = this.getNbPaths();
		int marketVectorDimension = this.marketProcess.getMarketVectorDimension();
		float[][] gridStates = new float[nbPaths][marketVectorDimension];
		float[][] identicalStartingStates = new float[nbPaths][marketVectorDimension];

		int k0 = 0;
		int k1 = 1; // Only one step
		for (int i = 0; i < nbPaths; i++) {
			float[][] sPathArray = this.processForGeneratingOnlyOneStep.generateBasicProcessesPaths();

			for (int j = 0; j < marketVectorDimension; j++) {
				gridStates[i][j] = sPathArray[j][k1];
				identicalStartingStates[i][j] = sPathArray[j][k0];
			}

			this.resetNextRandomSubstream();
		}

		this.updateStatePointsList(gridStates, identicalStartingStates);
	}

	private void assignIdenticalUnconditionalDensityToAllPeriods(MeshWeights meshWeights) {
		int nbPaths = this.getNbPaths();
		int nbSteps = this.getNbSteps();
		double[][] inverseUncondDensity = meshWeights.getInverseUnconditionalDensity();
		for (int k = 1; k < nbSteps; k++) {
			for (int j = 0; j < nbPaths; j++) {
				inverseUncondDensity[k][j] = inverseUncondDensity[0][j];
			}
		}
	}

	private void assignIdenticalWeightsToAllPeriods(MeshWeights meshWeights) {
		int nbSteps = this.getNbSteps();
		IntArrayList[][] weightsDestinationIndexList = meshWeights
				.getWeightsDestinationIndexListForStepAndOriginIndex();
		FloatArrayList[][] weightsList = meshWeights.getWeightsListForStepAndOriginIndex();
		int[][] nbWeights = meshWeights.getNumberOfWeightsForStepAndSimulationIndex();
		for (int k = nbSteps - 2; k > 0; k--) {
			weightsDestinationIndexList[k] = weightsDestinationIndexList[nbSteps - 1];
			weightsList[k] = weightsList[nbSteps - 1];
			nbWeights[k] = nbWeights[nbSteps - 1];
		}
	}

	private void updateStatePointsList(float[][] gridStates, float[][] identicalStartingStates) {
		List<float[][]> statePointsList = this.getStatePointsList();
		statePointsList.clear();
		statePointsList.add(0, identicalStartingStates);
		for (int kk = 1; kk <= this.getNbSteps(); kk++) {
			statePointsList.add(kk, gridStates);
		}
	}

	@Override
	protected MarketVectorProcess getMarketProcessForPricingUncondStateDensity() {
		return this.processForGeneratingOnlyOneStepForPricing;
	}

	@Override
	protected MarketVectorProcess getMarketProcessForUncondStateDensity() {
		return this.processForGeneratingOnlyOneStep;
	}

	@Override
	protected void normalizeAllWeights(MeshWeights weights) {
		int nbStepsToConsider = Math.min(2, this.getNbSteps());
		int nbPaths = this.getNbPaths();
		for (int k = 0; k < nbStepsToConsider; k++) {
			for (int i = 0; i < nbPaths; i++) {
				this.normalizeAllWeightsForStepAndIndex(weights, k, i);
			}
		}
	}
}
