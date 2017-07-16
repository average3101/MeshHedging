package meshmethods;

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import cern.colt.list.FloatArrayList;
import cern.colt.list.IntArrayList;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.rng.RandomStream;

public class MeshWeights {
	private long cardinalityOfActualGrid;
	private long cardinalityOfFullGrid;
	private double cutoffFactor;
	private double fractionOfWeightsUsed;
	protected double[][] inverseUncondDensity;
	private boolean isUsingRussianRoulette = false;
	private int nbCrossSectionWeights;
	private int[] nbWeightsForOuterMarketInfoPath;
	private int[][] numberOfWeightsForStepAndSimulationIndex;
	private MarketVector[] outerMarketVectorPath;
	private MarketVectorProcess processForDensity;
	private MarketVectorProcess processForUnconditionalDensity;
	private boolean removeTransitionsToIdenticalIndex = false;
	private double rouletteThreshold;
	private double rouletteThresholdInverse;
	private StochasticMesh stochasticMesh;
	private float[] tempWeights;
	private boolean usingCVTransformation = false;
	private boolean usingIdenticalWeightsAtFirstStep;
	private boolean usingNormalizedWeights = true;
	private RandomStream uWeightRoulette;
	private int[][] weightDestinationIndicesForOuterMarketInfoPath;
	private float[] weightsCrossSection;
	private int[] weightsCrossSectionDestinationIndices;
	private IntArrayList[][] weightsDestinationIndexListForStepAndOriginIndex; // [nbObs][nbPaths]
	private float[][] weightsForOuterMarketInfoPath;
	private FloatArrayList[][] weightsListForStepAndOriginIndex; // [nbObs][nbPaths]

	public MeshWeights() {
		this.setUsingIdenticalWeightsAtFirstStep(true);
	}

	public void addWeightAndDestinationIndexToLists(float weight, int destinationIndex, int timeStep, int originIndex) {
		this.weightsDestinationIndexListForStepAndOriginIndex[timeStep][originIndex].add(destinationIndex);
		this.weightsListForStepAndOriginIndex[timeStep][originIndex].add(weight);
		this.numberOfWeightsForStepAndSimulationIndex[timeStep][originIndex]++;
	}

	public double applyRussianRoulette(double weight) {
		if (weight <= this.rouletteThreshold) {
			double p = weight * this.rouletteThresholdInverse;
			if (p > this.uWeightRoulette.nextDouble()) {
				weight = this.rouletteThreshold;
			} else {
				weight = 0;
			}
		}
		return weight;
	}

	public void computeAllWeightsForMarketVectorPath(MarketVector[] marketInfoPath) {
		int nbPaths = this.getStochasticMesh().getNbPaths();
		int nbSteps = this.getStochasticMesh().getNbSteps();

		this.weightsForOuterMarketInfoPath = new float[nbSteps][nbPaths];
		this.weightDestinationIndicesForOuterMarketInfoPath = new int[nbSteps][nbPaths];
		this.nbWeightsForOuterMarketInfoPath = new int[nbSteps];
		this.outerMarketVectorPath = new MarketVector[nbSteps];

		for (int k = 0; k < nbSteps; k++) {
			this.computeWeightsCrossSectionForOuterPaths(k, marketInfoPath[k].getValues());

			for (int i = 0; i < nbPaths; i++) {
				this.weightsForOuterMarketInfoPath[k][i] = this.weightsCrossSection[i];
				this.weightDestinationIndicesForOuterMarketInfoPath[k][i] = this.weightsCrossSectionDestinationIndices[i];
			}
			this.nbWeightsForOuterMarketInfoPath[k] = this.nbCrossSectionWeights;
			this.outerMarketVectorPath[k] = marketInfoPath[k].clone();
		}
	}

	public double computeConditionalDensity(int k, int destinationIndex, float[] sk) {
		this.copyStatePointToTempMarketVectorDestination(k + 1, destinationIndex);
		float[] tempMarketVectorDestination = this.stochasticMesh.getTempMarketVectorDestination();

		double f = this.processForDensity.density(sk, tempMarketVectorDestination, k);
		return f;
	}

	public double computeFractionOfWeightsUsed() {
		int nbPaths = this.getStochasticMesh().getNbPaths();
		int nbSteps = this.getStochasticMesh().getNbSteps();

		long a = nbPaths * nbPaths;
		a *= (nbSteps - 1);
		a += nbPaths;
		this.cardinalityOfFullGrid = a;
		// = (nbObs-1) * nbPaths * nbPaths + nbPaths;
		// 1st step has a fixed starting point
		int[][] numberOfWeightsForStepAndSimulationIndex = this.getNumberOfWeightsForStepAndSimulationIndex();
		this.cardinalityOfActualGrid = numberOfWeightsForStepAndSimulationIndex[0][0];
		for (int k = 1; k < nbSteps; k++) {
			for (int i = 0; i < nbPaths; i++) {
				this.cardinalityOfActualGrid += numberOfWeightsForStepAndSimulationIndex[k][i];
			}
		}
		this.fractionOfWeightsUsed = ((double) this.cardinalityOfActualGrid) / this.cardinalityOfFullGrid;
		return this.fractionOfWeightsUsed;
	}

	public double computeLikelihoodRatio(int k, int destinationIndex) {
		float[] tempMarketVectorOrigin = this.stochasticMesh.getTempMarketVectorOrigin();

		return this.computeLikelihoodRatio(k, destinationIndex, tempMarketVectorOrigin);
	}

	public double computeLikelihoodRatio(int k, int destinationIndex, float[] sk) {
		double f = this.computeConditionalDensity(k, destinationIndex, sk);
		double g = this.inverseUncondDensity[k][destinationIndex];
		double w = f * g;

		return w;
	}

	public float[] computeWeightsCrossSectionForOuterPaths(int k, float[] sk) {
		this.nbCrossSectionWeights = 0;
		HashSet<Integer>[] availableStateIndices = this.stochasticMesh.getAttainableStateIndices();
		Iterator<Integer> availableIndicesIterator = availableStateIndices[k + 1].iterator();
		while (availableIndicesIterator.hasNext()) {
			int destinationIndex = availableIndicesIterator.next().intValue();
			double w = this.computeLikelihoodRatio(k, destinationIndex, sk);

			this.weightsCrossSection[this.nbCrossSectionWeights] = (float) w;
			this.weightsCrossSectionDestinationIndices[this.nbCrossSectionWeights] = destinationIndex;
			this.nbCrossSectionWeights++;
		}

		// Always normalize for outer paths (as opposed to inside mesh, where it
		// depends on the method)
		this.normalize(this.weightsCrossSection, this.nbCrossSectionWeights);

		return this.weightsCrossSection;
	}

	public void computeWeightsForAllDestinations(int timeStep, int originIndex) {
		int nbPaths = this.stochasticMesh.getNbPaths();

		for (int destinationIndex = 0; destinationIndex < nbPaths; destinationIndex++) {
			float weight = (float) this.computeLikelihoodRatio(timeStep, destinationIndex);

			if (this.isRemoveTransitionsToIdenticalIndex() && originIndex == destinationIndex && timeStep > 0) {
				weight = 0;
			}

			this.updateWeightInformation(weight, destinationIndex, timeStep, originIndex);
		}
	}

	public void computeWeightsForIntermediateStage(int timeStep) {
		int nbPaths = this.stochasticMesh.getNbPaths();
		for (int originIndex = 0; originIndex < nbPaths; originIndex++) {
			this.initializeWeightsListForStepAndOriginIndex(timeStep, originIndex);
			this.copyStatePointToTempMarketVectorOrigin(timeStep, originIndex);
			this.computeWeightsForAllDestinations(timeStep, originIndex);
		}
	}

	public void copyStatePointToTempMarketVectorDestination(int timeStep, int destinationIndex) {
		float[] tempMarketVectorDestination = this.stochasticMesh.getTempMarketVectorDestination();
		this.copyStatePointToTempMarketVector(timeStep, destinationIndex, tempMarketVectorDestination);
	}

	public void copyStatePointToTempMarketVectorOrigin(int timeStep, int originIndex) {
		float[] tempMarketVectorOrigin = this.stochasticMesh.getTempMarketVectorOrigin();
		this.copyStatePointToTempMarketVector(timeStep, originIndex, tempMarketVectorOrigin);
	}

	public String getDescription() {
		return "LR";
	}

	public double getFractionOfWeightsUsed() {
		this.computeFractionOfWeightsUsed();
		return this.fractionOfWeightsUsed;
	}

	public double[][] getInverseUnconditionalDensity() {
		return this.inverseUncondDensity;
	}

	public int getNbCrossSectionWeights() {
		return this.nbCrossSectionWeights;
	}

	public int[] getNbWeightsForMarketInfoPath() {
		return this.nbWeightsForOuterMarketInfoPath;
	}

	public int[][] getNbWeightsPerList() {
		return this.getNumberOfWeightsForStepAndSimulationIndex();
	}

	public int getNbWeightsPerList(int k, int i) {
		return this.getNumberOfWeightsForStepAndSimulationIndex()[k][i];
	}

	public int[][] getNumberOfWeightsForStepAndSimulationIndex() {
		return this.numberOfWeightsForStepAndSimulationIndex;
	}

	public MarketVector[] getOuterMarketVectorPath() {
		return this.outerMarketVectorPath;
	}

	public MarketVectorProcess getProcessForDensity() {
		return this.processForDensity;
	}

	public MarketVectorProcess getProcessForUnconditionalDensity() {
		return this.processForUnconditionalDensity;
	}

	public double getRouletteThreshold() {
		return this.rouletteThreshold;
	}

	public StochasticMesh getStochasticMesh() {
		return this.stochasticMesh;
	}

	public float[] getTempWeights() {
		return this.tempWeights;
	}

	public int[][] getWeightDestinationIndicesForMarketInfoPath() {
		return this.weightDestinationIndicesForOuterMarketInfoPath;
	}

	public int[] getWeightsCrossSectionIndices() {
		return this.weightsCrossSectionDestinationIndices;
	}

	public IntArrayList[][] getWeightsDestinationIndexListForStepAndOriginIndex() {
		return this.weightsDestinationIndexListForStepAndOriginIndex;
	}

	public float[][] getWeightsForOuterMarketInfoPath() {
		return this.weightsForOuterMarketInfoPath;
	}

	public IntArrayList[][] getWeightsIndexList() {
		return this.getWeightsDestinationIndexListForStepAndOriginIndex();
	}

	public IntArrayList getWeightsIndexList(int k, int i) {
		return this.getWeightsDestinationIndexListForStepAndOriginIndex()[k][i];
	}

	public FloatArrayList[][] getWeightsList() {
		return this.getWeightsListForStepAndOriginIndex();
	}

	public FloatArrayList getWeightsList(int k, int i) {
		return this.getWeightsListForStepAndOriginIndex()[k][i];
	}

	public FloatArrayList[][] getWeightsListForStepAndOriginIndex() {
		return this.weightsListForStepAndOriginIndex;
	}

	public void init() {
		int nbPaths = this.stochasticMesh.getNbPaths();
		int nbSteps = this.stochasticMesh.getNbSteps();

		this.rouletteThreshold = this.cutoffFactor / nbPaths;
		if (this.cutoffFactor == 0.0) {
			this.rouletteThresholdInverse = 0.0;
		} else {
			this.rouletteThresholdInverse = 1 / this.rouletteThreshold;
		}

		this.setTempWeights(new float[nbPaths]);
		this.weightsCrossSection = new float[nbPaths];
		this.weightsCrossSectionDestinationIndices = new int[nbPaths];

		this.weightsListForStepAndOriginIndex = new FloatArrayList[nbSteps][nbPaths];
		this.weightsDestinationIndexListForStepAndOriginIndex = new IntArrayList[nbSteps][nbPaths];
		this.numberOfWeightsForStepAndSimulationIndex = new int[nbSteps][nbPaths];
		this.inverseUncondDensity = new double[nbSteps][nbPaths];
	}

	public void initializeWeightsListForStepAndOriginIndex(int k, int i) {
		this.weightsListForStepAndOriginIndex[k][i] = new FloatArrayList();
		this.weightsDestinationIndexListForStepAndOriginIndex[k][i] = new IntArrayList();
		this.numberOfWeightsForStepAndSimulationIndex[k][i] = 0;
	}

	public boolean isRemoveTransitionsToIdenticalIndex() {
		return this.removeTransitionsToIdenticalIndex;
	}

	public boolean isUseIdenticalWeightsAtFirstStep() {
		return this.usingIdenticalWeightsAtFirstStep;
	}

	public boolean isUsingCVTransformation() {
		return this.usingCVTransformation;
	}

	public boolean isUsingNormalizedWeights() {
		return this.usingNormalizedWeights;
	}

	public boolean isUsingRussianRoulette() {
		return this.isUsingRussianRoulette;
	}

	public float[] normalize(float[] weightsCrossSection, int nbWeights) {
		double sum = 0.0;
		for (int i = 0; i < nbWeights; i++) {
			sum += weightsCrossSection[i];
		}
		for (int i = 0; i < nbWeights; i++) {
			weightsCrossSection[i] = (float) (weightsCrossSection[i] / sum);
		}
		return weightsCrossSection;
	}

	public FloatArrayList normalize(FloatArrayList weights) {
		double sum = 0.0;
		for (int i = 0; i < weights.size(); i++) {
			sum += weights.getQuick(i);
		}

		if (sum == 0.0) {
			return weights;
		} else {
			for (int i = 0; i < weights.size(); i++) {
				float value = (float) (weights.getQuick(i) / sum);
				weights.setQuick(i, value);
			}
		}

		return weights;
	}

	public void precomputeInverseUnconditionalDensity(int k, MarketVectorProcess process) {
		int nbPaths = this.getStochasticMesh().getNbPaths();
		int nbProcess = this.getStochasticMesh().getNbProcess();
		List<float[][]> statePointsList = this.getStochasticMesh().getStatePointsList();
		float[][] meshCS0 = statePointsList.get(k);
		float[][] meshCS1 = statePointsList.get(k + 1);
		float[] tempMarketVectorOrigin = this.getStochasticMesh().getTempMarketVectorOrigin();
		float[] tempMarketVectorDestination = this.getStochasticMesh().getTempMarketVectorDestination();
		for (int j = 0; j < nbPaths; j++) {
			for (int l = 0; l < nbProcess; l++) {
				tempMarketVectorOrigin[l] = meshCS0[j][l];
			}
			for (int l = 0; l < nbProcess; l++) {
				tempMarketVectorDestination[l] = meshCS1[j][l];
			}

			double g = process.density(tempMarketVectorOrigin, tempMarketVectorDestination, k);

			double inverse = 1 / (nbPaths * g);
			if (!Double.isFinite(inverse)) {
				inverse = 0.0;
			}
			this.inverseUncondDensity[k][j] = (float) inverse;
		}
	}

	public void setInverseUncondDensity(double[][] g) {
		this.inverseUncondDensity = g;
	}

	public void setProcessForDensity(MarketVectorProcess processForDensity) {
		this.processForDensity = processForDensity;
	}

	public void setProcessForUnconditionalDensity(MarketVectorProcess processForUnconditionalDensity) {
		this.processForUnconditionalDensity = processForUnconditionalDensity;
	}

	public void setRemoveTransitionsToIdenticalIndex(boolean removeTransitionsToIdenticalIndex) {
		this.removeTransitionsToIdenticalIndex = removeTransitionsToIdenticalIndex;
	}

	public void setStochasticMesh(StochasticMesh stochasticMesh) {
		this.stochasticMesh = stochasticMesh;
	}

	public void setTempWeights(float[] tempWeights) {
		this.tempWeights = tempWeights;
	}

	public void setUsingCVTransformation(boolean usingCVTransformation) {
		this.usingCVTransformation = usingCVTransformation;
	}

	public void setUsingIdenticalWeightsAtFirstStep(boolean usingIdenticalWeightsAtFirstStep) {
		this.usingIdenticalWeightsAtFirstStep = usingIdenticalWeightsAtFirstStep;
	}

	public void setUsingNormalizedWeights(boolean usingNormalizedWeights) {
		this.usingNormalizedWeights = usingNormalizedWeights;
	}

	public void updateWeightInformation(float weight, int destinationIndex, int timeStep, int originIndex) {
		if (this.isUsingRussianRoulette) {
			weight = (float) this.applyRussianRoulette(weight);
		}
		if (weight > 0) {
			this.addWeightAndDestinationIndexToLists(weight, destinationIndex, timeStep, originIndex);
		}
	}

	public void useRussianRoulette(double cutoffFactor, RandomStream uWeightRoulette) {
		this.cutoffFactor = cutoffFactor;
		this.uWeightRoulette = uWeightRoulette;
		this.isUsingRussianRoulette = true;
	}

	private void assignIdenticalWeightsToAllDestinations(int nbPaths, int timeStep, int originIndex) {
		for (int destinationIndex = 0; destinationIndex < nbPaths; destinationIndex++) {
			float weight = (float) (1.0 / nbPaths);
			this.addWeightAndDestinationIndexToLists(weight, destinationIndex, timeStep, originIndex);
		}
	}

	private void copyStatePointToTempMarketVector(int timeStep, int originIndex, float[] tempMarketVectorOrigin) {
		List<float[][]> statePointsList = this.stochasticMesh.getStatePointsList();
		float[][] meshCrossSection = statePointsList.get(timeStep);
		for (int processIndex = 0; processIndex < tempMarketVectorOrigin.length; processIndex++) {
			tempMarketVectorOrigin[processIndex] = meshCrossSection[originIndex][processIndex];
		}
	}

}
