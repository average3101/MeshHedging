package meshmethods;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import cern.colt.list.FloatArrayList;
import cern.colt.list.IntArrayList;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.rng.RandomStream;

public class StochasticMesh {
	private HashSet<Integer>[] attainableStateIndices;
	private String basicLabel;
	private boolean isUsingQMC = false;
	protected MarketVectorProcess marketProcess;
	protected MarketVectorProcess marketProcessForPricing;
	private MeshWeights meshWeights;
	private MeshWeights meshWeightsForPricing;
	private int nbPaths;
	private int nbProcessDim;
	private int nbSteps;
	private RandomStream randomStreamForMeshStateGeneration;
	private List<float[][]> statePointsList = new ArrayList<float[][]>(); // [nbPaths][nbProcess]
	private float[] tempMarketVectorDestination;
	private MarketVector tempMarketVectorGeneral;
	private float[] tempMarketVectorOrigin;
	private boolean usingDifferentProcessForPricing = false;
	protected double fractionAttainableStates;

	public StochasticMesh(MarketVectorProcess multiSp) {
		this.setMarketProcess(multiSp);
		this.nbProcessDim = multiSp.getDimension();
		this.setBasicLabel("SM");

		this.defineDefaultMeshWeights();
	}

	public void applyAdditionalTransformations(MeshWeights weights) {
		if (this.meshWeights.isUsingNormalizedWeights()) {
			this.normalizeAllWeights(weights);
		}
	}

	public void computeAllWeights(MeshWeights weights, MarketVectorProcess processForUncondDensity) {
		weights.init();

		for (int timeStep = this.getNbSteps() - 1; timeStep >= 0; timeStep--) {
			weights.precomputeInverseUnconditionalDensity(timeStep, processForUncondDensity);
			weights.computeWeightsForIntermediateStage(timeStep);
		}

		this.applyAdditionalTransformations(weights);
	}

	public void defineDefaultMeshWeights() {
		MeshWeights weights = this.newDefaultMeshWeights();
		this.setMeshWeights(weights);
	}

	public void generateMesh() {
		this.init();
		this.simulatePaths();

		this.computeAllWeights(this.meshWeights, this.getMarketProcessForUncondStateDensity());
		if (this.isUsingDifferentProcessForPricing()) {
			this.computeAllWeights(this.getMeshWeightsForPricing(),
					this.getMarketProcessForPricingUncondStateDensity());
		}

		this.identifyAttainableMeshStatesByForwardRecursion();
	}

	public HashSet<Integer>[] getAttainableStateIndices() {
		return this.attainableStateIndices;
	}

	public String getBasicLabel() {
		return this.basicLabel;
	}

	public String getDescription() {
		String s = this.getBasicLabel() + "-" + this.meshWeights.getDescription();
		if (this.isUsingQMC()) {
			s += "-Q";
		}
		if (this.meshWeights.isUsingRussianRoulette()) {
			s += "-R";
		}
		return s;
	}

	public double getFractionAttainableStates() {
		return this.fractionAttainableStates;
	}

	public MarketVectorProcess getMarketProcessForPricing() {
		return this.marketProcessForPricing;
	}

	public String getMarketVectorFullString(int k, int i) {
		String fullString = k + ";" + i + ";";
		fullString += this.getMarketVectorValuesString(k, i);
		return fullString;
	}

	public float getMarketVectorPrice(int k, int i) {
		float[][] meshCrossSection = this.statePointsList.get(k);
		float[] y = meshCrossSection[i];
		float s = y[MarketVector.PRICE];
		return s;
	}

	public MarketVectorProcess getMarketVectorProcess() {
		return this.marketProcess;
	}

	public float[] getMarketVectorValues(int k, int i) {
		float[][] meshCrossSection = this.statePointsList.get(k);
		float[] y = meshCrossSection[i];
		return y;
	}

	public String getMarketVectorValuesString(int k, int i) {
		float[] values = this.getMarketVectorValues(k, i);
		String vectorString = "";
		for (int j = 0; j < values.length; j++) {
			vectorString += values[j] + ";";
		}

		return vectorString;
	}

	public MeshWeights getMeshWeights() {
		return this.meshWeights;
	}

	public MeshWeights getMeshWeightsForPricing() {
		return this.meshWeightsForPricing;
	}

	public int getNbPaths() {
		return this.nbPaths;
	}

	public int getNbProcess() {
		return this.nbProcessDim;
	}

	public int getNbSteps() {
		return this.nbSteps;
	}

	public RandomStream getRandomStreamForMeshStateGeneration() {
		return this.randomStreamForMeshStateGeneration;
	}

	public List<float[][]> getStatePointsList() {
		return this.statePointsList;
	}

	public MarketVector getTempMarketVector(int k, int i) {
		float[] values = this.getMarketVectorValues(k, i);
		this.tempMarketVectorGeneral.setValues(values);
		this.tempMarketVectorGeneral.setSimulationIndex(i);
		return this.tempMarketVectorGeneral;
	}

	public float[] getTempMarketVectorDestination() {
		return this.tempMarketVectorDestination;
	}

	public MarketVector getTempMarketVectorGeneral() {
		return this.tempMarketVectorGeneral;
	}

	public float[] getTempMarketVectorOrigin() {
		return this.tempMarketVectorOrigin;
	}

	public void identifyAttainableMeshStatesByForwardRecursion() {
		this.attainableStateIndices = new HashSet[this.nbSteps + 1];
		this.attainableStateIndices[0] = new HashSet<Integer>();
		int totalNbStates = this.nbPaths * (this.nbSteps + 1);
		int attainableNbStates = 0;
		for (int i = 0; i < this.nbPaths; i++) {
			this.attainableStateIndices[0].add(i);
			attainableNbStates++;
		}

		for (int k = 1; k <= this.nbSteps; k++) {
			this.attainableStateIndices[k] = new HashSet<Integer>();
			Iterator<Integer> startingPointIterator = this.attainableStateIndices[k - 1].iterator();
			this.setAvailableMeshStateIndicesAtStep(k, startingPointIterator);
			attainableNbStates += this.attainableStateIndices[k].size();
		}
		this.fractionAttainableStates = ((double) attainableNbStates) / (double) totalNbStates;
	}

	public void init() {
		this.assertParameterNonZero(this.nbPaths, "nbPaths");
		this.assertParameterNonZero(this.nbSteps, "nbSteps");
		this.assertParameterNonZero(this.nbProcessDim, "nbProcess");

		int marketVectorDim = this.marketProcess.getMarketVectorDimension();
		this.tempMarketVectorOrigin = new float[marketVectorDim];
		this.tempMarketVectorDestination = new float[marketVectorDim];
		this.tempMarketVectorGeneral = this.marketProcess.newTempMarketVector();
	}

	public boolean isUsingDifferentProcessForPricing() {
		return this.usingDifferentProcessForPricing;
	}

	public boolean isUsingQMC() {
		return this.isUsingQMC;
	}

	public MeshWeights newDefaultMeshWeights() {
		MeshWeights weights = new MeshWeights();
		return weights;
	}

	public void resetNextRandomSubstream() {
		RandomStream stream = this.marketProcess.getStream();
		this.setRandomStreamForMeshStateGeneration(stream);
		stream.resetNextSubstream();
	}

	public void setBasicLabel(String label) {
		this.basicLabel = label;
	}

	public void setMarketProcess(MarketVectorProcess marketProcess) {
		this.marketProcess = marketProcess;
	}

	public void setMarketVectorProcessForPricing(MarketVectorProcess processForPricing) {
		if (!processForPricing.getClass().equals(this.marketProcess.getClass())) {
			this.setUsingDifferentProcessForPricing(true);
			this.marketProcessForPricing = processForPricing;
			MeshWeights weights = this.newDefaultMeshWeights();
			this.setMeshWeightsForDerivative(weights);
		}
	}

	public void setMeshWeights(MeshWeights weights) {
		this.meshWeights = weights;
		this.meshWeights.setStochasticMesh(this);
		this.meshWeights.setProcessForDensity(this.getMarketVectorProcess());
		this.meshWeights.setProcessForUnconditionalDensity(this.getMarketProcessForUncondStateDensity());
	}

	public void setMeshWeightsForDerivative(MeshWeights weights) {
		this.setMeshWeightsForPricing(weights);
		this.getMeshWeightsForPricing().setStochasticMesh(this);
		this.getMeshWeightsForPricing().setProcessForDensity(this.getMarketProcessForPricing());
		this.getMeshWeightsForPricing()
				.setProcessForUnconditionalDensity(this.getMarketProcessForPricingUncondStateDensity());
	}

	public void setMeshWeightsForPricing(MeshWeights meshWeightsForPricing) {
		this.meshWeightsForPricing = meshWeightsForPricing;
	}

	public void setNbPaths(int n) {
		this.nbPaths = n;
	}

	public void setNbSteps(int n) {
		this.nbSteps = n;
	}

	public void setObservationTimes(double[] obsTimes, int nbSteps) {
		this.nbSteps = nbSteps;
		this.marketProcess.setObservationTimes(obsTimes, this.nbSteps);

		if (this.isUsingDifferentProcessForPricing()) {
			this.getMarketProcessForPricing().setObservationTimes(obsTimes, this.nbSteps);
		}
	}

	public void setRandomStreamForMeshStateGeneration(RandomStream randomStreamForMeshStateGeneration) {
		this.randomStreamForMeshStateGeneration = randomStreamForMeshStateGeneration;
	}

	public void setUsingDifferentProcessForPricing(boolean usingDifferentProcessForPricing) {
		this.usingDifferentProcessForPricing = usingDifferentProcessForPricing;
	}

	public void setUsingQMC(boolean isUsingQMC) {
		this.isUsingQMC = isUsingQMC;
	}

	public void simulatePaths() {
		this.resetNextRandomSubstream();
		this.statePointsList.clear();
		int marketVectorDim = this.marketProcess.getMarketVectorDimension();
		for (int k = 0; k <= this.nbSteps; k++) {
			this.statePointsList.add(k, new float[this.nbPaths][marketVectorDim]);
		}

		float[][] sPathArray;
		float[][] stateCrossSection;
		for (int i = 0; i < this.nbPaths; i++) {
			sPathArray = this.marketProcess.generateBasicProcessesPaths();
			for (int k = 0; k <= this.nbSteps; k++) {
				stateCrossSection = this.statePointsList.get(k);
				for (int j = 0; j < marketVectorDim; j++) {
					stateCrossSection[i][j] = sPathArray[j][k];
				}
			}
			this.resetNextRandomSubstream();
		}
	}

	private void assertParameterNonZero(int value, String paramName) {
		try {
			if (value <= 0) {
				throw new Exception("Invalid parameter value : " + paramName + " = " + value);
			}
		} catch (Exception e) {

			e.printStackTrace();
		}
	}

	private void setAvailableMeshStateIndicesAtStep(int k, Iterator<Integer> startingPointIterator) {
		IntArrayList availableIndices;
		while (startingPointIterator.hasNext()) {
			int i = startingPointIterator.next().intValue();
			availableIndices = this.meshWeights.getWeightsIndexList(k - 1, i);
			for (int jj = 0; jj < availableIndices.size(); jj++) {
				int j = availableIndices.getQuick(jj);
				this.attainableStateIndices[k].add(j);
			}
		}
	}

	protected MarketVectorProcess getMarketProcessForPricingUncondStateDensity() {
		return this.getMarketProcessForPricing();
	}

	protected MarketVectorProcess getMarketProcessForUncondStateDensity() {
		return this.marketProcess;
	}

	protected void normalizeAllWeights(MeshWeights weights) {
		for (int k = 0; k < this.nbSteps; k++) {
			for (int i = 0; i < this.nbPaths; i++) {
				this.normalizeAllWeightsForStepAndIndex(weights, k, i);
			}
		}
	}

	protected void normalizeAllWeightsForStepAndIndex(MeshWeights weights, int k, int i) {
		FloatArrayList[][] weightsListForStepAndOriginIndex = weights.getWeightsListForStepAndOriginIndex();
		FloatArrayList weightsList = weightsListForStepAndOriginIndex[k][i];
		FloatArrayList normalizedWeights = weights.normalize(weightsList);
		weightsListForStepAndOriginIndex[k][i] = normalizedWeights;
	}

}
