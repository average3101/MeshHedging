package meshmethods;

import cern.colt.list.FloatArrayList;
import cern.colt.list.IntArrayList;
import numerics.OneDimFunction;
import stochprocess.MarketVector;

public class MeshFunctionExpectationEstimator extends MeshFullEstimator {

	private OneDimFunction meshBasedFunction;
	private double[][] valueArray;
	private double averageEstimate;

	public MeshFunctionExpectationEstimator(OneDimFunction f) {
		this.meshBasedFunction = f;
	}

	@Override
	public void computeForState(int timeStep, int index) {
		MarketVector y = this.getTempMarketVectorOrigin(timeStep, index);
		double value = this.evaluate(timeStep, y);
		int i = y.simulationIndex;
		this.valueArray[timeStep][i] = value;
	}

	@Override
	public void computeForStateOverLastStep(int i) {
		int stepNo = this.getNbSteps();
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		double s = stochasticMesh.getMarketVectorPrice(stepNo, i);
		double value = this.meshBasedFunction.evaluate(s);
		this.valueArray[stepNo][i] = value;
	}

	public double evaluate(int k, MarketVector y) {
		double weightedSum = 0.0;
		int nextStep = k + 1;
		int i = y.simulationIndex;

		StochasticMesh stochasticMesh = this.getStochasticMesh();
		MeshWeights meshWeights = stochasticMesh.getMeshWeights();
		FloatArrayList availableWeights = meshWeights.getWeightsList(k, i);
		IntArrayList availableIndices = meshWeights.getWeightsIndexList(k, i);
		int currentNbWeights = meshWeights.getNbWeightsPerList(k, i);
		// String startMarketVectorString =
		// stochasticMesh.getMarketVectorFullString(k, i);
		for (int jj = 0; jj < currentNbWeights; jj++) {
			int j = availableIndices.get(jj);
			double w = availableWeights.get(jj);
			double value = this.valueArray[nextStep][j];

			// String marketVectorFullString =
			// stochasticMesh.getMarketVectorFullString(nextStep, j);
			// String infoToPrint = meshWeights.getDescription() + ";" +
			// startMarketVectorString + ";"
			// + marketVectorFullString + ";" + w + ";" + value;
			// HedgingComparisonsOutput.writeLineToFile(infoToPrint,
			// "meshValidations20170307.csv", true);
			weightedSum += w * value;
		}

		return weightedSum;
	}

	public double evaluateAtInitialState() {
		int stepNo = 0;
		double sum = 0;
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		int nbPaths = stochasticMesh.getNbPaths();
		for (int i = 0; i < nbPaths; i++) {
			MarketVector tempMarketVector = this.getTempMarketVectorOrigin(stepNo, i);
			sum += this.evaluate(stepNo, tempMarketVector);
		}
		double avg = sum / nbPaths;
		return avg;
	}

	public double getAverageEstimate() {
		return this.averageEstimate;
	}

	public double[][] getValueArray() {
		return this.valueArray;
	}

	@Override
	public void initialize() {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		int nbPaths = stochasticMesh.getNbPaths();
		int nbSteps = stochasticMesh.getNbSteps();
		this.setNbSteps(nbSteps);
		this.valueArray = new double[nbSteps + 1][nbPaths];
	}

	@Override
	protected void computeAverageEstimate() {
		double sum = 0;
		int timeStep = 0;
		int nbPaths = this.getStochasticMesh().getNbPaths();
		for (int i = 0; i < nbPaths; i++) {
			sum += this.valueArray[timeStep][i];
		}
		this.averageEstimate = sum / nbPaths;
	}

}
