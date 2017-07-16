package derivatives;

import meshmethods.StochasticMesh;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;

public class DerivativePricedOnMesh extends Derivative {
	protected int currentNbWeights;
	protected int[] currentWeightIndices;
	protected float[] currentWeights;
	protected StochasticMesh stochasticMesh;
	private int[][] weightIndicesPath;
	private int[] weightNbPath;
	private float[][] weightPathForDerivativePath;
	protected Derivative basicDerivative;

	public DerivativePricedOnMesh(Derivative basicDerivative) {
		this.setBasicDerivative(basicDerivative);
	}

	@Override
	public Object clone() {
		DerivativePricedOnMesh derivativeCopy = new DerivativePricedOnMesh(this.getBasicDerivative());
		derivativeCopy.setStochasticMesh(this.getStochasticMesh());
		return derivativeCopy;
	}

	@Override
	public double delta(double s, double sigma, int k) {
		return this.basicDerivative.delta(s, sigma, k);
	}

	@Override
	public double evaluate(int k, float[] marketVector) {
		double weightedSum = 0.0;
		int nextStep = k + 1;
		double[][] derivativeValueArray = this.getDerivativeValueArray();
		for (int jj = 0; jj < this.currentNbWeights; jj++) {
			int j = this.currentWeightIndices[jj];
			double w = this.currentWeights[jj];
			double derivValueAtNextStep = derivativeValueArray[nextStep][j];

			weightedSum += w * (derivValueAtNextStep);
		}

		return weightedSum;
	}

	@Override
	public double evaluateAtT0(MarketVector y) {
		double sum = 0;
		int nbPaths = this.stochasticMesh.getNbPaths();
		for (int i = 0; i < nbPaths; i++) {
			sum += this.derivativeValueArray[0][i];
		}
		double average = sum / nbPaths;
		this.setInitDerivativeValue(average);

		return average;
	}

	@Override
	public double[] evaluatePath(MarketVector[] yArray) {
		double[] valuePath = this.getValuePath();
		valuePath[0] = this.getDerivativeValueArray()[0][0];
		int nbObs = this.getNbObs();
		for (int k = 1; k < nbObs; k++) {
			this.setCurrentWeights(this.weightPathForDerivativePath[k], this.weightIndicesPath[k],
					this.weightNbPath[k]);

			valuePath[k] = this.evaluate(k, yArray[k].values);
		}
		valuePath[nbObs] = this.payoff(yArray[nbObs].values[MarketVector.PRICE]);
		return valuePath;
	}

	@Override
	public double gamma(int k, float[] values) {
		return this.basicDerivative.gamma(k, values);
	}

	public Derivative getBasicDerivative() {
		return this.basicDerivative;
	}

	public int getCurrentNbWeights() {
		return this.currentNbWeights;
	}

	public int[] getCurrentWeightIndices() {
		return this.currentWeightIndices;
	}

	public float[] getCurrentWeights() {
		return this.currentWeights;
	}

	@Override
	public String getDescription() {
		return this.getBasicDerivative().getDescription() + " SM";
	}

	@Override
	public double getMaturity() {
		return this.getBasicDerivative().getMaturity();
	}

	@Override
	public double getMaxDeltaForHedge() {
		return this.getBasicDerivative().getMaxDeltaForHedge();
	}

	@Override
	public double getMinDeltaForHedge() {
		return this.getBasicDerivative().getMinDeltaForHedge();
	}

	public StochasticMesh getStochasticMesh() {
		return this.stochasticMesh;
	}

	@Override
	public double getStrike() {
		return this.basicDerivative.getStrike();
	}

	@Override
	public double payoff(float s) {
		return this.getBasicDerivative().payoff(s);
	}

	@Override
	public double payoff(float[] s) {
		return this.getBasicDerivative().payoff(s);
	}

	public void printDerivativeValuesAtStep(int k) {
		System.out.println();
		System.out.println("Derivative values at step " + k);
		System.out.println("i;hk;hk_cv");
		int nbStates = this.getDerivativeValueArray()[0].length;
		for (int i = 0; i < nbStates; i++) {
			String output = i + "; " + this.getDerivativeValueArray()[k][i] + "; " + this.getDerivativeCVArray()[k][i];
			System.out.println(output);
		}
	}

	public void setBasicDerivative(Derivative basicDerivative) {
		this.basicDerivative = basicDerivative;
	}

	public void setCurrentNbWeights(int currentNbWeights) {
		this.currentNbWeights = currentNbWeights;
	}

	public void setCurrentWeightIndices(int[] currentWeightIndices) {
		this.currentWeightIndices = currentWeightIndices;
	}

	public void setCurrentWeights(float[] currentWeights) {
		this.currentWeights = currentWeights;
	}

	public void setCurrentWeights(float[] weightPathForDerivativePath2, int[] currentWeightIndices,
			int currentNbWeights) {
		this.setCurrentNbWeights(currentNbWeights);
		this.currentWeights = weightPathForDerivativePath2;
		this.setCurrentWeightIndices(currentWeightIndices);
	}

	@Override
	public void setMaturity(double t) {
		this.getBasicDerivative().setMaturity(t);
	}

	@Override
	public void setObsTimes(double[] obsTimes) {
		super.setObsTimes(obsTimes);
		this.getBasicDerivative().setObsTimes(obsTimes);
	}

	public void setStochasticMesh(StochasticMesh sm) {
		this.stochasticMesh = sm;

		MarketVectorProcess process = sm.getMarketVectorProcess();
		if (process.requiresDifferentProcessForPricing()) {
			MarketVectorProcess processForPricing = process.newProcessForPricing();
			sm.setMarketVectorProcessForPricing(processForPricing);

		}
	}

	public void setWeightPathForDerivativePath(float[][] weightPathForDerivativePath, int[][] weightIndicesPath,
			int[] weightNbPath) {
		this.weightPathForDerivativePath = weightPathForDerivativePath;
		this.weightIndicesPath = weightIndicesPath;
		this.weightNbPath = weightNbPath;
	}

}
