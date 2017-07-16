package derivatives;

import java.util.List;

import stochprocess.MarketVector;

public class DerivativePricedOnMeshWithCVAnalytic extends DerivativePricedOnMesh {

	protected String evaluationMethod;
	protected DerivativePriceControlVariate control;
	protected DerivativePriceControlVariate.Label cvLabel;

	public DerivativePricedOnMeshWithCVAnalytic(Derivative basicDerivative,
			DerivativePriceControlVariate.Label cvLabel) {
		super(basicDerivative);

		this.evaluationMethod = "SM-CV-AN";

		this.cvLabel = cvLabel;
		this.control = DerivativePriceControlVariate.createControl(cvLabel, basicDerivative);
	}

	@Override
	public Object clone() {
		DerivativePricedOnMeshWithCVAnalytic derivativeCopy = new DerivativePricedOnMeshWithCVAnalytic(
				this.basicDerivative, this.cvLabel);
		derivativeCopy.setStochasticMesh(this.stochasticMesh);
		return derivativeCopy;
	}

	@Override
	public double evaluate(int k, float[] marketVector) {
		double weightedSum = 0.0;
		int nextStep = k + 1;

		List<float[][]> statePointList = this.stochasticMesh.getStatePointsList();
		float[][] meshCrossSection = statePointList.get(nextStep);

		double cv0 = this.control.getExpectedValue(k, marketVector);
		float sigma = marketVector[MarketVector.VOL];
		float s = marketVector[MarketVector.PRICE];

		double cvEffect = 0.0;
		double sumOfWeights = 0.0;
		for (int jj = 0; jj < this.currentNbWeights; jj++) {
			int j = this.currentWeightIndices[jj];
			double w = this.currentWeights[jj];
			float[] y = meshCrossSection[j];
			double h1 = this.derivativeValueArray[nextStep][j];

			double cv = this.control.evaluate(k, marketVector, nextStep, y) - cv0;
			cvEffect += w * cv;
			weightedSum += w * h1;
			sumOfWeights += w;
		}
		double beta = this.control.getMultiplicativeFactor(s, sigma, k);
		sumOfWeights = 1;
		double cvEstimate = (weightedSum - beta * cvEffect) / sumOfWeights;
		return cvEstimate;
	}

	@Override
	public String getDescription() {
		String cvLabelString = this.cvLabel.toString();
		String label = this.getBasicDerivative().getDescription() + " " + this.evaluationMethod + " (" + cvLabelString
				+ ")";
		return label;
	}

}
