package derivatives;

import java.util.List;

import derivatives.DerivativePriceControlVariate.Label;
import meshmethods.StochasticMesh;

public class DerivativePricedOnMeshWithCVEstimated extends DerivativePricedOnMeshWithCVAnalytic {
	protected int nbOfPilotRuns;
	protected int minNbOfPilotRuns = 20;
	protected double proportionOfPilotRuns = 0.10;

	public DerivativePricedOnMeshWithCVEstimated(Derivative basicDerivative, Label cvLabel) {
		super(basicDerivative, cvLabel);

		this.evaluationMethod = "SM-CV-ES";
	}

	@Override
	public Object clone() {
		DerivativePricedOnMeshWithCVEstimated derivativeCopy = new DerivativePricedOnMeshWithCVEstimated(
				this.basicDerivative, this.cvLabel);
		derivativeCopy.setStochasticMesh(this.stochasticMesh);
		return derivativeCopy;
	}

	// public double estimateBeta(double s, double sigma, int k) {
	// double a = 0;
	// double b = 0;
	//
	// int nextStep = k + 1;
	// List<float[][]> statePointList =
	// this.stochasticMesh.getStatePointsList();
	// float[][] meshCrossSection = statePointList.get(nextStep);
	//
	// double cv0 = this.control.getExpectedValue(s, sigma, k);
	// for (int jj = 0; jj < this.nbOfPilotRuns; jj++) {
	// int j = this.currentWeightIndices[jj];
	// double w = this.currentWeights[jj];
	//
	// float[] y = meshCrossSection[j];
	// float s1 = y[MarketVector.PRICE];
	// double h1 = this.derivativeValueArray[nextStep][j];
	//
	// double cv = this.control.evaluate(s1, sigma, nextStep) - cv0;
	//
	// a += w * h1 * cv;
	// b += w * cv * cv;
	// }
	//
	// double betaEstimate = a / b;
	// return betaEstimate;
	// }

	@Override
	public double evaluate(int k, float[] marketVector) {
		double a = 0;
		double b = 0;
		double c = 0;
		double d = 0;
		double e = 0;
		int nextStep = k + 1;

		List<float[][]> statePointList = this.stochasticMesh.getStatePointsList();
		float[][] meshCrossSection = statePointList.get(nextStep);

		double cv0 = this.control.getExpectedValue(k, marketVector);
		for (int jj = 0; jj < this.currentNbWeights; jj++) {
			int j = this.currentWeightIndices[jj];
			double w = this.currentWeights[jj];
			float[] y = meshCrossSection[j];
			double h1 = this.derivativeValueArray[nextStep][j];
			double cv = this.control.evaluate(k, marketVector, nextStep, y) - cv0;

			a += w;
			b += w * h1;
			c += w * cv;
			d += w * h1 * cv;
			e += w * cv * cv;
		}
		double betaEstimate = (a * d - c * b) / (a * e - c * c);
		double cvEstimate = (b - betaEstimate * c) / a;

		return cvEstimate;
	}

	@Override
	public void initializeAfterMeshComputation(int nbObs, int nbStatesPerStep, StochasticMesh stochasticMesh) {
		super.initializeAfterMeshComputation(nbObs, nbStatesPerStep, stochasticMesh);
		this.setNumberOfPilotRuns();
	}

	private void setNumberOfPilotRuns() {
		int targetNbOfPilotRuns = (int) (this.proportionOfPilotRuns * this.stochasticMesh.getNbPaths());
		this.nbOfPilotRuns = Math.max(this.minNbOfPilotRuns, targetNbOfPilotRuns);
	}

}
