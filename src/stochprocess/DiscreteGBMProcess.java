package stochprocess;

import org.apache.commons.math3.util.FastMath;

import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stochprocess.MultivariateBrownianMotion;
import umontreal.ssj.stochprocess.MultivariateBrownianMotionBridge;

public class DiscreteGBMProcess extends MarketVectorProcess {
	protected double densityNormalizationFactor = 1;
	protected double[] mu;
	protected double[][] multivariateBrownianMotionPath;
	protected double rootTwoPi;
	protected double twoNrootPi;
	protected double twoPi;

	public DiscreteGBMProcess(float[] y0Array, MultivariateBrownianMotion multiBM) {
		this.y0Array = y0Array;
		this.mu = multiBM.getMu();
		this.setMultivariateBrownianMotion(multiBM, multiBM.getStream());
		this.defineProcessLabel(multiBM);
	}

	@Override
	public Object clone() {
		DiscreteGBMProcess clone = new DiscreteGBMProcess(this.y0Array, this.multiDimensionalBrownianMotion);
		return clone;
	}

	@Override
	public double density(float[] x1, float[] x2, int t1) {
		// Valid only for ONE step transitions, i.e. t2 = t1 + 1
		// Assumes equal time steps

		double dT = this.dt[t1];
		double sigma = this.y0Array[1];
		double sigma2 = sigma * sigma;
		double sigma2DT = sigma2 * dT;

		this.densityXVector[0] = (Math.log(x2[0] / x1[0]) - (this.mu[0] * dT - 0.5 * sigma2 * dT));
		double distSquared = this.densityXVector[0] * this.densityXVector[0] / sigma2DT;
		double result = this.densityNormalizationFactor * FastMath.exp(-0.5 * (distSquared)) / x2[0];

		return result;
	}

	@Override
	public float[][] generateBasicProcessesPaths() {
		double[] bmPath = this.multiDimensionalBrownianMotion.generatePath();
		for (int i = 1; i <= this.d; i++) {
			double dW = bmPath[i] - bmPath[i - 1]; // with variance = dt, not 1
			double periodReturn = 0;
			float sigma = this.basicProcessesPaths[1][i - 1];
			double sigma2 = sigma * sigma;

			periodReturn += sigma * dW - 0.5 * sigma2 * this.dt[i - 1];

			this.basicProcessesPaths[0][i] = (float) (this.basicProcessesPaths[0][i - 1] * FastMath.exp(periodReturn));
			this.basicProcessesPaths[1][i] = sigma;

			this.path[i] = this.basicProcessesPaths[0][i];
		}
		return this.basicProcessesPaths;
	}

	@Override
	public double getExpectationAtStep(int k) {
		double s0 = this.y0Array[0];
		double expectation = s0;

		return expectation;
	}

	@Override
	public double getVarianceAtStep(int k) {
		// Variance of the value of the process starting at 0, observed at k
		double dt = this.t[k] - this.t[0];
		double sigma = this.basicProcessesPaths[1][0];
		double sigma2 = sigma * sigma;
		double x2 = FastMath.exp(sigma2 * dt) - 1.0;
		double s0 = this.y0Array[0];
		double variance = x2 * s0 * s0;
		return variance;
	}

	@Override
	public double integrateDensity(float[] x1, int t1, int nbPts) {
		double result = 0.0;

		float[] x2 = new float[1];

		float x;
		double minPt = -x1[0] + 0.001;
		double maxPt = 5 * this.sqrtdt[t1] * x1[0];
		double dP = (maxPt - minPt) / (nbPts - 1);
		double weight = dP;

		for (int i = 0; i < nbPts; i++) {
			x = (float) (x1[0] + minPt + i * dP);
			x2[0] = x;

			result += weight * this.density(x1, x2, t1);
		}

		return result;
	}

	@Override
	public MarketVector newMarketVector(float[] y) {
		MarketVectorGBM marketVector = this.getNewMarketVectorGBM(y);
		return marketVector;
	}

	private void defineProcessLabel(MultivariateBrownianMotion multiBM) {
		this.processLabel = "GBM";
		if (multiBM instanceof MultivariateBrownianMotionBridge) {
			this.processLabel += "_B";
		}
	}

	private MarketVectorGBM getNewMarketVectorGBM(float[] y) {
		MarketVectorGBM marketVector = new MarketVectorGBM(y);
		return marketVector;
	}

	private void setMultivariateBrownianMotion(MultivariateBrownianMotion multiBM, RandomStream stream) {
		this.uStream = stream;
		this.multiDimensionalBrownianMotion = multiBM;

		this.nbProcessDim = multiBM.getDimension();
	}

	protected void computeConstants() {
		this.twoPi = 2 * FastMath.PI;
		this.rootTwoPi = FastMath.sqrt(this.twoPi);
		this.twoNrootPi = FastMath.pow(this.twoPi, 0.5 * this.nbProcessDim);
	}

	protected void computeDensityNormalizationFactor() {
		double V1 = this.y0Array[1] * this.y0Array[1];
		double V1squaredDT = V1 * this.dt[0];
		this.densityNormalizationFactor = 1 / FastMath.sqrt(this.twoPi * V1squaredDT);
	}

	@Override
	protected void init() {
		super.init();

		this.initializeTimeSteps();
		this.computeConstants();
		this.computeDensityNormalizationFactor();
		this.initializePath();
	}

	protected void initializePath() {
		int basicProcessPathsDim = this.getMarketVectorDimension();
		this.basicProcessesPaths = new float[basicProcessPathsDim][this.d + 1];
		this.path = new double[this.d + 1];
		this.marketVectorPath = new MarketVector[this.d + 1];
		this.multivariateBrownianMotionPath = new double[this.nbProcessDim][this.d + 1];
		this.densityXVector = new double[this.nbProcessDim];

		for (int k = 0; k <= this.d; k++) {
			float[] y = new float[basicProcessPathsDim];
			this.basicProcessesPaths[0][k] = this.y0Array[0];
			this.path[k] = this.y0Array[0];
			y[0] = this.y0Array[0];
			y[1] = this.basicProcessesPaths[1][k] = this.y0Array[1];

			this.marketVectorPath[k] = this.newMarketVector(y);
		}
	}

	protected void initializeTimeSteps() {
		this.dt = new double[this.d];
		this.sqrtdt = new double[this.d];

		for (int j = 0; j < this.d; j++) {
			double dt2 = this.t[j + 1] - this.t[j];
			this.dt[j] = dt2;
			this.sqrtdt[j] = FastMath.sqrt(dt2);
		}
	}

}
