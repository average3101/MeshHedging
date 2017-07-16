package stochprocess;

import org.apache.commons.math3.util.FastMath;

import umontreal.ssj.stochprocess.MultivariateBrownianMotion;
import umontreal.ssj.stochprocess.MultivariateBrownianMotionBridge;

public class DiscreteExpOULagProcess extends DiscreteGBMProcess {

	protected double[] aVec;
	protected double[] bVec;
	protected double[] expFactor;
	protected double[] expFactorSp;
	protected double logLongTermVolatility;
	protected double[] logVolatilityPath;
	protected double[] sigVec;
	protected double[] volatilityPath;
	protected double volOfLogVolFactor;
	protected double priceVolCorrelation;
	protected double correlationScalingFactor;

	public DiscreteExpOULagProcess(float[] y0Array, double[] mu, double[] aVec, double[] bVec, double[] sigma,
			double priceVolCorrelation, MultivariateBrownianMotion multiBM) {
		super(y0Array, multiBM);

		this.setParameters(aVec, bVec, sigma, priceVolCorrelation);

		if (multiBM instanceof MultivariateBrownianMotionBridge) {
			this.processLabel += "_B";
		}
	}

	@Override
	public Object clone() {
		DiscreteExpOULagProcess clone = new DiscreteExpOULagProcess(this.y0Array, this.mu, this.aVec, this.bVec,
				this.sigVec, this.priceVolCorrelation, this.multiDimensionalBrownianMotion);
		return clone;
	}

	@Override
	public double density(float[] x1, float[] x2, int t1) { // , int t2) {
		// Valid only for ONE step transitions, i.e. t2 = t1 + 1
		// Assumes: time steps (t2-t1) are constants
		// dim=2 : S, V
		double result = 0;
		try {
			double V1 = x1[1]; // volatility @ t = 0
			double V2 = x2[1]; // volatility @ t = 1

			double dT = this.dt[t1];
			double returnVol = this.getReturnVolatilityForDensity(V1, V2);
			double volSqrtdT = returnVol * this.sqrtdt[t1];

			double z1 = (Math.log(V2)
					- (this.logLongTermVolatility * (1 - this.expFactor[t1]) + Math.log(V1) * this.expFactor[t1]))
					/ this.volOfLogVolFactor;
			this.densityXVector[1] = z1;

			double drift = (this.mu[0] - 0.5 * returnVol * returnVol) * dT;

			double z0 = (Math.log(x2[0] / x1[0]) - drift) / volSqrtdT;
			double w0 = (z0 - this.priceVolCorrelation * z1) / this.correlationScalingFactor;
			this.densityXVector[0] = w0;

			double innerProduct = this.computeInnerProduct(this.densityXVector);
			result = this.densityNormalizationFactor * FastMath.exp(-0.5 * innerProduct) / (x2[0] * x2[1] * volSqrtdT);
		} catch (Exception e) {
			System.out.println("1");
		}

		return result;
	}

	@Override
	public float[][] generateBasicProcessesPaths() {
		double[] bmPath = this.multiDimensionalBrownianMotion.generatePath();
		for (int i = 1; i <= this.d; i++) {
			double w0SqrtDt = bmPath[this.nbProcessDim * (i) + (0)] - bmPath[this.nbProcessDim * (i - 1) + (0)];
			double w1SqrtDt = bmPath[this.nbProcessDim * (i) + (1)] - bmPath[this.nbProcessDim * (i - 1) + (1)];

			double z0SqrtDt = this.correlationScalingFactor * w0SqrtDt + this.priceVolCorrelation * w1SqrtDt;
			double z1SqrtDt = w1SqrtDt;

			double lnVolj = this.logVolatilityPath[i - 1];

			double z1 = z1SqrtDt / this.sqrtdt[i - 1];
			this.logVolatilityPath[i] = this.logLongTermVolatility * (1 - this.expFactor[i - 1])
					+ lnVolj * this.expFactor[i - 1] + this.volOfLogVolFactor * z1;

			double volatility = Math.exp(this.logVolatilityPath[i]);
			this.volatilityPath[i] = volatility;

			double returnVolatility = this.getReturnVolatilityForPaths(i);
			double drift = this.mu[0] * this.dt[i - 1] - 0.5 * returnVolatility * returnVolatility * this.dt[i - 1];
			double periodReturn = drift + returnVolatility * z0SqrtDt;

			this.basicProcessesPaths[0][i] = (float) (this.basicProcessesPaths[0][i - 1] * Math.exp(periodReturn));
			this.basicProcessesPaths[1][i] = (float) volatility;

			this.path[i] = this.basicProcessesPaths[0][i];
		}
		return this.basicProcessesPaths;
	}

	@Override
	public double integrateDensity(float[] x1, int t1, int nbPts) {
		double result = 0.0;

		float[] x2 = new float[2];

		float x, y;
		double minPt = -x1[0] + 10e-10;
		double maxPt = 10 * this.sqrtdt[t1] * x1[0] * FastMath.sqrt(x1[1]);
		double minPtV = -x1[1] + 10e-10;
		double maxPtV = 10 * this.sqrtdt[t1] * FastMath.sqrt(x1[1]) * this.sigVec[0];
		double dP = (maxPt - minPt) / (nbPts - 1);
		double dV = (maxPtV - minPtV) / (nbPts - 1);
		double weight = dP * dV;

		for (int i = 0; i < nbPts; i++) {
			x = (float) (x1[0] + minPt + i * dP);
			for (int j = 0; j < nbPts; j++) {
				y = (float) (x1[1] + minPtV + j * dV);

				x2[0] = x;
				x2[1] = y;

				result += weight * this.density(x1, x2, t1);
			}
		}

		return result;
	}

	@Override
	public MarketVector newMarketVector(float[] y) {
		return new MarketVectorExpOU(y);
	}

	@Override
	public void setInitialVolatilityToLongRunAverage() {
		double initialVol = this.y0Array[MarketVector.VOL];
		double z = this.sigVec[0] * this.sigVec[0] / (2 * this.bVec[0]);
		double factor = FastMath.exp(z);

		double longRunVolatility = initialVol * factor;
		this.y0Array[1] = (float) longRunVolatility;
	}

	@Override
	protected void computeConstants() {
		super.computeConstants();

		this.expFactor = new double[this.d];
		this.expFactorSp = new double[this.d];

		for (int j = 0; j < this.d; j++) {
			this.expFactor[j] = FastMath.exp(-this.dt[j] * this.bVec[0]);
			this.expFactorSp[j] = FastMath.sqrt((1 - this.expFactor[j] * this.expFactor[j]) / (2 * this.bVec[0]));
		}

		this.volOfLogVolFactor = this.sigVec[0] * this.expFactorSp[0];
	}

	@Override
	protected void computeDensityNormalizationFactor() {
		this.densityNormalizationFactor = 1.0
				/ (this.twoNrootPi * this.volOfLogVolFactor * this.correlationScalingFactor);
	}

	protected double computeInnerProduct(double[] vector) {
		double innerProduct = 0.0;
		for (int i = 0; i < this.nbProcessDim; i++) {
			innerProduct += vector[i] * vector[i];
		}
		return innerProduct;
	}

	protected double getReturnVolatilityForDensity(double laggedVolatility, double synchedVolatility) {
		return laggedVolatility;
	}

	protected double getReturnVolatilityForPaths(int i) {
		return this.volatilityPath[i - 1];
	}

	@Override
	protected void initializePath() {
		super.initializePath();

		this.volatilityPath = new double[this.d + 1];
		this.logVolatilityPath = new double[this.d + 1];

		for (int k = 0; k <= this.d; k++) {
			this.volatilityPath[k] = this.y0Array[1];
			this.logVolatilityPath[k] = Math.log(this.volatilityPath[k]);
		}
	}

	protected void setParameters(double[] aVec, double[] bVec, double[] sigVec, double priceVolCorrelation) {
		this.nbProcessDim = 2;
		this.processLabel = "expOU-lag";

		this.aVec = aVec;
		this.bVec = bVec;
		this.sigVec = sigVec;
		this.logLongTermVolatility = Math.log(aVec[0]);
		this.priceVolCorrelation = priceVolCorrelation;
		this.correlationScalingFactor = FastMath.sqrt(1 - priceVolCorrelation * priceVolCorrelation);
	}

}
