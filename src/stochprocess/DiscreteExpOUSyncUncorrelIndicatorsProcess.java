package stochprocess;

import org.apache.commons.math3.util.FastMath;

import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stochprocess.MultivariateBrownianMotion;

public class DiscreteExpOUSyncUncorrelIndicatorsProcess extends DiscreteExpOUSyncUncorrelProcess {

	public static double getInnovationLoading(double indicatorLoading) {
		double a = FastMath.sqrt(1 - indicatorLoading * indicatorLoading);
		return a;
	}

	protected RandomStream streamForInitialIndicatorValues;
	private double[] indicatorLoadings;
	private float[][] indicatorPaths;
	private double[] innovationLoadings;

	private int nbIndicators = 2;

	public DiscreteExpOUSyncUncorrelIndicatorsProcess(float[] y0Array, double[] mu, double[] aVec, double[] bVec,
			double[] sigVec, MultivariateBrownianMotion multiBM, double[] indicatorLoadings,
			RandomStream streamForInitialIndicatorValues) {
		super(y0Array, mu, aVec, bVec, sigVec, multiBM);
		this.indicatorLoadings = indicatorLoadings;
		this.streamForInitialIndicatorValues = streamForInitialIndicatorValues;

		this.initializeLoadings();
	}

	@Override
	public Object clone() {
		DiscreteExpOUSyncUncorrelIndicatorsProcess clone = new DiscreteExpOUSyncUncorrelIndicatorsProcess(this.y0Array,
				this.mu, this.aVec, this.bVec, this.sigVec, this.multiDimensionalBrownianMotion, this.indicatorLoadings,
				this.streamForInitialIndicatorValues);

		return clone;
	}

	@Override
	public double density(float[] x1, float[] x2, int t1) {
		// Valid only for ONE step transitions, i.e. t2 = t1 + 1
		// Assumes: time steps (t2-t1) are constants
		// dim=4 : S, V, I1, I2
		double V1 = x1[1]; // vol @ t = 0
		double V2 = x2[1]; // vol @ t = 1
		double indicator11 = x1[2];
		double indicator12 = x1[3];
		double indicator21 = x2[2];
		double indicator22 = x2[3];

		double dT = this.dt[t1];
		double returnVol = this.getReturnVolatilityForDensity(V1, V2);
		double volSqrtdT = returnVol * this.sqrtdt[t1];

		double indicatorTerm2 = this.volOfLogVolFactor * this.indicatorLoadings[1] * indicator12;

		double temp = Math.log(V2)
				- (this.logLongTermVolatility * (1 - this.expFactor[t1]) + Math.log(V1) * this.expFactor[t1]);
		this.densityXVector[1] = (temp - indicatorTerm2) / (this.volOfLogVolFactor * this.innovationLoadings[1]);

		double drift = (this.mu[0] - 0.5 * returnVol * returnVol) * dT;
		double indicatorTerm1 = volSqrtdT * this.indicatorLoadings[0] * indicator11;
		this.densityXVector[0] = ((Math.log(x2[0] / x1[0])) - drift - indicatorTerm1)
				/ (volSqrtdT * this.innovationLoadings[0]);

		this.densityXVector[2] = indicator21;
		this.densityXVector[3] = indicator22;

		double innerProduct = this.computeInnerProduct(this.densityXVector);
		double result = this.densityNormalizationFactor * FastMath.exp(-0.5 * innerProduct)
				/ (x2[0] * x2[1] * volSqrtdT);

		return result;
	}

	public double densityUncond(float[] x1, float[] x2, int t1) {
		double value = super.density(x1, x2, t1);
		return value;
	}

	@Override
	public float expectedPriceConditionalOnMarketVector(int step0, float[] values0, int step1, float[] values1) {
		float volatility = values0[MarketVector.VOL];
		float indicator0 = values0[MarketVector.INDIC0];
		float s = values0[MarketVector.PRICE];
		double a = volatility * this.sqrtdt[step0] * this.indicatorLoadings[0];
		double b = -0.5 * a * a;
		double c = a * indicator0 + b;
		return (float) (s * FastMath.exp(c));
	}

	@Override
	public float[][] generateBasicProcessesPaths() {
		double[] bmPath = this.multiDimensionalBrownianMotion.generatePath();
		this.generateInitialIndicatorValues();

		for (int i = 1; i <= this.d; i++) {

			double indicator00 = this.indicatorPaths[0][i - 1] * this.sqrtdt[i - 1];
			double indicator01 = this.indicatorPaths[1][i - 1] * this.sqrtdt[i - 1];

			double indicator10 = bmPath[this.nbProcessDim * (i) + (0)] - bmPath[this.nbProcessDim * (i - 1) + (0)];
			double innovation00 = bmPath[this.nbProcessDim * (i) + (1)] - bmPath[this.nbProcessDim * (i - 1) + (1)];
			double indicator11 = bmPath[this.nbProcessDim * (i) + (2)] - bmPath[this.nbProcessDim * (i - 1) + (2)];
			double innovation01 = bmPath[this.nbProcessDim * (i) + (3)] - bmPath[this.nbProcessDim * (i - 1) + (3)];

			double z0SqrtDt = this.indicatorLoadings[0] * indicator00 + this.innovationLoadings[0] * innovation00;
			double z1SqrtDt = this.indicatorLoadings[1] * indicator01 + this.innovationLoadings[1] * innovation01;

			double lnVolj = this.logVolatilityPath[i - 1];
			double z1 = z1SqrtDt / this.sqrtdt[i - 1];
			this.logVolatilityPath[i] = this.logLongTermVolatility * (1 - this.expFactor[i - 1])
					+ lnVolj * this.expFactor[i - 1] + this.sigVec[0] * this.expFactorSp[i - 1] * z1;

			double volatility = FastMath.exp(this.logVolatilityPath[i]);
			this.volatilityPath[i] = volatility;

			double returnVolatility = this.getReturnVolatilityForPaths(i);
			double drift = this.mu[0] * this.dt[i - 1] - 0.5 * returnVolatility * returnVolatility * this.dt[i - 1];
			double periodReturn = drift + returnVolatility * z0SqrtDt;

			this.basicProcessesPaths[0][i] = (float) (this.basicProcessesPaths[0][i - 1] * FastMath.exp(periodReturn));
			this.basicProcessesPaths[1][i] = (float) volatility;

			this.path[i] = this.basicProcessesPaths[0][i];

			this.indicatorPaths[0][i] = (float) (indicator10 / this.sqrtdt[i - 1]);
			this.indicatorPaths[1][i] = (float) (indicator11 / this.sqrtdt[i - 1]);

			this.basicProcessesPaths[2][i] = this.indicatorPaths[0][i];
			this.basicProcessesPaths[3][i] = this.indicatorPaths[1][i];
		}

		return this.basicProcessesPaths;
	}

	@Override
	public double getVarianceAtStep(int i) {
		// TODO
		return Double.NaN;
	}

	@Override
	public MarketVector newMarketVector(float[] y) {
		return new MarketVectorExpOUWithIndicators(y);
	}

	private void defineInnovationLoadings() throws Exception {
		this.innovationLoadings = new double[this.nbIndicators];
		for (int i = 0; i < this.nbIndicators; i++) {
			double b = this.indicatorLoadings[i];
			double a = getInnovationLoading(b);
			this.innovationLoadings[i] = a;
			if (a < 0 || b < 0) {
				String message = "(DiscreteExpOUWithIndicatorsProcess) Invalid loading for indicator (#" + i + ") : "
						+ b;
				throw new Exception(message);
			}
		}
	}

	private void initializeLoadings() {
		try {
			this.defineInnovationLoadings();
		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
		}
	}

	@Override
	protected void computeDensityNormalizationFactor() {
		this.densityNormalizationFactor = 1
				/ (this.twoNrootPi * this.volOfLogVolFactor * this.innovationLoadings[0] * this.innovationLoadings[1]);
	}

	protected void generateInitialIndicatorValues() {
		this.indicatorPaths[0][0] = (float) NormalGen.nextDouble(this.streamForInitialIndicatorValues, 0, 1);
		this.indicatorPaths[1][0] = (float) NormalGen.nextDouble(this.streamForInitialIndicatorValues, 0, 1);

		this.basicProcessesPaths[2][0] = this.indicatorPaths[0][0];
		this.basicProcessesPaths[3][0] = this.indicatorPaths[1][0];
	}

	@Override
	protected void initializePath() {
		super.initializePath();

		this.indicatorPaths = new float[this.nbIndicators][this.d + 1];
	}

	@Override
	protected void setParameters(double[] aVec, double[] bVec, double[] sigVec, double priceVolCorrelation) {
		super.setParameters(aVec, bVec, sigVec, priceVolCorrelation);
		this.nbProcessDim = 4;
		this.processLabel = "expOU-wIndic";
	}

}
