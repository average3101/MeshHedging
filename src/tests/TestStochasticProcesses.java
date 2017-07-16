package tests;

import experiments.ProblemDefinition;
import stochprocess.DiscreteExpOULagIndicatorsProcess;
import stochprocess.DiscreteExpOULagProcess;
import stochprocess.DiscreteExpOUSyncUncorrelProcess;
import stochprocess.DiscreteGBMProcess;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.MultivariateBrownianMotion;
import umontreal.ssj.stochprocess.MultivariateBrownianMotionBridge;

public class TestStochasticProcesses extends TestBase {

	private double[][] correlationMatrix;
	private double horizon = 0.5;
	private double indicatorLoading1 = 0.9; // 0.5;
	private double indicatorLoading2 = 0; // 0.5;
	private double meanReversion = 2.6;
	private double meanVol = 0.20;
	private int nbPaths = 100000;
	private int nbPointsForDensityIntegral = 1000;
	private double priceVolCorrelation = -0.5;
	private double rfRate = 0.0;
	private double s0 = 10.0;
	private double sigma = 0.4; // 0.40;
	private double[] unitVec;
	private double volOfLogVar = 0.6; // 1.2; // 2.5; // 0.001; // 2.5;
	private float[] y0Array;
	private double[] zeroVec;

	public static TestBase GetInstance() {
		TestStochasticProcesses testStochasticProcesses = new TestStochasticProcesses();

		return testStochasticProcesses;
	}

	public static void main(String[] args) {
		TestBase tests = GetInstance();
		tests.execute();
	}

	@Override
	public void defineTestLabel() {
		this.setTestLabel("stochprocess");
	}

	public DiscreteExpOULagProcess getDiscreteExpOULagProcess() {
		return this.getDiscreteExpOUProcess(false, false, false);
	}

	public DiscreteExpOULagProcess getDiscreteExpOULagProcess(boolean usingBrownianBridge) {
		return this.getDiscreteExpOUProcess(usingBrownianBridge, false, false);
	}

	public DiscreteExpOULagProcess getDiscreteExpOUProcess(boolean usingBrownianBridge, boolean usingIndicators,
			boolean usingSynchronousVol) {
		this.generateParameterArrays(usingIndicators);

		int dimensionVol = 1;
		double[] muVec = new double[dimensionVol];
		double[] aVec = new double[dimensionVol];
		double[] bVec = new double[dimensionVol];
		double[] sigVec = new double[dimensionVol];

		muVec[0] = this.rfRate;
		aVec[0] = this.meanVol; // simParams.sigma1;
		bVec[0] = this.meanReversion;
		sigVec[0] = this.volOfLogVar;

		DiscreteExpOULagProcess process;
		if (usingIndicators) {
			double[] indicatorLoadings = new double[2];
			indicatorLoadings[0] = this.indicatorLoading1;
			indicatorLoadings[1] = this.indicatorLoading2;
			int processDim = 4;
			MultivariateBrownianMotion multiBM = this.getBrownianMotionProcess(usingBrownianBridge, processDim);
			process = new DiscreteExpOULagIndicatorsProcess(this.y0Array, muVec, aVec, bVec, sigVec,
					this.priceVolCorrelation, multiBM, indicatorLoadings, new MRG32k3a());
		} else {
			int processDim = 2;
			MultivariateBrownianMotion multiBM = this.getBrownianMotionProcess(usingBrownianBridge, processDim);
			if (usingSynchronousVol) {
				process = new DiscreteExpOUSyncUncorrelProcess(this.y0Array, muVec, aVec, bVec, sigVec, multiBM);
			} else {
				process = new DiscreteExpOULagProcess(this.y0Array, muVec, aVec, bVec, sigVec, this.priceVolCorrelation,
						multiBM);
			}
		}

		return process;
	}

	public DiscreteExpOUSyncUncorrelProcess getDiscreteExpOUSyncProcess() {
		return (DiscreteExpOUSyncUncorrelProcess) this.getDiscreteExpOUProcess(false, false, true);
	}

	public DiscreteExpOUSyncUncorrelProcess getDiscreteExpOUSyncProcess(boolean usingBrownianBridge) {
		return (DiscreteExpOUSyncUncorrelProcess) this.getDiscreteExpOUProcess(usingBrownianBridge, false, true);
	}

	public DiscreteExpOULagIndicatorsProcess getDiscreteExpOUWithIndicatorsProcess() {
		return (DiscreteExpOULagIndicatorsProcess) this.getDiscreteExpOUProcess(false, true, true);
	}

	public DiscreteExpOULagIndicatorsProcess getDiscreteExpOUWithIndicatorsProcess(boolean usingBrownianBridge) {
		return (DiscreteExpOULagIndicatorsProcess) this.getDiscreteExpOUProcess(usingBrownianBridge, true, true);
	}

	public DiscreteGBMProcess getDiscreteGBMProcess() {
		return this.getDiscreteGBMProcess(false);
	}

	public DiscreteGBMProcess getDiscreteGBMProcess(boolean usingBrownianBridge) {
		this.generateParameterArrays(false);

		RandomStream uStream = new MRG32k3a();
		NormalDist normDist = new NormalDist();
		NormalGen normGen = new NormalGen(uStream, normDist);

		int dimension = 1;

		MultivariateBrownianMotion multiBM = this.getBrownianMotionProcess(usingBrownianBridge, normGen, dimension);
		DiscreteGBMProcess process = new DiscreteGBMProcess(this.y0Array, multiBM);

		return process;
	}

	@Override
	public void run() {
		int nbSteps = 4;
		DiscreteGBMProcess discreteBSMProcess = this.getDiscreteGBMProcess();
		this.setProcessNbSteps(nbSteps, discreteBSMProcess);
		this.testMarketVectorProcess(discreteBSMProcess);

		nbSteps = 1;
		DiscreteExpOULagProcess discreteExpOULagProcess = this.getDiscreteExpOULagProcess();
		this.setProcessNbSteps(nbSteps, discreteExpOULagProcess);
		this.testMarketVectorProcess(discreteExpOULagProcess);

		nbSteps = 1;
		DiscreteExpOUSyncUncorrelProcess discreteExpOUSyncProcess = this.getDiscreteExpOUSyncProcess();
		this.setProcessNbSteps(nbSteps, discreteExpOUSyncProcess);
		this.testMarketVectorProcess(discreteExpOUSyncProcess);

		nbSteps = 1;
		DiscreteExpOULagIndicatorsProcess discreteExpOUWithIndicatorsProcess = this
				.getDiscreteExpOUWithIndicatorsProcess();
		this.setProcessNbSteps(nbSteps, discreteExpOUWithIndicatorsProcess);
		this.testMarketVectorProcess(discreteExpOUWithIndicatorsProcess);
	}

	public void setProcessNbSteps(int nbSteps, MarketVectorProcess process) {
		double dt = this.horizon / nbSteps;
		process.setObservationTimes(dt, nbSteps);
	}

	private void generateParameterArrays(boolean usingIndicators) {
		int processDim;
		if (usingIndicators) {
			processDim = 4;
			this.correlationMatrix = ProblemDefinition.defineCorrelationMatrixForExpOUWithIndicators(
					this.priceVolCorrelation, this.indicatorLoading1, this.indicatorLoading2);
		} else {
			processDim = 2;
			this.correlationMatrix = ProblemDefinition.defineCorrelationMatrixForExpOU(this.priceVolCorrelation);
		}

		this.zeroVec = new double[processDim];
		this.unitVec = new double[processDim];
		this.y0Array = new float[processDim];

		for (int i = 0; i < processDim; i++) {
			this.zeroVec[i] = 0.0;
			this.unitVec[i] = 1.0;
			this.y0Array[i] = 0;
		}

		this.y0Array[0] = (float) this.s0;
		this.y0Array[1] = (float) this.sigma;
	}

	private MultivariateBrownianMotion getBrownianMotionProcess(boolean usingBrownianBridge, int dimensionPriceAndVol) {
		RandomStream uStream = new MRG32k3a();
		NormalDist normDist = new NormalDist();
		NormalGen normGen = new NormalGen(uStream, normDist);
		MultivariateBrownianMotion multiBM = this.getBrownianMotionProcess(usingBrownianBridge, normGen,
				dimensionPriceAndVol);
		return multiBM;
	}

	private MultivariateBrownianMotion getBrownianMotionProcess(boolean usingBrownianBridge, NormalGen normGen,
			int dimension) {
		MultivariateBrownianMotion multiBM;
		if (usingBrownianBridge) {
			multiBM = new MultivariateBrownianMotionBridge(dimension, this.zeroVec, this.zeroVec, this.unitVec,
					this.correlationMatrix, normGen);
		} else {
			multiBM = new MultivariateBrownianMotion(dimension, this.zeroVec, this.zeroVec, this.unitVec,
					this.correlationMatrix, normGen);
		}
		return multiBM;
	}

	private boolean testFirstMoment(MarketVectorProcess process, Tally tally, int nbSteps) {
		String processType = process.getClass().toString();
		String source = processType + " average";
		double sampleAverage = tally.average();
		double expectation = process.getExpectationAtStep(nbSteps);
		boolean resultExpectation = this.assertEqualRelative(sampleAverage, expectation, source);
		return resultExpectation;
	}

	private void testMarketVectorProcess(MarketVectorProcess marketVectorProcess) {
		if (marketVectorProcess.getDimension() <= 2) {
			this.testMarketVectorProcessDensities(marketVectorProcess);
		}
		this.testMarketVectorProcessMoments(marketVectorProcess);
	}

	private boolean testMarketVectorProcessDensities(MarketVectorProcess process) {
		int t0 = 0;

		String processType = process.getClass().toString();
		String source = processType + " density integral";
		float[] x = new float[2];
		x[0] = (float) this.s0;
		x[1] = (float) this.sigma;
		double densityIntegralApprox = process.integrateDensity(x, t0, this.nbPointsForDensityIntegral);
		double densityIntegralExact = 1.0;
		boolean resultDensityIntegral = this.assertEqualRelative(densityIntegralApprox, densityIntegralExact, source);

		return resultDensityIntegral;
	}

	private boolean testMarketVectorProcessMoments(MarketVectorProcess process) {
		Tally tally = new Tally();
		int nbSteps = process.getNbObservationTimes();
		for (int i = 0; i < this.nbPaths; i++) {
			float[][] basicPaths = process.generateBasicProcessesPaths();
			float sk = basicPaths[0][nbSteps];
			tally.add(sk);
		}

		boolean resultExpectation = this.testFirstMoment(process, tally, nbSteps);

		boolean resultVariance = true;
		if (!(process instanceof DiscreteExpOUSyncUncorrelProcess)) {
			resultVariance = this.testSecondMoment(process, tally, nbSteps);
		}

		return resultExpectation && resultVariance;
	}

	private boolean testSecondMoment(MarketVectorProcess process, Tally tally, int nbSteps) {
		String processType = process.getClass().toString();
		String source = processType + " variance";
		double sampleVariance = tally.variance();
		double variance = process.getVarianceAtStep(nbSteps);
		boolean resultVariance = this.assertEqualRelative(sampleVariance, variance, source);
		return resultVariance;
	}

}
