package tests;

import java.util.List;

import derivatives.Derivative;
import derivatives.DerivativePayoffFunction;
import experiments.ProblemDefinition;
import experiments.ProblemDefinition.PriceModel;
import meshmethods.MeshFunctionExpectationEstimator;
import meshmethods.MeshWeights;
import meshmethods.StochasticMesh;
import meshmethods.StochasticMeshAverageDensity;
import meshmethods.StochasticMeshAverageDensityQMC;
import meshmethods.StochasticMeshSingleGrid;
import meshmethods.StochasticMeshSingleGridQMC;
import numerics.ConstantFunction;
import numerics.LinearFunction;
import numerics.OneDimFunction;
import stochprocess.DiscreteExpOULagProcess;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.stat.Tally;

public class TestMeshMethods extends TestBase {

	private static double relativeToleranceForMesh = 0; // 0.05; // 0.10;

	public static TestBase GetInstance() {
		TestMeshMethods test = new TestMeshMethods();
		return test;
	}

	public static void main(String[] args) {
		TestBase tests = GetInstance();

		tests.setRelativeTolerance(relativeToleranceForMesh);
		tests.execute();
	}

	private int nbMeshPaths = 128; // 512;// 256;// 1024;
	private int nbRepsForConfInt = 20;
	private int nbPathsForPriceEstimates = 100000;
	private MarketVectorProcess processForPriceEstimates;
	private double rouletteThreshold = 0.01;

	private TestStochasticProcesses testProcesses;

	protected int nbDataPointsForTesting = 21;

	@Override
	public void defineTestLabel() {
		String label = "mesh methods (N=" + this.nbMeshPaths + ")";
		this.setTestLabel(label);
	}

	@Override
	public void run() {
		this.testProcesses = new TestStochasticProcesses();

		// this.runForNbSteps(1, ProblemDefinition.PriceModel.GBM);
		// this.runForNbSteps(4, ProblemDefinition.PriceModel.GBM);
		// this.runForNbSteps(1, ProblemDefinition.PriceModel.EXP_OU);
		this.runForNbSteps(4, ProblemDefinition.PriceModel.EXP_OU);
		// this.runForNbSteps(1, ProblemDefinition.PriceModel.EXP_OU_SYNC);
		// this.runForNbSteps(2, ProblemDefinition.PriceModel.EXP_OU_SYNC);
		// this.runForNbSteps(1, ProblemDefinition.PriceModel.EXP_OU_INDIC);
		this.runForNbSteps(4, ProblemDefinition.PriceModel.EXP_OU_INDIC);
	}

	private MarketVectorProcess createProcess(PriceModel priceModel) {
		return this.createProcesses(priceModel, false);
	}

	private MarketVectorProcess createProcesses(PriceModel priceModel, boolean usingBrownianBridge) {
		MarketVectorProcess process = this.createSingleProcess(priceModel, usingBrownianBridge);
		this.processForPriceEstimates = this.createSingleProcess(priceModel, false);
		return process;
	}

	private MarketVectorProcess createSingleProcess(PriceModel priceModel, boolean usingBrownianBridge) {
		MarketVectorProcess process = null;
		switch (priceModel) {
		case GBM:
			process = this.testProcesses.getDiscreteGBMProcess(usingBrownianBridge);
			break;
		case EXP_OU:
			process = this.testProcesses.getDiscreteExpOULagProcess(usingBrownianBridge);
			break;
		case EXP_OU_SYNC:
			process = this.testProcesses.getDiscreteExpOUSyncProcess(usingBrownianBridge);
			break;
		case EXP_OU_INDIC:
			process = this.testProcesses.getDiscreteExpOUWithIndicatorsProcess(usingBrownianBridge);
			break;
		default:
		}
		return process;
	}

	private StochasticMesh initializeStochasticMesh(StochasticMesh stochasticMesh, MarketVectorProcess process) {
		double[] obsTimes = process.getObservationTimes();
		int nbSteps = process.getNbObservationTimes();
		stochasticMesh.setObservationTimes(obsTimes, nbSteps);
		stochasticMesh.setNbPaths(this.nbMeshPaths);
		stochasticMesh.generateMesh();
		return stochasticMesh;
	}

	private OneDimFunction newIdentityFunction() {
		double slope = 1;
		double intercept = 0;
		OneDimFunction linearFunction = new LinearFunction(slope, intercept);
		return linearFunction;
	}

	private void runForNbSteps(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		// this.runMeshLikelihoodRatio(nbSteps, priceModel);
		this.runMeshAverageDensity(nbSteps, priceModel);
		// this.runMeshAverageDensityRoulette(nbSteps, priceModel);
		// this.runMeshAverageDensityQMC(nbSteps, priceModel);
		// this.runGridLikelihoodRatio(nbSteps, priceModel);
		// this.runGridLikelihoodRatioQMC(nbSteps, priceModel);
		// this.runGridLikelihoodRatioRouletteQMC(nbSteps, priceModel);
	}

	private void runForNbStepsAndMesh(int nbSteps, StochasticMesh stochasticMesh, MarketVectorProcess process) {
		stochasticMesh.setMarketProcess(process);
		this.testProcesses.setProcessNbSteps(nbSteps, process);
		this.testProcesses.setProcessNbSteps(nbSteps, this.processForPriceEstimates);
		this.initializeStochasticMesh(stochasticMesh, process);

		this.testConstantsEstimation(nbSteps, stochasticMesh);
		this.testFirstMomentEstimation(nbSteps, process, stochasticMesh);
		if (nbSteps > 1) {
			this.testPriceExpectationsAcrossAllStates(nbSteps, stochasticMesh);
		}
		this.testOptionPricing(nbSteps, process, stochasticMesh);
	}

	private void runGridLikelihoodRatio(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		MarketVectorProcess process = this.createProcess(priceModel);
		StochasticMesh stochmeshSG = new StochasticMeshSingleGrid(process);
		this.runForNbStepsAndMesh(nbSteps, stochmeshSG, process);
	}

	private void runGridLikelihoodRatioQMC(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		MarketVectorProcess process = this.createProcess(priceModel);
		PointSet qmcPointSetSG = new SobolSequence(this.nbMeshPaths, process.getDimension());
		StochasticMesh stochmeshSGQ = new StochasticMeshSingleGridQMC(process, qmcPointSetSG);
		this.runForNbStepsAndMesh(nbSteps, stochmeshSGQ, process);
	}

	private void runGridLikelihoodRatioRouletteQMC(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		MarketVectorProcess process = this.createProcess(priceModel);
		PointSet qmcPointSetSGR = new SobolSequence(this.nbMeshPaths, process.getDimension());
		StochasticMesh stochmeshSGRQ = new StochasticMeshSingleGridQMC(process, qmcPointSetSGR);
		MeshWeights weightsSGRQ = stochmeshSGRQ.getMeshWeights();
		weightsSGRQ.useRussianRoulette(this.rouletteThreshold, new MRG32k3a());
		this.runForNbStepsAndMesh(nbSteps, stochmeshSGRQ, process);
	}

	private void runMeshAverageDensity(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		MarketVectorProcess process = this.createProcess(priceModel);
		StochasticMesh stochmeshAD = new StochasticMeshAverageDensity(process);
		this.runForNbStepsAndMesh(nbSteps, stochmeshAD, process);
	}

	private void runMeshAverageDensityQMC(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		boolean usingBrownianBridge = true;
		MarketVectorProcess process = this.createProcesses(priceModel, usingBrownianBridge);
		PointSet qmcPointSetAD = new SobolSequence(this.nbMeshPaths, nbSteps * process.getDimension());
		StochasticMesh stochmeshADQ = new StochasticMeshAverageDensityQMC(process, qmcPointSetAD);
		this.runForNbStepsAndMesh(nbSteps, stochmeshADQ, process);
	}

	private void runMeshAverageDensityRoulette(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		MarketVectorProcess process = this.createProcess(priceModel);
		StochasticMesh stochmeshADR = new StochasticMeshAverageDensity(process);
		MeshWeights weightsADR = stochmeshADR.getMeshWeights();
		weightsADR.useRussianRoulette(this.rouletteThreshold, new MRG32k3a());
		this.runForNbStepsAndMesh(nbSteps, stochmeshADR, process);
	}

	private void runMeshLikelihoodRatio(int nbSteps, ProblemDefinition.PriceModel priceModel) {
		MarketVectorProcess process = this.createProcess(priceModel);
		StochasticMesh stochmeshLR = new StochasticMesh(process);
		this.runForNbStepsAndMesh(nbSteps, stochmeshLR, process);
		stochmeshLR = null;
	}

	private void testConstantsEstimation(int nbSteps, StochasticMesh stochasticMesh) {
		double value = 1;
		OneDimFunction constantFunction = new ConstantFunction(value);

		this.testFunctionExpectationEstimation(nbSteps, stochasticMesh, constantFunction, value);
	}

	private void testFirstMomentEstimation(int nbSteps, MarketVectorProcess process, StochasticMesh stochasticMesh) {
		MarketVector y = process.getInitialMarketVector();
		double expectation = y.values[MarketVector.PRICE];

		OneDimFunction identityFunction = this.newIdentityFunction();

		this.testFunctionExpectationEstimation(nbSteps, stochasticMesh, identityFunction, expectation);
	}

	private void testFunctionExpectationEstimation(int nbSteps, StochasticMesh stochasticMesh, OneDimFunction f,
			double expectation) {
		String meshMethod = stochasticMesh.getDescription();
		String model = stochasticMesh.getMarketVectorProcess().getProcessLabel();
		String source = f.getName() + " expectation estimation, K=" + nbSteps + ", Model:" + model + ", Method:"
				+ meshMethod;

		MeshFunctionExpectationEstimator meshEstimator = new MeshFunctionExpectationEstimator(f);
		meshEstimator.setStochasticMesh(stochasticMesh);

		Tally obsTally = new Tally();
		for (int i = 0; i < this.nbRepsForConfInt; i++) {
			stochasticMesh.generateMesh();
			meshEstimator.computeOverAllSteps();
			double observed = meshEstimator.evaluateAtInitialState();
			obsTally.add(observed);
		}

		this.assertEqualRelative(obsTally, expectation, source);
	}

	private void testOptionPricing(int nbSteps, MarketVectorProcess process, StochasticMesh stochasticMesh) {
		TestDerivatives testDerivatives = new TestDerivatives();
		Derivative callOption = testDerivatives.getCallOption(process);
		// Derivative callOption =
		// testDerivatives.getCallOptionEvaluatedByStochasticMesh(process,
		// stochasticMesh, false);
		OneDimFunction optionPayoff = new DerivativePayoffFunction(callOption);

		MarketVector y0 = process.getInitialMarketVector();
		double expectation = callOption.evaluateAtT0(y0);

		if (process instanceof DiscreteExpOULagProcess) {

			double simulationEstimate = callOption.evaluateBySimulation(this.processForPriceEstimates,
					this.nbPathsForPriceEstimates);
			this.testFunctionExpectationEstimation(nbSteps, stochasticMesh, optionPayoff, simulationEstimate);
		} else {
			this.testFunctionExpectationEstimation(nbSteps, stochasticMesh, optionPayoff, expectation);
		}
	}

	private void testPriceExpectationsAcrossAllStates(int nbSteps, StochasticMesh stochasticMesh) {
		OneDimFunction identityFunction = this.newIdentityFunction();
		String meshMethod = stochasticMesh.getDescription();
		MarketVectorProcess process = stochasticMesh.getMarketVectorProcess();
		String model = process.getProcessLabel();
		String source = "Average of estimated price relative error, K=" + nbSteps + ", Model:" + model + ", Method:"
				+ meshMethod;

		MeshFunctionExpectationEstimator meshEstimator = new MeshFunctionExpectationEstimator(identityFunction);

		meshEstimator.setStochasticMesh(stochasticMesh);
		meshEstimator.computeOverAllSteps();
		double[][] meshEstimates = meshEstimator.getValueArray();
		List<float[][]> statePointsList = stochasticMesh.getStatePointsList();
		double sumAbsoluteErrors = 0;
		for (int k = 1; k < nbSteps; k++) {
			float[][] statePoints = statePointsList.get(k);
			for (int i = 0; i < this.nbMeshPaths; i++) {
				double meshPriceEstimate = meshEstimates[k][i];
				float[] marketVectorValues = statePoints[i];
				// double expectedPrice =
				// marketVectorValues[MarketVector.PRICE];
				double expectedPrice = process.expectedPriceConditionalOnMarketVector(k, marketVectorValues, -1, null);
				double error = meshPriceEstimate - expectedPrice;
				sumAbsoluteErrors += Math.abs(error);
			}
		}
		double targetError = 0;
		double referencePrice = statePointsList.get(0)[0][MarketVector.PRICE];
		double meanAbsoluteError = sumAbsoluteErrors / (referencePrice * this.nbMeshPaths * (nbSteps - 1));
		this.assertEqualAbsolute(meanAbsoluteError, targetError, source);
	}

}
