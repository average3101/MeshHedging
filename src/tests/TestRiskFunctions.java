package tests;

import derivatives.Derivative;
import derivatives.PutOption;
import dynprog.DynamicHedgingProgram;
import dynprog.DynamicHedgingProgramUsingLQApproximation;
import dynprog.ExponentialLossFunction;
import dynprog.HedgingBoundary;
import dynprog.LossFunction;
import dynprog.ProportionalCosts;
import dynprog.QFunction;
import dynprog.RiskFunctionApproximation;
import dynprog.RiskFunctionAtFirstStep;
import dynprog.State;
import dynprog.TransactionCosts;
import experiments.ProblemDefinition;
import experiments.SimulationParameters;
import meshmethods.MeshWeights;
import meshmethods.StochasticMesh;
import meshmethods.StochasticMeshAverageDensity;
import numerics.OneDimFunction;
import policies.HedgingPortfolio;
import stochprocess.DiscreteGBMProcess;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.rng.RandomStream;

public class TestRiskFunctions extends TestBase {

	public static TestBase GetInstance() {
		TestRiskFunctions test = new TestRiskFunctions();

		return test;
	}

	public static void main(String[] args) {
		TestBase tests = GetInstance();
		tests.execute();
	}

	private static double slope(OneDimFunction f, double x0, double x1) {
		double v0 = f.evaluate(x0);
		double v1 = f.evaluate(x1);
		double m = (v1 - v0) / (x1 - x0);
		return m;
	}

	private double costParamA = 0.0;
	private double costParamB = 0.02;
	private double divYield = 0;
	private double maturity;
	private double maxLossForTest = 100.0; // 1.0;
	private int nbDataPointsForTesting = 21;
	private int nbMeshPaths = 100;
	private double[] obsTimes;
	private double rfRate = 0;
	private double riskAversion = 1.0;
	private double strike;
	private double toleranceForSlope = 10e-4;
	private double u0 = -1.0;

	private double u1 = 1.0;

	private RandomStream uWeightRoulette = null;

	private double weightCutoff = 0.0;

	public boolean assertConvexity(OneDimFunction f, double u0, double u1, int nbPoints, String source) {
		double du = (u1 - u0) / (nbPoints - 1);
		double lastSlope = -Double.MAX_VALUE;
		boolean isConvex = true;
		double notConvexNearPoint = -Double.MAX_VALUE;
		for (int i = 1; i < nbPoints; i++) {
			double x0 = u0 + (i - 1) * du;
			double x1 = x0 + du;

			double slope = slope(f, x0, x1);

			if (lastSlope >= slope + this.toleranceForSlope) {
				isConvex = false;
				notConvexNearPoint = x0;
				break;
			}

			lastSlope = slope;
		}

		String message = "Function is non-convex near point x0:" + notConvexNearPoint;
		return this.assertTrue(isConvex, message, source);
	}

	@Override
	public void defineTestLabel() {
		this.setTestLabel("risk functions");
	}

	public DynamicHedgingProgram getInitializedDynamicProgram(int nbSteps) {
		TestStochasticProcesses testProcesses = new TestStochasticProcesses();
		TestDerivatives testDerivatives = new TestDerivatives();

		DiscreteGBMProcess discreteGBMProcess = testProcesses.getDiscreteGBMProcess();
		testProcesses.setProcessNbSteps(nbSteps, discreteGBMProcess);
		Derivative callOption = testDerivatives.getCallOption(discreteGBMProcess);

		DynamicHedgingProgram program = this.defineDynamicHedgingProgramUsingLQApproximation(callOption,
				discreteGBMProcess);
		program.initializeQFunctions();
		program.initializeStructuresDependentOnMeshConstruction();

		return program;
	}

	public Derivative getPutOption(MarketVectorProcess process) {
		this.getOptionParametersFromMarketProcess(process);

		Derivative putOption = new PutOption(this.strike, this.maturity, this.rfRate, this.divYield);
		putOption.setObsTimes(this.obsTimes);

		return putOption;
	}

	@Override
	public void run() {
		this.runForNbSteps(1);
		this.runForNbSteps(4);
	}

	private void convexityTests(DynamicHedgingProgram program, int nbSteps, int k, int i) {
		QFunction qFunction = program.getQFunction();

		String source = "Loss function";
		this.testConvexityInInterval(program.getLossFunction(), this.u0, this.u1, source);

		source = "Q function (k=" + k + "/" + nbSteps + ")";
		this.testConvexityInInterval(qFunction, this.u0, this.u1, source);
		// this.graphFunction(qFunction, this.u0, this.u1, source);

		source = "Exact risk function (k=" + k + "/" + nbSteps + ")";
		RiskFunctionAtFirstStep jFunction = program.getRiskFunctionAtFirstStep();
		this.testConvexityInInterval(jFunction, this.u0, this.u1, source);
		// this.graphFunction(jFunction, this.u0, this.u1, source);

		source = "Risk function approx. (k=" + k + "/" + nbSteps + ")";
		RiskFunctionApproximation jFunctionApprox = program.jStarArray[k][i];
		this.testConvexityInInterval(jFunctionApprox, this.u0, this.u1, source);
	}

	private DynamicHedgingProgram defineDynamicHedgingProgramUsingLQApproximation(Derivative derivative,
			MarketVectorProcess process) {
		ProblemDefinition problem = this.defineProblem(derivative, process);

		DynamicHedgingProgram program = new DynamicHedgingProgramUsingLQApproximation();
		program.setProblemDefinition(problem);
		program.setObservationTimes(process.getObservationTimes());
		program.setStochasticMesh(problem.stochasticMesh);

		SimulationParameters simParams = new SimulationParameters();
		simParams.derivative = derivative;
		program.setSimulationParameters(simParams);

		return program;
	}

	private ProblemDefinition defineProblem(Derivative derivative, MarketVectorProcess process) {
		StochasticMesh stochasticMesh = this.defineStochasticMesh(process);
		LossFunction lossFunction = new ExponentialLossFunction(-this.riskAversion);
		TransactionCosts transactionCosts = new ProportionalCosts(this.costParamA, this.costParamB);
		HedgingPortfolio hedgingPortfolio = new HedgingPortfolio(derivative, transactionCosts, this.rfRate);

		ProblemDefinition problem = new ProblemDefinition();
		problem.stochasticMesh = stochasticMesh;
		problem.lossFunction = lossFunction;
		problem.tcStructure = transactionCosts;
		problem.hedgedPtf = hedgingPortfolio;
		problem.derivative = derivative;

		return problem;
	}

	private StochasticMesh defineStochasticMesh(MarketVectorProcess process) {
		StochasticMesh stochasticMesh = new StochasticMeshAverageDensity(process);
		MeshWeights weights = stochasticMesh.getMeshWeights();
		weights.useRussianRoulette(this.weightCutoff, this.uWeightRoulette);
		double[] obsTimes = process.getObservationTimes();
		int nbObs = obsTimes.length - 1;
		stochasticMesh.setObservationTimes(obsTimes, nbObs);
		stochasticMesh.setNbPaths(this.nbMeshPaths);
		stochasticMesh.generateMesh();
		return stochasticMesh;
	}

	private void getOptionParametersFromMarketProcess(MarketVectorProcess process) {
		MarketVector y = process.getInitialMarketVector();
		this.strike = y.values[MarketVector.PRICE];
		this.obsTimes = process.getObservationTimes();
		this.maturity = this.obsTimes[this.obsTimes.length - 1];
	}

	private void graphFunction(OneDimFunction f, double u0, double u1, String source) {
		System.out.println(source);
		System.out.println("u;f");

		double du = (u1 - u0) / (this.nbDataPointsForTesting - 1);
		for (int i = 0; i < this.nbDataPointsForTesting; i++) {
			double u = u0 + i * du;

			double value = f.evaluate(u);
			System.out.println(u + ";" + value);
		}
	}

	private void graphFunctions(DynamicHedgingProgram program, int nbSteps, int k, int i) {

		// DynamicHedgingProgramDiagnostics diagnostics = new
		// DynamicHedgingProgramDiagnostics(program,
		// this.nbDataPointsForTesting);
		// QFunction qFunction = program.getQFunction();
		// QFunction qFunctionForLowerBoundary =
		// program.getQFunctionForLowerBoundary();
		// QFunction qFunctionForUpperBoundary =
		// program.getQFunctionForUpperBoundary();

	}

	private void runForNbSteps(int nbSteps) {
		DynamicHedgingProgram program = this.getInitializedDynamicProgram(nbSteps);
		HedgingPortfolio hedgingPortfolio = program.getHedgingPortfolio();
		hedgingPortfolio.setMaxAllowedOneStepLoss(this.maxLossForTest);

		program.computeOverAllSteps();

		int k = 0;
		int i = 0;

		this.graphFunctions(program, nbSteps, k, i);
		this.convexityTests(program, nbSteps, k, i);

		RiskFunctionAtFirstStep jFunction = program.getRiskFunctionAtFirstStep();
		RiskFunctionApproximation jFunctionApprox = program.jStarArray[k][i];

		this.testBreakPoints(nbSteps, k, jFunction, jFunctionApprox);
		this.testFunctionEqualities(program, nbSteps, k, jFunction, jFunctionApprox);
	}

	private void testBreakPoints(int nbSteps, int k, RiskFunctionAtFirstStep jFunction,
			RiskFunctionApproximation jFunctionApprox) {
		String source;
		source = "J value at lower hedging boundary (k=" + k + "/" + nbSteps + ")";
		HedgingBoundary lowerBoundary = jFunctionApprox.getLowerBoundary();
		this.testHedgingBoundary(lowerBoundary, jFunction, source);

		source = "J value at upper hedging boundary (k=" + k + "/" + nbSteps + ")";
		HedgingBoundary upperBoundary = jFunctionApprox.getUpperBoundary();
		this.testHedgingBoundary(upperBoundary, jFunction, source);
	}

	private void testConvexityInInterval(OneDimFunction function, double x0, double x1, String source) {
		boolean isConvex = this.assertConvexity(function, x0, x1, this.nbDataPointsForTesting, source);
		if (!isConvex) {
			this.graphFunction(function, this.u0, this.u1, source);
		}
	}

	private void testFunctionEqualities(DynamicHedgingProgram program, int nbSteps, int k,
			RiskFunctionAtFirstStep jFunction, RiskFunctionApproximation jFunctionApprox) {
		QFunction qFunction = program.getQFunction();
		QFunction qFunctionForLowerBoundary = program.getQFunctionForLowerBoundary();
		QFunction qFunctionForUpperBoundary = program.getQFunctionForUpperBoundary();

		double uk = 0;
		int i = 0;

		State x = new State();
		x.cashAccount = 0;
		x.stockQuantity = uk;
		x.marketVector = program.getNewMarketVector(k, i);
		qFunction.setState(x);
		qFunctionForLowerBoundary.setState(x);
		qFunctionForUpperBoundary.setState(x);

		String source;
		source = "Q-function for trade boundary (-), K=" + nbSteps;
		this.testFunctionEqualityInInterval(qFunctionForLowerBoundary, qFunction, uk, this.u1, source);
		source = "Q-function for trade boundary (+), K=" + nbSteps;
		this.testFunctionEqualityInInterval(qFunctionForUpperBoundary, qFunction, this.u0, uk, source);
		source = "Sum of errors of risk function approx - exact (k=" + k + "/" + nbSteps + ")";
		this.testFunctionEqualityInInterval(jFunctionApprox, jFunction, this.u0, this.u1, source);
	}

	private boolean testFunctionEqualityInInterval(OneDimFunction f1, OneDimFunction f2, double u0, double u1,
			String source) {
		double expectedError = 0;
		double sumOfSquaredErrors = 0;

		double du = (u1 - u0) / (this.nbDataPointsForTesting - 1);
		for (int i = 0; i < this.nbDataPointsForTesting; i++) {
			double u = u0 + i * du;

			double v2 = f2.evaluate(u);
			double v1 = f1.evaluate(u);

			double error = v1 - v2;

			sumOfSquaredErrors += error * error;
		}

		return this.assertEqualAbsolute(sumOfSquaredErrors, expectedError, source);
	}

	private boolean testHedgingBoundary(HedgingBoundary boundary, RiskFunctionAtFirstStep jFunction, String source) {
		double x = boundary.boundaryValue;
		double y0 = boundary.riskFunctionValueAtBoundary;

		double y1 = jFunction.evaluate(x);

		source += " (x=" + x + ")";

		return this.assertEqualAbsolute(y1, y0, source);
	}

}
