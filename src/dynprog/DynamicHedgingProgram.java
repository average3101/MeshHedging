package dynprog;

import java.util.HashSet;
import java.util.Iterator;

import derivatives.Derivative;
import derivatives.DerivativePricedOnMesh;
import experiments.ProblemDefinition;
import experiments.SimulationParameters;
import meshmethods.MeshFullEstimator;
import meshmethods.MeshWeights;
import meshmethods.StochasticMesh;
import numerics.OneDimSolverResults;
import policies.BSMHedging;
import policies.HedgingPolicy;
import policies.HedgingPolicy.HedgingMethod;
import policies.HedgingPortfolio;
import policies.RiskFunctionBasedPolicy;
import stochprocess.MarketVector;

public abstract class DynamicHedgingProgram extends MeshFullEstimator {

	private double boundaryTol = 1E-8;
	private double[][] cummulApproxError;
	private Derivative derivative;
	private HedgingPortfolio hedgingPortfolio;
	private HedgingPolicy heuristicForInitialApprox;
	public RiskFunctionApproximation[][] jStarArray;
	private LossFunction lossFunction;
	private int nbPtsForErrorCalc;
	private int nbStatesPerStep;
	private double[] obsTimes;
	private HedgingMethod policyHedgingMethod;
	private QFunction qFunction;
	private QFunctionForTradeBoundaries qFunctionForLowerBoundary;
	private QFunctionForTradeBoundaries qFunctionForUpperBoundary;
	private double qFunctionTol = 1e-8;
	private State tempStateDestination;
	private State tempStateOrigin;
	private TransactionCosts transactionCosts;
	private boolean usingDerivativeEvaluatedByStochMesh = false;

	public DynamicHedgingProgram() {
		this.defineHeuristicForOptimizationStartPoint();
	}

	public void computeCummulApproxErrorForAllSteps() {
		int nbSteps = this.getNbSteps();
		double minDelta = this.qFunction.getMinDeltaForOptimization();
		double maxDelta = this.qFunction.getMaxDeltaForOptimization();
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		HashSet<Integer>[] availableStateIndices = stochasticMesh.getAttainableStateIndices();
		if (nbSteps > 1) {
			for (int k = nbSteps - 1; k > 0; k--) {
				Iterator<Integer> availableIndicesIterator = availableStateIndices[k].iterator();
				while (availableIndicesIterator.hasNext()) {
					int i = availableIndicesIterator.next().intValue();
					this.computeNewCummulativeApproximationError(this.jStarArray[k][i], minDelta, maxDelta);
				}
			}
		}
		int k = 0;
		this.computeNewCummulativeApproximationError(this.jStarArray[k][0], minDelta, maxDelta);
	}

	@Override
	public void computeForStateOverLastStep(int i) {
		int k = this.getNbSteps();
		State stateForApprox = this.getTempStateOrigin(k, i);
		float sK = stateForApprox.marketVector.values[MarketVector.PRICE];
		double hK = this.derivative.payoff(sK);
		this.jStarArray[k][i] = new RiskFunctionApproxLastStep(this, stateForApprox, hK, this.hedgingPortfolio);
	}

	public void computeNewCummulativeApproximationError(RiskFunctionApproximation riskFunction, double uMin,
			double uMax) {
		int k = riskFunction.approximationOriginState.timeStep;
		int i = riskFunction.approximationOriginState.marketVector.simulationIndex;

		double uBackupValue = riskFunction.approximationOriginState.stockQuantity;

		this.setMeshRelatedValuesWithPrecomputedDerivativeValue(k, i);

		QFunction qFunction = this.getQFunction();
		int currentNbWeights = qFunction.getCurrentNbWeights();
		int[] currentWeightIndicies = qFunction.getCurrentWeightIndices();
		float[] currentWeights = qFunction.getCurrentWeights();

		double avgCummulApproxErrorAtNextStep = 0;
		for (int jj = 0; jj < currentNbWeights; jj++) {
			int j = currentWeightIndicies[jj];
			double w = currentWeights[jj];
			double errorForState = this.cummulApproxError[j][k + 1];
			avgCummulApproxErrorAtNextStep += w * errorForState;
		}

		double newOneStepErrorTerm = this.computeNewOneStepErrorTerm(riskFunction, uMin, uMax);

		riskFunction.approximationOriginState.stockQuantity = uBackupValue;

		this.cummulApproxError[i][k] = avgCummulApproxErrorAtNextStep + newOneStepErrorTerm;
	}

	public double computeNewOneStepErrorTerm(RiskFunctionApproximation riskFunction, double uMin, double uMax) {
		double errorSup = Double.MIN_VALUE;
		double oneStepErrorTerm = 0;
		int nbOfPtsForErrorCalc = this.getNbPtsForErrorCalc();
		double du = (uMax - uMin) / (nbOfPtsForErrorCalc - 1.0);
		for (int j = 0; j < nbOfPtsForErrorCalc; j++) {
			double u = uMin + j * du;
			riskFunction.approximationOriginState.stockQuantity = u;
			OneDimSolverResults results = this.evaluate(riskFunction.approximationOriginState);
			double exactValue = results.optim_f;

			double newC = 0; // c + u * s; // To account for fact cost function
								// approx is based on c=0, u=0
			double approxValue = riskFunction.evaluate(newC, u);

			oneStepErrorTerm = Math.abs(exactValue - approxValue);

			if (errorSup < oneStepErrorTerm) {
				errorSup = oneStepErrorTerm;
			}
		}
		return oneStepErrorTerm;
	}

	public OneDimSolverResults evaluate(State x) {

		return this.evaluateForQFunction(x, this.qFunction);
	}

	public double evaluateAtLastStep(State x) {
		double pnl = this.hedgingPortfolio.evaluate(x) - this.hedgingPortfolio.getInitialValue();
		return this.lossFunction.evaluate(pnl);
	}

	public OneDimSolverResults evaluateForQFunction(State x, QFunction qf) {
		qf.setState(x);

		// DynamicHedgingProgramDiagnostics diag = new
		// DynamicHedgingProgramDiagnostics(this, 21);//
		// diag.printQFunctionValuesForTradeBoundariesToScreen(x.timeStep,
		// x.marketVector.simulationIndex, -1, 0);

		double bsmDeltaAsInitialApproximation = this.heuristicForInitialApprox.evaluate(x);
		OneDimSolverResults results = qf.computeCurrentMinimum(bsmDeltaAsInitialApproximation);

		return results;
	}

	public double getBoundaryTol() {
		return this.boundaryTol;
	}

	public RiskFunctionApproximation[][] getCostFunctionApproximationArray() {
		return this.jStarArray;
	}

	public double getCummulativeApproximationError() {
		return this.cummulApproxError[0][0];
	}

	public Derivative getDerivative() {
		return this.derivative;
	}

	public double getFractionAttainableStates() {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		return stochasticMesh.getFractionAttainableStates();
	}

	public double getFractionOfWeightsUsed() {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		MeshWeights weights = stochasticMesh.getMeshWeights();
		return weights.getFractionOfWeightsUsed();
	}

	public HedgingPortfolio getHedgingPortfolio() {
		return this.hedgingPortfolio;
	}

	public LossFunction getLossFunction() {
		return this.lossFunction;
	}

	public int getNbPtsForErrorCalc() {
		return this.nbPtsForErrorCalc;
	}

	public int getNbStatesPerStep() {
		return this.nbStatesPerStep;
	}

	public HedgingPolicy getNewHedgingPolicy() {
		return new RiskFunctionBasedPolicy(this);
	}

	public double[] getObservationTimes() {
		return this.obsTimes;
	}

	public HedgingPolicy.HedgingMethod getPolicyHedgingMethod() {
		return this.policyHedgingMethod;
	}

	public QFunction getQFunction() {
		return this.qFunction;
	}

	public QFunctionForTradeBoundaries getQFunctionForLowerBoundary() {
		return this.qFunctionForLowerBoundary;
	}

	public QFunctionForTradeBoundaries getQFunctionForUpperBoundary() {
		return this.qFunctionForUpperBoundary;
	}

	public RiskFunctionAtFirstStep getRiskFunctionAtFirstStep() {
		RiskFunctionAtFirstStep jFunction = new RiskFunctionAtFirstStep(this);
		return jFunction;
	}

	public State getTempStateDestination(int k, int i) {
		MarketVector tempMarketVectorDestination = this.getTempMarketVectorDestination(k, i);
		this.initializeTempStateValues(this.tempStateDestination, k, i, tempMarketVectorDestination);
		return this.tempStateDestination;
	}

	public State getTempStateOrigin(int k, int i) {
		MarketVector tempMarketVectorOrigin = this.getTempMarketVectorOrigin(k, i);
		this.initializeTempStateValues(this.tempStateOrigin, k, i, tempMarketVectorOrigin);
		return this.tempStateOrigin;
	}

	public TransactionCosts getTransactionCosts() {
		return this.transactionCosts;
	}

	public void init() {
		this.heuristicForInitialApprox.setObservationTimes(this.obsTimes);
		this.tempStateOrigin = this.createGenericState();
		this.tempStateDestination = this.createGenericState();
	}

	public void initializeDerivativePayoffOnMesh() {
		int k = this.getNbSteps();
		for (int i = 0; i < this.nbStatesPerStep; i++) {
			State stateForApprox = this.getTempStateOrigin(k, i);
			float sK = stateForApprox.marketVector.values[MarketVector.PRICE];
			double hK = this.derivative.payoff(sK);
			this.derivative.setDerivativeValueArrayValue(k, i, hK);
		}
	}

	public void initializeQFunctions() {
		this.qFunction = new QFunction(this);
		this.qFunctionForLowerBoundary = new QFunctionForTradeBoundaries(this);
		this.qFunctionForUpperBoundary = new QFunctionForTradeBoundaries(this);

		this.qFunctionForLowerBoundary.setSignOfTradeNegative();
		this.qFunctionForUpperBoundary.setSignOfTradePositive();

		this.qFunction.initializeSolver(this.qFunctionTol, this.derivative.getMinDeltaForHedge(),
				this.derivative.getMaxDeltaForHedge());
		this.qFunctionForLowerBoundary.initializeSolver(this.qFunctionTol, this.derivative.getMinDeltaForHedge(),
				this.derivative.getMaxDeltaForHedge());
		this.qFunctionForUpperBoundary.initializeSolver(this.qFunctionTol, this.derivative.getMinDeltaForHedge(),
				this.derivative.getMaxDeltaForHedge());
	}

	public void initializeRiskFunctionArray() {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		this.nbStatesPerStep = stochasticMesh.getNbPaths();
		int nbSteps = this.getNbSteps();
		this.jStarArray = new RiskFunctionApproximation[nbSteps + 1][this.nbStatesPerStep];
		this.cummulApproxError = new double[this.nbStatesPerStep][nbSteps + 1];
	}

	public void initializeStructuresDependentOnMeshConstruction() {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		this.nbStatesPerStep = stochasticMesh.getNbPaths();

		this.qFunction.setCurrentWeights(new float[this.nbStatesPerStep]);
		this.qFunction.setCurrentWeightsForPricing(new float[this.nbStatesPerStep]);
		this.qFunction.setStochasticMesh(stochasticMesh);
		this.qFunctionForLowerBoundary.setCurrentWeights(new float[this.nbStatesPerStep]);
		this.qFunctionForUpperBoundary.setCurrentWeights(new float[this.nbStatesPerStep]);

		this.qFunction.setCurrentWeightIndices(new int[this.nbStatesPerStep]);
		this.qFunction.setCurrentWeightIndicesForPricing(new int[this.nbStatesPerStep]);
		this.qFunctionForLowerBoundary.setCurrentWeightIndices(new int[this.nbStatesPerStep]);
		this.qFunctionForUpperBoundary.setCurrentWeightIndices(new int[this.nbStatesPerStep]);

		this.derivative.initializeAfterMeshComputation(this.getNbSteps(), this.nbStatesPerStep, stochasticMesh);

		this.initializeRiskFunctionArray();
		this.initializeDerivativePayoffOnMesh();
	}

	public void setCostFunctionApproximationForPredefinedState(int timeStep, int simulationIndex,
			RiskFunctionApproximation approximation) {
		this.jStarArray[timeStep][simulationIndex] = approximation;
	}

	public void setMeshRelatedValuesIncludingDerivativeValueComputation(int k, int i) {
		this.setMeshRelatedValuesWithPrecomputedDerivativeValue(k, i);
		MarketVector y = this.getTempMarketVectorOrigin(k, i);
		double h = this.derivative.evaluate(k, y.values);
		this.derivative.setDerivativeValueArrayValue(k, i, h);

		this.setCurrentDerivativeValue(h);
	}

	public void setMeshRelatedValuesWithPrecomputedDerivativeValue(int k, int i) {
		this.qFunction.setCurrentWeights(k, i);

		float[] currentWeights = this.qFunction.getCurrentWeights();
		int[] currentWeightIndicies = this.qFunction.getCurrentWeightIndices();
		int currentNbWeights = this.qFunction.getCurrentNbWeights();

		this.qFunctionForLowerBoundary.setCurrentRiskFunctionWeights(currentWeights, currentWeightIndicies,
				currentNbWeights);
		this.qFunctionForUpperBoundary.setCurrentRiskFunctionWeights(currentWeights, currentWeightIndicies,
				currentNbWeights);

		if (this.usingDerivativeEvaluatedByStochMesh) {
			float[] currentWeightsForPricing = currentWeights;
			int[] currentWeightIndiciesForPricing = currentWeightIndicies;
			int currentNbWeightsForPricing = currentNbWeights;

			if (this.getStochasticMesh().isUsingDifferentProcessForPricing()) {
				currentWeightsForPricing = this.qFunction.getCurrentWeightsForPricing();
				currentWeightIndiciesForPricing = this.qFunction.getCurrentWeightIndicesForPricing();
				currentNbWeightsForPricing = this.qFunction.getCurrentNbWeightsForPricing();
			}

			((DerivativePricedOnMesh) this.derivative).setCurrentWeights(currentWeightsForPricing,
					currentWeightIndiciesForPricing, currentNbWeightsForPricing);
		}

		double h = this.derivative.getDerivativeValueArray()[k][i];
		this.setCurrentDerivativeValue(h);
	}

	public void setObservationTimes(double[] t) {
		this.obsTimes = t;
		this.setNbSteps(t.length - 1);
		this.init();
	}

	public void setPolicyHedgingMethod(HedgingMethod hedgingMethod) {
		this.policyHedgingMethod = hedgingMethod;
	}

	public void setProblemDefinition(ProblemDefinition problem) {
		this.setStochasticMesh(problem.stochasticMesh);

		this.lossFunction = problem.lossFunction;
		this.hedgingPortfolio = problem.hedgedPtf;
		this.transactionCosts = problem.tcStructure;
		this.derivative = problem.derivative;
	}

	public void setSimulationParameters(SimulationParameters simParams) {
		this.nbPtsForErrorCalc = simParams.nbPtsForErrorCalc;

		this.derivative = simParams.derivative;
		this.heuristicForInitialApprox.setParameters(simParams);

		if (this.derivative instanceof DerivativePricedOnMesh) {
			this.usingDerivativeEvaluatedByStochMesh = true;
		}
	}

	public void setWeightsForConditionalExpectationsBasedOnOuterPaths(int k) {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		MeshWeights meshWeights = stochasticMesh.getMeshWeights();
		float[] weights = meshWeights.getWeightsForOuterMarketInfoPath()[k];
		int[] indices = meshWeights.getWeightDestinationIndicesForMarketInfoPath()[k];
		int nbIndices = meshWeights.getNbWeightsForMarketInfoPath()[k];

		this.qFunction.setCurrentRiskFunctionWeights(weights, indices, nbIndices);
		if (this.usingDerivativeEvaluatedByStochMesh) {

			if (this.usingDerivativeEvaluatedByStochMesh) {
				float[] weightsForPricing = weights;
				int[] indicesPricing = indices;
				int nbIndicesPricing = nbIndices;

				if (this.getStochasticMesh().isUsingDifferentProcessForPricing()) {
					MeshWeights meshWeightsForPricing = stochasticMesh.getMeshWeightsForPricing();

					weightsForPricing = meshWeightsForPricing.getWeightsForOuterMarketInfoPath()[k];
					indicesPricing = meshWeightsForPricing.getWeightDestinationIndicesForMarketInfoPath()[k];
					nbIndicesPricing = meshWeightsForPricing.getNbWeightsForMarketInfoPath()[k];
				}

				((DerivativePricedOnMesh) this.derivative).setCurrentWeights(weightsForPricing, indicesPricing,
						nbIndicesPricing);
			}
		}

		MarketVector[] outerMarketVectorPath = meshWeights.getOuterMarketVectorPath();
		MarketVector y = outerMarketVectorPath[k];
		double h = this.derivative.evaluate(k, y.values);

		this.qFunction.setCurrentDerivativeValue(h);
	}

	private State createGenericState() {
		State stateForApprox = new State();
		stateForApprox.cashAccount = 0.0;
		stateForApprox.stockQuantity = 0.0;
		return stateForApprox;
	}

	private void defineHeuristicForOptimizationStartPoint() {
		this.heuristicForInitialApprox = new BSMHedging();

	}

	private void initializeTempStateValues(State tempState, int k, int i, MarketVector marketVector) {
		tempState.timeStep = k;
		tempState.cashAccount = 0;
		tempState.stockQuantity = 0;
		tempState.marketVector = marketVector;
	}

	private void setCurrentDerivativeValue(double h) {
		this.qFunction.setCurrentDerivativeValue(h);
		this.qFunctionForLowerBoundary.setCurrentDerivativeValue(h);
		this.qFunctionForUpperBoundary.setCurrentDerivativeValue(h);
	}

}