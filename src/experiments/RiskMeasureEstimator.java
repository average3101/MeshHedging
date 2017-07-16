package experiments;

import derivatives.DerivativePricedOnMesh;
import dynprog.LossFunction;
import dynprog.State;
import meshmethods.MeshWeights;
import meshmethods.StochasticMesh;
import policies.HedgingPolicy;
import policies.HedgingPortfolio;

public class RiskMeasureEstimator {
	private State initialState;
	private LossFunction lossFunctionForOuterSim;
	protected ProblemDefinition problem;
	private TestedPolicies testedPolicies;

	public RiskMeasureEstimator(ProblemDefinition problem, TestedPolicies testPolicies) {
		this.problem = problem;
		this.testedPolicies = testPolicies;
	}

	public SimulationResults estimatePolicyCosts(SimulationParameters simParams) {
		this.initialize(simParams);

		SimulationResults results = new SimulationResults();
		results.initialize(simParams, this.testedPolicies.getNbPolicies());

		for (int i = 0; i < simParams.nbReps; i++) {
			this.testedPolicies.computeMesh(simParams, this.problem, results);
			this.testedPolicies.computePolicies(simParams, this.problem, results);
			this.testedPolicies.computeInitialValues(simParams, this.problem, results);

			results.initializeOuterSimulationRiskEstimates();
			results = this.estimateAverageRiskForPoliciesOverAllPaths(simParams, results);

			this.storeRiskMeasures(simParams, results);
			this.storeInitialPolicyDecisions(results);
			this.storeMeshLowBiasedEstimator(simParams, results);
		}
		return results;
	}

	private void initialize(SimulationParameters simParams) {
		this.problem.initializeProcesses(simParams);
		this.problem.initializePortfolioHedgingProblem(simParams);
		this.problem.defineObservationTimes(simParams, this.testedPolicies);
		this.problem.defineStochasticMesh(simParams);
		this.problem.assignStochasticMeshToDerivative();
		this.testedPolicies.initializePolicySettings(simParams, this.problem);

		this.lossFunctionForOuterSim = this.problem.lossFunction;

		this.initialState = new State();
		this.initialState.cashAccount = this.problem.c0Start;
		this.initialState.stockQuantity = simParams.u0;
		this.initialState.timeStep = 0;
		this.initialState.marketVector = this.problem.mktInfoPath[0];
	}

	private void storeInitialPolicyDecisions(SimulationResults results) {
		int initialStepIndex = 0;
		for (int k = 0; k < this.testedPolicies.getNbPolicies(); k++) {
			results.decisions[k].add(this.testedPolicies.getPolicyArray()[k].getPastDecisions()[initialStepIndex]);
		}
	}

	private void storeMeshLowBiasedEstimator(SimulationParameters simParams, SimulationResults results) {
		double lbEstimate = Double.NaN;
		double errorEstimate = Double.NaN;

		if (this.testedPolicies.includesMeshBasedEstimator()) {
			int index = this.testedPolicies.getMeshBasedEstimatorIndex();
			HedgingPolicy[] policyArray = this.testedPolicies.getPolicyArray();
			double riskEstimate = policyArray[index].getRiskEstimateAtInitialState(this.problem.c0Start, simParams.u0);
			lbEstimate = this.problem.lossFunction.applyStandardScaling(riskEstimate);
			errorEstimate = policyArray[index].getErrorEstimateAtInitialeState();
		}
		results.lbTally.add(lbEstimate);
		results.approxErrorTally.add(errorEstimate);
		results.lbWithErrorTally.add(lbEstimate - errorEstimate);
	}

	private void storeRiskMeasures(SimulationParameters simParams, SimulationResults results) {
		for (int k = 0; k < this.testedPolicies.getNbPolicies(); k++) {
			double averageRiskUnscaled = this.getAverageRiskUnscaled(simParams, results, k);
			double averageRisk = this.lossFunctionForOuterSim.applyStandardScaling(averageRiskUnscaled);
			results.riskMeasureValues[k].add(averageRisk);
			results.ptfPnLVarTally[k].add(results.ptfPnLForVarTally[k].variance());

			double averageRiskWithoutCV = this.lossFunctionForOuterSim
					.applyStandardScaling(results.outerSimAverageRiskWithoutCV[k] / simParams.nbOuterSims);
			results.riskMeasureValuesWithoutCV[k].add(averageRiskWithoutCV);
		}
	}

	protected SimulationResults computeCurrentPathLossFunctionForTestPolicies(SimulationParameters simParams,
			SimulationResults results, int simIndex) {
		HedgingPolicy[] policyArray = this.testedPolicies.getPolicyArray();
		HedgingPortfolio portfolio = this.problem.hedgedPtf;
		for (int k = 0; k < this.testedPolicies.getNbPolicies(); k++) {
			results.initCompTime = results.timer.getSeconds();

			double ptfValueAtMaturity = portfolio.terminalPtfValueWithPolicy(this.initialState,
					this.problem.mktInfoPath, policyArray[k]);
			double ptfInitialValue = portfolio.getInitialValue();
			double pnl = (ptfValueAtMaturity - ptfInitialValue);

			double lossPnl = this.lossFunctionForOuterSim.evaluate(pnl);
			double transactionCosts = portfolio.getLastSimulatedPathTransactionCosts();

			results.outerSimRiskEstimate[k][simIndex] = lossPnl;

			// CV : pnl without costs has expectation 0
			results.outerSimCVEstimate[k][simIndex] = pnl - transactionCosts;

			results.outerSimAverageRisk[k] += lossPnl;
			results.outerSimAverageRiskWithoutCV[k] += lossPnl;
			results.ptfPnLTally[k].add(pnl);
			results.ptfPnLForVarTally[k].add(pnl);
			results.policyEvalTime[k] += results.timer.getSeconds() - results.initCompTime;
		}
		return results;
	}

	protected void computeWeightsAlongPath(SimulationParameters simParams) {
		StochasticMesh sm = this.problem.stochasticMesh;
		MeshWeights meshWeights = sm.getMeshWeights();
		meshWeights.computeAllWeightsForMarketVectorPath(this.problem.mktInfoPath);

		if (simParams.derivative instanceof DerivativePricedOnMesh) {
			MeshWeights meshWeightsForPricing = meshWeights;
			if (sm.isUsingDifferentProcessForPricing()) {
				meshWeightsForPricing = sm.getMeshWeightsForPricing();
				meshWeightsForPricing.computeAllWeightsForMarketVectorPath(this.problem.mktInfoPath);
			}

			DerivativePricedOnMesh d = ((DerivativePricedOnMesh) simParams.derivative);
			d.setWeightPathForDerivativePath(meshWeightsForPricing.getWeightsForOuterMarketInfoPath(),
					meshWeightsForPricing.getWeightDestinationIndicesForMarketInfoPath(),
					meshWeightsForPricing.getNbWeightsForMarketInfoPath());
		}
	}

	protected SimulationResults estimateAverageRiskForPoliciesOverAllPaths(SimulationParameters simParams,
			SimulationResults results) {
		double derivAvg = 0.0;
		for (int ii = 0; ii < simParams.nbOuterSims; ii++) {
			this.simulatePath();
			this.computeWeightsAlongPath(simParams);
			results = this.computeCurrentPathLossFunctionForTestPolicies(simParams, results, ii);

			double[] path = this.problem.stochprocessForExperiments.getPath();
			float sk = (float) path[simParams.nbObs];
			double derivPayoff = simParams.derivative.payoff(sk);
			derivAvg += derivPayoff;
		}

		double hOuterValueEstimate = derivAvg / simParams.nbOuterSims;
		results.hOuterValue.add(hOuterValueEstimate);
		return results;
	}

	protected double getAverageRiskUnscaled(SimulationParameters simParams, SimulationResults results, int policyNo) {
		return results.outerSimAverageRisk[policyNo] / simParams.nbOuterSims;
	}

	protected void simulatePath() {
		this.problem.stochprocessForExperiments.generateBasicProcessesPaths();
		this.problem.mktInfoPath = this.problem.stochprocessForExperiments.getMarketVectorPath();
		this.problem.stochprocessForExperiments.resetNextSubstream();
	}
}
