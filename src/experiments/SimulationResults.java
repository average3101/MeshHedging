package experiments;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

public class SimulationResults {

	Tally approxErrorTally;
	Tally[] decisions;
	double[] dpCompTime;
	String errorMessage;
	Tally fracAttainablePoints;
	Tally fracWeightsUsed;
	boolean generatesError;
	Tally hInitValue;
	Tally hOuterValue;
	double initCompTime;
	double[][] initDec;
	double initializationTime = 0.0;
	Tally lbTally;
	Tally lbWithErrorTally;
	double[] meshCompTime;
	double[] outerSimAverageRisk;
	double[] outerSimAverageRiskWithoutCV;
	double[][] pathDecisions;
	double[] policyEvalTime;
	Tally[] ptfPnLForVarTally;
	Tally[] ptfPnLTally;
	Tally[] ptfPnLVarTally;
	Tally[] riskMeasureValues;
	Tally[] riskMeasureValuesWithoutCV;
	Chrono timer;
	double[][] outerSimRiskEstimate;
	double[][] outerSimCVEstimate;

	public SimulationResults() {
	}

	public void initialize(SimulationParameters simulationParameters, int nbPolicies) {
		this.defineStatisticalAccumulators(nbPolicies, simulationParameters);
		this.defineTimer();
		this.initializeStatisticalAccumulators(nbPolicies);
	}

	public void initializeOuterSimulationRiskEstimates() {
		int nbPolicies = this.outerSimAverageRisk.length;
		for (int k = 0; k < nbPolicies; k++) {
			this.outerSimAverageRisk[k] = 0.0;
			this.outerSimAverageRiskWithoutCV[k] = 0.0;
		}
	}

	private void defineStatisticalAccumulators(int nbPolicies, SimulationParameters simParams) {
		this.riskMeasureValues = new Tally[nbPolicies];
		this.riskMeasureValuesWithoutCV = new Tally[nbPolicies];
		this.decisions = new Tally[nbPolicies];
		this.ptfPnLTally = new Tally[nbPolicies];
		this.ptfPnLVarTally = new Tally[nbPolicies];
		this.ptfPnLForVarTally = new Tally[nbPolicies];

		this.lbTally = new Tally();
		this.approxErrorTally = new Tally();
		this.lbWithErrorTally = new Tally();
		this.policyEvalTime = new double[nbPolicies];
		this.meshCompTime = new double[nbPolicies];
		this.dpCompTime = new double[nbPolicies];

		for (int k = 0; k < nbPolicies; k++) {
			this.riskMeasureValues[k] = new Tally();
			this.riskMeasureValuesWithoutCV[k] = new Tally();
			this.decisions[k] = new Tally();
			this.ptfPnLTally[k] = new Tally();
			this.ptfPnLVarTally[k] = new Tally();
			this.ptfPnLForVarTally[k] = new Tally();
		}
		this.hInitValue = new Tally();
		this.hOuterValue = new Tally();
		this.fracWeightsUsed = new Tally();
		this.fracAttainablePoints = new Tally();

		this.outerSimAverageRisk = new double[nbPolicies];
		this.outerSimAverageRiskWithoutCV = new double[nbPolicies];

		this.outerSimRiskEstimate = new double[nbPolicies][simParams.nbOuterSims];
		this.outerSimCVEstimate = new double[nbPolicies][simParams.nbOuterSims];
	}

	private void defineTimer() {
		this.timer = new Chrono();
		this.timer.init();
	}

	private void initializeStatisticalAccumulators(int nbPolicies) {
		for (int k = 0; k < nbPolicies; k++) {
			this.policyEvalTime[k] = 0.0;
			this.meshCompTime[k] = 0.0;
			this.dpCompTime[k] = 0.0;
		}

		for (int k = 0; k < nbPolicies; k++) {
			this.riskMeasureValues[k].init();
			this.riskMeasureValuesWithoutCV[k].init();
			this.decisions[k].init();
		}

		this.lbTally.init();
		this.approxErrorTally.init();
		this.lbWithErrorTally.init();

		this.hInitValue.init();
		this.hOuterValue.init();
		this.fracAttainablePoints.init();
		this.fracWeightsUsed.init();

		for (int k = 0; k < nbPolicies; k++) {
			this.ptfPnLForVarTally[k].init();
		}
	}
}
