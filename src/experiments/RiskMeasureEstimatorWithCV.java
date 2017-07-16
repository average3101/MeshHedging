package experiments;

public class RiskMeasureEstimatorWithCV extends RiskMeasureEstimator {

	private int nbPilotRuns;
	private double basicRiskEstimate;
	private double avgCV;

	public RiskMeasureEstimatorWithCV(ProblemDefinition problem, TestedPolicies testPolicies, int nbPilotRuns) {
		super(problem, testPolicies);

		this.nbPilotRuns = nbPilotRuns;
	}

	private void calculateBasicEstimate(SimulationParameters simParams, SimulationResults results, int policyNo) {
		this.basicRiskEstimate = 0;
		this.avgCV = 0;
		for (int i = this.nbPilotRuns; i < simParams.nbOuterSims; i++) {
			this.basicRiskEstimate += results.outerSimRiskEstimate[policyNo][i];
			this.avgCV += results.outerSimCVEstimate[policyNo][i];
		}

		int nbPaths = simParams.nbOuterSims - this.nbPilotRuns;
		this.basicRiskEstimate /= nbPaths;
		this.avgCV /= nbPaths;
	}

	private double estimateBeta(SimulationResults results, int policyNo) {
		double avgCV = 0;
		double avgRisk = 0;
		for (int i = 0; i < this.nbPilotRuns; i++) {
			avgCV += results.outerSimCVEstimate[policyNo][i];
			avgRisk += results.outerSimRiskEstimate[policyNo][i];
		}
		avgCV /= this.nbPilotRuns;
		avgRisk /= this.nbPilotRuns;

		double covarEstim = 0;
		double varEstim = 0;

		for (int i = 0; i < this.nbPilotRuns; i++) {
			double x = (results.outerSimCVEstimate[policyNo][i] - avgCV);
			double y = (results.outerSimRiskEstimate[policyNo][i] - avgRisk);
			// System.out
			// .println(results.outerSimRiskEstimate[policyNo][i] + ";" +
			// results.outerSimCVEstimate[policyNo][i]);
			covarEstim += x * y;
			varEstim += y * y;
		}

		double betaEstim = covarEstim / varEstim;

		System.out.println(policyNo + ";" + betaEstim + ";" + avgCV + ";" + avgRisk);

		return betaEstim;
	}

	@Override
	protected double getAverageRiskUnscaled(SimulationParameters simParams, SimulationResults results, int policyNo) {
		double betaEstimate = this.estimateBeta(results, policyNo);

		this.calculateBasicEstimate(simParams, results, policyNo);
		double controlledEstimate = this.basicRiskEstimate - betaEstimate * this.avgCV;

		return controlledEstimate;
	}

}
