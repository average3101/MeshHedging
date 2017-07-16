package policies;

import dynprog.DynamicHedgingProgram;
import dynprog.RiskFunctionApproximation;
import dynprog.State;
import numerics.OneDimSolverResults;

public class RiskFunctionBasedPolicy extends HedgingPolicy {

	private DynamicHedgingProgram program;

	public RiskFunctionBasedPolicy(DynamicHedgingProgram program) {
		this.program = program;
		this.setPolicyHedgingMethod(program.getPolicyHedgingMethod());
	}

	@Override
	public void computeCummulApproxError() {
		this.program.computeCummulApproxErrorForAllSteps();
	}

	@Override
	public void computePolicy() {
		this.program.computeOverAllSteps();
	}

	@Override
	public double evaluate(State x) {
		this.program.setWeightsForConditionalExpectationsBasedOnOuterPaths(x.timeStep);

		OneDimSolverResults results = this.program.evaluate(x);
		double decision = results.optim_v;

		this.getPastDecisions()[x.timeStep] = decision;
		return decision;
	}

	public DynamicHedgingProgram getDynamicHedgingProgram() {
		return this.program;
	}

	@Override
	public double getErrorEstimateAtInitialeState() {
		double cummulError = this.program.getCummulativeApproximationError();
		return cummulError;
	}

	public double getFractionAttainableStates() {
		return this.program.getFractionAttainableStates();
	}

	public double getFractionOfWeightsUsed() {
		return this.program.getFractionOfWeightsUsed();
	}

	public double getRiskEstimate(State x) {
		int k = x.timeStep;
		int i = x.marketVector.simulationIndex;

		RiskFunctionApproximation approximation = this.program.jStarArray[k][i];
		return approximation.evaluate(x.cashAccount, x.stockQuantity);
	}

	@Override
	public double getRiskEstimateAtInitialState(double c0, double u0) {
		double sum = 0;
		int nbPaths = this.program.getNbStatesPerStep();
		for (int i = 0; i < nbPaths; i++) {
			RiskFunctionApproximation approximation = this.program.jStarArray[0][0];
			double result = approximation.evaluate(c0, u0);
			sum += result;
		}

		double average = sum / nbPaths;
		return average;
	}

	@Override
	public void init() {
		super.init();
		this.program.setObservationTimes(this.getObservationTimes());
	}

	public void initializeAfterMeshComputation() {
		this.program.initializeQFunctions();
		this.program.initializeStructuresDependentOnMeshConstruction();
	}

	public void setProgram(DynamicHedgingProgram program) {
		this.program = program;
	}

}
