package dynprog;

import numerics.OneDimSolverResults;
import stochprocess.MarketVector;

public class RiskFunctionAtFirstStep extends RiskFunctionApproximation {

	private State initialState;
	private DynamicHedgingProgram program;

	public RiskFunctionAtFirstStep(DynamicHedgingProgram program) {
		this.name = "RiskFunctionAtFirstStep";
		this.program = program;

		this.defineInitialState();
	}

	@Override
	public double gradient(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double evaluate(double v) {
		this.initialState.stockQuantity = v;
		OneDimSolverResults results = this.program.evaluate(this.initialState);
		return results.optim_f;
	}

	@Override
	public double evaluate(double c, double u) {
		this.initialState.stockQuantity = u;
		this.initialState.cashAccount = c;
		OneDimSolverResults results = this.program.evaluate(this.initialState);
		return results.optim_f;
	}

	@Override
	public double evaluateIncludingCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double evaluateOnlyCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double hessian(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	private void defineInitialState() {
		this.initialState = new State();
		int k = 0;
		int cashAccount = 0;
		int stateIndex = 0;
		MarketVector marketVector = this.program.getNewMarketVector(k, stateIndex);

		this.initialState.timeStep = k;
		this.initialState.cashAccount = cashAccount;
		this.initialState.marketVector = marketVector;
	}

}
