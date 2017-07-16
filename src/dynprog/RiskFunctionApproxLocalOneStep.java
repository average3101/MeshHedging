package dynprog;

import stochprocess.MarketVector;

public class RiskFunctionApproxLocalOneStep extends RiskFunctionApproxLinQuad {

	public RiskFunctionApproxLocalOneStep(DynamicHedgingProgram program, State x) {
		super(program, x, null, null);

		this.name = "RiskFunctionApproxLocalOneStep";

	}

	@Override
	public double evaluate(double c, double u) {
		return Double.NaN; // DO NOT USE
	}

	@Override
	public double evaluateIncludingCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts) {
		return this.evaluateOnlyCurrentCost(currentState, v, currentDerivValue, transactionCosts);
	}

	@Override
	public double evaluateOnlyCurrentCost(State currentState, double v, double currentDerivativeValue,
			double transactionCosts) {
		int k = currentState.timeStep;
		MarketVector destinationMarketVector = this.approximationOriginState.marketVector;
		float s1 = this.approximationOriginState.marketVector.values[MarketVector.PRICE];
		int nextStateIndex = destinationMarketVector.simulationIndex;

		double dh = this.derivative.derivativeValueArray[k + 1][nextStateIndex] - currentDerivativeValue;
		double ptfVariation = this.hedgingPortfolio.calculatePortfolioVariation(currentState, s1, v, dh,
				transactionCosts);

		return this.lf.evaluate(ptfVariation);
	}

	@Override
	protected void computeLQParameters() {

	}

}
