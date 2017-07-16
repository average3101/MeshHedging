package dynprog;

import policies.HedgingPortfolio;
import stochprocess.MarketVector;

public class RiskFunctionApproxLastStep extends RiskFunctionApproxLocalOneStep {

	double sK, hK, initPtfValue;

	public RiskFunctionApproxLastStep(DynamicHedgingProgram program, State x, double hK,
			HedgingPortfolio hedgingPortfolio) {
		super(program, x);
		this.name = "RiskFunctionApproxLastStep";

		this.sK = x.marketVector.values[MarketVector.PRICE];
		this.hK = hK;
		this.hedgingPortfolio = hedgingPortfolio;

	}

	@Override
	public double evaluate(double c, double u) {
		return 1.0;
	}

	@Override
	public double evaluateIncludingCurrentCost(State currentState, double v, double currentDerivativeValue,
			double transactionCosts) {
		double currentCost = this.evaluateOnlyCurrentCost(currentState, v, currentDerivativeValue, transactionCosts);
		return currentCost;
	}

	@Override
	public double evaluateOnlyCurrentCost(State currentState, double v, double currentDerivativeValue,
			double transactionCosts) {
		int k = currentState.timeStep;
		MarketVector destinationMarketVector = this.approximationOriginState.marketVector;
		float s1 = this.approximationOriginState.marketVector.values[MarketVector.PRICE];
		int nextStateIndex = destinationMarketVector.simulationIndex;

		double dh = this.derivative.derivativeValueArray[k + 1][nextStateIndex] - currentDerivativeValue;
		// double dh = this.derivative.evaluateIncrement(currentDerivValue, k +
		// 1, nextStateIndex);
		double ptfVariation = this.hedgingPortfolio.calculatePortfolioVariation(currentState, s1, v, dh,
				transactionCosts);

		double currentCost = this.lf.evaluate(ptfVariation);
		return currentCost;
	}
}
