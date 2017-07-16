package dynprog;

import org.apache.commons.math3.util.FastMath;

import numerics.OneDimSolverResults;
import stochprocess.MarketVector;

public class HedgingBoundary {

	public double boundaryValue;
	public double riskFunctionValueAtBoundary;

	private double gamma;
	private DynamicHedgingProgram program;
	private QFunctionForTradeBoundaries qFunctionForBoundary;
	private double signOfTrade;
	private TransactionCosts transactionCosts;

	public HedgingBoundary(DynamicHedgingProgram program, double riskAversionParameter,
			QFunctionForTradeBoundaries qFunctionForBoundary) {
		this.program = program;
		this.gamma = riskAversionParameter;
		this.qFunctionForBoundary = qFunctionForBoundary;
		this.signOfTrade = qFunctionForBoundary.getSignOfTradeBoundary();
		this.transactionCosts = program.getTransactionCosts();
	}

	public void estimate(State x) {
		OneDimSolverResults results = this.program.evaluateForQFunction(x, this.qFunctionForBoundary);

		this.boundaryValue = results.optim_v;

		double s0 = x.marketVector.values[MarketVector.PRICE];
		double costs = this.transactionCosts.evaluatePerUnitStock(s0);
		double z = -this.signOfTrade * this.gamma * costs * this.boundaryValue;
		this.riskFunctionValueAtBoundary = results.optim_f * FastMath.exp(z);
	}

	public double evaluate(double u) {
		return this.boundaryValue;
	}
}
