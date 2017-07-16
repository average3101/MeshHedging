package dynprog;

import derivatives.Derivative;
import numerics.OneDimFunction;
import policies.HedgingPortfolio;
import stochprocess.MarketVector;

public abstract class RiskFunctionApproximation implements OneDimFunction {

	protected State approximationOriginState;
	protected MarketVector approximationMarketVector;
	protected Derivative derivative;
	protected double[][] derivativeValueArray;
	protected float[] marketVectorValues;
	protected HedgingPortfolio hedgingPortfolio;
	protected LossFunction lf;
	protected HedgingBoundary lowerBoundary;
	protected String name = "undefined";
	protected DynamicHedgingProgram program;
	protected HedgingBoundary upperBoundary;

	@Override
	public double evaluate(double u) {
		double cashAccount = 0;

		return this.evaluate(cashAccount, u);
	}

	public abstract double evaluate(double c, double u);

	public abstract double evaluateIncludingCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts);

	public abstract double evaluateOnlyCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts);

	public double getFactorForDerivative() {
		// TODO Auto-generated method stub
		return 0;
	}

	public HedgingBoundary getLowerBoundary() {
		return this.lowerBoundary;
	}

	@Override
	public String getName() {
		return this.name;
	}

	public HedgingBoundary getUpperBoundary() {
		return this.upperBoundary;
	}

	public void printBoundaryValues() {
	}

	public void printValuesToScreen(double x1, double x2, int nbPts) {
		double c = 0;
		double dx = (x2 - x1) / nbPts;
		System.out.println("*** Fct values for " + this);
		for (int i = 0; i <= nbPts; i++) {
			double x = x1 + dx * i;
			double f = this.evaluate(c, x);
			System.out.println(x + "; " + f);

		}
	}

	protected void initialize() {

	}
}