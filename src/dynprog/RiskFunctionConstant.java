package dynprog;

public class RiskFunctionConstant extends RiskFunctionApproximation {

	private double value;

	public RiskFunctionConstant(double value) {
		this.value = value;
	}

	@Override
	public double gradient(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double evaluate(double c, double u) {
		return this.value;
	}

	@Override
	public double evaluateIncludingCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts) {
		return this.value;
	}

	@Override
	public double evaluateOnlyCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts) {
		return 1;
	}

	@Override
	public double hessian(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

}
