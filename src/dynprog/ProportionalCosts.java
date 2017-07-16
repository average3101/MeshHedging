package dynprog;

import org.apache.commons.math3.util.FastMath;

public class ProportionalCosts implements TransactionCosts {
	double a, b;

	public ProportionalCosts(double a, double b) {
		this.a = a;
		this.b = b;
	}

	@Override
	public double evaluate(double v, double s) {
		return FastMath.abs(v) * (this.a + this.b * s);
	}

	@Override
	public double evaluatePerUnitStock(double s) {
		return (this.a + this.b * s);
	}

	@Override
	public double getCostsParameterA() {
		return this.a;
	}

	@Override
	public double getCostsParameterB() {
		return this.b;
	}
}