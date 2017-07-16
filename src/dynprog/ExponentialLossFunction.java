package dynprog;

import org.apache.commons.math3.util.FastMath;

public class ExponentialLossFunction extends LossFunction {

	double b;

	public ExponentialLossFunction(double b) {
		this.b = b;
	}

	@Override
	public double applyStandardScaling(double y) {
		return (1.0 - y) / this.b;
	}

	@Override
	public double evaluate(double v) {
		return FastMath.exp(this.b * v);
	}

	@Override
	public double evaluateLog(double v) {
		return this.b * v;
	}

	@Override
	public String getName() {
		return "exp";
	}

	public double getParameterB() {
		return this.b;
	}

	@Override
	public double getRiskAversionParameter() {
		return this.b;
	}

	@Override
	public double gradient(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double hessian(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

}
