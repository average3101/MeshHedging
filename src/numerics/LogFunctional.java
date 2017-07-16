package numerics;

import org.apache.commons.math3.util.FastMath;

public class LogFunctional implements OneDimFunction {

	public static double inverse(double y) {
		return FastMath.exp(y);
	}

	private OneDimFunction basicFunction;

	public LogFunctional(OneDimFunction basicFunction) {
		this.basicFunction = basicFunction;
	}

	@Override
	public double evaluate(double x) {
		return Math.log(this.basicFunction.evaluate(x));
	}

	@Override
	public String getName() {
		return "Log(" + this.basicFunction.getName() + ")";
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
