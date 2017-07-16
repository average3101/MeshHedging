package numerics;

public class ConstantFunction implements OneDimFunction {

	private double cValue;

	public ConstantFunction(double c) {
		this.cValue = c;
	}

	@Override
	public double evaluate(double v) {
		return this.cValue;
	}

	@Override
	public String getName() {
		return "ConstantFunction";
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
