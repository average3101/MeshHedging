package numerics;

public class LinearFunction implements OneDimFunction {

	private double intercept;
	private double slope;

	public LinearFunction(double slope, double intercept) {
		this.slope = slope;
		this.intercept = intercept;
	}

	@Override
	public double gradient(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double evaluate(double x) {
		double v = this.intercept + this.slope * x;
		return v;
	}

	@Override
	public String getName() {
		return "LinearFunction";
	}

	@Override
	public double hessian(double x) {
		// TODO Auto-generated method stub
		return 0;
	}

}
