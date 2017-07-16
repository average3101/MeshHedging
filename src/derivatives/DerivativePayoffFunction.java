package derivatives;

import numerics.OneDimFunction;

public class DerivativePayoffFunction implements OneDimFunction {
	private Derivative derivative;

	public DerivativePayoffFunction(Derivative derivative) {
		this.derivative = derivative;
	}

	@Override
	public double evaluate(double v) {
		double value = this.derivative.payoff((float) v);
		return value;
	}

	@Override
	public String getName() {
		String label = this.derivative.getDescription() + " payoff";
		return label;
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
