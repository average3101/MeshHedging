package numerics;

import umontreal.ssj.functions.MathFunction;

public interface OneDimFunction extends MathFunction {

	@Override
	public double evaluate(double v);

	public String getName();

	public abstract double gradient(double x);

	public abstract double hessian(double x);

}
