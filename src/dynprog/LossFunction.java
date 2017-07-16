package dynprog;

import numerics.OneDimFunction;

public abstract class LossFunction implements OneDimFunction {

	public abstract double applyStandardScaling(double y);

	@Override
	public abstract double evaluate(double x);

	public abstract double evaluateLog(double x);

	@Override
	public abstract String getName();

	public abstract double getRiskAversionParameter();

}
