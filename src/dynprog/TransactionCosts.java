package dynprog;

public interface TransactionCosts {
	public double evaluate(double v, double s);

	public double evaluatePerUnitStock(double s);

	public double getCostsParameterA();

	public double getCostsParameterB();
}