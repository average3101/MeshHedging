package dynprog;

import meshmethods.MeshPointEstimator;
import numerics.BrentSolver;
import numerics.OneDimFunction;
import numerics.OneDimSolver;
import numerics.OneDimSolverResults;
import stochprocess.MarketVector;

public class QFunction extends MeshPointEstimator implements OneDimFunction {
	public final static int FUNCTION = 0;
	public final static int GRADIENT = 1;
	public final static int HESSIAN = 2;

	private double currentDerivativeValue;
	private State currentState = new State();
	private double interestRatePerPeriod;
	private double maxDelta;
	private double minDelta;
	private State nextState = new State();
	protected OneDimSolver oneDimSolver;
	private final DynamicHedgingProgram program;
	private double signOfTradeBoundary = Double.NaN;
	private TransactionCosts transactionCosts;

	// returns v[0] = fct value, v[1] = gradient, v[2] = hessian
	protected double[] valueArray = new double[3];

	public QFunction(DynamicHedgingProgram program) {
		this.program = program;

		this.transactionCosts = program.getTransactionCosts();
	}

	public OneDimSolverResults computeCurrentMinimum(double v0) {
		// NOTE: Must have |maxDelta - minDelta| >> tol ...
		// Otherwise we run into problems

		double uL = this.minDelta;
		double uH = this.maxDelta;
		double uLast = v0;

		OneDimSolverResults results = null;
		results = this.oneDimSolver.minimize(uL, uLast, uH);

		return results;
	}

	public double computeFee(double du, double s0) {
		return -this.transactionCosts.evaluate(du, s0) * (1 + this.interestRatePerPeriod);
	}

	public double computeFeeDerivative(double du, double s0) {
		return -this.transactionCosts.evaluatePerUnitStock(s0) * (1 + this.interestRatePerPeriod);
	}

	@Override
	public double evaluate(double v) {
		// NOTES: 1) Weights need to be computed first !
		// 2) Depends on the current value of 'u' ! (Must be updated...)

		int nextStep = this.currentState.timeStep + 1;
		double du = v - this.currentState.stockQuantity;
		double s0 = this.currentState.marketVector.values[MarketVector.PRICE];
		double fee = this.computeFee(du, s0);

		double weightedSum = 0.0;
		RiskFunctionApproximation[][] jStarArray = this.getProgram().jStarArray;
		// double lastGradient = 0;
		// double lastHessian = 0;
		// double gamma =
		// this.getProgram().getLossFunction().getRiskAversionParameter();
		for (int jj = 0; jj < this.currentNbWeights; jj++) {
			int j = this.currentWeightIndices[jj];
			double w = this.currentWeights[jj];

			RiskFunctionApproximation approximation = jStarArray[nextStep][j];

			double costFctAtNextStep = approximation.evaluateIncludingCurrentCost(this.currentState, v,
					this.currentDerivativeValue, fee);
			weightedSum += w * costFctAtNextStep;

			// double s1 =
			// approximation.approximationMarketVector.values[MarketVector.PRICE];
			// double d1 = gamma * (this.signOfTradeBoundary *
			// this.transactionCosts.evaluatePerUnitStock(s0) + (s1 - s0));
			// double d2 = approximation.getFactorForDerivative();
			// double d = d1 + d2;
			// lastGradient += w * d * costFctAtNextStep;
			// lastHessian += w * d * d * costFctAtNextStep;
		}
		//
		// this.valueArray[QFunction.FUNCTION] = weightedSum;
		// this.valueArray[QFunction.GRADIENT] = lastGradient;
		// this.valueArray[QFunction.HESSIAN] = lastHessian;

		// System.out.println("*;" + v + ";" + weightedSum + ";" +
		// lastGradient);

		return weightedSum;
	}

	public double evaluateOnlyCurrentCost(double v) {
		// NOTES: 1) Weights need to be computed first !
		// 2) Depends on the current value of 'u' ! (Must be updated...)

		int nextStep = this.currentState.timeStep + 1;
		double du = v - this.currentState.stockQuantity;
		double s0 = this.currentState.marketVector.values[MarketVector.PRICE];
		double fee = this.computeFee(du, s0);

		double weightedSum = 0.0;
		for (int jj = 0; jj < this.currentNbWeights; jj++) {
			int j = this.currentWeightIndices[jj];
			double w = this.currentWeights[jj];

			RiskFunctionApproximation approximation = this.getProgram().jStarArray[nextStep][j];
			double costFctAtNextStep = approximation.evaluateOnlyCurrentCost(this.currentState, v,
					this.currentDerivativeValue, fee);
			weightedSum += w * costFctAtNextStep;
		}

		return weightedSum;
	}

	public double evaluateWithoutCurrentCost(double v) {
		double weightedSum = 0.0;
		int nextStep = this.currentState.timeStep + 1;
		double du = v - this.currentState.stockQuantity;
		int k = this.currentState.timeStep;

		this.nextState.timeStep = k + 1;
		this.nextState.stockQuantity = this.currentState.stockQuantity + du;

		for (int jj = 0; jj < this.currentNbWeights; jj++) {
			int j = this.currentWeightIndices[jj];
			double w = this.currentWeights[jj];

			RiskFunctionApproximation approximation = this.getProgram().jStarArray[nextStep][j];
			double costFctAtNextStep = approximation.evaluate(0, v);

			weightedSum += w * costFctAtNextStep;
		}

		return weightedSum;
	}

	public State getCurrentState() {
		return this.currentState;
	}

	public double getMaxDeltaForOptimization() {
		return this.maxDelta;
	}

	public double getMinDeltaForOptimization() {
		return this.minDelta;
	}

	@Override
	public String getName() {
		return "Q-Function";
	}

	public DynamicHedgingProgram getProgram() {
		return this.program;
	}

	public double getSignOfTradeBoundary() {
		return this.signOfTradeBoundary;
	}

	public double[] getValueArray() {
		return this.valueArray;
	}

	@Override
	public double gradient(double x) {
		return this.valueArray[QFunction.GRADIENT];
	}

	@Override
	public double hessian(double x) {
		return this.valueArray[QFunction.HESSIAN];
	}

	public void initializeSolver(double tol, double minDelta, double maxDelta) {
		this.minDelta = minDelta;
		this.maxDelta = maxDelta;

		this.defineOneDimSolver();
		this.oneDimSolver.setTolerance(tol);
	}

	public void setCurrentDerivativeValue(double h) {
		this.currentDerivativeValue = h;
	}

	public void setSignOfTradeBoundary(double signOfTradeBoundary) {
		this.signOfTradeBoundary = signOfTradeBoundary;
	}

	public void setSignOfTradeNegative() {
		this.setSignOfTradeBoundary(-1);
	}

	public void setSignOfTradePositive() {
		this.setSignOfTradeBoundary(+1);
	}

	public void setState(State x) {
		this.currentState.cashAccount = x.cashAccount;
		this.currentState.timeStep = x.timeStep;
		this.currentState.marketVector = x.marketVector;
		this.currentState.stockQuantity = x.stockQuantity;
	}

	protected void defineOneDimSolver() {
		this.oneDimSolver = new BrentSolver(this);
	}
}