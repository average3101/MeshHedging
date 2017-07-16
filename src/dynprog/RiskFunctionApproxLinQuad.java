package dynprog;

import org.apache.commons.math3.util.FastMath;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import stochprocess.MarketVector;

public class RiskFunctionApproxLinQuad extends RiskFunctionApproximation {

	double A, B, C;
	Algebra alg = new Algebra();
	double log_y1, log_y2; // params for linear parts in the wings
	double m;
	double x1, x2; // upper and lower bounds for no trade region
	double x3, log_y3; // mid point used to compute quadratic part
	double riskAversion;
	double factorForDerivative = 0;

	protected double aTC;
	protected double bTC;
	protected double gamma;
	protected double s0;

	public RiskFunctionApproxLinQuad(DynamicHedgingProgram program, State x, HedgingBoundary lowerBoundary,
			HedgingBoundary upperBoundary) {
		this.name = "RiskFunctionApproxLinQuad";

		this.approximationOriginState = x.clone();
		this.approximationMarketVector = this.approximationOriginState.marketVector;
		this.program = program;
		this.hedgingPortfolio = this.program.getHedgingPortfolio();
		this.lf = program.getLossFunction();
		this.derivative = this.program.getDerivative();

		this.riskAversion = this.lf.getRiskAversionParameter();

		this.lowerBoundary = lowerBoundary;
		this.upperBoundary = upperBoundary;

		this.initialize();
		this.computeLQParameters();
	}

	public DoubleMatrix1D computeQuadraticFunctionParameters(double x1, double x2, double x3, double y1, double y2,
			double y3) {
		double[][] data = { { x1 * x1, x1, 1 }, { x2 * x2, x2, 1 }, { x3 * x3, x3, 1 } };
		double[] yVec = { y1, y2, y3 };
		DoubleMatrix1D paramVector = null;
		DoubleMatrix2D xDataMatrix = DoubleFactory2D.sparse.make(data);
		DoubleMatrix1D yVector = DoubleFactory1D.sparse.make(yVec);
		paramVector = this.alg.mult(this.alg.inverse(xDataMatrix), yVector);

		return paramVector;
	}

	@Override
	public double evaluate(double cash, double u) {
		double z = this.evaluateLog(cash, u);
		return FastMath.exp(z);
	}

	@Override
	public double evaluateIncludingCurrentCost(State currentState, double v, double currentDerivValue,
			double transactionCosts) {
		double currentCostLog = this.evaluateOnlyCurrentCostLog(currentState, v, currentDerivValue, transactionCosts);
		double nextStepCostFunctionValueLog = currentCostLog + this.evaluateLog(0, v);

		return FastMath.exp(nextStepCostFunctionValueLog);
	}

	public double evaluateLog(double cash, double u) {
		double d1 = u - this.x1;
		double d2 = u - this.x2;
		if (d1 < 0) {
			this.factorForDerivative = -this.m;
			return -this.m * d1 + this.log_y1;
		} else if (d2 > 0) {
			this.factorForDerivative = this.m;
			return this.m * d2 + this.log_y2;
		} else {
			this.factorForDerivative = 2 * this.A * u + this.B;
			return this.A * u * u + this.B * u + this.C;
		}
	}

	@Override
	public double evaluateOnlyCurrentCost(State currentState, double v, double currentDerivativeValue,
			double transactionCosts) {
		double currentCostLog = this.evaluateOnlyCurrentCostLog(currentState, v, currentDerivativeValue,
				transactionCosts);
		return FastMath.exp(currentCostLog);
	}

	public double evaluateOnlyCurrentCostLog(State currentState, double v, double currentDerivativeValue,
			double transactionCosts) {
		int k = currentState.timeStep;
		float s1 = this.approximationMarketVector.values[MarketVector.PRICE];
		int nextStateIndex = this.approximationMarketVector.simulationIndex;
		double dh = this.derivative.derivativeValueArray[k + 1][nextStateIndex] - currentDerivativeValue;
		double ptfVariation = this.hedgingPortfolio.calculatePortfolioVariation(currentState, s1, v, dh,
				transactionCosts);

		return -this.gamma * ptfVariation;
	}

	@Override
	public double getFactorForDerivative() {
		return this.factorForDerivative;
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

	@Override
	public void initialize() {
		super.initialize();

		TransactionCosts tcStructure = this.hedgingPortfolio.getTCStructure();
		this.aTC = tcStructure.getCostsParameterA();
		this.bTC = tcStructure.getCostsParameterB();
		this.gamma = -((ExponentialLossFunction) this.lf).getParameterB();
		this.s0 = this.approximationOriginState.marketVector.values[MarketVector.PRICE];
	}

	@Override
	public void printBoundaryValues() {
		System.out.println(this.x1 + ";" + this.x2);
	}

	private void computePointsForBoundaryAndCenter() {
		this.m = this.getSlopeParameter();

		double y1 = this.lowerBoundary.riskFunctionValueAtBoundary;
		this.x1 = this.lowerBoundary.boundaryValue;
		double y2 = this.upperBoundary.riskFunctionValueAtBoundary;
		this.x2 = this.upperBoundary.boundaryValue;

		this.x3 = (this.x1 + this.x2) / 2;

		this.approximationOriginState.cashAccount = 0;
		this.approximationOriginState.stockQuantity = this.x3;

		// OneDimSolver.OneDimSolverResults solverResults =
		// this.program.evaluate(this.approximationOriginState);
		// double y3 = solverResults.optim_f;

		QFunction qFunction = this.program.getQFunction();
		qFunction.setState(this.approximationOriginState);
		double y3 = qFunction.evaluate(this.x3);
		// System.out.println(yy + ";" + y3);

		this.log_y1 = Math.log(y1);
		this.log_y2 = Math.log(y2);
		this.log_y3 = Math.log(y3);
	}

	private void evaluateQuadraticApproximationParameters() {
		boolean boundariesAreSeparate = (Math.abs(this.x2 - this.x1) > this.program.getBoundaryTol());
		if (boundariesAreSeparate) {
			try {
				DoubleMatrix1D results = this.computeQuadraticFunctionParameters(this.x1, this.x2, this.x3, this.log_y1,
						this.log_y2, this.log_y3);
				this.A = results.getQuick(0);
				this.B = results.getQuick(1);
				this.C = results.getQuick(2);
			} catch (Exception e) {
				this.setQuadraticApproximationToConstant();
			}
		} else {
			this.setQuadraticApproximationToConstant();
		}
	}

	private double getSlopeParameter() {
		double s0 = this.approximationMarketVector.values[MarketVector.PRICE];
		double slope = this.gamma * (this.aTC + this.bTC * s0);
		return slope;
	}

	private void setQuadraticApproximationToConstant() {
		this.A = 0;
		this.B = 0;
		this.C = this.log_y3;
	}

	protected void computeLQParameters() {
		this.computePointsForBoundaryAndCenter();
		this.evaluateQuadraticApproximationParameters();
	}
}