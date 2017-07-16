package dynprog;

import policies.HedgingPortfolio;

public class QFunctionForTradeBoundaries extends QFunction {

	private HedgingPortfolio hedgingPortfolio;
	private TransactionCosts transactionCosts;

	public QFunctionForTradeBoundaries(DynamicHedgingProgram program) {
		super(program);
		this.transactionCosts = program.getTransactionCosts();
		this.hedgingPortfolio = program.getHedgingPortfolio();
	}

	@Override
	public double computeFee(double du, double s0) {
		double signOfTradeBoundary = this.getSignOfTradeBoundary();
		double fee = signOfTradeBoundary * du * this.transactionCosts.evaluatePerUnitStock(s0)
				* (1 + this.hedgingPortfolio.getInterestRatePerPeriod());
		return fee;
	}

	@Override
	public void setSignOfTradeNegative() {
		this.setSignOfTradeBoundary(-1);
	}

	@Override
	public void setSignOfTradePositive() {
		this.setSignOfTradeBoundary(+1);
	}

	// @Override
	// protected void defineOneDimSolver() {
	// this.oneDimSolver = new BrentSolver(this);
	// // this.oneDimSolver = new BrentSolverWithNewtonStart(this);
	// // this.oneDimSolver = new SecantQuadSolver(this); // n
	// // this.oneDimSolver = new SecantDiffSolver(this); // new
	// // this.oneDimSolver = new SecantSolver(this); // new
	// // this.oneDimSolver = new QuadraticSeqSolver(this); // new
	// // this.oneDimSolver = new QuadraticSolver(this); // new
	// // this.oneDimSolver = new NewtonSolverWithGradientCheck(this); // new
	// // BrentSolverWithGradientCheck(this);
	// // this.oneDimSolver = new BisectionSolver(this);
	// // this.oneDimSolver = new BrentSolverWithGradientCheck(this);
	// // this.oneDimSolver = new UncminSolver(this);
	// }
}