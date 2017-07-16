package policies;

import org.apache.commons.math3.util.FastMath;

import derivatives.Derivative;
import dynprog.State;
import dynprog.TransactionCosts;
import stochprocess.MarketVector;

public class HedgingPortfolio {

	private State currentState = new State();
	private Derivative derivative;

	private double initialValue;
	private double interestRatePerPeriod;
	private double[] lastDecisionPath;
	private double[] lastPortfolioValuePath;
	private double lastPortfolioValuePathTransactionCosts;
	private int nbObs;
	private double newCashValueAfterStateChange;
	private double rfRate;
	private TransactionCosts tcStructure;

	public HedgingPortfolio(Derivative derivative, TransactionCosts tcStructure, double rfRate) {
		this.derivative = derivative;
		this.tcStructure = tcStructure;
		this.rfRate = rfRate;
	}

	public double calculateInitialValue(MarketVector marketVector, double initialCash, double initialPosition) {
		return this.calculateInitialValue(marketVector, initialCash, initialPosition,
				this.derivative.evaluate(0, marketVector.values));
	}

	public double calculateInitialValue(MarketVector mktInfo, double initialCash, double initialPosition,
			double derivativeInitialValue) {
		this.initialValue = initialCash + derivativeInitialValue;
		double s = mktInfo.values[MarketVector.PRICE];
		this.initialValue += initialPosition * s;
		return this.initialValue;
	}

	public double calculatePortfolioVariation(State x0, double s1, double u1, double derivativeValueIncrement,
			double transactionCosts) {
		double s0 = x0.marketVector.values[MarketVector.PRICE];
		double c0 = x0.cashAccount;
		double u0 = x0.stockQuantity;

		double c1 = (1 + this.interestRatePerPeriod) * (c0 - (u1 - u0) * s0 + transactionCosts);
		double dV = (c1 + u1 * s1) - (c0 + u0 * s0) + derivativeValueIncrement;

		// double correctionForLossLimit = 0.0;
		// boolean maxLossReached = (dV < -this.maxAllowedOneStepLoss);
		// if (maxLossReached) {
		// correctionForLossLimit = -this.maxAllowedOneStepLoss - dV;
		// }

		this.newCashValueAfterStateChange = c1;

		return dV;
	}

	public double calculateStockAndCashValue(State state) {
		return state.cashAccount + state.stockQuantity * state.marketVector.values[MarketVector.PRICE];
	}

	public double computeDiscreteTimeRate(double[] obsTimes, int nbObs, double rfRate) {
		double dt = obsTimes[nbObs] - obsTimes[0];
		double f = FastMath.exp(rfRate * dt);
		return FastMath.pow(f, 1.0 / nbObs) - 1.0;
	}

	public double evaluate(State x) {
		MarketVector marketVector = x.marketVector;
		double derivativeValue = this.derivative.evaluate(x.timeStep, marketVector.values);
		float s = marketVector.values[MarketVector.PRICE];
		double value = x.cashAccount + x.stockQuantity * s + derivativeValue;

		return value;
	}

	public Derivative getDerivative() {
		return this.derivative;
	}

	public double getInitialValue() {
		return this.initialValue;
	}

	public double getInterestRatePerPeriod() {
		return this.interestRatePerPeriod;
	}

	public double[] getLastDecisionPath() {
		return this.lastDecisionPath;
	}

	public double[] getLastSimulatedPath() {
		return this.lastPortfolioValuePath;
	}

	public double getLastSimulatedPathTransactionCosts() {
		return this.lastPortfolioValuePathTransactionCosts;
	}

	public TransactionCosts getTCStructure() {
		return this.tcStructure;
	}

	public void setInitialValue(double v) {
		this.initialValue = v;
	}

	public void setInterestRatePerPeriod(double interestRatePerPeriod) {
		this.interestRatePerPeriod = interestRatePerPeriod;
	}

	public void setMaxAllowedOneStepLoss(double maxLoss) {
	}

	public void setObservationTimes(double[] obsTimes) {
		this.nbObs = obsTimes.length - 1;

		this.lastPortfolioValuePath = new double[this.nbObs + 1];
		this.lastDecisionPath = new double[this.nbObs + 1];

		this.derivative.setObsTimes(obsTimes);

		this.setInterestRatePerPeriod(this.computeDiscreteTimeRate(obsTimes, this.nbObs, this.rfRate));
	}

	public double terminalPtfValueWithPolicy(State startingState, MarketVector[] mktInfoPath, HedgingPolicy policy) {
		int tIndex = startingState.timeStep;
		double[] derivativePath = this.derivative.evaluatePath(mktInfoPath);

		this.lastPortfolioValuePathTransactionCosts = 0;
		this.lastPortfolioValuePath[tIndex] = this.calculateStockAndCashValue(startingState) + derivativePath[tIndex];

		for (int i = tIndex; i < this.nbObs; i++) {
			MarketVector y1 = mktInfoPath[i + 1];
			double u1 = policy.evaluate(startingState);
			this.lastDecisionPath[i] = u1;

			double transactionCosts = this.calculateTransactionCosts(startingState, u1);
			this.calculatePortfolioVariation(startingState, y1.values[MarketVector.PRICE], u1,
					derivativePath[i + 1] - derivativePath[i], transactionCosts);

			double previousTotalCosts = this.lastPortfolioValuePathTransactionCosts;
			this.lastPortfolioValuePathTransactionCosts = (previousTotalCosts + transactionCosts)
					* (1 + this.interestRatePerPeriod);

			startingState = this.currentState;
			this.currentState = this.defineNewCurrentState(i, y1, u1);

			this.lastPortfolioValuePath[i + 1] = this.calculateStockAndCashValue(this.currentState)
					+ derivativePath[i + 1];
		}

		this.lastDecisionPath[this.nbObs] = this.lastDecisionPath[this.nbObs - 1];

		return this.lastPortfolioValuePath[this.nbObs];
	}

	private double calculateTransactionCosts(State startingState, double u1) {
		double u0 = startingState.stockQuantity;
		double s0 = startingState.marketVector.values[MarketVector.PRICE];
		double transactionCosts = -this.tcStructure.evaluate(u1 - u0, s0);
		return transactionCosts;
	}

	private State defineNewCurrentState(int i, MarketVector y1, double u1) {
		this.currentState.cashAccount = this.newCashValueAfterStateChange;
		this.currentState.marketVector = y1;
		this.currentState.stockQuantity = u1;
		this.currentState.timeStep = i + 1;
		return this.currentState;
	}
}
