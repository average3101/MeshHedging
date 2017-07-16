package derivatives;

import org.apache.commons.math3.util.FastMath;

import meshmethods.StochasticMesh;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;

public abstract class Derivative {

	private double[][] derivativeCVArray;
	public double[][] derivativeValueArray;
	private double gammaFiniteDiffPrecision = 0.001;
	private double initialValue = Double.NaN;
	private double minDelta, maxDelta;
	private int nbObs;
	private double[] valuePath;
	protected double divYield;
	protected double[] exp_qT;
	protected double[] exp_rT;
	protected double maturity;
	protected double[] obsTimes;
	protected double rfRate;
	protected double[] sqrtT;
	protected double strike;
	protected MarketVectorProcess process;

	@Override
	public abstract Object clone();

	public abstract double delta(double s, double sigma, int k);

	public abstract double evaluate(int k, float[] values);

	public double evaluateAtT0(MarketVector y) {
		this.initialValue = this.evaluate(0, y.values);
		return this.initialValue;
	}

	public double evaluateBySimulation(MarketVectorProcess sp, int nbPaths) {
		double value = 0;
		int nbObs = this.getNbObs();
		for (int i = 0; i < nbPaths; i++) {
			sp.resetNextSubstream();
			float[][] basicPaths = sp.generateBasicProcessesPaths();
			float sk = basicPaths[MarketVector.PRICE][nbObs];
			double p = this.payoff(sk);
			value += p;
		}
		value /= nbPaths;

		return value;
	}

	public double[] evaluatePath(MarketVector[] marketVectorPath) {
		double[] valuePath = this.getValuePath();
		valuePath[0] = this.initialValue;
		int nbObs = this.getNbObs();
		for (int k = 1; k < nbObs; k++) {
			valuePath[k] = this.evaluate(k, marketVectorPath[k].values);
		}
		float sk = marketVectorPath[nbObs].values[MarketVector.PRICE];
		valuePath[nbObs] = this.payoff(sk);

		return valuePath;
	}

	public double gamma(int k, float[] values) {
		return this.gammaFD(k, values);
	}

	public double gammaFD(int k, float[] values) {
		float s = values[MarketVector.PRICE];
		float s_l = (float) (s * (1.0 + this.gammaFiniteDiffPrecision));
		float s_h = (float) (s * (1.0 - this.gammaFiniteDiffPrecision));

		values[MarketVector.PRICE] = s_l;
		double p_low = this.evaluate(k, values);
		values[MarketVector.PRICE] = s_h;
		double p_high = this.evaluate(k, values);
		values[MarketVector.PRICE] = s;
		double p_mid = this.evaluate(k, values);

		return (p_low + p_high - 2.0 * p_mid) / FastMath.pow(this.gammaFiniteDiffPrecision * s, 2.0);
	}

	public double[][] getDerivativeCVArray() {
		return this.derivativeCVArray;
	}

	public double[][] getDerivativeValueArray() {
		return this.derivativeValueArray;
	}

	public abstract String getDescription();

	public double getDivYield() {
		return this.divYield;
	}

	public double[] getExp_qT() {
		return this.exp_qT;
	}

	public double[] getExp_rT() {
		return this.exp_rT;
	}

	public double getInitDerivativeValue() {
		return this.initialValue;
	}

	public double getMaturity() {
		return this.maturity;
	}

	public double getMaxDeltaForHedge() {
		return this.maxDelta;
	}

	public double getMinDeltaForHedge() {
		return this.minDelta;
	}

	public int getNbObs() {
		return this.nbObs;
	}

	public double[] getObsTimes() {
		return this.obsTimes;
	}

	public double getRfRate() {
		return this.rfRate;
	}

	public double[] getSqrtT() {
		return this.sqrtT;
	}

	public double getStrike() {
		return this.strike;
	}

	public double[] getValuePath() {
		return this.valuePath;
	}

	public void init() {
		this.setExp_qT(new double[this.getNbObs() + 1]);
		this.setExp_rT(new double[this.getNbObs() + 1]);
		this.setSqrtT(new double[this.getNbObs() + 1]);

		for (int k = 0; k <= this.getNbObs(); k++) {
			double dT = FastMath.max(this.maturity - this.getObsTimes()[k], 0.0);
			this.getExp_rT()[k] = FastMath.exp(-this.getRfRate() * dT);
			this.getExp_qT()[k] = FastMath.exp(-this.getDivYield() * dT);
			this.getSqrtT()[k] = FastMath.sqrt(dT);
		}
	}

	public void initializeAfterMeshComputation(int nbObs, int nbStatesPerStep, StochasticMesh stochasticMesh) {
		this.derivativeValueArray = new double[nbObs + 1][nbStatesPerStep];

		if (stochasticMesh.isUsingDifferentProcessForPricing()) {
			this.process = stochasticMesh.getMarketProcessForPricing();
		} else {
			this.process = stochasticMesh.getMarketVectorProcess();
		}

		if (this instanceof DerivativePricedOnMeshWithCVAnalytic) {
			((DerivativePricedOnMeshWithCVAnalytic) this).control.marketProcess = this.process;
		}
	}

	public abstract double payoff(float s);

	public abstract double payoff(float[] sPath);

	public void setDerivativeCVArray(double[][] derivativeCVArray) {
		this.derivativeCVArray = derivativeCVArray;
	}

	public void setDerivativeValueArrayValue(int k, int i, double h) {
		this.derivativeValueArray[k][i] = h;
	}

	public void setDivYield(double divYield) {
		this.divYield = divYield;
	}

	public void setExp_qT(double[] exp_qT) {
		this.exp_qT = exp_qT;
	}

	public void setExp_rT(double[] exp_rT) {
		this.exp_rT = exp_rT;
	}

	public void setInitDerivativeValue(double v) {
		this.initialValue = v;
	}

	public void setMaturity(double t) {
		this.maturity = t;
	}

	public void setMaxDeltaForHedge(double v) {
		this.maxDelta = v;
	}

	public void setMinDeltaForHedge(double v) {
		this.minDelta = v;
	}

	public void setNbObs(int nbObs) {
		this.nbObs = nbObs;
	}

	public void setObsTimes(double[] obsTimes) {
		this.obsTimes = obsTimes;
		this.setNbObs(obsTimes.length - 1);

		this.setValuePath(new double[this.getNbObs() + 1]);

		this.init();
	}

	public void setRfRate(double rfRate) {
		this.rfRate = rfRate;
	}

	public void setSqrtT(double[] sqrtT) {
		this.sqrtT = sqrtT;
	}

	public void setStrike(double strike) {
		this.strike = strike;
	}

	public void setValuePath(double[] valuePath) {
		this.valuePath = valuePath;
	}

}