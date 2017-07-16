package derivatives;

import stochprocess.MarketVector;

public class DerivativeBasket extends Derivative {

	Derivative[] derivativeArray;
	String description;
	double[] positionArray;

	public DerivativeBasket(String description, Derivative[] derivativeArray, double[] positionArray,
			double minDeltaForHedge, double maxDeltaForHedge) {
		this.derivativeArray = derivativeArray;
		this.positionArray = positionArray;
		this.description = description;
		this.setMinDeltaForHedge(minDeltaForHedge);
		this.setMaxDeltaForHedge(maxDeltaForHedge);
		// Assuming all maturities are the same...
		this.setMaturity(derivativeArray[0].getMaturity());
	}

	@Override
	public Object clone() {
		Derivative[] newDerivativeArray = new Derivative[this.derivativeArray.length];
		double[] newPositionArray = new double[this.positionArray.length];
		int n = this.derivativeArray.length;
		for (int i = 0; i < n; i++) {
			newPositionArray[i] = this.positionArray[i];
			newDerivativeArray[i] = (Derivative) this.derivativeArray[i].clone();
		}
		return new DerivativeBasket(this.description, newDerivativeArray, newPositionArray, this.getMinDeltaForHedge(),
				this.getMaxDeltaForHedge());
	}

	@Override
	public double delta(double s, double sigma, int k) {
		int n = this.derivativeArray.length;
		double delta = 0;
		for (int i = 0; i < n; i++) {
			delta += this.positionArray[i] * this.derivativeArray[i].delta(s, sigma, k);
		}
		return delta;
	}

	@Override
	public double evaluate(int k, float[] values) {
		int n = this.derivativeArray.length;
		double delta = 0;
		for (int i = 0; i < n; i++) {
			delta += this.positionArray[i] * this.derivativeArray[i].evaluate(k, values);
		}
		return delta;
	}

	@Override
	public double evaluateAtT0(MarketVector y) {
		int n = this.derivativeArray.length;
		double p0 = 0;
		for (int i = 0; i < n; i++) {
			p0 += this.positionArray[i] * this.derivativeArray[i].evaluateAtT0(y);
		}

		this.setInitDerivativeValue(p0);
		return p0;
	}

	@Override
	public double gamma(int k, float[] values) {
		int n = this.derivativeArray.length;
		double delta = 0;
		for (int i = 0; i < n; i++) {
			delta += this.positionArray[i] * this.derivativeArray[i].gamma(k, values);
		}
		return delta;
	}

	@Override
	public String getDescription() {
		return this.description;
	}

	@Override
	public double getStrike() {
		return this.derivativeArray[0].getStrike();
	}

	@Override
	public double payoff(float s) {
		int n = this.derivativeArray.length;
		double delta = 0;
		for (int i = 0; i < n; i++) {
			delta += this.positionArray[i] * this.derivativeArray[i].payoff(s);
		}
		return delta;
	}

	@Override
	public double payoff(float[] s) {
		int n = this.derivativeArray.length;
		double delta = 0;
		for (int i = 0; i < n; i++) {
			delta += this.positionArray[i] * this.derivativeArray[i].payoff(s);
		}
		return delta;
	}

	@Override
	public void setMaturity(double t) {
		super.setMaturity(t);

		int n = this.derivativeArray.length;
		for (int i = 0; i < n; i++) {
			this.derivativeArray[i].setMaturity(t);
		}
	}

	@Override
	public void setObsTimes(double[] obsTimes) {
		super.setObsTimes(obsTimes);

		int n = this.derivativeArray.length;
		for (int i = 0; i < n; i++) {
			this.derivativeArray[i].setObsTimes(obsTimes);
		}
	}

}
