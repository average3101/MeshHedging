package derivatives;

import org.apache.commons.math3.util.FastMath;

import stochprocess.MarketVector;

public class PutOption extends Derivative {

	public PutOption(double strike, double maturity, double rfRate, double divYield) {
		this.setStrike(strike);
		this.setMaturity(maturity);
		this.setRfRate(rfRate);
		this.setDivYield(divYield);
	}

	@Override
	public Object clone() {
		return new PutOption(this.strike, this.maturity, this.rfRate, this.divYield);
	}

	@Override
	public double delta(double s, double sigma, int k) {
		double dT = this.maturity - this.obsTimes[k];
		double sigSqrT = sigma * this.sqrtT[k];

		return OptionPricing.deltaPutBSFast(s, this.strike, dT, sigma, this.rfRate, this.divYield, sigSqrT,
				this.exp_rT[k], this.exp_qT[k]);
	}

	@Override
	public double evaluate(int k, float[] values) {
		double dT = this.maturity - this.obsTimes[k];
		double s = values[MarketVector.PRICE];
		double sigma = values[MarketVector.VOL];
		double sigSqrT = sigma * this.getSqrtT()[k];
		double p1 = OptionPricing.putBSFast(s, this.strike, dT, sigma, this.rfRate, this.divYield, sigSqrT,
				this.exp_rT[k], this.exp_qT[k]);

		return p1;
	}

	@Override
	public double gamma(int k, float[] values) {
		double dT = this.maturity - this.obsTimes[k];
		double s = values[MarketVector.PRICE];
		double sigma = values[MarketVector.VOL];
		return OptionPricing.gammaBS(s, this.strike, dT, sigma, this.rfRate, this.divYield);
	}

	@Override
	public String getDescription() {
		return "Put";
	}

	@Override
	public double getMaxDeltaForHedge() {
		return 1.0;
	}

	@Override
	public double getMinDeltaForHedge() {
		return 0.0;
	}

	@Override
	public double payoff(float s) {
		return FastMath.max(0, this.getStrike() - s);
	}

	@Override
	public double payoff(float[] s) {
		return FastMath.max(0, this.getStrike() - s[s.length - 1]);
	}

}