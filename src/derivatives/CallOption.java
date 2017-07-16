package derivatives;

import org.apache.commons.math3.util.FastMath;

import stochprocess.MarketVector;

public class CallOption extends PutOption {

	public CallOption(double strike, double maturity, double rfRate, double divYield) {
		super(strike, maturity, rfRate, divYield);
	}

	@Override
	public Object clone() {
		return new CallOption(this.strike, this.maturity, this.rfRate, this.divYield);
	}

	@Override
	public double delta(double s, double sigma, int k) {
		double dT = this.maturity - this.obsTimes[k];
		double sigSqrT = sigma * this.sqrtT[k];
		return OptionPricing.deltaCallBSFast(s, this.strike, dT, sigma, this.rfRate, this.divYield, sigSqrT,
				this.exp_rT[k], this.exp_qT[k]);
	}

	@Override
	public double evaluate(int k, float[] values) {
		double dT = this.maturity - this.obsTimes[k];
		double s = values[MarketVector.PRICE];
		double sigma = values[MarketVector.VOL];
		double sigSqrT = sigma * this.sqrtT[k];
		double p1 = OptionPricing.callBSFast(s, this.strike, dT, sigma, this.rfRate, this.divYield, sigSqrT,
				this.exp_rT[k], this.exp_qT[k]);

		return p1;
	}

	@Override
	public String getDescription() {
		return "Call";
	}

	@Override
	public double getMaxDeltaForHedge() {
		return 0.0;
	}

	@Override
	public double getMinDeltaForHedge() {
		return -1.0;
	}

	@Override
	public double payoff(float s) {
		return FastMath.max(0, s - this.strike);
	}

	@Override
	public double payoff(float[] s) {
		return FastMath.max(0, s[s.length - 1] - this.strike);
	}

}