package derivatives;

import stochprocess.MarketVector;

public class BinaryCallOption extends CallOption {

	public BinaryCallOption(double strike, double maturity, double rfRate, double divYield) {
		super(strike, maturity, rfRate, divYield);
	}

	@Override
	public Object clone() {
		return new CallOption(this.strike, this.maturity, this.rfRate, this.divYield);
	}

	@Override
	public double delta(double s, double sigma, int k) {
		double dT = this.maturity - this.obsTimes[k];
		return OptionPricing.deltaBinaryCallBS(s, this.strike, dT, sigma, this.rfRate, this.divYield);
	}

	@Override
	public double evaluate(int k, float[] values) {
		double dT = this.maturity - this.obsTimes[k];
		double s = values[MarketVector.PRICE];
		double sigma = values[MarketVector.VOL];
		double sigSqrT = sigma * this.getSqrtT()[k];
		double p1 = OptionPricing.binaryCallBSFast(s, this.strike, dT, sigma, this.rfRate, this.divYield, sigSqrT,
				this.exp_rT[k], this.exp_qT[k]);

		return p1;
	}

	@Override
	public String getDescription() {
		return "BinaryCall";
	}

	@Override
	public double payoff(float s) {
		if (s > this.getStrike()) {
			return 1.0;
		} else {
			return 0.0;
		}
	}

	@Override
	public double payoff(float[] s) {
		float terminalPrice = s[s.length - 1];
		return this.payoff(terminalPrice);
	}

}
