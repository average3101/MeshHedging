package derivatives;

public class ShortCallOption extends CallOption {

	public ShortCallOption(double strike, double maturity, double rfRate, double divYield) {
		super(strike, maturity, rfRate, divYield);
	}

	@Override
	public Object clone() {
		return new ShortCallOption(this.getStrike(), this.getMaturity(), this.getRfRate(), this.getDivYield());
	}

	@Override
	public double delta(double s, double sigma, int k) {
		return -super.delta(s, sigma, k);
	}

	@Override
	public double evaluate(int k, float[] values) {
		return -super.evaluate(k, values);
	}

	@Override
	public double gamma(int k, float[] values) {
		return -super.gamma(k, values);
	}

	@Override
	public String getDescription() {
		return "ShortCall";
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
		return -super.payoff(s);
	}

	@Override
	public double payoff(float[] s) {
		return -super.payoff(s);
	}

}