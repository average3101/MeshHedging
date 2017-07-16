package derivatives;

public class ShortPutOption extends PutOption {

	public ShortPutOption(double strike, double maturity, double rfRate, double divYield) {
		super(strike, maturity, rfRate, divYield);
	}

	@Override
	public Object clone() {
		return new ShortPutOption(this.getStrike(), this.getMaturity(), this.getRfRate(), this.getDivYield());
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
		return "ShortPut";
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
		return -super.payoff(s);
	}

	@Override
	public double payoff(float[] s) {
		return -super.payoff(s);
	}

}