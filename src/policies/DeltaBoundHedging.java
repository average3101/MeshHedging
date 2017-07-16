package policies;

import dynprog.State;
import experiments.SimulationParameters;

public abstract class DeltaBoundHedging extends HedgingPolicy {

	private double aTC;
	private double bTC;
	private double maturity;
	private double rfRate;
	private double riskAversionCoeff;

	public abstract double computeDelta(State x);

	public abstract double computeDeltaBound(State x);

	@Override
	public double evaluate(State x) {
		int tIndex = x.timeStep;
		double u0 = x.stockQuantity;

		this.precomputeCommonQuantities(x);
		double delta = this.computeDelta(x);
		double bound = this.computeDeltaBound(x);

		double decision = u0; // By default
		if (u0 > (delta + bound)) {
			decision = delta + bound;
		} else if (u0 < (delta - bound)) {
			decision = delta - bound;
		}
		this.getPastDecisions()[tIndex] = decision;

		return decision;
	}

	public double getaTC() {
		return this.aTC;
	}

	public double getbTC() {
		return this.bTC;
	}

	public double getMaturity() {
		return this.maturity;
	}

	public double getRfRate() {
		return this.rfRate;
	}

	public double getRiskAversionCoeff() {
		return this.riskAversionCoeff;
	}

	public abstract void precomputeCommonQuantities(State x);

	public void setaTC(double aTC) {
		this.aTC = aTC;
	}

	public void setbTC(double bTC) {
		this.bTC = bTC;
	}

	public void setMaturity(double maturity) {
		this.maturity = maturity;
	}

	@Override
	public void setParameters(SimulationParameters simParams) {
		super.setParameters(simParams);
		this.setMaturity(this.getDerivative().getMaturity());
		this.setRfRate(simParams.rfRate);
		this.setRiskAversionCoeff(simParams.riskAversion);
		this.setaTC(simParams.aTC);
		this.setbTC(simParams.bTC);
	}

	public void setRfRate(double rfRate) {
		this.rfRate = rfRate;
	}

	public void setRiskAversionCoeff(double riskAversionCoeff) {
		this.riskAversionCoeff = riskAversionCoeff;
	}

}
