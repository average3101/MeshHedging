package policies;

import derivatives.Derivative;
import dynprog.State;
import stochprocess.MarketVector;

public class BSMHedging extends HedgingPolicy {

	public BSMHedging() {
		this.setPolicyHedgingMethod(HedgingMethod.DH);
	}

	@Override
	public double evaluate(State x) {
		Derivative derivative = this.getDerivative();
		float[] values = x.marketVector.getValues();
		double decision = -derivative.delta(values[MarketVector.PRICE], values[MarketVector.VOL], x.timeStep);
		this.getPastDecisions()[x.timeStep] = decision;

		return decision;
	}

}