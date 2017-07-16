package policies;

import dynprog.State;

public class NoHedging extends HedgingPolicy {

	public NoHedging() {
		setPolicyHedgingMethod(HedgingMethod.NH);
	}

	@Override
	public double evaluate(State x) {
		double decision = x.stockQuantity;
		getPastDecisions()[x.timeStep] = decision;
		return decision;
	}

}
