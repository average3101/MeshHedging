package dynprog;

import stochprocess.MarketVector;

public class State {
	public static void copy(State originalState, State destinationState) {
		destinationState.timeStep = originalState.timeStep;
		destinationState.cashAccount = originalState.cashAccount;
		destinationState.stockQuantity = originalState.stockQuantity;
		MarketVector.copy(originalState.marketVector, destinationState.marketVector);
	}

	public double cashAccount;
	public MarketVector marketVector;
	public double stockQuantity;

	public int timeStep;

	public State() {
	}

	public State(int k, double c, double u, MarketVector y) {
		this.timeStep = k;
		this.cashAccount = c;
		this.stockQuantity = u;
		this.marketVector = y;
	}

	@Override
	public State clone() {
		State newState = new State(this.timeStep, this.cashAccount, this.stockQuantity, this.marketVector.clone());
		return newState;
	}

}