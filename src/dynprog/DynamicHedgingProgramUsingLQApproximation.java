package dynprog;

import policies.HedgingPolicy;

public class DynamicHedgingProgramUsingLQApproximation extends DynamicHedgingProgram {

	public DynamicHedgingProgramUsingLQApproximation() {
		super();

		this.setPolicyHedgingMethod(HedgingPolicy.HedgingMethod.SM_LQ);
	}

	@Override
	public void computeForState(int k, int i) {
		State stateForApprox = this.getTempStateOrigin(k, i);
		this.setMeshRelatedValuesIncludingDerivativeValueComputation(k, i);

		LossFunction lossFunction = this.getLossFunction();
		double riskAversionParameter = lossFunction.getRiskAversionParameter();

		HedgingBoundary upperBoundary = new HedgingBoundary(this, riskAversionParameter,
				this.getQFunctionForUpperBoundary());
		HedgingBoundary lowerBoundary = new HedgingBoundary(this, riskAversionParameter,
				this.getQFunctionForLowerBoundary());

		lowerBoundary.estimate(stateForApprox);
		upperBoundary.estimate(stateForApprox);

		RiskFunctionApproximation approximation = new RiskFunctionApproxLinQuad(this, stateForApprox, lowerBoundary,
				upperBoundary);

		this.setCostFunctionApproximationForPredefinedState(k, i, approximation);
	}

}
