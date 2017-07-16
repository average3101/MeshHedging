package dynprog;

import policies.HedgingPolicy;

public class DynamicHedgingProgramApproxLocalOneStep extends DynamicHedgingProgram {

	public DynamicHedgingProgramApproxLocalOneStep() {
		super();

		setPolicyHedgingMethod(HedgingPolicy.HedgingMethod.SM_LOC_1S);
	}

	@Override
	public void computeForState(int k, int i) {
		State stateForApprox = getTempStateOrigin(k, i);
		setMeshRelatedValuesIncludingDerivativeValueComputation(k, i);

		RiskFunctionApproximation approximation = new RiskFunctionApproxLocalOneStep(this, stateForApprox);
		setCostFunctionApproximationForPredefinedState(k, i, approximation);
	}

}
