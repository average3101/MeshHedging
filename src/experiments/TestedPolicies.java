package experiments;

import dynprog.DynamicHedgingProgram;
import dynprog.DynamicHedgingProgramApproxLocalOneStep;
import dynprog.DynamicHedgingProgramUsingLQApproximation;
import dynprog.ExponentialLossFunction;
import dynprog.ProportionalCosts;
import policies.BSMHedging;
import policies.HedgingPolicy;
import policies.HedgingPolicy.HedgingMethod;
import policies.Local1SApproxHedging;
import policies.NoHedging;
import policies.RiskFunctionBasedPolicy;
import policies.WWHedging;
import policies.ZQF2006Hedging;

public class TestedPolicies {

	private boolean includesMeshBasedEstimator = false;
	private int meshBasedEstimatorIndex = -1;
	private int nbDecisionPts = 101;
	private int nbPolicies;
	private HedgingPolicy[] policyArray;
	private double[][] policyDecisionsArray;

	public TestedPolicies(boolean excludeLocalHedging) {
		this.definePolicies(excludeLocalHedging);
	}

	public void computeInitialValues(SimulationParameters simParams, ProblemDefinition problem,
			SimulationResults results) {
		double initialDerivativeValue = simParams.derivative.evaluateAtT0(problem.mktInfoPath[0]);
		results.hInitValue.add(initialDerivativeValue);
		problem.hedgedPtf.calculateInitialValue(problem.mktInfoPath[0], problem.c0Start, simParams.u0,
				initialDerivativeValue);

	}

	public void computeMesh(SimulationParameters simParams, ProblemDefinition problem, SimulationResults results) {
		double t0 = results.timer.getSeconds();

		problem.computeMesh(simParams);

		double t1 = results.timer.getSeconds();
		results.meshCompTime[0] += t1 - t0;
	}

	public void computePolicies(SimulationParameters simParams, ProblemDefinition problem, SimulationResults results) {
		this.setPolicyDecisionsArray(new double[this.getNbPolicies() + 1][this.getNbDecisionPts()]);

		for (int policyNo = 0; policyNo < this.getNbPolicies(); policyNo++) {
			double initialTime = results.timer.getSeconds();

			HedgingPolicy currentPolicy = this.getPolicyArray()[policyNo];
			HedgingMethod currentHedgingMethod = currentPolicy.getPolicyHedgingMethod();
			currentPolicy.setObservationTimes(problem.obsTimes);

			if (currentHedgingMethod == HedgingPolicy.HedgingMethod.SM_LQ) {
				double fp = ((RiskFunctionBasedPolicy) currentPolicy).getFractionAttainableStates();
				double fw = ((RiskFunctionBasedPolicy) currentPolicy).getFractionOfWeightsUsed();
				results.fracAttainablePoints.add(fp);
				results.fracWeightsUsed.add(fw);
			}

			if (currentHedgingMethod == HedgingPolicy.HedgingMethod.SM_LQ
					|| currentHedgingMethod == HedgingPolicy.HedgingMethod.SM_LOC_1S) {
				((RiskFunctionBasedPolicy) currentPolicy).initializeAfterMeshComputation();

				currentPolicy.computePolicy();

				// FunctionPlotter.printOptimalDecisionForQFunction(0,
				// (RiskFunctionBasedPolicy) currentPolicy, 21);

				// FunctionPlotter.printLinearizedQFunctions(0,
				// (RiskFunctionBasedPolicy) currentPolicy);

				if (currentHedgingMethod == HedgingPolicy.HedgingMethod.SM_LQ) {
					currentPolicy.computeCummulApproxError();
				}
			}

			results.dpCompTime[policyNo] += (results.timer.getSeconds() - initialTime);
		}
	}

	public void definePolicies(boolean excludeLocalHedging) {
		if (excludeLocalHedging) {
			this.nbPolicies = 6;
		} else {
			this.nbPolicies = 7;
		}

		this.policyArray = new HedgingPolicy[this.nbPolicies];

		this.policyArray[0] = new BSMHedging();
		this.policyArray[1] = new WWHedging();
		this.policyArray[2] = new ZQF2006Hedging();
		this.policyArray[3] = new NoHedging();

		DynamicHedgingProgram jFunction = new DynamicHedgingProgramUsingLQApproximation();
		this.policyArray[4] = new RiskFunctionBasedPolicy(jFunction);
		this.meshBasedEstimatorIndex = 4;
		this.includesMeshBasedEstimator = true;

		if (excludeLocalHedging) {
			this.getPolicyArray()[5] = new Local1SApproxHedging();
		} else {
			DynamicHedgingProgram jFunctionLocalOneStep = new DynamicHedgingProgramApproxLocalOneStep();
			this.getPolicyArray()[5] = new RiskFunctionBasedPolicy(jFunctionLocalOneStep);
			this.getPolicyArray()[6] = new Local1SApproxHedging();
		}
	}

	public int getMeshBasedEstimatorIndex() {
		return this.meshBasedEstimatorIndex;
	}

	public int getNbDecisionPts() {
		return this.nbDecisionPts;
	}

	public int getNbPolicies() {
		return this.nbPolicies;
	}

	public HedgingPolicy[] getPolicyArray() {
		return this.policyArray;
	}

	public double[][] getPolicyDecisionsArray() {
		return this.policyDecisionsArray;
	}

	public boolean includesMeshBasedEstimator() {
		return this.includesMeshBasedEstimator;
	}

	public void initializePolicySettings(SimulationParameters simParams, ProblemDefinition problem) {
		problem.lossFunction = new ExponentialLossFunction(-simParams.riskAversion);
		problem.tcStructure = new ProportionalCosts(simParams.aTC, simParams.bTC);

		for (int i = 0; i < this.getNbPolicies(); i++) {
			HedgingPolicy p = this.getPolicyArray()[i];
			p.setParameters(simParams);
			if (p instanceof RiskFunctionBasedPolicy) {
				DynamicHedgingProgram program = ((RiskFunctionBasedPolicy) p).getDynamicHedgingProgram();
				program.setProblemDefinition(problem);
				program.setSimulationParameters(simParams);
			}
		}
	}

	public void setNbDecisionPts(int nbDecisionPts) {
		this.nbDecisionPts = nbDecisionPts;
	}

	public void setPolicyDecisionsArray(double[][] policyDecisionsArray) {
		this.policyDecisionsArray = policyDecisionsArray;
	}

}
