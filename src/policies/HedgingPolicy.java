package policies;

import derivatives.Derivative;
import dynprog.DynamicHedgingProgram;
import dynprog.State;
import experiments.SimulationParameters;
import stochprocess.MarketVectorProcess;

public abstract class HedgingPolicy {

	public enum HedgingMethod {
		DH, LOC_1SA, NH, OSL, SM_LOC_1S, SM_LQ, WW, Z, ZQF
	};

	private Derivative derivative;
	private State initialState;
	private String label;
	private int nbObs;
	private double[] obsTimes;
	private double[] pastDecisions;
	private HedgingMethod policyHedgingMethod = null;

	public void computeCummulApproxError() {
	}

	public void computeMesh() {
	}

	public void computePolicy() {
	}

	public abstract double evaluate(State x);

	public Derivative getDerivative() {
		return this.derivative;
	}

	public double getErrorEstimateAtInitialeState() {
		return 0;
	}

	public State getInitialState() {
		return this.initialState;
	}

	public String getLabel() {
		return this.label;
	}

	public double[] getObservationTimes() {
		return this.obsTimes;
	}

	public double[] getPastDecisions() {
		return this.pastDecisions;
	}

	public HedgingMethod getPolicyHedgingMethod() {
		return this.policyHedgingMethod;
	}

	public double[] getPolicyValuesAsFunctionOfS0(int nbPts, double uMin, double uMax, int k, double c0, double u0,
			double sigma1, double sigma2, MarketVectorProcess sp) {

		double[] uData = this.defineEvaluationPoints(nbPts, uMin, uMax);
		State stateForGraph = this.defineStateForGraph(k, c0, u0);

		float[] y = new float[3];
		y[1] = (float) (sigma1 * sigma1);
		y[2] = (float) (sigma2 * sigma2);

		double[] decisionData = new double[nbPts];
		for (int j = 0; j < nbPts; j++) {
			float u = (float) uData[j]; //

			y[0] = u;
			stateForGraph.marketVector = sp.newMarketVector(y);

			double v_s = this.evaluate(stateForGraph);
			decisionData[j] = v_s;
		}
		return decisionData;
	}

	public double[] getPolicyValuesAsFunctionOfU0(int nbPts, double uMin, double uMax, int k, double c0, double s0,
			double sigma1, double sigma2, MarketVectorProcess sp) {

		double[] uData = this.defineEvaluationPoints(nbPts, uMin, uMax);
		double u0 = 0;
		State stateForGraph = this.defineStateForGraph(k, c0, u0);

		float[] y = new float[3];
		y[0] = (float) s0;
		y[1] = (float) (sigma1 * sigma1);
		y[2] = (float) (sigma2 * sigma2);

		stateForGraph.marketVector = sp.newMarketVector(y);
		double[] decisionData = new double[nbPts];
		for (int j = 0; j < nbPts; j++) {
			double u = uData[j]; //
			stateForGraph.stockQuantity = u;

			double v_s = this.evaluate(stateForGraph);
			decisionData[j] = v_s;
		}
		return decisionData;
	}

	public double getRiskEstimateAtInitialState(double u0, double c0) {
		return 0.0;
	}

	public void init() { // Policy dependent
		this.pastDecisions = new double[this.nbObs];
	}

	public String outputCode() {
		return this.getPolicyHedgingMethod().toString();
	}

	public void printPolicyValuesToScreen(int nbPts, double uMin, double uMax, int k, double c0, double s0,
			double sigma1, double sigma2, MarketVectorProcess sp) {
		double[] uData = this.defineEvaluationPoints(nbPts, uMin, uMax);
		double u0 = 0;
		State stateForGraph = this.defineStateForGraph(k, c0, u0);

		float[] y = new float[3];
		y[0] = (float) s0;
		y[1] = (float) (sigma1 * sigma1);
		y[2] = (float) (sigma2 * sigma2);

		stateForGraph.marketVector = sp.newMarketVector(y);

		System.out.println(this.getPolicyHedgingMethod() + ";*********");
		System.out.println("index;u0;v*");
		for (int j = 0; j < nbPts; j++) {
			double u = uData[j]; //
			stateForGraph.stockQuantity = u;

			double v_s = this.evaluate(stateForGraph);
			System.out.println(j + ";" + u + ";" + v_s);
		}
	}

	public void setInitialState(State initialState) {
		this.initialState = initialState;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public void setObservationTimes(double[] t) {
		this.obsTimes = t;
		this.nbObs = t.length - 1;
		this.init();
	}

	public void setParameters(SimulationParameters simParams) {
		this.derivative = simParams.derivative;

		if (this instanceof RiskFunctionBasedPolicy) {
			DynamicHedgingProgram program = ((RiskFunctionBasedPolicy) this).getDynamicHedgingProgram();
			program.setSimulationParameters(simParams);
		}
	}

	public void setPolicyHedgingMethod(HedgingMethod policyHedgingMethod) {
		this.policyHedgingMethod = policyHedgingMethod;
	}

	private double[] defineEvaluationPoints(int nbPts, double uMin, double uMax) {
		double[] uData = new double[nbPts];
		double du = (uMax - uMin) / (nbPts - 1.0);
		for (int j = 0; j < nbPts; j++) {
			uData[j] = uMin + j * du;
		}
		return uData;
	}

	private State defineStateForGraph(int k, double c0, double u0) {
		State stateForGraph = new State();
		stateForGraph.cashAccount = c0;
		stateForGraph.timeStep = k;
		stateForGraph.stockQuantity = u0;
		return stateForGraph;
	}

}