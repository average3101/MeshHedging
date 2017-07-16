package dynprog;

import cern.colt.list.FloatArrayList;
import derivatives.Derivative;
import meshmethods.MeshWeights;
import numerics.OneDimSolverResults;
import stochprocess.MarketVector;

public class DynamicHedgingProgramDiagnostics {

	private double[] jData;
	private int nbDataPointsForTesting;
	private DynamicHedgingProgram program;
	private double[] uData;

	public DynamicHedgingProgramDiagnostics(DynamicHedgingProgram dynamicHedgingProgram, int nbDataPointsForTesting) {
		this.program = dynamicHedgingProgram;
		this.nbDataPointsForTesting = nbDataPointsForTesting;
	}

	public void printApproxCostFunctionValuesToScreen(int k, int i) {
		this.initializeValuesForTests();
		// Note: Approx based on zero initial cash
		double c = 0;
		// Compute 'J*' (objective) values corresponding to initial 'u' values
		System.out.println("Approx cost function values for step " + k + ", state " + i);
		System.out.println("pointNo;u;j_approx");
		RiskFunctionApproximation approximation = this.program.jStarArray[k][i];
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			double u = this.uData[j];
			double j_s = approximation.evaluate(c, u);
			System.out.println(j + ";" + u + ";" + j_s);
		}
	}

	public void printDerivativeValuesAtStep(int k) {
		System.out.println();
		System.out.println("Derivative values at step " + k);
		System.out.println("i;sk;hk;hk_cv");

		for (int i = 0; i < this.program.getNbStatesPerStep(); i++) {
			MarketVector marketVector = this.program.getTempMarketVectorDestination(k, i);
			String output = i + "; " + marketVector.values[MarketVector.PRICE] + "; "
					+ this.program.getDerivative().getDerivativeValueArray()[k][i];
			System.out.println(output);
		}
	}

	public void printExactCostFunctionValuesForTradeBoundariesToScreen(int k, int i, double sign) {
		this.initializeValuesForTests();
		State stateForGraph = new State();
		stateForGraph.cashAccount = 0.0;

		stateForGraph.timeStep = k;
		stateForGraph.marketVector = this.program.getTempMarketVectorDestination(k, i);
		this.program.setMeshRelatedValuesWithPrecomputedDerivativeValue(k, i);

		// Note: Approx based on zero initial cash
		// Compute 'J*' (objective) values corresponding to initial 'u' values

		System.out.println("+/- Cost function values for step " + k + ", state " + i + ", sign " + sign);
		QFunction qFunctionForTradeBoundaries;
		if (sign < 0) {
			qFunctionForTradeBoundaries = this.program.getQFunctionForLowerBoundary();
		} else {
			qFunctionForTradeBoundaries = this.program.getQFunctionForUpperBoundary();
		}
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			double u = this.uData[j]; //
			stateForGraph.stockQuantity = u;

			OneDimSolverResults results = this.program.evaluateForQFunction(stateForGraph,
					qFunctionForTradeBoundaries);

			double j_s = results.optim_f;
			double v_s = results.optim_v;
			System.out.println(j + ";" + u + ";" + j_s + ";" + v_s);
		}
	}

	public void printExactCostFunctionValuesToScreen(int k, int i) {
		this.initializeValuesForTests();
		State stateForGraph = new State();
		stateForGraph.cashAccount = 0.0;

		stateForGraph.timeStep = k;
		stateForGraph.marketVector = this.program.getTempMarketVectorDestination(k, i);
		this.program.setMeshRelatedValuesWithPrecomputedDerivativeValue(k, i);

		// Note: Approx based on zero initial cash
		// Compute 'J*' (objective) values corresponding to initial 'u' values
		System.out.println();
		double s = stateForGraph.marketVector.values[MarketVector.PRICE];
		double sigma = stateForGraph.marketVector.values[MarketVector.VOL];
		System.out.println("Exact cost function values for step " + k + ", state " + i + "; sk=" + s + "; vk=" + sigma
				+ "; hk=" + this.program.getDerivative().getDerivativeValueArray()[k][i]);
		System.out.println("index;u0;J*;v*");

		this.initializeValuesForTests();
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			double u = this.uData[j]; //
			stateForGraph.stockQuantity = u;
			OneDimSolverResults results = this.program.evaluate(stateForGraph);
			double j_s = results.optim_f;
			double v_s = results.optim_v;

			System.out.println(j + "; " + u + "; " + j_s + "; " + v_s);
		}
	}

	public void printInitialCostFunctionValuesToScreen() {
		this.printApproxCostFunctionValuesToScreen(0, 0);
	}

	public void printQFunctionValuesForIntervalToScreen(double x1, double x2, int nbPts) {
		double dx = (x2 - x1) / nbPts;
		System.out.println("*** Fct values for " + this);
		QFunction qFunction = this.program.getQFunction();
		for (int i = 0; i <= nbPts; i++) {
			double x = x1 + dx * i;
			double f = qFunction.evaluate(x);
			System.out.println(x + "; " + f);

		}
	}

	public void printQFunctionValuesForTradeBoundariesToScreen(int k, int i, double sign, double u0) {
		this.initializeValuesForTests();
		State stateForGraph = new State();
		stateForGraph.cashAccount = 0.0;
		stateForGraph.timeStep = k;
		stateForGraph.stockQuantity = u0;
		stateForGraph.marketVector = this.program.getTempMarketVectorDestination(k, i);
		this.program.setMeshRelatedValuesWithPrecomputedDerivativeValue(k, i);

		// Note: Approx based on zero initial cash
		// Compute 'J*' (objective) values corresponding to initial 'u' values

		System.out.println("+/- Q function values for step " + k + ", state " + i + ", sign " + sign + ", u0 " + u0);
		QFunction qFunctionForTradeBoundaries;
		if (sign < 0) {
			qFunctionForTradeBoundaries = this.program.getQFunctionForLowerBoundary();
		} else {
			qFunctionForTradeBoundaries = this.program.getQFunctionForUpperBoundary();
		}

		qFunctionForTradeBoundaries.setState(stateForGraph);
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			double v = this.uData[j]; //

			double q = qFunctionForTradeBoundaries.evaluate(v);

			System.out.println(j + ";" + v + ";" + q);
		}

	}

	public void printQFunctionValuesOnlyCurrentCostToScreen() {
		this.initializeValuesForTests();
		// Note: Approx based on zero initial cash
		// Compute 'J*' (objective) values corresponding to initial 'u'
		// values
		System.out.println("Q-function values (ONLY current cost)");
		QFunction qFunction = this.program.getQFunction();
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			double v = this.uData[j];
			double q_s = qFunction.evaluateOnlyCurrentCost(v);
			System.out.println(j + ";" + v + ";" + q_s);
		}
	}

	public void printQFunctionValuesToScreen() {
		QFunction qFunction = this.program.getQFunction();
		this.printQFunctionValuesToScreen(qFunction);
	}

	public void printQFunctionValuesToScreen(QFunction qFunction) {
		Derivative derivative = this.program.getDerivative();
		double u0 = derivative.getMinDeltaForHedge();
		double u1 = derivative.getMaxDeltaForHedge();
		this.printQFunctionValuesToScreen(qFunction, u0, u1);
	}

	public void printQFunctionValuesToScreen(QFunction qFunction, double u0, double u1) {
		this.initializeValuesForTests(u0, u1);
		// Note: Approx based on zero initial cash
		// Compute 'J*' (objective) values corresponding to initial 'u'
		// values
		System.out.println("Q-function values");
		System.out.println("#;v;f;g;df;h");
		double dv = 1e-6;
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			double v = this.uData[j];
			double q_s = qFunction.evaluate(v);
			double deriv = qFunction.gradient(v);
			double q2 = qFunction.evaluate(v + dv);
			double dq = (q2 - q_s) / dv;
			double h = qFunction.hessian(v);
			System.out.println(j + ";" + v + ";" + q_s + ";" + deriv + ";" + dq + ";" + h);
		}
	}

	public void printQFunctionValuesWithoutCurrentCostToScreen() {
		this.initializeValuesForTests();
		// Note: Approx based on zero initial cash
		// Compute 'J*' (objective) values corresponding to initial 'u'
		// values
		System.out.println("Q-function values (without current cost)");
		QFunction qFunction = this.program.getQFunction();
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			double v = this.uData[j];
			double q_s = qFunction.evaluateWithoutCurrentCost(v);
			System.out.println(j + ";" + v + ";" + q_s);
		}
	}

	public void printSumOfWeightsToScreen(int k, int i) {
		double sum = 0.0;
		MeshWeights weights = this.program.getStochasticMesh().getMeshWeights();
		FloatArrayList wList = weights.getWeightsList(k, i);
		for (int j = 0; j < this.program.getNbStatesPerStep(); j++) {
			sum += wList.getQuick(j);
		}
		System.out.println(k + ";" + i + ";" + sum);
	}

	public void showCostFunctionApproxCrossSectionToScreen(int k) {
		System.out.println("***** k = " + k);
		int n = this.program.getQFunction().getCurrentNbWeights();
		if (k == 0) {
			n = 1; // Special case
		}
		for (int j = 0; j < n; j++) {
			MarketVector marketVector = this.program.getTempMarketVectorDestination(k, j);
			double s = marketVector.values[MarketVector.PRICE];
			RiskFunctionApproximation approximation = this.program.jStarArray[k][j];
			double jApprox = approximation.evaluate(0, 0);
			System.out.println(j + "; " + s + "; " + jApprox);
		}
	}

	public void showCostFunctionApproxResultsToScreen(State x) {
		this.initializeValuesForTests();
		int k = x.timeStep;
		int i = x.marketVector.simulationIndex;
		System.out.println("***** " + k);

		RiskFunctionApproximation approximation = this.program.jStarArray[k][i];
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			x.stockQuantity = this.uData[j];
			System.out.println(j + "; " + x.stockQuantity + "; " + this.jData[j] + "; "
					+ approximation.evaluate(x.cashAccount, x.stockQuantity));
			this.jData[j] = this.program.evaluateAtLastStep(x);
		}
	}

	void initializeValuesForTests() {
		Derivative derivative = this.program.getDerivative();
		this.initializeValuesForTests(derivative.getMinDeltaForHedge(), derivative.getMaxDeltaForHedge());
	}

	void initializeValuesForTests(double u0, double u1) {
		double du = (u1 - u0) / (this.nbDataPointsForTesting - 1);

		this.uData = new double[this.nbDataPointsForTesting];
		this.jData = new double[this.nbDataPointsForTesting];
		for (int j = 0; j < this.nbDataPointsForTesting; j++) {
			this.uData[j] = u0 + j * du;
		}
	}

}