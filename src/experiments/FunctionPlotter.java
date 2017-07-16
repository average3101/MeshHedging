package experiments;

import derivatives.Derivative;
import dynprog.DynamicHedgingProgram;
import dynprog.DynamicHedgingProgramDiagnostics;
import policies.RiskFunctionBasedPolicy;
import stochprocess.MarketVector;

public class FunctionPlotter {

	public static void printDerivativeValueAndPayoffToScreen(Derivative derivative, double s_min, double s_max,
			int nbPoints, MarketVector y0) {
		MarketVector y_clone = y0.clone();
		float[] values = y_clone.getValues();
		float sigma = values[MarketVector.VOL];
		double ds = (s_max - s_min) / nbPoints;

		System.out.println("S;InitValue;Payoff;Delta");
		for (int i = 0; i <= nbPoints; i++) {
			double s = s_min + i * ds;
			values[MarketVector.PRICE] = (float) s;
			double h0 = derivative.evaluateAtT0(y_clone);
			double hK = derivative.payoff((float) s);

			double delta = derivative.delta(s, sigma, 0);

			System.out.println(s + ";" + h0 + ";" + hK + ";" + delta);
		}
	}

	public static void printLinearizedQFunctions(int k, RiskFunctionBasedPolicy costFunctionBasedPolicy) {
		int nbDataPointsForTesting = 10;
		DynamicHedgingProgram cf = costFunctionBasedPolicy.getDynamicHedgingProgram();
		DynamicHedgingProgramDiagnostics diagnostics = new DynamicHedgingProgramDiagnostics(cf, nbDataPointsForTesting);

		diagnostics.printExactCostFunctionValuesForTradeBoundariesToScreen(k, 0, 1);
		diagnostics.printExactCostFunctionValuesForTradeBoundariesToScreen(k, 0, -1);

		diagnostics.printQFunctionValuesForTradeBoundariesToScreen(k, 0, 1, -1.0);
		diagnostics.printQFunctionValuesForTradeBoundariesToScreen(k, 0, -1, -1.0);
		diagnostics.printQFunctionValuesForTradeBoundariesToScreen(k, 0, 1, -.5);
		diagnostics.printQFunctionValuesForTradeBoundariesToScreen(k, 0, -1, -.5);
		diagnostics.printQFunctionValuesForTradeBoundariesToScreen(k, 0, 1, 0.0);
		diagnostics.printQFunctionValuesForTradeBoundariesToScreen(k, 0, -1, 0.0);
	}

	public static void printOptimalDecisionForQFunction(int k, RiskFunctionBasedPolicy riskFunctionBasedPolicy,
			int nbDataPointsForTesting) {
		DynamicHedgingProgram cf = riskFunctionBasedPolicy.getDynamicHedgingProgram();
		DynamicHedgingProgramDiagnostics diagnostics = new DynamicHedgingProgramDiagnostics(cf, nbDataPointsForTesting);

		diagnostics.printExactCostFunctionValuesToScreen(k, 0);
	}

	public static void printPolicyDecisionsToScreen(TestedPolicies testPolicies) {
		String output = "";
		for (int policyNo = 0; policyNo < testPolicies.getNbPolicies(); policyNo++) {
			output += testPolicies.getPolicyArray()[policyNo].getPolicyHedgingMethod() + ";";
		}
		System.out.println(output);

		for (int j = 0; j < testPolicies.getNbDecisionPts(); j++) {
			output = "";
			for (int policyNo = 0; policyNo < testPolicies.getNbPolicies(); policyNo++) {
				output += testPolicies.getPolicyDecisionsArray()[policyNo][j] + ";";
			}
			System.out.println(output);
		}
	}

	public static void printSamplePathsToFile(SimulationParameters simParams, ProblemDefinition problem, int k) {
		problem.hedgedPtf.getLastSimulatedPath();
		problem.hedgedPtf.getLastDecisionPath();

		// 1) Print Sk path
		String outputString = "";
		outputString = "*** Print Sk path ***";
		HedgingComparisonsOutput.writeLineToFile(outputString, "SamplePath.txt", false);

		for (int kk = 0; kk <= simParams.nbObs; kk++) {
			outputString = kk + ";" + problem.mktInfoPath[kk].getValues()[0] + "; "
					+ Math.sqrt(problem.mktInfoPath[kk].getValues()[1]) + "; ";
			HedgingComparisonsOutput.writeLineToFile(outputString, "SamplePath.txt", true);
		}
	}

	public static void printSamplePathsToScreen(SimulationParameters simParams, ProblemDefinition problemDefinition,
			TestedPolicies testPolicies, int k) {
		double[] ptfPath = problemDefinition.hedgedPtf.getLastSimulatedPath();
		double[] decisionPath = problemDefinition.hedgedPtf.getLastDecisionPath();

		// 1) Print Sk path
		System.out.println("*** Print Sk path ***");
		String outputCode = testPolicies.getPolicyArray()[k].outputCode();
		String outputString = outputCode + " " + k + " ; ";
		for (int kk = 0; kk <= simParams.nbObs; kk++) {
			outputString += problemDefinition.mktInfoPath[kk].getValues()[0] + "; ";
		}
		System.out.println(outputString);

		// Vol, if available
		if (problemDefinition.mktInfoPath[0].getValues().length > 1) {
			outputString = outputCode + " vol " + k + " ; ";
			for (int kk = 0; kk <= simParams.nbObs; kk++) {
				outputString += Math.sqrt(problemDefinition.mktInfoPath[kk].getValues()[1]) + "; ";
			}
			System.out.println(outputString);
		}
		// 2) Print Ptf return path
		outputString = "" + outputCode + " ptf ; ";
		for (int kk = 0; kk <= simParams.nbObs; kk++) {
			// outputString += +((ptfPath[kk]-ptfPath[0])/ptfPath[0])+"; ";
			outputString += +((ptfPath[kk] - ptfPath[0])) + "; ";
		}
		System.out.println(outputString);
		// 3) Print Decision Path
		outputString = outputCode + " u";
		for (int kk = 0; kk <= simParams.nbObs; kk++) {
			outputString += // +ptfPath[kk]+";";
					"; " + decisionPath[kk];
		}
		System.out.println(outputString);
	}

	double[] x;

	double[][] y;

}
