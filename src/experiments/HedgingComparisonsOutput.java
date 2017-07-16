package experiments;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import policies.HedgingPolicy;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.stat.Tally;

public class HedgingComparisonsOutput {

	public static final String DATE_FORMAT_NOW = "yyyy-MM-dd HH:mm:ss";

	public static void writeLineToFile(String line, String filename, boolean appendLine) {
		BufferedWriter bufferedWriter = null;
		try {
			bufferedWriter = new BufferedWriter(new FileWriter(filename, appendLine));
			bufferedWriter.write(line);
			bufferedWriter.newLine();
		} catch (Exception ex) {
			ex.printStackTrace();
		} finally {
			try {
				if (bufferedWriter != null) {
					bufferedWriter.flush();
					bufferedWriter.close();
				}
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}

	private String dataFileName;
	private boolean firstTimeWriting = false; // Always append to file
	private String hostName = "";
	private InputParameters inputParameters;
	private DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");

	public HedgingComparisonsOutput(InputParameters inputParameters, SimulationParameters simulationParameters) {
		this.inputParameters = inputParameters;
		this.dataFileName = this.getDataFileName(simulationParameters);
		this.getHostName();
	}

	public void writeFinalInformation() {
		this.writeLineToFile("*** Experiment ended; " + this.getDateTime() + ";FileName:" + this.dataFileName,
				this.dataFileName);
	}

	public void writeInitialInformation() {
		this.writeLineToFile("*** Experiment started; HostName:" + this.hostName.replace('@', ';') + ";"
				+ this.getDateTime() + ";FileName:" + this.dataFileName, this.dataFileName);

	}

	public void writeResults(SimulationParameters simParams, ProblemDefinition problem, TestedPolicies testedPolicies,
			SimulationResults results) {
		MarketVectorProcess spExpt = problem.stochprocessForExperiments;
		String output = simParams.toString();
		output += "; " + spExpt.getProcessLabel();

		if (results.generatesError) {
			output += "; " + results.errorMessage;
		} else {
			output = this.addDerivativePriceInformation(results, output);
			output = this.addFractionOfMeshPointsUsed(results, output);
			output = this.addMeshApproximationErrorEstimates(results, output);
			output = this.addMeshBasedRiskFunctionEstimates(results, output, testedPolicies);
			output = this.addOutOfSampleRiskResults(results, output, testedPolicies);
			output = this.addPolicyDecisions(results, output, testedPolicies);
			output = this.addPnLStatistics(results, output, testedPolicies);
			output = this.addComputationTimes(results, output, testedPolicies);
			output = this.addExtraInformation(simParams, problem, results, output, testedPolicies);
		}

		this.writeLineToFile(output, this.dataFileName);
	}

	public void writeResultsHeader(SimulationParameters simulationParameters, TestedPolicies testedPolicies) {

		String headerLine = SimulationParameters.headerString;

		headerLine = this.addHeaderDetails(headerLine);
		headerLine = this.addPolicySpecificHeaderDetails(headerLine, testedPolicies);
		headerLine = this.addRemainingHeaderDetails(headerLine, testedPolicies);

		this.writeLineToFile(headerLine, this.dataFileName);
	}

	public void writeThrowableInformation(Throwable e) {
		String message = e.toString();
		this.writeLineToFile("*** Exception;" + message, this.dataFileName);

		StackTraceElement[] stacktrace = e.getStackTrace();
		for (int i = 0; i < stacktrace.length; i++) {
			this.writeLineToFile("*** StackTrace[" + i + "];" + stacktrace[i].toString(), this.dataFileName);
		}
	}

	private String addComputationTimes(SimulationResults results, String output, TestedPolicies testedPolicies) {
		for (int k = 0; k < testedPolicies.getNbPolicies(); k++) {
			output += "; " + results.policyEvalTime[k] + "; " + results.dpCompTime[k];
		}
		output += "; " + results.meshCompTime[0];
		output += "; " + results.timer.getSeconds();
		output += "; " + this.getDateTime();
		return output;
	}

	private String addDerivativePriceInformation(SimulationResults results, String output) {
		Tally hInitValue = results.hInitValue;
		Tally hOuterValue = results.hOuterValue;
		output += "; " + hInitValue.average() + "; " + (hInitValue.standardDeviation());
		output += "; " + hOuterValue.average() + "; " + (hOuterValue.standardDeviation());
		return output;
	}

	private String addExtraInformation(SimulationParameters simParams, ProblemDefinition problem,
			SimulationResults results, String output, TestedPolicies testedPolicies) {
		output += "; " + simParams.usingBridgeForOuterSims;

		String cvLabel = "";
		if (simParams.usingCVForOuterSims) {
			cvLabel = "-CV";
		}
		output += "; " + (problem.lossFunction.getName() + cvLabel);
		output += "; " + simParams.maxAllowedLoss;
		output += "; " + simParams.nbPilotRunsForOuterSims;

		return output;
	}

	private String addFractionOfMeshPointsUsed(SimulationResults results, String output) {
		output += "; " + results.fracAttainablePoints.average() + "; "
				+ (results.fracAttainablePoints.standardDeviation());
		output += "; " + results.fracWeightsUsed.average() + "; " + (results.fracWeightsUsed.standardDeviation());
		return output;
	}

	private String addHeaderDetails(String headerLine) {
		headerLine += "; Process";
		headerLine += "; h-in avg; h-in sdev";
		headerLine += "; h-out avg; h-out sdev";
		headerLine += "; fracPoints avg; fracPoints sdev";
		headerLine += "; fracWeights avg; fracWeights sdev";
		headerLine += "; Err avg; Err stdv; Err sup";
		return headerLine;
	}

	private String addMeshApproximationErrorEstimates(SimulationResults results, String output) {
		output += "; " + results.approxErrorTally.average() + "; " + results.approxErrorTally.standardDeviation() + "; "
				+ results.approxErrorTally.max();
		return output;
	}

	private String addMeshBasedRiskFunctionEstimates(SimulationResults results, String output,
			TestedPolicies testedPolicies) {
		HedgingPolicy[] policyArray = testedPolicies.getPolicyArray();
		for (int k = 0; k < testedPolicies.getNbPolicies(); k++) {
			if (policyArray[k].getPolicyHedgingMethod() == HedgingPolicy.HedgingMethod.SM_LQ) {
				output += "; " + results.lbTally.average() + "; " + results.lbTally.standardDeviation() + "; "
						+ results.lbWithErrorTally.average() + "; " + results.lbWithErrorTally.standardDeviation();
			}
		}
		return output;
	}

	private String addOutOfSampleRiskResults(SimulationResults results, String output, TestedPolicies testedPolicies) {
		for (int k = 0; k < testedPolicies.getNbPolicies(); k++) {
			if (results.riskMeasureValues[k].numberObs() > 0) {
				output += "; " + results.riskMeasureValues[k].average() + "; "
						+ (results.riskMeasureValues[k].standardDeviation());
			} else {
				output += "; ; ";
			}
		}
		return output;
	}

	private String addPnLStatistics(SimulationResults results, String output, TestedPolicies testedPolicies) {
		for (int k = 0; k < testedPolicies.getNbPolicies(); k++) {
			output += "; " + results.ptfPnLTally[k].average() + "; " + results.ptfPnLTally[k].standardDeviation() + "; "
					+ results.ptfPnLVarTally[k].standardDeviation();
		}
		return output;
	}

	private String addPolicyDecisions(SimulationResults results, String output, TestedPolicies testedPolicies) {
		for (int k = 0; k < testedPolicies.getNbPolicies(); k++) {
			output += "; " + results.decisions[k].average() + "; " + results.decisions[k].standardDeviation();
		}
		return output;
	}

	private String addPolicySpecificHeaderDetails(String headerLine, TestedPolicies testedPolicies) {
		int nbPolicies = testedPolicies.getNbPolicies();
		HedgingPolicy[] policyArray = testedPolicies.getPolicyArray();
		for (int k = 0; k < nbPolicies; k++) {
			if (testedPolicies.getPolicyArray()[k].getPolicyHedgingMethod() == HedgingPolicy.HedgingMethod.SM_LQ) {
				headerLine += "; " + "LB " + policyArray[k].outputCode() + "; " + "Sdev LB "
						+ policyArray[k].outputCode() + "; " + "LBwErr " + policyArray[k].outputCode() + "; "
						+ "Sdev LBwErr " + policyArray[k].outputCode();
			}
		}

		String exptCode = "";
		for (int k = 0; k < nbPolicies; k++) {
			exptCode = policyArray[k].outputCode();
			headerLine += "; " + exptCode + " avg; " + exptCode + " sdev";
		}

		for (int k = 0; k < nbPolicies; k++) {
			headerLine += "; " + policyArray[k].outputCode() + " u1" + "; " + policyArray[k].outputCode() + " u1 sdv";
		}

		for (int k = 0; k < nbPolicies; k++) {
			headerLine += "; " + policyArray[k].outputCode() + " pnl" + "; " + policyArray[k].outputCode() + " pnl sdv"
					+ "; " + policyArray[k].outputCode() + " pnl var sdv";
		}

		for (int k = 0; k < nbPolicies; k++) {
			headerLine += "; " + policyArray[k].outputCode() + " Eval T." + "; " + policyArray[k].outputCode()
					+ " CT-DP";
		}
		return headerLine;
	}

	private String addRemainingHeaderDetails(String headerLine, TestedPolicies testedPolicies) {
		headerLine += "; " + "CT-Mesh";
		headerLine += "; " + "CT-Total";
		headerLine += "; DateStamp";
		headerLine += "; usingOuterBB";
		headerLine += "; lossFunction";
		headerLine += "; maxAllowedLoss";
		headerLine += "; nbPilotRunsForOuterCV";

		return headerLine;
	}

	private String getDataFileName(SimulationParameters simulationParameters) {
		String dataFileName = this.inputParameters.fileName + " " + simulationParameters.derivative.getDescription()
				+ " " + simulationParameters.priceModel + " K" + this.inputParameters.k0 + " N"
				+ simulationParameters.nbReps + "n" + simulationParameters.nbOuterSims + " "
				+ simulationParameters.meshMethod + " " + Math.round(simulationParameters.t * 12) + "m" + " v"
				+ InputParameters.getRangeString(this.inputParameters.sigma0, this.inputParameters.sigma1, 100) + " s"
				+ Math.round(this.inputParameters.s0) + " u" + Math.round(100 * this.inputParameters.u0) + " b"
				+ Math.round(100 * this.inputParameters.bTC0) + ".csv";
		return dataFileName;
	}

	private String getDateTime() {
		Date date = new Date();
		return this.dateFormat.format(date);
	}

	private void getHostName() {
		try {
			this.hostName = ManagementFactory.getRuntimeMXBean().getName();
		} catch (Exception e) {
			this.hostName = "Unkown Host";
		}
	}

	private void writeLineToFile(String line, String filename) {
		HedgingComparisonsOutput.writeLineToFile(line, filename, !this.firstTimeWriting);
		this.firstTimeWriting = false;
	}

}
