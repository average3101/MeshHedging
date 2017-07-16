package experiments;

import java.util.List;

public class HedgingComparisons {

	private static final int FILENAME = 0;
	private static final int COLUMN_INDEX = 1;
	private static final int EXCLUDE_LOCAL_HEDGING = 2;

	public static void main(String[] args) {
		String separator = ";";
		String filename = args[FILENAME];
		int columnIndex = getColumnIndexToUse(args);

		boolean excludeLocalHedging = getExcludeLocalHedgingFlag(args);

		if (columnIndex > 0) {
			runComparison(separator, filename, columnIndex, excludeLocalHedging);
		} else {
			runAllComparisons(separator, filename, excludeLocalHedging);
		}
	}

	private static int getColumnIndexToUse(String[] args) {
		int defaultIndex = -1;
		int columnIndex = defaultIndex;
		if (args.length > COLUMN_INDEX) {
			columnIndex = Integer.parseInt(args[COLUMN_INDEX]);
		}
		return columnIndex;
	}

	private static boolean getExcludeLocalHedgingFlag(String[] args) {
		boolean excludeLocalHedging = false;
		if (args.length > EXCLUDE_LOCAL_HEDGING && args[EXCLUDE_LOCAL_HEDGING].equalsIgnoreCase("nolocal")) {
			excludeLocalHedging = true;
		}

		return excludeLocalHedging;
	}

	private static void runAllComparisons(String separator, String filename, boolean excludeLocalHedging) {
		int nbSettings = InputParameters.getNumberOfParameterSettingsFromFile(filename, separator);
		for (int i = 0; i < nbSettings; i++) {
			int columnIndex = i + 1;
			runComparison(separator, filename, columnIndex, excludeLocalHedging);
		}
	}

	private static void runComparison(String separator, String filename, int columnIndex, boolean excludeLocalHedging) {
		InputParameters inputParameters = InputParameters.loadInputParametersFromFile(filename, separator, columnIndex);
		new HedgingComparisons(inputParameters, excludeLocalHedging);
	}

	HedgingComparisonsOutput output;
	ProblemDefinition problem;
	RiskMeasureEstimator riskMeasureEstimator;
	TestedPolicies testPolicies;

	public HedgingComparisons(InputParameters inputParameters, boolean excludeLocalHedging) {
		List<SimulationParameters> simParamList = SimulationParameters
				.initializeSimulationParametersList(inputParameters);

		SimulationParameters simParamExample = simParamList.get(0);
		this.defineHelperClasses(inputParameters, simParamExample, excludeLocalHedging);
		this.writeStaticInformationToOutput(simParamExample);
		try {
			for (int i = 0; i < simParamList.size(); i++) {
				SimulationParameters simulationParameters = simParamList.get(i);
				this.defineRiskMeasureEstimator(simulationParameters);
				SimulationResults results = this.riskMeasureEstimator.estimatePolicyCosts(simulationParameters);
				this.output.writeResults(simulationParameters, this.problem, this.testPolicies, results);
			}
		} catch (Throwable e) {
			this.output.writeThrowableInformation(e);
		} finally {
			this.output.writeFinalInformation();
		}
	}

	private void defineHelperClasses(InputParameters inputParameters, SimulationParameters simParamExample,
			boolean excludeLocalHedging) {
		this.problem = new ProblemDefinition();
		this.testPolicies = new TestedPolicies(excludeLocalHedging);
		this.output = new HedgingComparisonsOutput(inputParameters, simParamExample);
	}

	private void defineRiskMeasureEstimator(SimulationParameters simParams) {
		if (simParams.usingCVForOuterSims) {
			this.riskMeasureEstimator = new RiskMeasureEstimatorWithCV(this.problem, this.testPolicies,
					simParams.nbPilotRunsForOuterSims);
		} else {
			this.riskMeasureEstimator = new RiskMeasureEstimator(this.problem, this.testPolicies);
		}
	}

	private void writeStaticInformationToOutput(SimulationParameters simParamExample) {
		this.output.writeInitialInformation();
		this.output.writeResultsHeader(simParamExample, this.testPolicies);
	}
}