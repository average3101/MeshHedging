package experiments;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class InputParameters {

	double bTC0;
	int cc0;
	double correl;
	String cvFlag1;
	String cvFlag2;
	String derivEvalByStochMesh;
	String derivType;
	double divYield;
	String fileName;
	double horizon;
	int k0;
	double maturity;
	double maxAllowedLoss;
	int maxMeshSize;
	double meanReversionRate;
	int meshMethodNo;
	int minMeshSize;
	int modelNo;
	int nbOuterSims;
	int nbPtsForErrorCalc;
	int nbReps;
	double RA0;
	double rfRate;
	double s0;
	double sigma0;
	double sigma1;
	double strike;
	double u0;
	double volOfLogVar;
	double indicatorLoading0;
	double indicatorLoading1;
	int nbPilotRuns;

	public InputParameters(String[] args) {
		this.fileName = args[0];
		this.sigma0 = Double.valueOf(args[1]);
		this.sigma1 = Double.valueOf(args[2]);
		this.horizon = Double.valueOf(args[3]);
		this.maturity = Double.valueOf(args[4]);
		this.strike = Double.valueOf(args[5]);
		this.correl = Double.valueOf(args[6]);
		this.minMeshSize = Integer.valueOf(args[7]);
		this.maxMeshSize = Integer.valueOf(args[8]);
		this.meshMethodNo = Integer.valueOf(args[9]);
		this.modelNo = Integer.valueOf(args[10]);
		this.cc0 = Integer.valueOf(args[11]);
		this.nbReps = Integer.valueOf(args[12]);
		this.nbOuterSims = Integer.valueOf(args[13]);
		this.nbPtsForErrorCalc = Integer.valueOf(args[14]);
		this.k0 = Integer.valueOf(args[15]);
		this.s0 = Double.valueOf(args[16]);
		this.u0 = Double.valueOf(args[17]);
		this.volOfLogVar = Double.valueOf(args[18]);
		this.bTC0 = Double.valueOf(args[19]);
		this.derivType = args[20];
		this.derivEvalByStochMesh = args[21];
		this.cvFlag1 = args[22]; // AN or ES for Derivative price CV
		this.cvFlag2 = args[23]; // Y or N for Outer Sim Risk CV
		this.meanReversionRate = Double.valueOf(args[24]);
		this.RA0 = Double.valueOf(args[25]);
		this.maxAllowedLoss = Double.valueOf(args[26]);
		this.indicatorLoading0 = Double.valueOf(args[27]);
		this.indicatorLoading1 = Double.valueOf(args[28]);
		this.nbPilotRuns = Integer.valueOf(args[29]);

		this.initializePreDefinedParameters();
	}

	public static int getNumberOfParameterSettingsFromFile(String filename, String fileSeparator) {
		int nbSettings = 0;
		List<String> lines = getLinesFromFile(filename);

		String[] parsedLine = lines.get(0).split(fileSeparator);
		nbSettings = parsedLine.length - 1;
		return nbSettings;
	}

	public static String getRangeString(double x0, double x1, double multiplier) {
		String rangeString = Math.round(x0 * multiplier) + "-" + Math.round(x1 * multiplier);
		return rangeString;
	}

	public static InputParameters loadInputParametersFromFile(String filename, String fileSeparator, int columnIndex) {
		List<String> lines = getLinesFromFile(filename);

		int nbLines = lines.size();
		String[] args = new String[nbLines - 1];
		for (int i = 1; i < nbLines; i++) {
			String[] parsedLine = lines.get(i).split(fileSeparator);
			args[i - 1] = parsedLine[columnIndex];
		}

		return new InputParameters(args);
	}

	private static List<String> getLinesFromFile(String filename) {
		Path path = Paths.get(filename);
		List<String> lines = null;

		try {
			lines = Files.readAllLines(path);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return lines;
	}

	public void initializePreDefinedParameters() {
		this.rfRate = 0.00;
		this.divYield = 0.0;
	}

}
