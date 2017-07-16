package experiments;

import java.util.ArrayList;
import java.util.List;

import derivatives.CallOption;
import derivatives.Derivative;
import derivatives.DerivativeBasket;
import derivatives.DerivativePriceControlVariate;
import derivatives.DerivativePricedOnMesh;
import derivatives.DerivativePricedOnMeshWithCVAnalytic;
import derivatives.DerivativePricedOnMeshWithCVEstimated;
import derivatives.PutOption;
import derivatives.ShortCallOption;
import derivatives.ShortPutOption;
import experiments.ProblemDefinition.MeshMethod;
import experiments.ProblemDefinition.PriceModel;
import numerics.NormalGenerationMethod;

public class SimulationParameters {

	static final String headerString = "H Desc.; s0; strike; u0; sigma0; sigma1; t; maturity; nbObs; nbInnerSims; nbOuterSims; nbReps; aTC; bTC; correl; Rever; VoLV; N Gen; riskAversion; nbPtsForError; MeshMethod; RouletteCutoff; IndicLoading1; IndicLoading2";
	static NormalGenerationMethod[] normalGenMethodArray = { NormalGenerationMethod.chol, NormalGenerationMethod.pca };

	public double aTC;
	public double bTC;
	public Derivative derivative;
	public double indicatorLoading1;
	public double indicatorLoading2;
	public double maturity;
	public double maxAllowedLoss;
	public double meanReversionRate;
	public MeshMethod meshMethod;
	public String meshMethodString;
	public int nbInnerSims;
	public int nbObs;
	public int nbOuterSims;
	public int nbPathsForInitialPrice = 100000;
	public int nbPtsForErrorCalc;
	public int nbReps;
	public NormalGenerationMethod normalGenMethod;
	public PriceModel priceModel;
	public double priceVolCorrelation;
	public double rfRate;
	public double riskAversion;
	public double s0;
	public double sigma0;
	public double sigma1;
	public double strike;
	public double t;
	public double u0;
	public boolean usingCVForOuterSims = false;
	public boolean usingStochMeshBasedPolicies = false;
	public double volOfLogVar;
	public double weightCutoff;
	public int nbPilotRunsForOuterSims;
	public boolean usingBridgeForOuterSims;

	public SimulationParameters() {
	}

	public SimulationParameters(MeshMethod meshMethod, double s0, double strike, double u0, double sigma0,
			double sigma1, double t, double maturity, int nbObs, int nbInnerSims, int nbOuterSims, int nbReps,
			double aTC, double bTC, double correl, double meanReversionRate, double volOfLogVar,
			NormalGenerationMethod normalGenMethod, double weightCutoff, double riskAversion, int nbPtsForErrorCalc,
			Derivative derivative, PriceModel priceModel, boolean usingCV1, boolean usingOuterCV, double maxAllowedLoss,
			double rfRate, double indicatorLoading1, double indicatorLoading2, int nbPilotRuns) {
		this.meshMethod = meshMethod;
		this.s0 = s0;
		this.strike = strike;
		this.u0 = u0;
		this.sigma0 = sigma0;
		this.sigma1 = sigma1;
		this.t = t;
		this.maturity = maturity;
		this.nbObs = nbObs;
		this.nbInnerSims = nbInnerSims;
		this.nbOuterSims = nbOuterSims;
		this.nbReps = nbReps;
		this.aTC = aTC;
		this.bTC = bTC;
		this.priceVolCorrelation = correl;
		this.meanReversionRate = meanReversionRate;
		this.volOfLogVar = volOfLogVar;
		this.normalGenMethod = normalGenMethod;
		this.weightCutoff = weightCutoff;
		this.riskAversion = riskAversion;
		this.nbPtsForErrorCalc = nbPtsForErrorCalc;
		this.derivative = derivative;
		this.priceModel = priceModel;
		// this.usingCV1 = usingCV1;
		this.usingCVForOuterSims = usingOuterCV;
		this.maxAllowedLoss = maxAllowedLoss;
		this.indicatorLoading1 = indicatorLoading1;
		this.indicatorLoading2 = indicatorLoading2;
		this.meshMethodString = meshMethod.toString();
		this.rfRate = rfRate;
		this.nbPilotRunsForOuterSims = nbPilotRuns;
	}

	public static List<SimulationParameters> initializeSimulationParametersList(InputParameters inputs) {
		boolean usingCV1 = false;
		boolean usingCV2 = false;
		if (inputs.cvFlag1.equals("Y")) {
			usingCV1 = true;
		}
		if (inputs.cvFlag2.equals("Y")) {
			usingCV2 = true;
		}

		List<SimulationParameters> simParamList = new ArrayList<SimulationParameters>();
		for (int meshSize = inputs.minMeshSize; meshSize <= inputs.maxMeshSize; meshSize *= 2) {
			int nbObs = inputs.k0;
			Derivative h = createDerivative(inputs);
			double aTC = 0;
			double bTC = inputs.bTC0;

			int cc = inputs.cc0;
			double weightCutoffFactor = Math.pow(10, -1 * cc);
			if (cc < 0) {
				weightCutoffFactor = 0;
			}

			MeshMethod meshMethod = MeshMethod.values()[inputs.meshMethodNo];
			PriceModel priceModel = PriceModel.values()[inputs.modelNo];

			simParamList.add(new SimulationParameters(meshMethod, inputs.s0, h.getStrike(), inputs.u0, inputs.sigma0,
					inputs.sigma1, inputs.horizon, inputs.maturity, nbObs, meshSize, inputs.nbOuterSims, inputs.nbReps,
					aTC, bTC, inputs.correl, inputs.meanReversionRate, inputs.volOfLogVar, normalGenMethodArray[0],
					weightCutoffFactor, inputs.RA0, inputs.nbPtsForErrorCalc, h, priceModel, usingCV1, usingCV2,
					inputs.maxAllowedLoss, inputs.rfRate, inputs.indicatorLoading0, inputs.indicatorLoading1,
					inputs.nbPilotRuns));
		}
		return simParamList;
	}

	private static Derivative createBasicDerivative(InputParameters inputs) {
		Derivative deriv = null;
		if (inputs.derivType.equals("CALL")) {
			deriv = new CallOption(inputs.strike, inputs.maturity, inputs.rfRate, inputs.divYield);
		} else if (inputs.derivType.equals("SCALL")) {
			deriv = new ShortCallOption(inputs.strike, inputs.maturity, inputs.rfRate, inputs.divYield);
		} else if (inputs.derivType.equals("PUT")) {
			deriv = new PutOption(inputs.strike, inputs.maturity, inputs.rfRate, inputs.divYield);
		} else if (inputs.derivType.equals("SPUT")) {
			deriv = new ShortPutOption(inputs.strike, inputs.maturity, inputs.rfRate, inputs.divYield);
		} else if (inputs.derivType.equals("CS")) {
			deriv = createCallSpread(inputs);
		} else if (inputs.derivType.equals("SS")) {
			deriv = createShortStrangle(inputs);
		}
		return deriv;
	}

	private static Derivative createCallSpread(InputParameters inputs) {
		Derivative deriv;
		double strike1 = 9;
		double strike2 = 11;

		Derivative[] basket = new Derivative[2];
		basket[0] = new CallOption(strike1, inputs.maturity, inputs.rfRate, inputs.divYield);
		basket[1] = new CallOption(strike2, inputs.maturity, inputs.rfRate, inputs.divYield);

		double[] positions = { -1, 1 };

		deriv = new DerivativeBasket("CallSpread", basket, positions, 0.0, 1.0);
		return deriv;
	}

	private static Derivative createDerivative(InputParameters inputs) {
		Derivative derivative = null;
		Derivative basicDerivative = createBasicDerivative(inputs);
		if (inputs.derivEvalByStochMesh.equals("Y")) {
			derivative = new DerivativePricedOnMesh(basicDerivative);
		} else if (inputs.derivEvalByStochMesh.equals("C-D")) {
			derivative = createDerivativePricedOnMesh(inputs, basicDerivative,
					DerivativePriceControlVariate.Label.DELTA);
		} else if (inputs.derivEvalByStochMesh.equals("C-CE")) {
			derivative = createDerivativePricedOnMesh(inputs, basicDerivative,
					DerivativePriceControlVariate.Label.COND_EXPECT);
		} else if (inputs.derivEvalByStochMesh.equals("N")) {
			derivative = basicDerivative;
		} else {
			derivative = null;
			System.out.println("Cannot interpret input for 'derivEvalByStochMesh' : " + inputs.derivEvalByStochMesh);
		}
		return derivative;
	}

	private static Derivative createDerivativePricedOnMesh(InputParameters inputs, Derivative basicDerivative,
			DerivativePriceControlVariate.Label label) {
		Derivative derivative = null;
		if (inputs.cvFlag1.equals("AN")) {
			derivative = new DerivativePricedOnMeshWithCVAnalytic(basicDerivative, label);
		} else if (inputs.cvFlag1.equals("ES")) {
			derivative = new DerivativePricedOnMeshWithCVEstimated(basicDerivative, label);
		} else {
			System.out.println("Cannot interpret input for 'cvFlag1' : " + inputs.cvFlag1);
		}
		return derivative;
	}

	private static Derivative createShortStrangle(InputParameters inputs) {
		Derivative deriv;
		double strike1 = 9;
		double strike2 = 11;

		Derivative[] basket = new Derivative[2];
		basket[0] = new PutOption(strike1, inputs.maturity, inputs.rfRate, inputs.divYield);
		basket[1] = new CallOption(strike2, inputs.maturity, inputs.rfRate, inputs.divYield);

		double[] positions = { -1, -1 };

		deriv = new DerivativeBasket("ShortStrangle", basket, positions, -1.0, 0.0);
		return deriv;
	}

	@Override
	public String toString() {
		String a = this.derivative.getDescription() + "; " + this.s0 + "; " + this.strike + "; " + this.u0 + "; "
				+ this.sigma0 + "; " + this.sigma1 + "; " + this.t + "; " + this.maturity + "; " + this.nbObs + "; "
				+ this.nbInnerSims + "; " + this.nbOuterSims + "; " + this.nbReps + "; " + this.aTC + "; " + this.bTC
				+ "; " + this.priceVolCorrelation + "; " + this.meanReversionRate + "; " + this.volOfLogVar + "; "
				+ this.normalGenMethod.getName() + "; " + this.riskAversion + "; " + this.nbPtsForErrorCalc + "; "
				+ this.meshMethodString + ";" + this.weightCutoff + ";" + this.indicatorLoading1 + ";"
				+ this.indicatorLoading2;
		return a;
	}

}