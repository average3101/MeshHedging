package tests;

import derivatives.BinaryCallOption;
import derivatives.CallOption;
import derivatives.Derivative;
import derivatives.DerivativePriceControlVariate;
import derivatives.DerivativePricedOnMesh;
import derivatives.DerivativePricedOnMeshWithCVAnalytic;
import derivatives.PutOption;
import meshmethods.StochasticMesh;
import stochprocess.DiscreteGBMProcess;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.stat.Tally;

public class TestDerivatives extends TestBase {

	public static TestBase GetInstance() {
		return new TestDerivatives();
	}

	public static void main(String[] args) {
		TestBase tests = GetInstance();
		tests.execute();
	}

	private double divYield = 0;
	private double maturity;
	private int nbPaths = 2000000;
	private double[] obsTimes;

	private double rfRate = 0;

	private double strike;

	@Override
	public void defineTestLabel() {
		this.setTestLabel("derivatives");
	}

	public Derivative getBinaryCallOption(MarketVectorProcess process) {
		this.getOptionParametersFromMarketProcess(process);

		Derivative binaryCallOption = new BinaryCallOption(this.strike, this.maturity, this.rfRate, this.divYield);
		binaryCallOption.setObsTimes(this.obsTimes);

		return binaryCallOption;
	}

	public Derivative getCallOption(MarketVectorProcess process) {
		this.getOptionParametersFromMarketProcess(process);

		Derivative callOption = new CallOption(this.strike, this.maturity, this.rfRate, this.divYield);
		callOption.setObsTimes(this.obsTimes);

		return callOption;
	}

	public DerivativePricedOnMesh getCallOptionEvaluatedByStochasticMesh(MarketVectorProcess process,
			StochasticMesh stochasticMesh, boolean usingControlVariate) {
		this.getOptionParametersFromMarketProcess(process);

		Derivative callOption = new CallOption(this.strike, this.maturity, this.rfRate, this.divYield);
		DerivativePricedOnMesh callOptionEvaluatedByStochasticMesh;
		if (usingControlVariate) {
			callOptionEvaluatedByStochasticMesh = new DerivativePricedOnMeshWithCVAnalytic(callOption,
					DerivativePriceControlVariate.Label.COND_EXPECT);
		} else {
			callOptionEvaluatedByStochasticMesh = new DerivativePricedOnMesh(callOption);
		}

		callOptionEvaluatedByStochasticMesh.setStochasticMesh(stochasticMesh);
		callOptionEvaluatedByStochasticMesh.setObsTimes(this.obsTimes);

		return callOptionEvaluatedByStochasticMesh;
	}

	public Derivative getPutOption(MarketVectorProcess process) {
		this.getOptionParametersFromMarketProcess(process);

		Derivative putOption = new PutOption(this.strike, this.maturity, this.rfRate, this.divYield);
		putOption.setObsTimes(this.obsTimes);

		return putOption;
	}

	@Override
	public void run() {
		int nbSteps = 4;
		TestStochasticProcesses testProcesses = new TestStochasticProcesses();
		DiscreteGBMProcess discreteBSMProcess = testProcesses.getDiscreteGBMProcess();
		testProcesses.setProcessNbSteps(nbSteps, discreteBSMProcess);

		Derivative callOption = this.getCallOption(discreteBSMProcess);
		this.testOptionPricing(discreteBSMProcess, callOption);

		Derivative putOption = this.getPutOption(discreteBSMProcess);
		this.testOptionPricing(discreteBSMProcess, putOption);

		Derivative binaryCallOption = this.getBinaryCallOption(discreteBSMProcess);
		this.testOptionPricing(discreteBSMProcess, binaryCallOption);

		this.testOptionGammaFiniteDifference(discreteBSMProcess, callOption);
	}

	private void getOptionParametersFromMarketProcess(MarketVectorProcess process) {
		MarketVector y = process.getInitialMarketVector();
		this.strike = y.values[MarketVector.PRICE];
		this.obsTimes = process.getObservationTimes();
		this.maturity = this.obsTimes[this.obsTimes.length - 1];
	}

	private boolean testOptionGammaFiniteDifference(MarketVectorProcess process, Derivative derivative) {
		String derivativeType = derivative.getClass().toString();
		String source = derivativeType + " gammaFD";

		MarketVector y0 = process.getInitialMarketVector();
		int initialIndex = 0;
		double gamma = derivative.gamma(initialIndex, y0.values);
		double gammaFD = derivative.gammaFD(initialIndex, y0.values);

		boolean result = this.assertEqualRelative(gammaFD, gamma, source);

		return result;
	}

	private boolean testOptionPricing(MarketVectorProcess process, Derivative derivative) {
		Tally tally = new Tally();
		process.getNbObservationTimes();
		for (int i = 0; i < this.nbPaths; i++) {
			float[][] basicPaths = process.generateBasicProcessesPaths();
			float[] sPath = basicPaths[0];
			double payoff = derivative.payoff(sPath);
			tally.add(payoff);
		}

		String derivativeType = derivative.getClass().toString();
		String source = derivativeType + " average";

		double sampleAverage = tally.average();

		MarketVector y0 = process.getInitialMarketVector();
		double expectation = derivative.evaluateAtT0(y0);

		boolean resultExpectation = this.assertEqualRelative(sampleAverage, expectation, source);

		return resultExpectation;
	}

}
