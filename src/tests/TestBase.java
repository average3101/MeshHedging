package tests;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;

public class TestBase {

	private double absoluteTolerance = 0.05;

	private int nbFails = 0;
	private int nbTests = 0;
	private double relativeTolerance = 0.05;
	private String testLabel;
	private Chrono timer = new Chrono();

	public static TestBase GetInstance() {
		return new TestBase();
	}

	public boolean assertEqualAbsolute(double observed, double expected, String source) {
		double absoluteError = Math.abs((observed - expected));
		boolean isEqual = (absoluteError < this.getAbsoluteTolerance());

		String absoluteErrorString = String.format("%s", absoluteError);
		String failureMessage = "  Observed = " + observed + ", Expected = " + expected + ", Absolute Error = "
				+ absoluteErrorString;
		return this.assertTrue(isEqual, failureMessage, source);
	}

	public boolean assertEqualRelative(double observed, double expected, String source) {
		double relativeError = Math.abs((observed - expected) / expected);
		boolean isEqual = (relativeError < this.getRelativeTolerance());
		int relativeErrorInPercent = (int) (100 * relativeError);
		String failureMessage = "  Observed = " + observed + ", Expected = " + expected + ", Relative Error = "
				+ relativeErrorInPercent + "%";
		return this.assertTrue(isEqual, failureMessage, source);
	}

	public boolean assertEqualRelative(Tally obsTally, double expectation, String source) {
		double mean = obsTally.average();
		boolean passedTest = this.assertEqualRelative(mean, expectation, source);
		if (!passedTest) {
			System.out.println(obsTally.reportAndCIStudent(0.95));
		}

		return passedTest;
	}

	public boolean assertTrue(boolean isTrue, String failureMessage, String source) {
		this.nbTests++;
		String message = null;
		if (isTrue) {
			message = ".";
		} else {
			this.nbFails++;
			message = "\nFailed test : " + source;
			message += "\n " + failureMessage + "\n";
		}
		System.out.print(message);
		return isTrue;
	}

	public void defineTestLabel() {
		this.setTestLabel("*undefined*");
	}

	public void execute() {
		this.startTimer();
		this.defineTestLabel();
		this.printInitialMessage();

		this.run();

		this.printFailSummary();
		this.printFinalMessages();
	}

	public double getRelativeTolerance() {
		return this.relativeTolerance;
	}

	public void run() {
	}

	public void setRelativeTolerance(double relativeTolerance) {
		this.relativeTolerance = relativeTolerance;
	}

	public void setTestLabel(String s) {
		this.testLabel = s;
	}

	private double getAbsoluteTolerance() {

		return this.absoluteTolerance;
	}

	private void printFailSummary() {
		if (this.nbFails > 0) {
			System.out.println("\nNumber of failed tests = " + this.nbFails + " (" + this.nbTests + " tests).");
		} else {
			System.out.println("\nNo failed tests (" + this.nbTests + " tests).");
		}
	}

	private void printFinalMessages() {
		System.out.println("End of tests for : " + this.testLabel);

		// double h = Math.round(this.timer.getHours());
		// double m = Math.round(this.timer.getMinutes());
		double s = Math.round(this.timer.getSeconds());
		String timerString = s + "s"; // h + "h" + m + "m" + s + "s";
		System.out.println("Total execution time : " + timerString);
	}

	private void printInitialMessage() {
		System.out.println("Start of tests for : " + this.testLabel);
	}

	private void startTimer() {
		this.timer.init();
	}

}
