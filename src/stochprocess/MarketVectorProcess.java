package stochprocess;

import umontreal.ssj.hups.PointSet;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stochprocess.MultivariateBrownianMotion;
import umontreal.ssj.stochprocess.StochasticProcess;

public abstract class MarketVectorProcess extends StochasticProcess {

	protected StochasticProcess[] basicProcesses;
	protected float[][] basicProcessesPaths;
	protected double[][] correlationMatrix;
	protected double[] densityXVector;
	protected double[] dt;
	protected MarketVector[] marketVectorPath;
	protected double[] mu;
	protected MultivariateBrownianMotion multiDimensionalBrownianMotion;
	protected int nbProcessDim;
	protected PointSet pointSet;
	protected String processLabel = "Undefined";
	protected double[] sqrtdt;
	protected boolean usingQMC = false;
	protected RandomStream uStream;
	protected RandomStream uStreamShift;
	protected float[] y0Array;

	// public float[][] applyPathTransformationForStochasticGridMethod(float[][]
	// sPathArray) {
	// return sPathArray;
	// }
	//
	// public void applySettingsForStochasticGridMethod() {
	// }

	@Override
	public Object clone() {
		return null;
	}

	public double density(float[] forDensity, float[] forDensity2, int k) {
		return this.density(forDensity, forDensity2, k, k + 1);
	}

	public double density(float[] x1, float[] x2, int t1, int t2) {
		return Double.NaN;
	}

	public float expectedPriceConditionalOnMarketVector(int step0, float[] marketVector0, int step1,
			float[] marketVector1) {
		return marketVector0[MarketVector.PRICE];
	}

	public float expectedVolatilityConditionalOnMarketVector(int k0, float[] marketVector0, int k1,
			float[] marketVector1) {
		return marketVector0[MarketVector.VOL];
	}

	public double[][] generateAntitheticPathArray() {
		return null;
	}

	public abstract float[][] generateBasicProcessesPaths();

	@Override
	public double[] generatePath() {
		return null;
	}

	public double[][] getCorrelationMatrix() {
		return null;
	}

	public int getDimension() {
		return this.nbProcessDim;
	}

	public abstract double getExpectationAtStep(int nbSteps);

	public MarketVector getInitialMarketVector() {
		return this.newMarketVector(this.y0Array);
	}

	public int getMarketVectorDimension() {
		int d = Math.max(2, this.nbProcessDim);
		return d;
	}

	public MarketVector[] getMarketVectorPath() {
		int nbVectorDim = this.getMarketVectorDimension();
		for (int i = 0; i <= this.d; i++) {
			for (int j = 0; j < nbVectorDim; j++) {
				this.marketVectorPath[i].setValue(j, this.basicProcessesPaths[j][i]);
			}
		}

		return this.marketVectorPath;
	}

	public float[][] getPathArray() {
		return this.basicProcessesPaths;
	}

	public MarketVectorProcess getProcessForGeneratingOnlyOneStep() {
		return ((MarketVectorProcess) this.clone());
	}

	public String getProcessLabel() {
		return this.processLabel;
	}

	@Override
	public RandomStream getStream() {
		return this.uStream;
	}

	public abstract double getVarianceAtStep(int nbSteps);

	public double[] getVolArray() {
		return null;
	}

	public float[] getX0Array() {
		return this.y0Array;
	}

	public double integrateDensity(float[] x1, int t1, int nbPts) {
		return Double.NaN;
	}

	public abstract MarketVector newMarketVector(float[] y);

	public MarketVectorProcess newProcessForPricing() {
		return (MarketVectorProcess) this.clone();
	}

	public MarketVector newTempMarketVector() {
		int d = this.getMarketVectorDimension();
		float[] y = new float[d];
		return this.newMarketVector(y);
	}

	public boolean requiresDifferentProcessForPricing() {
		return false;
	}

	public void resetNextSubstream() {
		this.uStream.resetNextSubstream();
	}

	public void resetStartStream() {
		if (this.usingQMC) {
			this.pointSet.randomize(this.uStreamShift);
			RandomStream qmcStream = this.pointSet.iterator();
			this.setStream(qmcStream);
		}
		this.uStream.resetStartStream();
	}

	public void resetStartSubtream() {
		this.uStream.resetStartSubstream();
	}

	public void setInitialVolatilityToLongRunAverage() {
	}

	public void setQMCPointSet(PointSet pointSet, RandomStream uStreamShift) {
		this.usingQMC = true;
		this.pointSet = pointSet;
		this.uStreamShift = uStreamShift;
	}

	@Override
	public void setStream(RandomStream stream) {
		this.uStream = stream;
		this.multiDimensionalBrownianMotion.setStream(stream);
	}

	public void setX0(float[] x0) {
		this.y0Array = x0;
	}

	public double squaredDistance(double[] forDensity, double[] forDensity2, int k, int i) {
		return 0;
	}

	@Override
	protected void init() {
		this.multiDimensionalBrownianMotion.setObservationTimes(this.t, this.d);
	}

}
