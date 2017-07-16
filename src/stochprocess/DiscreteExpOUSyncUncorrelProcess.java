package stochprocess;

import umontreal.ssj.stochprocess.MultivariateBrownianMotion;

public class DiscreteExpOUSyncUncorrelProcess extends DiscreteExpOULagProcess {

	public DiscreteExpOUSyncUncorrelProcess(float[] y0Array, double[] mu, double[] aVec, double[] bVec, double[] sigma,
			MultivariateBrownianMotion multiBM) {
		super(y0Array, mu, aVec, bVec, sigma, 0, multiBM);
	}

	@Override
	public Object clone() {
		DiscreteExpOUSyncUncorrelProcess clone = new DiscreteExpOUSyncUncorrelProcess(this.y0Array, this.mu, this.aVec, this.bVec,
				this.sigVec, this.multiDimensionalBrownianMotion);
		return clone;
	}

	@Override
	public double getVarianceAtStep(int i) {
		// TODO
		return Double.NaN;
	}

	@Override
	protected double getReturnVolatilityForDensity(double laggedVolatility, double synchedVolatility) {
		return synchedVolatility;
	}

	@Override
	protected double getReturnVolatilityForPaths(int i) {
		return this.volatilityPath[i];
	}

	@Override
	protected void setParameters(double[] aVec, double[] bVec, double[] sigVec, double priceVolCorrelation) {
		super.setParameters(aVec, bVec, sigVec, priceVolCorrelation);
		this.processLabel = "expOU-sync";
	}
}
