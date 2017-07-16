package meshmethods;

import cern.colt.list.FloatArrayList;
import cern.colt.list.IntArrayList;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;

public abstract class MeshPointEstimator {

	protected int currentNbWeights;
	protected int[] currentWeightIndices;
	protected float[] currentWeights;
	protected int currentNbWeightsForPricing;
	private int[] currentWeightIndicesForPricing;
	private float[] currentWeightsForPricing;
	private StochasticMesh stochasticMesh;
	private MarketVector tempMarketVectorDestination;
	private MarketVector tempMarketVectorOrigin;

	public int getCurrentNbWeights() {
		return this.currentNbWeights;
	}

	public int getCurrentNbWeightsForPricing() {
		return this.currentNbWeightsForPricing;
	}

	public int[] getCurrentWeightIndices() {
		return this.currentWeightIndices;
	}

	public int[] getCurrentWeightIndicesForPricing() {
		return this.currentWeightIndicesForPricing;
	}

	public float[] getCurrentWeights() {
		return this.currentWeights;
	}

	public float[] getCurrentWeightsForPricing() {
		return this.currentWeightsForPricing;
	}

	public MarketVector getNewMarketVector(int stepNo, int index) {
		MarketVectorProcess process = this.stochasticMesh.marketProcess;

		float[] y0 = this.stochasticMesh.getMarketVectorValues(stepNo, index);
		MarketVector newMarketVector = process.newMarketVector(y0);

		return newMarketVector;
	}

	public StochasticMesh getStochasticMesh() {
		return this.stochasticMesh;
	}

	public MarketVector getTempMarketVectorDestination() {
		return this.tempMarketVectorDestination;
	}

	public MarketVector getTempMarketVectorDestination(int stepNo, int index) {
		float[] y0 = this.stochasticMesh.getMarketVectorValues(stepNo, index);
		this.tempMarketVectorDestination.copyValues(index, y0);

		return this.tempMarketVectorDestination;
	}

	public MarketVector getTempMarketVectorOrigin() {
		return this.tempMarketVectorOrigin;
	}

	public MarketVector getTempMarketVectorOrigin(int stepNo, int index) {
		float[] y0 = this.stochasticMesh.getMarketVectorValues(stepNo, index);
		this.tempMarketVectorOrigin.copyValues(index, y0);

		return this.tempMarketVectorOrigin;
	}

	public void setCurrentNbWeights(int currentNbWeights) {
		this.currentNbWeights = currentNbWeights;
	}

	public void setCurrentRiskFunctionWeights(float[] currentWeights, int[] currentWeightIndices,
			int currentNbWeights) {
		this.currentNbWeights = currentNbWeights;
		this.currentWeights = currentWeights;
		this.currentWeightIndices = currentWeightIndices;
	}

	public void setCurrentWeightIndices(int[] currentWeightIndices) {
		this.currentWeightIndices = currentWeightIndices;
	}

	public void setCurrentWeights(float[] currentWeights) {
		this.currentWeights = currentWeights;
	}

	public void setCurrentWeights(int k, int i) {
		MeshWeights meshWeights = this.stochasticMesh.getMeshWeights();
		FloatArrayList availableWeights = meshWeights.getWeightsList(k, i);
		IntArrayList availableIndices = meshWeights.getWeightsIndexList(k, i);

		this.currentNbWeights = meshWeights.getNbWeightsPerList(k, i);
		for (int ii = 0; ii < this.currentNbWeights; ii++) {
			int j = availableIndices.getQuick(ii);
			this.currentWeightIndices[ii] = j;
			this.currentWeights[ii] = availableWeights.getQuick(ii);
		}

		if (this.stochasticMesh.isUsingDifferentProcessForPricing()) {
			MeshWeights meshWeightsForPricing = this.stochasticMesh.getMeshWeightsForPricing();
			FloatArrayList availableWeightsForPricing = meshWeightsForPricing.getWeightsList(k, i);
			IntArrayList availableIndicesForPricing = meshWeightsForPricing.getWeightsIndexList(k, i);

			this.currentNbWeightsForPricing = meshWeightsForPricing.getNbWeightsPerList(k, i);
			for (int ii = 0; ii < this.currentNbWeightsForPricing; ii++) {
				int j = availableIndicesForPricing.getQuick(ii);
				this.getCurrentWeightIndicesForPricing()[ii] = j;
				this.getCurrentWeightsForPricing()[ii] = availableWeightsForPricing.getQuick(ii);
			}
		} else {
			this.currentNbWeightsForPricing = this.currentNbWeights;
			this.setCurrentWeightIndicesForPricing(this.currentWeightIndices);
			this.setCurrentWeightsForPricing(this.currentWeights);
		}
	}

	public void setStochasticMesh(StochasticMesh stochasticMesh) {
		this.stochasticMesh = stochasticMesh;
		MarketVectorProcess process = stochasticMesh.getMarketVectorProcess();
		this.tempMarketVectorOrigin = process.newTempMarketVector();
		this.tempMarketVectorDestination = process.newTempMarketVector();
	}

	public void setCurrentWeightsForPricing(float[] currentWeightsForPricing) {
		this.currentWeightsForPricing = currentWeightsForPricing;
	}

	public void setCurrentWeightIndicesForPricing(int[] currentWeightIndicesForPricing) {
		this.currentWeightIndicesForPricing = currentWeightIndicesForPricing;
	}

}
