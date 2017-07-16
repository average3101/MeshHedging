package meshmethods;

import stochprocess.MarketVectorProcess;

public class StochasticMeshAverageDensity extends StochasticMesh {

	public StochasticMeshAverageDensity(MarketVectorProcess multiSp) {
		super(multiSp);
	}

	@Override
	public MeshWeights newDefaultMeshWeights() {
		MeshWeights weights = new MeshWeightsAverageDensity();
		return weights;
	}

}
