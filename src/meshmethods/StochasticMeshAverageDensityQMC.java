package meshmethods;

import stochprocess.MarketVectorProcess;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;

public class StochasticMeshAverageDensityQMC extends StochasticMeshAverageDensity {

	RandomStream meshRandomizingStream = new MRG32k3a();
	PointSet qmcPointSet;

	public StochasticMeshAverageDensityQMC(MarketVectorProcess multiSp, PointSet qmcPointSet) {
		super(multiSp);
		this.qmcPointSet = qmcPointSet;
		this.setUsingQMC(true);
	}

	@Override
	public void resetNextRandomSubstream() {
		this.qmcPointSet.randomize(this.meshRandomizingStream);
		RandomStream stream = this.qmcPointSet.iterator();
		this.setRandomStreamForMeshStateGeneration(stream);
		stream.resetNextSubstream();

		this.marketProcess.setStream(stream);
	}
}