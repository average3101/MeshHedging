package meshmethods;

import java.util.HashSet;
import java.util.Iterator;

public abstract class MeshFullEstimator extends MeshPointEstimator {

	private int nbSteps;

	public abstract void computeForState(int k, int i);

	public abstract void computeForStateOverLastStep(int i);

	public void computeOverAllSteps() {
		this.initialize();

		this.computeOverLastStep();
		// if (this.nbSteps > 1) {
		this.computeOverIntermediateSteps();
		// }
		this.computeAverageEstimate();
	}

	public void computeOverLastStep() {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		HashSet<Integer>[] availableStateIndices = stochasticMesh.getAttainableStateIndices();
		Iterator<Integer> availableIndicesIterator = availableStateIndices[this.nbSteps].iterator();

		while (availableIndicesIterator.hasNext()) {
			int i = availableIndicesIterator.next().intValue();

			this.computeForStateOverLastStep(i);
		}
	}

	public int getNbSteps() {
		return this.nbSteps;
	}

	public void initialize() {
	}

	public void setNbSteps(int nbSteps) {
		this.nbSteps = nbSteps;
	}

	private void computeOverIntermediateSteps() {
		StochasticMesh stochasticMesh = this.getStochasticMesh();
		HashSet<Integer>[] availableStateIndices = stochasticMesh.getAttainableStateIndices();
		for (int k = this.nbSteps - 1; k >= 0; k--) {
			Iterator<Integer> availableIndicesIterator = availableStateIndices[k].iterator();
			while (availableIndicesIterator.hasNext()) {
				int i = availableIndicesIterator.next().intValue();

				this.computeForState(k, i);
			}
		}
	}

	protected void computeAverageEstimate() {
		// int stepNo = 0;
		// int index = 0;
		// this.computeForState(stepNo, index);
	}

}
