package meshmethods;

import java.text.DecimalFormat;

import dynprog.State;

public class StateInformationPrinter {

	private StochasticMesh stochasticMesh;

	public StateInformationPrinter(StochasticMesh stochasticMesh) {
		this.stochasticMesh = stochasticMesh;
	}

	public void printListOfStatesAtStep(int k) {
		DecimalFormat df = new DecimalFormat("0.00");

		float[][] meshCS0 = this.stochasticMesh.getStatePointsList().get(k);

		System.out.println();
		System.out.println("State info at step " + k);
		for (int i = 0; i < this.stochasticMesh.getNbPaths(); i++) {
			String outputLine = Integer.toString(i);

			for (int l = 0; l < this.stochasticMesh.getNbProcess(); l++) {
				outputLine += "; " + df.format(meshCS0[i][l]);
			}

			System.out.println(outputLine);
		}
	}

	public void printListOfWeightsAtState(int k, int index) {
		DecimalFormat df = new DecimalFormat("0.000000");

		System.out.println();
		System.out.println("Weights from state k=" + k + ", i=" + index);
		double sum = 0;
		MeshWeights meshWeights = this.stochasticMesh.getMeshWeights();
		for (int i = 0; i < this.stochasticMesh.getNbPaths(); i++) {
			double w = meshWeights.getWeightsListForStepAndOriginIndex()[k][index].getQuick(i);
			String outputLine = Integer.toString(i) + "; " + df.format(w);
			sum += w;

			System.out.println(outputLine);
		}
		System.out.println("Sum; " + sum);
	}

	public void printListOfWeightsAtState(State s) {
		this.printListOfWeightsAtState(s.timeStep, s.marketVector.simulationIndex);
	}

	public void printMemory() {
		long mem1 = Runtime.getRuntime().totalMemory();
		long mem2 = Runtime.getRuntime().freeMemory();
		long mem3 = mem1 - mem2;
		System.out.println("Tot mem :" + mem1 + ", Free : " + mem2 + ", Used: " + mem3);
	}
}
