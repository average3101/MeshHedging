package numerics;

public final class OneDimSolverResults {
	public int nbSteps;
	public double optim_f; // objective min
	public double optim_v; // argmin
	public int terminationCode;

	public OneDimSolverResults() {
	}
}