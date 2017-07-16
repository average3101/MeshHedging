package numerics;

public abstract class OneDimSolver {

	protected OneDimFunction f;
	protected OneDimSolverResults results = new OneDimSolverResults();
	protected double tol;

	public abstract OneDimSolverResults minimize(double ax, double bx, double cx);

	public void setTolerance(double tol) {
		this.tol = tol;
	}
}
