package experiments;

import derivatives.Derivative;
import derivatives.DerivativePricedOnMesh;
import dynprog.LossFunction;
import dynprog.ProportionalCosts;
import dynprog.TransactionCosts;
import meshmethods.MeshWeights;
import meshmethods.StochasticMesh;
import meshmethods.StochasticMeshAverageDensity;
import meshmethods.StochasticMeshAverageDensityQMC;
import meshmethods.StochasticMeshSingleGrid;
import meshmethods.StochasticMeshSingleGridQMC;
import policies.HedgingPortfolio;
import stochprocess.DiscreteExpOULagIndicatorsProcess;
import stochprocess.DiscreteExpOULagProcess;
import stochprocess.DiscreteExpOUSyncUncorrelProcess;
import stochprocess.DiscreteGBMProcess;
import stochprocess.MarketVector;
import stochprocess.MarketVectorProcess;
import umontreal.ssj.hups.PointSet;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stochprocess.MultivariateBrownianMotion;
import umontreal.ssj.stochprocess.MultivariateBrownianMotionBridge;

public class ProblemDefinition {

	public enum MeshMethod {
		AD, ADQ, SG, SGQ
	};

	public enum PriceModel {
		GBM, GBM_B, EXP_OU, EXP_OU_B, EXP_OU_INDIC, EXP_OU_INDIC_B, EXP_OU_SYNC, EXP_OU_SYNC_B,
	};

	public double c0Start = 0.0;
	public Derivative derivative;
	public HedgingPortfolio hedgedPtf;
	public LossFunction lossFunction;
	public MarketVector[] mktInfoPath;
	public double[] obsTimes;
	public StochasticMesh stochasticMesh;
	public MarketVectorProcess stochprocessForExperiments;
	public TransactionCosts tcStructure;
	private double[] aVec = { Double.NaN }; // mean-vol
	private double[] bVec = { Double.NaN }; // mean-reversion strength
	private double[][] correlationMatrixExpOU;
	private double[][] correlationMatrixExpOUWithIndicators;
	private double[][] correlationMatrixGBM;
	private MultivariateBrownianMotion multiBM_expt = null;
	private MultivariateBrownianMotion multiBM_mesh = null;
	private double[] muVec = { Double.NaN };
	private PointSet psetAD;
	private PointSet psetSG;
	private double[] sigVec = { Double.NaN }; // vol of Log Variance
	private MarketVectorProcess spMesh;
	private RandomStream uStreamExpt;
	private RandomStream uStreamMesh;
	private float[] y0Array;
	private double priceVolCorrelation;
	private double[] indicatorLoadings;

	public ProblemDefinition() {
	}

	public static double[][] defineCorrelationMatrixForExpOU(double priceVolCorrelation) {
		int d = 2;
		double[][] squareMatrix = newUnitSquareMatrix(d);

		return squareMatrix;
	}

	public static double[][] defineCorrelationMatrixForExpOUWithIndicators(double priceVolCorrelation,
			double indicatorLoading1, double indicatorLoading2) {
		int d = 4;
		double[][] squareMatrix = newUnitSquareMatrix(d);

		return squareMatrix;
	}

	public static double[][] defineCorrelationMatrixForGBM() {
		int d = 1;
		double[][] squareMatrix = newUnitSquareMatrix(d);

		return squareMatrix;
	}

	private static double[][] newUnitSquareMatrix(int d) {
		double[][] unitSquareMatrix = new double[d][d];
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				if (i == j) {
					unitSquareMatrix[i][j] = 1.0;
				}
			}
		}
		return unitSquareMatrix;
	}

	public void assignStochasticMeshToDerivative() {
		boolean derivativeRequiresStochasticMesh = (this.derivative instanceof DerivativePricedOnMesh);
		if (derivativeRequiresStochasticMesh) {
			((DerivativePricedOnMesh) this.derivative).setStochasticMesh(this.stochasticMesh);
		}
	}

	public void computeMesh(SimulationParameters simParams) {
		this.stochasticMesh.setObservationTimes(this.obsTimes, simParams.nbObs);
		this.stochasticMesh.setNbPaths(simParams.nbInnerSims);
		this.stochasticMesh.generateMesh();

	}

	public void defineObservationTimes(SimulationParameters simulationParameters, TestedPolicies testPolicies) {
		double dt = simulationParameters.t / simulationParameters.nbObs;
		this.stochprocessForExperiments.setObservationTimes(dt, simulationParameters.nbObs);

		this.obsTimes = this.stochprocessForExperiments.getObservationTimes();
		this.hedgedPtf.setObservationTimes(this.obsTimes);
		for (int k = 0; k < testPolicies.getNbPolicies(); k++) {
			testPolicies.getPolicyArray()[k].setObservationTimes(this.obsTimes);
		}

		this.mktInfoPath = this.stochprocessForExperiments.getMarketVectorPath();
	}

	public void defineStochasticMesh(SimulationParameters simParams) {
		int nbProcessDimensions = this.spMesh.getDimension();
		int nbOfPoints = simParams.nbInnerSims;
		int nbSteps = simParams.nbObs;
		switch (simParams.meshMethod) {
			case AD:
				this.stochasticMesh = new StochasticMeshAverageDensity(this.spMesh);
				break;
			case ADQ:
				this.psetAD = new SobolSequence(nbOfPoints, nbSteps * nbProcessDimensions);
				this.stochasticMesh = new StochasticMeshAverageDensityQMC(this.spMesh, this.psetAD);
				break;
			case SG:
				this.stochasticMesh = new StochasticMeshSingleGrid(this.spMesh);
				break;
			case SGQ:
				this.psetSG = new SobolSequence(nbOfPoints, nbProcessDimensions);
				this.stochasticMesh = new StochasticMeshSingleGridQMC(this.spMesh, this.psetSG);
				break;
			default:
				this.stochasticMesh = null;
		}

		if (this.stochasticMesh != null) {
			MeshWeights meshWeights = this.stochasticMesh.getMeshWeights();
			meshWeights.useRussianRoulette(simParams.weightCutoff, new MRG32k3a());
		}
	}

	public void initializePortfolioHedgingProblem(SimulationParameters simParams) {
		this.tcStructure = new ProportionalCosts(simParams.aTC, simParams.bTC);
		this.hedgedPtf = new HedgingPortfolio(simParams.derivative, this.tcStructure, simParams.rfRate);
		this.hedgedPtf.setMaxAllowedOneStepLoss(simParams.maxAllowedLoss);
		this.derivative = simParams.derivative;
	}

	public void initializeProcesses(SimulationParameters simParams) {
		this.initializeProcessParameters(simParams);
		this.defineRandomStreams();
		this.defineBrownianMotions(simParams);
		this.defineProcessesAccordingToModel(simParams.priceModel);
		this.setInitialValueOfProcessForExperiments();
	}

	private void createBrownianMotions(int nbDim, NormalGen normGen_expt, NormalGen normGen_mesh,
			double[][] correlationMatrix) {
		double[] zeroVec = new double[nbDim];
		double[] unitVec = this.getOneVector(nbDim);
		this.multiBM_expt = new MultivariateBrownianMotion(nbDim, zeroVec, zeroVec, unitVec, correlationMatrix,
				normGen_expt);
		this.multiBM_mesh = new MultivariateBrownianMotion(nbDim, zeroVec, zeroVec, unitVec, correlationMatrix,
				normGen_mesh);
	}

	private void createBrownianMotionsBridge(int nbDim, NormalGen normGen_expt, NormalGen normGen_mesh,
			double[][] correlationMatrix) {
		double[] zeroVec = new double[nbDim];
		double[] unitVec = this.getOneVector(nbDim);
		this.multiBM_expt = new MultivariateBrownianMotionBridge(nbDim, zeroVec, zeroVec, unitVec, correlationMatrix,
				normGen_expt);
		this.multiBM_mesh = new MultivariateBrownianMotionBridge(nbDim, zeroVec, zeroVec, unitVec, correlationMatrix,
				normGen_mesh);
	}

	private void defineBrownianMotions(SimulationParameters simParams) {
		this.defineCorrelationMatrix(simParams);

		NormalGen normGen_expt = new NormalGen(this.uStreamExpt, new NormalDist());
		NormalGen normGen_mesh = new NormalGen(this.uStreamMesh, new NormalDist());

		switch (simParams.priceModel) {
			case GBM:
				this.createBrownianMotions(1, normGen_expt, normGen_mesh, this.correlationMatrixGBM);
				break;
			case EXP_OU:
			case EXP_OU_SYNC:
				this.createBrownianMotions(2, normGen_expt, normGen_mesh, this.correlationMatrixExpOU);
				break;
			case GBM_B:
				this.createBrownianMotionsBridge(1, normGen_expt, normGen_mesh, this.correlationMatrixGBM);
				break;
			case EXP_OU_B:
			case EXP_OU_SYNC_B:
				this.createBrownianMotionsBridge(2, normGen_expt, normGen_mesh, this.correlationMatrixExpOU);
				break;
			case EXP_OU_INDIC:
				this.createBrownianMotions(4, normGen_expt, normGen_mesh, this.correlationMatrixExpOUWithIndicators);
				break;
			case EXP_OU_INDIC_B:
				this.createBrownianMotionsBridge(4, normGen_expt, normGen_mesh,
						this.correlationMatrixExpOUWithIndicators);
				break;
			default:
		}
	}

	private void defineCorrelationMatrix(SimulationParameters simParams) {
		this.correlationMatrixGBM = ProblemDefinition.defineCorrelationMatrixForGBM();
		this.correlationMatrixExpOU = defineCorrelationMatrixForExpOU(simParams.priceVolCorrelation);
		this.correlationMatrixExpOUWithIndicators = defineCorrelationMatrixForExpOUWithIndicators(
				simParams.priceVolCorrelation, simParams.indicatorLoading1, simParams.indicatorLoading2);
	}

	private void defineProcessesAccordingToModel(PriceModel model) {
		switch (model) {
			case GBM:
			case GBM_B:
				this.stochprocessForExperiments = new DiscreteGBMProcess(this.y0Array, this.multiBM_expt);
				this.spMesh = new DiscreteGBMProcess(this.y0Array, this.multiBM_mesh);
				break;

			case EXP_OU:
			case EXP_OU_B:
				this.stochprocessForExperiments = new DiscreteExpOULagProcess(this.y0Array, this.muVec, this.aVec,
						this.bVec, this.sigVec, this.priceVolCorrelation, this.multiBM_expt);
				this.spMesh = new DiscreteExpOULagProcess(this.y0Array, this.muVec, this.aVec, this.bVec, this.sigVec,
						this.priceVolCorrelation, this.multiBM_mesh);
				break;

			case EXP_OU_SYNC:
			case EXP_OU_SYNC_B:
				this.stochprocessForExperiments = new DiscreteExpOUSyncUncorrelProcess(this.y0Array, this.muVec,
						this.aVec, this.bVec, this.sigVec, this.multiBM_expt);
				this.spMesh = new DiscreteExpOUSyncUncorrelProcess(this.y0Array, this.muVec, this.aVec, this.bVec,
						this.sigVec, this.multiBM_mesh);
				break;

			case EXP_OU_INDIC:
			case EXP_OU_INDIC_B:
				this.stochprocessForExperiments = new DiscreteExpOULagIndicatorsProcess(this.y0Array, this.muVec,
						this.aVec, this.bVec, this.sigVec, this.priceVolCorrelation, this.multiBM_expt,
						this.indicatorLoadings, new MRG32k3a());
				this.spMesh = new DiscreteExpOULagIndicatorsProcess(this.y0Array, this.muVec, this.aVec, this.bVec,
						this.sigVec, this.priceVolCorrelation, this.multiBM_mesh, this.indicatorLoadings,
						new MRG32k3a());
				break;

			default:
		}
	}

	private void defineRandomStreams() {
		this.uStreamExpt = new MRG32k3a();
		this.uStreamMesh = new MRG32k3a();
	}

	private double[] getOneVector(int nbDim) {
		double[] oneVec = new double[nbDim];
		for (int i = 0; i < nbDim; i++) {
			oneVec[i] = 1.0;
		}
		return oneVec;
	}

	private void initializeProcessParameters(SimulationParameters simParams) {
		this.y0Array = new float[2];
		this.y0Array[0] = (float) simParams.s0;
		this.y0Array[1] = (float) simParams.sigma0;

		this.aVec[0] = simParams.sigma1;
		this.bVec[0] = simParams.meanReversionRate;
		this.sigVec[0] = simParams.volOfLogVar;
		this.priceVolCorrelation = simParams.priceVolCorrelation;

		this.indicatorLoadings = new double[2];
		this.indicatorLoadings[0] = simParams.indicatorLoading1;
		this.indicatorLoadings[1] = simParams.indicatorLoading2;
	}

	private void setInitialValueOfProcessForExperiments() {
		this.stochprocessForExperiments.setX0(this.y0Array);
	}

}
