package numerics;

public enum NormalGenerationMethod {
	chol(0, "CHOL"), pca(1, "PCA");
	private final String name;
	private final int value;

	private NormalGenerationMethod(int value, String name) {
		this.value = value;
		this.name = name;
	}

	public String getName() {
		return this.name;
	}

	public int getValue() {
		return this.value;
	}
};
