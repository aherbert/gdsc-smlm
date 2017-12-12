package gdsc.smlm.filters;

public abstract class IntDataFilter
{
	/** The name. */
	final String name;

	/** Flag to indicate that the filter can support non-integer box sizes. */
	final boolean isInterpolated;

	/** The min box size. */
	final int minBoxSize;

	/**
	 * Instantiates a new data filter.
	 *
	 * @param name
	 *            the name
	 * @param isInterpolated
	 *            Flag to indicate that the filter can support non-integer box sizes
	 */
	public IntDataFilter(String name, boolean isInterpolated)
	{
		this(name, isInterpolated, 0);
	}

	/**
	 * Instantiates a new data filter.
	 *
	 * @param name
	 *            the name
	 * @param isInterpolated
	 *            Flag to indicate that the filter can support non-integer box sizes
	 * @param minBoxSize
	 *            the min box size
	 */
	public IntDataFilter(String name, boolean isInterpolated, int minBoxSize)
	{
		this.name = name;
		this.isInterpolated = isInterpolated;
		this.minBoxSize = minBoxSize;
	}

	public abstract void filter(int[] data, int width, int height, int boxSize);

	public abstract void filterInternal(int[] data, int width, int height, int boxSize);

	public abstract void setWeights(float[] weights, int width, int height);
}
