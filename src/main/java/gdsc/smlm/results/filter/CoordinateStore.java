package gdsc.smlm.results.filter;

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows checking for duplicates.
 */
public interface CoordinateStore
{
	/**
	 * Gets the resolution of the store.
	 *
	 * @return the resolution
	 */
	public double getResolution();

	/**
	 * Queue a coordinate to the store.
	 * <p>
	 * It is not added to the store until flush is called.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 */
	public void addToQueue(double x, double y);

	/**
	 * Flush the queue to the store
	 */
	public void flush();

	/**
	 * Add a coordinate to the store. Assumes that the coordinates are within the size of the grid otherwise they will be ignored.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 */
	public void add(double x, double y);

	/**
	 * Clear to the store.
	 */
	public void clear();

	/**
	 * Check if the store contains the coordinates within the configured resolution.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return true, if the store contains another coordinate closer than the resolution
	 */
	public boolean contains(double x, double y);

	/**
	 * Find the closest coordinate within the configured resolution.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return the coordinate closer than the resolution (or null)
	 */
	public double[] find(double x, double y);
}