package gdsc.smlm.results;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specifies a peak fitting result
 */
public class PeakResult implements Comparable<PeakResult>
{
	/** Index of the background in the parameters array */
	public static final int BACKGROUND = 0;
	/** Index of the intensity in the parameters array */
	public static final int INTENSITY = 1;
	/** Index of the x-position in the parameters array */
	public static final int X = 2;
	/** Index of the y-position in the parameters array */
	public static final int Y = 3;
	/** Index of the z-position in the parameters array */
	public static final int Z = 4;
	/** Number of standard parameters */
	public static final int STANDARD_PARAMETERS = 5;

	private int frame;
	public int origX;
	public int origY;
	public float origValue;
	public double error;
	public float noise;

	// TODO - make this private
	float[] params;
	// TODO - make this private
	float[] paramsStdDev;

	/**
	 * Instantiates a new peak result.
	 *
	 * @param frame
	 *            the frame
	 * @param origX
	 *            the original X position
	 * @param origY
	 *            the original Y position
	 * @param origValue
	 *            the original value
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param params
	 *            the params (must not be null and must have at least {@value #STANDARD_PARAMETERS} parameters)
	 * @param paramsStdDev
	 *            the params standard deviations (if not null must match the length of the {@link #params} array)
	 * @throws IllegalArgumentException
	 *             the illegal argument exception if the parameters are invalid
	 */
	public PeakResult(int frame, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev) throws IllegalArgumentException
	{
		if (params == null)
			throw new IllegalArgumentException("Parameters must not be null");
		if (params.length < STANDARD_PARAMETERS)
			throw new IllegalArgumentException("Parameters must contain all standard parameters");
		if (paramsStdDev != null && paramsStdDev.length != params.length)
			throw new IllegalArgumentException("Parameters must contain all standard parameters");
		this.frame = frame;
		this.origX = origX;
		this.origY = origY;
		this.origValue = origValue;
		this.error = error;
		this.noise = noise;
		this.params = params;
		this.paramsStdDev = paramsStdDev;
	}

	/**
	 * Simple constructor to create a result with frame, location, and intensity.
	 *
	 * @param frame
	 *            the frame
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param intensity
	 *            the intensity
	 */
	public PeakResult(int frame, float x, float y, float intensity)
	{
		setFrame(frame);
		origX = (int) x;
		origY = (int) y;
		params = new float[5];
		params[X] = x;
		params[Y] = y;
		params[INTENSITY] = intensity;
	}

	/**
	 * Simple constructor to create a result with location and intensity.
	 *
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param intensity
	 *            the intensity
	 */
	public PeakResult(float x, float y, float intensity)
	{
		this(0, x, y, intensity);
	}

	/**
	 * Creates the params array for a peak result.
	 *
	 * @param background
	 *            the background
	 * @param intensity
	 *            the intensity
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @return the params array
	 */
	public static float[] createParams(float background, float intensity, float x, float y, float z)
	{
		return new float[] { background, intensity, x, y, z };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(PeakResult o)
	{
		// Sort by peak number: Ascending
		if (frame < o.frame)
			return -1;
		if (frame > o.frame)
			return 1;
		// Sort by peak height: Descending
		if (params[INTENSITY] > o.params[INTENSITY])
			return -1;
		if (params[INTENSITY] < o.params[INTENSITY])
			return 1;
		return 0;
	}

	/**
	 * @return The background for the first peak
	 */
	public float getBackground()
	{
		return params[BACKGROUND];
	}

	/**
	 * Sets the background.
	 *
	 * @param b
	 *            the new background
	 */
	public void setBackground(float b)
	{
		params[BACKGROUND] = b;
	}

	/**
	 * Get the signal strength (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @return The signal of the first peak
	 */
	public float getSignal()
	{
		return params[INTENSITY];
	}

	/**
	 * Sets the signal.
	 *
	 * @param s
	 *            the new signal
	 */
	public void setSignal(float s)
	{
		params[INTENSITY] = s;
	}

	/**
	 * @return The x position for the first peak
	 */
	public float getXPosition()
	{
		return params[X];
	}

	/**
	 * Sets the x position.
	 *
	 * @param x
	 *            the new x position
	 */
	public void setXPosition(float x)
	{
		params[X] = x;
	}

	/**
	 * @return The y position for the first peak
	 */
	public float getYPosition()
	{
		return params[Y];
	}

	/**
	 * Sets the y position.
	 *
	 * @param y
	 *            the new y position
	 */
	public void setYPosition(float y)
	{
		params[Y] = y;
	}

	/**
	 * @return The z position for the first peak
	 */
	public float getZPosition()
	{
		return params[Z];
	}

	/**
	 * Sets the z position.
	 *
	 * @param z
	 *            the new z position
	 */
	public void setZPosition(float z)
	{
		params[Z] = z;
	}

	/**
	 * Gets the frame.
	 *
	 * @return The time frame that this result corresponds to
	 */
	public int getFrame()
	{
		return frame;
	}

	/**
	 * Sets the frame.
	 *
	 * @param frame
	 *            The time frame that this result corresponds to
	 */
	public void setFrame(int frame)
	{
		this.frame = frame;
	}

	/**
	 * Checks for end frame. Derived classes can override this.
	 *
	 * @return true, if successful
	 */
	public boolean hasEndFrame()
	{
		return false;
	}

	/**
	 * Gets the end frame. Default = {@link #getFrame()}. Derived classes can override this.
	 *
	 * @return The last time frame that this result corresponds to
	 */
	public int getEndFrame()
	{
		return frame;
	}

	/**
	 * Checks for id. Derived classes can override this.
	 *
	 * @return true, if successful
	 */
	public boolean hasId()
	{
		return false;
	}

	/**
	 * Gets the id. Default = 0. Derived classes can override this.
	 *
	 * @return The results identifier
	 */
	public int getId()
	{
		return 0;
	}

	/**
	 * Checks for precision. Derived classes can override this.
	 *
	 * @return true, if successful
	 */
	public boolean hasPrecision()
	{
		return false;
	}

	/**
	 * Gets the localisation precision. Default = Double.NaN. Derived classes can override this.
	 * <p>
	 * This is provided so that results can be loaded from external sources.
	 *
	 * @return the precision (in nm)
	 */
	public double getPrecision()
	{
		return Double.NaN;
	}

	/**
	 * Return the true positive score for use in classification analysis
	 * 
	 * @return The true positive score
	 */
	public double getTruePositiveScore()
	{
		return (origValue != 0) ? 1 : 0;
	}

	/**
	 * Return the false positive score for use in classification analysis
	 * 
	 * @return The false positive score
	 */
	public double getFalsePositiveScore()
	{
		return 1 - getTruePositiveScore();
	}

	/**
	 * Return the true negative score for use in classification analysis
	 * 
	 * @return The true negative score
	 */
	public double getTrueNegativeScore()
	{
		return (origValue != 0) ? 0 : 1;
	}

	/**
	 * Return the false negative score for use in classification analysis
	 * 
	 * @return The false negative score
	 */
	public double getFalseNegativeScore()
	{
		return 1 - getTrueNegativeScore();
	}

	/**
	 * Return the squared distance to the other peak result
	 * 
	 * @param r
	 *            The result
	 * @return The squared distance
	 */
	public double distance2(PeakResult r)
	{
		final double dx = getXPosition() - r.getXPosition();
		final double dy = getYPosition() - r.getYPosition();
		return dx * dx + dy * dy;
	}

	/**
	 * Return the distance to the other peak result
	 * 
	 * @param r
	 *            The result
	 * @return The distance
	 */
	public double distance(PeakResult r)
	{
		return Math.sqrt(distance2(r));
	}

	/**
	 * Return the squared distance to the other coordinate.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return The squared distance
	 */
	public double distance2(final double x, final double y)
	{
		final double dx = getXPosition() - x;
		final double dy = getYPosition() - y;
		return dx * dx + dy * dy;
	}

	/**
	 * Return the distance to the other coordinate.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return The distance
	 */
	public double distance(final double x, final double y)
	{
		return Math.sqrt(distance2(x, y));
	}

	/**
	 * This methods return the x-position. To allow filters to use the actual shift requires either off-setting the
	 * position with the initial fit position, or extending this class so the shift can be stored.
	 */
	public float getXShift()
	{
		return getXPosition();
	}

	/**
	 * This methods return the y-position. To allow filters to use the actual shift requires either off-setting the
	 * position with the initial fit position, or extending this class so the shift can be stored.
	 */
	public float getYShift()
	{
		return getYPosition();
	}

	/**
	 * Gets the noise.
	 *
	 * @return the noise
	 */
	public float getNoise()
	{
		return noise;
	}

	/**
	 * Gets the parameters. This is a direct reference to the instance parameter array so use with caution.
	 *
	 * @return the parameters
	 */
	public float[] getParameters()
	{
		return params;
	}

	/**
	 * Gets the parameter deviations. This is a direct reference to the instance parameter array so use with caution.
	 *
	 * @return the parameter deviations
	 */
	public float[] getParameterDeviations()
	{
		return paramsStdDev;
	}

	/**
	 * Gets the number of parameters.
	 *
	 * @return the number of parameters
	 */
	public int getNumberOfParameters()
	{
		return params.length;
	}

	/**
	 * Gets the parameter for the given index.
	 *
	 * @param i
	 *            the index
	 * @return the parameter
	 */
	public float getParameter(int i)
	{
		return params[i];
	}

	/**
	 * Gets the parameter deviation for the given index.
	 *
	 * @param i
	 *            the index
	 * @return the parameter deviation
	 */
	public float getParameterDeviation(int i)
	{
		return paramsStdDev[i];
	}
}
