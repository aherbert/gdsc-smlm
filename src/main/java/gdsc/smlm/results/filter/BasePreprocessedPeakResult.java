package gdsc.smlm.results.filter;

import gdsc.core.match.FractionalAssignment;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.results.PeakResult;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specifies a peak fitting result for use in filtering.
 */
public class BasePreprocessedPeakResult implements AssignablePreprocessedPeakResult
{
	public enum ResultType
	{
		NEW, EXISTING, CANDIDATE
	}

	private final int frame;
	private final float signal;
	private final float photons;
	private final float snr;
	private final float noise;
	private final float sd;
	private final float b;
	private final float amp;
	private final float angle;
	private final float x;
	private final float y;
	private final float xshift2;
	private final float yshift2;
	private final float xsd;
	private final float ysd;
	private final float xwf;
	private final float ywf;
	private final double variance;
	private final boolean existingResult;
	private final boolean newResult;

	private ResultAssignment[] assignments;

	//@formatter:off
	/**
	 * Create a new BasePreprocessedPeakResult
	 * @param frame The frame
	 * @param signal The signal
	 * @param photons The signal calibrated as photons
	 * @param noise the noise estimate
	 * @param b The background level
	 * @param angle The angle of the fit
	 * @param x The x-position
	 * @param y The y-position
	 * @param x0 The initial x-position 
	 * @param y0 The initial y-position
	 * @param xsd The x standard deviation
	 * @param ysd The y standard deviation
	 * @param xsd0 The initial x standard deviation
	 * @param ysd0 The initial y standard deviation
	 * @param variance The estimate of the localisation variance
	 * @param resultType The type of result
	 */
	public BasePreprocessedPeakResult(
			int frame,
			float signal,
			float photons,
			float noise,
			float b,
			float angle,
			float x,
			float y,
			float x0,
			float y0,
			float xsd,
			float ysd,
			float xsd0,
			float ysd0,
			double variance,
			ResultType resultType
			)
	{
		//@formatter:on
		this.frame = frame;
		this.signal = signal;
		this.photons = photons;
		this.snr = signal / noise;
		this.noise = noise;
		this.sd = PeakResult.getSD(xsd, ysd);
		this.b = b;
		this.amp = (float) (signal / (2 * Math.PI * xsd * ysd));
		this.angle = angle;
		this.x = x;
		this.y = y;
		this.xshift2 = squared((x - x0) / xsd0);
		this.yshift2 = squared((y - y0) / ysd0);
		this.xsd = xsd;
		this.ysd = ysd;
		this.xwf = xsd / xsd0;
		this.ywf = ysd / ysd0;
		this.variance = variance;
		this.existingResult = resultType == ResultType.EXISTING;
		this.newResult = resultType == ResultType.NEW;
	}

	private static float squared(float f)
	{
		return f * f;
	}

	public int getFrame()
	{
		return frame;
	}

	public float getSignal()
	{
		return signal;
	}

	public float getPhotons()
	{
		return photons;
	}

	public float getSNR()
	{
		return snr;
	}

	public float getNoise()
	{
		return noise;
	}

	public double getLocationVariance()
	{
		return variance;
	}

	public float getSD()
	{
		return sd;
	}

	public float getBackground()
	{
		return b;
	}

	public float getAmplitude()
	{
		return amp;
	}

	public float getAngle()
	{
		return angle;
	}

	public float getX()
	{
		return x;
	}

	public float getY()
	{
		return y;
	}

	public float getXRelativeShift2()
	{
		return xshift2;
	}

	public float getYRelativeShift2()
	{
		return yshift2;
	}

	public float getXSD()
	{
		return xsd;
	}

	public float getYSD()
	{
		return ysd;
	}

	public float getXSDFactor()
	{
		return xwf;
	}

	public float getYSDFactor()
	{
		return ywf;
	}

	public boolean isExistingResult()
	{
		return existingResult;
	}

	public boolean isNewResult()
	{
		return newResult;
	}

	/**
	 * Returns a new array and so is thread-safe (unless another thread updates the assignments concurrently). It should
	 * be thread safe for use in scoring of the result using a multi-path filter.
	 * 
	 * @see gdsc.smlm.results.filter.PreprocessedPeakResult#getAssignments(int)
	 */
	public FractionalAssignment[] getAssignments(final int predictedId)
	{
		if (assignments == null || assignments.length == 0)
			return null;
		// Create a new set of assignments. Since this will be new and all other members are final the class is thread-safe.  
		final FractionalAssignment[] out = new FractionalAssignment[assignments.length];
		for (int i = 0; i < out.length; i++)
			out[i] = assignments[i].toFractionalAssignment(predictedId);
		return out;
	}

	public void setAssignments(ResultAssignment[] assignments)
	{
		this.assignments = assignments;
	}

	/**
	 * Convert this to the parameters for a Gaussian2DFunction
	 * 
	 * @return the parameters
	 */
	public double[] toGaussian2DParameters()
	{
		final double[] p = new double[7];
		p[Gaussian2DFunction.BACKGROUND] = b;
		p[Gaussian2DFunction.SIGNAL] = signal;
		p[Gaussian2DFunction.ANGLE] = angle;
		p[Gaussian2DFunction.X_POSITION] = x;
		p[Gaussian2DFunction.Y_POSITION] = y;
		p[Gaussian2DFunction.X_SD] = xsd;
		p[Gaussian2DFunction.Y_SD] = ysd;
		return p;
	}
}
