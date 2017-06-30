package gdsc.smlm.function.gaussian.erf;

import gdsc.smlm.function.ValueProcedure;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The function will calculate the value of the Gaussian and evaluate the gradient of a set of parameters. The class can
 * specify which of the following parameters the function will evaluate:<br/>
 * background, signal, z-depth, position0, position1, sd0, sd1
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class MultiErfGaussian2DFunction extends ErfGaussian2DFunction
{
	protected final int nPeaks;
	protected final int[] gradientIndices;

	// Required for the PSF
	protected final double[] tI;

	/**
	 * Instantiates a new multi-peak erf gaussian 2D function.
	 *
	 * @param nPeaks
	 *            The number of peaks
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public MultiErfGaussian2DFunction(int nPeaks, int maxx, int maxy)
	{
		super(nPeaks, maxx, maxy);
		this.nPeaks = nPeaks;
		this.gradientIndices = createGradientIndices();
		tI = new double[nPeaks];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.function.GaussianFunction#getNPeaks()
	 */
	@Override
	public int getNPeaks()
	{
		return nPeaks;
	}

	/**
	 * Creates the gradient indices.
	 *
	 * @return the the gradient indices
	 */
	abstract protected int[] createGradientIndices();

	/**
	 * Replicate the gradient indices from a single peak for the configured number of peaks.
	 *
	 * @param singleGradientIndices
	 *            the single gradient indices
	 */
	protected int[] replicateGradientIndices(int[] singleGradientIndices)
	{
		final int start = (evaluatesBackground() ? 1 : 0);
		final int m = singleGradientIndices.length;
		final int[] indices = new int[start + nPeaks * (m - start)];
		int p = 0;
		if (evaluatesBackground())
			indices[p++] = 0;
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			for (int j = start; j < m; j++)
				indices[p++] = i + singleGradientIndices[j];
		}
		return indices;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return gradientIndices;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#getNumberOfGradients()
	 */
	public int getNumberOfGradients()
	{
		return gradientIndices.length;
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 * 
	 * @param i
	 *            Input predictor
	 * @return The Gaussian value
	 * 
	 * @see gdsc.fitting.function.NonLinearFunction#eval(int)
	 */
	public double eval(final int i)
	{
		// Unpack the predictor into the dimensions
		int yy = i / maxx;
		int xx = i % maxx;

		double I = tB;
		for (int n = 0; n < nPeaks; n++, xx += maxx, yy += maxy)
			I += tI[n] * deltaEx[xx] * deltaEy[yy];
		return I;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.GradientFunction.ValueProcedure)
	 */
	public void forEach(ValueProcedure procedure)
	{
		// Unroll for the number of peaks
		if (nPeaks == 2)
		{
			if (tB == 0)
			{
				for (int y = 0; y < maxy; y++)
				{
					// Pre-compute
					final double tI_deltaEy0 = tI[0] * deltaEy[y];
					final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

					for (int x = 0; x < maxx; x++)
					{
						procedure.execute(tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
					}
				}
			}
			else
			{
				for (int y = 0; y < maxy; y++)
				{
					// Pre-compute
					final double tI_deltaEy0 = tI[0] * deltaEy[y];
					final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

					for (int x = 0; x < maxx; x++)
					{
						procedure.execute(tB + tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
					}
				}
			}
		}
		else
		{
			final double[] tI_deltaEy = new double[nPeaks];
			for (int y = 0; y < maxy; y++)
			{
				// Pre-compute
				for (int n = 0, yy = y; n < nPeaks; n++, yy += maxy)
					tI_deltaEy[n] = tI[n] * deltaEy[yy];

				for (int x = 0; x < maxx; x++)
				{
					double I = tB;
					for (int n = 0, xx = x; n < nPeaks; n++, xx += maxx)
						I += tI_deltaEy[n] * deltaEx[xx];
					procedure.execute(I);
				}

				// No pre-compute
				//			for (int x = 0; x < maxx; x++)
				//			{
				//				double I = tB;
				//				for (int n = 0, xx = x, yy = y; n < nPeaks; n++, xx += maxx, yy += maxy)
				//					I += tI[n] * deltaEx[xx] * deltaEy[yy];
				//				procedure.execute(I);
				//			}
			}
		}
	}
}
