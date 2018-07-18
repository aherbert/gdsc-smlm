/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;

/**
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The function will calculate the value of the Gaussian and evaluate the gradient of a set of parameters. The class can
 * specify which of the following parameters the function will evaluate:<br/>
 * background, signal, z-depth, position0, position1, sd0, sd1
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class SingleErfGaussian2DFunction extends ErfGaussian2DFunction
{
	// Required for the PSF

	/** The intensity. */
	protected double tI;

	/**
	 * Instantiates a new erf gaussian 2D function.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public SingleErfGaussian2DFunction(int maxx, int maxy)
	{
		super(1, maxx, maxy);
	}

	@Override
	public int getNPeaks()
	{
		return 1;
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 *
	 * @param i
	 *            Input predictor
	 * @return The Gaussian value
	 */
	@Override
	public double eval(final int i)
	{
		// Unpack the predictor into the dimensions
		final int y = i / maxx;
		final int x = i % maxx;

		return tB + tI * deltaEx[x] * deltaEy[y];
	}

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 *
	 * @param i
	 *            Input predictor
	 * @param duda
	 *            Partial gradient of function with respect to each coefficient
	 * @return The predicted value
	 *
	 * @see uk.ac.sussex.gdsc.smlm.function.NonLinearFunction#eval(int, double[])
	 */
	@Override
	public abstract double eval(final int i, final double[] duda);

	/**
	 * Evaluates an 2-dimensional Gaussian function for a single peak.
	 *
	 * @param i
	 *            Input predictor
	 * @param duda
	 *            Partial first gradient of function with respect to each coefficient
	 * @param d2uda2
	 *            Partial second gradient of function with respect to each coefficient
	 * @return The predicted value
	 */
	@Override
	public abstract double eval(final int i, final double[] duda, final double[] d2uda2);

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.function.GradientFunction#forEach(uk.ac.sussex.gdsc.smlm.function.GradientFunction.ValueProcedure)
	 */
	@Override
	public void forEach(ValueProcedure procedure)
	{
		if (tB == 0)
			for (int y = 0; y < maxy; y++)
			{
				final double tI_deltaEy = tI * deltaEy[y];
				for (int x = 0; x < maxx; x++)
					procedure.execute(tI_deltaEy * deltaEx[x]);
			}
		else
			for (int y = 0; y < maxy; y++)
			{
				final double tI_deltaEy = tI * deltaEy[y];
				for (int x = 0; x < maxx; x++)
					procedure.execute(tB + tI_deltaEy * deltaEx[x]);
			}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.function.ExtendedNonLinearFunction#computeValues(double[])
	 */
	@Override
	public double[] computeValues(double[] variables)
	{
		initialise0(variables);
		final double[] values = new double[size()];
		if (tB == 0)
			for (int y = 0, i = 0; y < maxy; y++)
			{
				final double tI_deltaEy = tI * deltaEy[y];
				for (int x = 0; x < maxx; x++)
					values[i++] = tI_deltaEy * deltaEx[x];
			}
		else
			for (int y = 0, i = 0; y < maxy; y++)
			{
				final double tI_deltaEy = tI * deltaEy[y];
				for (int x = 0; x < maxx; x++)
					values[i++] = tB + tI_deltaEy * deltaEx[x];
			}
		return values;
	}

	// Force implementation
	@Override
	public abstract int getNumberOfGradients();

	// Force implementation
	@Override
	public abstract double integral(double[] a);
}
