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
package uk.ac.sussex.gdsc.smlm.function;

/**
 * Class for evaluating a function
 */
public class StandardGradient2Procedure implements Gradient2Procedure
{
	private int i;

	/** The values from the last call to {@link #getValues(Gradient2Function, double[])}. */
	public double[] values;
	/** The gradients from the last call to {@link #getValues(Gradient2Function, double[])}. */
	public double[][] dyda;
	/** The second order gradients from the last call to {@link #getValues(Gradient2Function, double[])}. */
	public double[][] d2yda2;

	/**
	 * Gets the values.
	 *
	 * @param f
	 *            the function
	 * @param a
	 *            the function coefficients
	 * @return the values
	 */
	public double[] getValues(Gradient2Function f, double[] a)
	{
		values = new double[f.size()];
		dyda = new double[values.length][];
		d2yda2 = new double[values.length][];
		i = 0;
		f.initialise2(a);
		f.forEach(this);
		return values;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	@Override
	public void execute(double value, double[] dy_da, double[] d2y_da2)
	{
		values[i] = value;
		dyda[i] = dy_da.clone();
		d2yda2[i] = d2y_da2.clone();
		i++;
	}
}