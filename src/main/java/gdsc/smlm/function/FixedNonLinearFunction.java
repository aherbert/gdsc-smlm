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
package gdsc.smlm.function;

/**
 * Implement a fixed value non-linear fitting function
 */
public class FixedNonLinearFunction implements NonLinearFunction
{
	/** The values. */
	final double[] values;

	/**
	 * Instantiates a new fixed non linear function.
	 *
	 * @param values
	 *            the values
	 */
	public FixedNonLinearFunction(double[] values)
	{
		this.values = values;
	}

	@Override
	public void initialise(double[] a)
	{
		// Do nothing
	}

	@Override
	public int[] gradientIndices()
	{
		return new int[0];
	}

	@Override
	public int getNumberOfGradients()
	{
		return 0;
	}

	@Override
	public double eval(int x, double[] dyda)
	{
		return values[x];
	}

	@Override
	public double eval(int x)
	{
		return values[x];
	}

	@Override
	public double eval(int x, double[] dyda, double[] w)
	{
		return values[x];
	}

	@Override
	public double evalw(int x, double[] w)
	{
		return values[x];
	}

	@Override
	public boolean canComputeWeights()
	{
		return false;
	}
}
