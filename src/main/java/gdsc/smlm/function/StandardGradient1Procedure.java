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
 * Class for evaluating a function
 */
public class StandardGradient1Procedure implements Gradient1Procedure
{
	private int i;

	/** The values from the last call to {@link #getValues(Gradient1Function, double[])}. */
	public double[] values;
	/** The gradients from the last call to {@link #getValues(Gradient1Function, double[])}. */
	public double[][] dyda;

	/**
	 * Gets the values.
	 *
	 * @param f
	 *            the function
	 * @param a
	 *            the function coefficients
	 * @return the values
	 */
	public double[] getValues(Gradient1Function f, double[] a)
	{
		values = new double[f.size()];
		dyda = new double[values.length][];
		i = 0;
		f.initialise1(a);
		f.forEach(this);
		return values;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	@Override
	public void execute(double value, double[] dy_da)
	{
		values[i] = value;
		dyda[i] = dy_da.clone();
		i++;
	}
}
