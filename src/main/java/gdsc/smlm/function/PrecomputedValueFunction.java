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
 * Wraps a set of function values to implement the forEach procedure
 */
public class PrecomputedValueFunction implements ValueFunction
{
	/** The values. */
	protected final double[] values;

	/**
	 * Instantiates a new pre-computed value function.
	 *
	 * @param values
	 *            the pre-computed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	public PrecomputedValueFunction(double[] values)
	{
		if (values == null)
			throw new IllegalArgumentException("Value is null");
		this.values = values;
	}

	/**
	 * Gets a reference to the first order gradients
	 *
	 * @return the values
	 */
	public double[] getValuesRef()
	{
		return values;
	}

	@Override
	public int size()
	{
		return values.length;
	}

	@Override
	public void initialise0(double[] a)
	{
		// Ignore
	}

	@Override
	public void forEach(ValueProcedure procedure)
	{
		for (int i = 0; i < values.length; i++)
			procedure.execute(values[i]);
	}
}
