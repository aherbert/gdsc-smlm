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
public class StandardFloatValueProcedure implements ValueProcedure
{
	private int i;
	/** The values from the last call to {@link #getValues(ValueFunction, float[])}. */
	public float[] values;

	/**
	 * Gets the values.
	 *
	 * @param f
	 *            the function
	 * @param a
	 *            the function coefficients
	 * @return the values
	 */
	public float[] getValues(ValueFunction f, double[] a)
	{
		values = new float[f.size()];
		i = 0;
		f.initialise0(a);
		f.forEach(this);
		return values;
	}

	/**
	 * Gets the values into the buffer.
	 *
	 * @param f
	 *            the function
	 * @param a
	 *            the function coefficients
	 * @param buffer
	 *            the buffer
	 * @param offset
	 *            the offset
	 */
	public void getValues(ValueFunction f, double[] a, float[] buffer, int offset)
	{
		if (buffer == null || buffer.length < offset + f.size())
			throw new IllegalArgumentException("Buffer is not large enough for the function values");
		values = buffer;
		i = offset;
		f.initialise0(a);
		f.forEach(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(float)
	 */
	public void execute(double value)
	{
		values[i++] = (float) value;
	}
}
