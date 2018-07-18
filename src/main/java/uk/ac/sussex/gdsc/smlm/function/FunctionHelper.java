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

import java.util.Arrays;

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.utils.Maths;

/**
 * Utility class for functions.
 */
public class FunctionHelper
{
	/**
	 * Gets the mean value using a fraction of the cumulative value, when values are sorted in descending order. All
	 * values must be positive. The input values are modified by sorting.
	 * <p>
	 * If fraction is <=0 then the max value is returned. If fraction is >=1 then the mean of the data is returned.
	 *
	 * @param values
	 *            the values
	 * @param fraction
	 *            the fraction
	 * @return the mean value
	 * @throws DataException
	 *             if the values are not positive.
	 */
	public static double getMeanValue(double[] values, double fraction) throws DataException
	{
		if (fraction <= 0)
			return Maths.max(values);
		double sum = 0;
		for (int i = 0; i < values.length; i++)
		{
			if (values[i] < 0)
				throw new DataException("Values must be positive");
			sum += values[i];
		}
		if (fraction >= 1)
			return sum / values.length;
		final double target = sum * fraction;
		sum = 0;
		Arrays.sort(values);
		for (int i = values.length; i-- > 0;)
		{
			sum += values[i];
			if (sum >= target)
			{
				// Interpolate the count X to obtain the target
				final int x = values.length - i;
				return target / Maths.interpolateX(x - 1, sum - values[i], x, sum, target);
			}
		}
		// Edge case
		return sum / values.length;
	}

	/**
	 * Gets the x-value corresponding to a fraction of the cumulative value, when values are sorted in descending order.
	 * All
	 * values must be positive. The input values are modified by sorting.
	 * <p>
	 * If fraction is <=0 then zero is returned. If fraction is >=1 then data.length is returned.
	 *
	 * @param values
	 *            the values
	 * @param fraction
	 *            the fraction
	 * @return the x-value
	 * @throws DataException
	 *             if the values are not positive.
	 */
	public static double getXValue(double[] values, double fraction) throws DataException
	{
		if (fraction <= 0)
			return 0;
		if (fraction >= 1)
			return values.length;
		double sum = 0;
		for (int i = 0; i < values.length; i++)
		{
			if (values[i] < 0)
				throw new DataException("Values must be positive");
			sum += values[i];
		}
		final double target = sum * fraction;
		sum = 0;
		Arrays.sort(values);
		for (int i = values.length; i-- > 0;)
		{
			sum += values[i];
			if (sum >= target)
			{
				// Interpolate the count X to obtain the target
				final int x = values.length - i;
				return Maths.interpolateX(x - 1, sum - values[i], x, sum, target);
			}
		}
		// Edge case
		return values.length;
	}
}
