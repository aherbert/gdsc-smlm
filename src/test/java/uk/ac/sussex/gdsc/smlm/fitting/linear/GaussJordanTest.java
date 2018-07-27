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
package uk.ac.sussex.gdsc.smlm.fitting.linear;

import java.util.Arrays;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;import uk.ac.sussex.gdsc.test.junit5.SeededTest;import uk.ac.sussex.gdsc.test.junit5.RandomSeed;import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class GaussJordanTest
{
	@Test
	public void canSolveLinearEquation()
	{
		final GaussJordan solver = new GaussJordan();

		// Solves (one) linear equation, a x = b, for x[n]

		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		final float[][] a = new float[][] { new float[] { 1, 2, 3 }, new float[] { 4, 5, 6 }, new float[] { 7, 8, 10 } };
		final float[] b = new float[] { 3, 3, 4 };
		final float[] expecteds = new float[] { -2, 1, 1 };

		final boolean result = solver.solve(a, b);

		Assertions.assertTrue(result);
		Assertions.assertArrayEquals(expecteds, b, 1e-4f);

		if (TestSettings.allow(LogLevel.INFO))
		{
			log("x = %s\n", Arrays.toString(b));
			for (int i = 0; i < b.length; i++)
				log("a[%d] = %s\n", i, Arrays.toString(a[i]));
		}
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
