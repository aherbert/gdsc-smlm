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
package gdsc.smlm.function.cspline;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.test.TestSettings;

public class CubicSplineDataTest
{
	@Test
	public void canExternaliseDoubleFunction() throws IOException
	{
		canExternaliseFunction(false);
	}

	@Test
	public void canExternaliseFloatFunction() throws IOException
	{
		canExternaliseFunction(true);
	}

	private void canExternaliseFunction(boolean singlePrecision) throws IOException
	{
		RandomGenerator r = TestSettings.getRandomGenerator();
		int x = 6, y = 5, z = 4;

		int size = x * y;
		CustomTricubicFunction[][] splines = new CustomTricubicFunction[z][x * y];
		double[] a = new double[64];
		for (int zz = 0; zz < z; zz++)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < 64; j++)
					a[j] = r.nextDouble();
				splines[zz][i] = CustomTricubicFunction.create(a.clone());
				if (singlePrecision)
					splines[zz][i] = splines[zz][i].toSinglePrecision();
			}
		}
		CubicSplineData f1 = new CubicSplineData(x, y, splines);

		ByteArrayOutputStream b = new ByteArrayOutputStream();
		f1.write(b);

		byte[] bytes = b.toByteArray();
		CubicSplineData f2 = CubicSplineData.read(new ByteArrayInputStream(bytes));

		for (int zz = 0; zz < z; zz++)
		{
			for (int i = 0; i < size; i++)
			{
				Assert.assertArrayEquals(f1.splines[zz][i].getA(), f2.splines[zz][i].getA(), 0);
			}
		}
	}
}
