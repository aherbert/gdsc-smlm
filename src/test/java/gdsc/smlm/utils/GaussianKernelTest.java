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
package gdsc.smlm.utils;

import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Maths;
import gdsc.test.TestAssert;

@SuppressWarnings({ "javadoc" })
public class GaussianKernelTest
{
	@Test
	public void canGetConversionFactor()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final double norm = 1.0 / (Math.sqrt(2 * Math.PI) * s);
			final double var2 = 2 * s * s;

			final GaussianKernel k = new GaussianKernel(s);
			final double[] o = k.getGaussianKernel(1, 5, false);
			final double f = k.getConversionFactor(o);
			for (int u = o.length / 2, x = u; x >= 0; x--)
			{
				final double e = norm * FastMath.exp(-(x - u) * (x - u) / var2);
				Assert.assertEquals(e, f * o[x], e * 1e-10);
			}
		}
	}

	@Test
	public void canComputeGaussianKernelIncScaleIncRange()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int scale = 1; scale < 16; scale *= 2)
				for (int range = 3; range < 5; range++)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
		}
	}

	@Test
	public void canComputeGaussianKernelDecScaleIncRange()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int scale = 16; scale > 1; scale /= 2)
				for (int range = 3; range < 5; range++)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
		}
	}

	@Test
	public void canComputeGaussianKernelIncRangeIncScale()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int range = 3; range < 5; range++)
				for (int scale = 1; scale < 16; scale *= 2)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
		}
	}

	@Test
	public void canComputeGaussianKernelIncRangeDecScale()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int range = 3; range < 5; range++)
				for (int scale = 16; scale > 1; scale /= 2)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s * scale, range, edgeCorrection);
						final double[] o = k.getGaussianKernel(scale, range, edgeCorrection);

						Assert.assertArrayEquals(e, o, 0.0);
					}
		}
	}

	@Test
	public void canComputeDownscaleGaussianKernelIncScaleIncRange()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int scale = 1; scale <= 5; scale++)
				for (int range = 3; range < 5; range++)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
						final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

						TestAssert.assertArrayEqualsRelative(e, o, Maths.isPow2(scale) ? 0 : 1e-10);
					}
		}
	}

	@Test
	public void canComputeDownscaleGaussianKernelDecScaleIncRange()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int scale = 5; scale >= 1; scale--)
				for (int range = 3; range < 5; range++)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
						final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

						TestAssert.assertArrayEqualsRelative(e, o, Maths.isPow2(scale) ? 0 : 1e-10);
					}
		}
	}

	@Test
	public void canComputeDownscaleGaussianKernelIncRangeIncScale()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int range = 3; range < 5; range++)
				for (int scale = 1; scale <= 5; scale++)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
						final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

						TestAssert.assertArrayEqualsRelative(e, o, Maths.isPow2(scale) ? 0 : 1e-10);
					}
		}
	}

	@Test
	public void canComputeDownscaleGaussianKernelIncRangeDecScale()
	{
		for (int i = 0; i < 5; i++)
		{
			final double s = 0.33 * (1 << i);

			final GaussianKernel k = new GaussianKernel(s);

			for (int range = 3; range < 5; range++)
				for (int scale = 5; scale >= 1; scale--)
					for (final boolean edgeCorrection : new boolean[] { true, false })
					{
						final double[] e = GaussianKernel.makeGaussianKernel(s / scale, range, edgeCorrection);
						final double[] o = k.getDownscaleGaussianKernel(scale, range, edgeCorrection);

						TestAssert.assertArrayEqualsRelative(e, o, Maths.isPow2(scale) ? 0 : 1e-10);
					}
		}
	}
}
