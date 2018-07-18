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
package gdsc.smlm.function.gaussian;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.test.LogLevel;
import gdsc.test.TestLog;
import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class QuadraticAstigmatismZModelTest
{
	protected DoubleEquality eq = new DoubleEquality(1e-5, 1e-7);

	// Compute as per Numerical Recipes 5.7.
	// Approximate error accuracy in single precision: Ef
	// Step size for derivatives:
	// h ~ (Ef)^(1/3) * xc
	// xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
	protected double h_ = 0.0001; //(double) (Math.pow(1e-3f, 1.0 / 3));

	@Test
	public void canStaticComputeGradient()
	{
		canStaticComputeGradient(1.2);
		canStaticComputeGradient(0.3);
	}

	private void canStaticComputeGradient(double zDepth)
	{
		final double[] ds_dz = new double[1];
		final double[] ds_dz2 = new double[2];
		final double[] ds_duz = new double[1];
		final double[] ds_dlz = new double[1];
		final boolean record = TestSettings.allow(LogLevel.INFO);
		for (double z = -0.5; z < 0.5; z += 0.01)
		{
			final double s0 = QuadraticAstigmatismZModel.getS(z, zDepth);
			final double s1 = QuadraticAstigmatismZModel.getS1(z, zDepth, ds_dz);
			final double s2 = QuadraticAstigmatismZModel.getS2(z, zDepth, ds_dz2);

			Assert.assertEquals(s0, s1, 0);
			Assert.assertEquals(s0, s2, 0);
			Assert.assertEquals(ds_dz[0], ds_dz2[0], 0);

			final double uz = z + h_;
			final double lz = z - h_;
			final double upper = QuadraticAstigmatismZModel.getS1(uz, zDepth, ds_duz);
			final double lower = QuadraticAstigmatismZModel.getS1(lz, zDepth, ds_dlz);

			final double e1 = (upper - lower) / (uz - lz);
			final double o1 = ds_dz[0];

			// Second gradient
			final double e2 = (ds_duz[0] - ds_dlz[0]) / (uz - lz);
			final double o2 = ds_dz2[1];

			if (record)
				TestLog.info("z=%f s=%f : ds_dz=%g  %g  (%g): d2s_dz2=%g   %g  (%g)\n", z, s0, e1, o1,
						DoubleEquality.relativeError(o1, e1), e2, o2, DoubleEquality.relativeError(o2, e2));

			//double error = DoubleEquality.relativeError(o, e);
			if (Math.abs(z) > 0.02)
				Assert.assertTrue(e1 + " sign != " + o1, (e1 * o1) >= 0);
			Assert.assertTrue(e1 + " != " + o1, eq.almostEqualRelativeOrAbsolute(e1, o1));

			if (Math.abs(z) > 0.02)
				Assert.assertTrue(e2 + " sign != " + o2, (e2 * o2) >= 0);
			Assert.assertTrue(e2 + " != " + o2, eq.almostEqualRelativeOrAbsolute(e2, o2));
		}
	}
}
