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
package org.apache.commons.math3.analysis.integration;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Maths;
import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class CustomSimpsonIntegratorTest
{
	private interface TestUnivariateFunction extends UnivariateFunction
	{
		public double sum(double a, double b);
	}

	private class LinearTestUnivariateFunction implements TestUnivariateFunction
	{
		@Override
		public double value(double x)
		{
			return x;
		}

		@Override
		public double sum(double a, double b)
		{
			// y=x => x^2/2
			return (b * b - a * a) / 2;
		}
	}

	private class QuadraticTestUnivariateFunction implements TestUnivariateFunction
	{
		@Override
		public double value(double x)
		{
			return x * x - x;
		}

		@Override
		public double sum(double a, double b)
		{
			// y=x^2 - x => x^3/3 - x^2 / 2
			double u = Maths.pow3(b) / 3 - Maths.pow2(b) / 2;
			double l = Maths.pow3(a) / 3 - Maths.pow2(a) / 2;
			return u - l;
		}
	}

	private class CubicTestUnivariateFunction implements TestUnivariateFunction
	{
		@Override
		public double value(double x)
		{
			return x * x * x - x * x;
		}

		@Override
		public double sum(double a, double b)
		{
			// y=x^3 - x^2 => x^4/4 - x^3 / 3
			double u = Maths.pow4(b) / 4 - Maths.pow3(b) / 3;
			double l = Maths.pow4(a) / 4 - Maths.pow3(a) / 3;
			return u - l;
		}
	}

	private class ReciprocalTestUnivariateFunction implements TestUnivariateFunction
	{
		@Override
		public double value(double x)
		{
			return 1.0 / x;
		}

		@Override
		public double sum(double a, double b)
		{
			// y=1/x => log(x)
			return Math.log(b) - Math.log(a);
		}
	}

	@Test
	public void canIntegrateFunction()
	{
		TestUnivariateFunction f;

		f = new LinearTestUnivariateFunction();
		canIntegrateFunction(f, 0.5, 2, 1);
		canIntegrateFunction(f, 0.5, 2, 2);
		canIntegrateFunction(f, 0.5, 2, 3);

		f = new QuadraticTestUnivariateFunction();
		canIntegrateFunction(f, 0.5, 2, 1);
		canIntegrateFunction(f, 0.5, 2, 2);
		canIntegrateFunction(f, 0.5, 2, 3);

		f = new CubicTestUnivariateFunction();
		canIntegrateFunction(f, 0.5, 2, 1);
		canIntegrateFunction(f, 0.5, 2, 2);
		canIntegrateFunction(f, 0.5, 2, 3);

		// Harder so use more iterations
		f = new ReciprocalTestUnivariateFunction();
		canIntegrateFunction(f, 0.5, 2, 5);
	}

	private void canIntegrateFunction(TestUnivariateFunction f, double a, double b, int c)
	{
		final double relativeAccuracy = 1e-4;
		// So it stops at the min iterations
		final double absoluteAccuracy = Double.POSITIVE_INFINITY;

		CustomSimpsonIntegrator in = new CustomSimpsonIntegrator(
				//new SimpsonIntegrator(
				relativeAccuracy, absoluteAccuracy, c, CustomSimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);

		double e = f.sum(a, b);
		double ee = simpson(f, a, b, c);
		double o = in.integrate(Integer.MAX_VALUE, f, a, b);

		TestSettings.info("%s c=%d  %g-%g  e=%g  ee=%g  o=%g\n", f.getClass().getSimpleName(), c, a, b, e, ee, o);

		double delta = Math.abs(e) * 1e-6;
		Assert.assertEquals(e, ee, delta);
		Assert.assertEquals(e, o, delta);

		// These should be the same
		Assert.assertEquals(ee, o, delta * 1e-6);
	}

	private double simpson(UnivariateFunction f, double a, double b, int c)
	{
		// Simple Simpson integration:
		// https://en.wikipedia.org/wiki/Simpson%27s_rule

		// Number of sub intervals
		int n = 1 << c + 1;
		double h = (b - a) / n; //sub-interval width 
		double sum2 = 0, sum4 = 0;
		for (int j = 1; j <= n / 2 - 1; j++)
			sum2 += f.value(a + (2 * j) * h);
		for (int j = 1; j <= n / 2; j++)
			sum4 += f.value(a + (2 * j - 1) * h);
		double sum = (h / 3) * (f.value(a) + 2 * sum2 + 4 * sum4 + f.value(b));
		return sum;
	}
}
