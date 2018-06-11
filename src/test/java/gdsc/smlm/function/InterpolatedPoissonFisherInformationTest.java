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

import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;

public class InterpolatedPoissonFisherInformationTest
{
	@Test
	public void canInterpolateFisherInformation()
	{
		for (int i = 0; i < 4; i++)
		{
			double s = (1 << i) * 0.25;
			canInterpolateFisherInformation(s);
		}
	}

	private void canInterpolateFisherInformation(double s)
	{
		canComputeFisherInformation(new PoissonGaussianApproximationFisherInformation(s));
	}

	private void canComputeFisherInformation(BasePoissonFisherInformation f)
	{
		// Build a range for the Fisher information
		int min = -100;
		int max = 20;
		double[] logU = new double[max - min + 1];
		double[] alpha = new double[logU.length];

		for (int exp = min, i = 0; exp <= max; exp++, i++)
		{
			logU[i] = exp;
			alpha[i] = f.getAlpha(FastMath.exp(exp));
		}

		InterpolatedPoissonFisherInformation fi = new InterpolatedPoissonFisherInformation(logU, alpha, false, f);

		// No check for beyond the range since a separate test does this.
		check(f, fi, min, 1e-2);

		// Within
		for (int exp = min; exp < max; exp++)
		{
			check(f, fi, exp + 0.5, 1e-2);
		}

		// Upper bound
		check(f, fi, max, 1e-2);
		
		// Beyond upper bound
		check(f, fi, max + 1, 1e-2);
	}

	private void check(BasePoissonFisherInformation f, InterpolatedPoissonFisherInformation fi, double logU, double tol)
	{
		double u = FastMath.exp(logU);
		double e = f.getAlpha(u);
		double o = fi.getAlpha(u);
		//System.out.printf("logU=%g  u=%g  e=%g  o=%g  error=%g\n", logU, u, e, o, DoubleEquality.relativeError(o, e));

		// Small numbers may have a large relative error but the absolute error is small
		Assert.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e, o, 5e-3, 1e-20));
	}

	@Test
	public void canInterpolateLowerFisherInformation()
	{
		for (int i = 0; i < 4; i++)
		{
			double s = (1 << i) * 0.25;
			canInterpolateLowerFisherInformation(s);
		}
	}

	private void canInterpolateLowerFisherInformation(double s)
	{
		canComputeLowerFisherInformation(new PoissonGaussianApproximationFisherInformation(s));
	}

	private void canComputeLowerFisherInformation(BasePoissonFisherInformation f)
	{
		// Build a range for the Fisher information where it should plateau
		int min = -500;
		int max = min + 4;
		double[] logU = new double[max - min + 1];
		double[] alpha = new double[logU.length];

		for (int exp = min, i = 0; exp <= max; exp++, i++)
		{
			logU[i] = exp;
			alpha[i] = f.getAlpha(FastMath.exp(exp));
		}

		// Lower fixed I
		InterpolatedPoissonFisherInformation fi = new InterpolatedPoissonFisherInformation(logU, alpha, true, f);
		
		final double I = f.getFisherInformation(FastMath.exp(min));
		BasePoissonFisherInformation fixedI = new BasePoissonFisherInformation()
		{
			public double getFisherInformation(double t) throws IllegalArgumentException
			{
				return I;
			}

			@Override
			public double getAlpha(double t)
			{
				return t * I;
			}

			@Override
			protected void postClone()
			{
			}
		};
		check(fixedI, fi, min - 1, 0);
		check(fixedI, fi, min - 2, 0);
		check(fixedI, fi, min - 20, 0);

		// Lower fixed alpha
		fi = new InterpolatedPoissonFisherInformation(logU, alpha, false, f);
		
		final double A = alpha[0]; 
		BasePoissonFisherInformation fixedA = new BasePoissonFisherInformation()
		{
			public double getFisherInformation(double t) throws IllegalArgumentException
			{
				return t / A;
			}

			@Override
			public double getAlpha(double t)
			{
				return A;
			}

			@Override
			protected void postClone()
			{
			}
		};

		check(fixedA, fi, min - 1, 0);
		check(fixedA, fi, min - 2, 0);
		check(fixedA, fi, min - 20, 0);
	}
}
