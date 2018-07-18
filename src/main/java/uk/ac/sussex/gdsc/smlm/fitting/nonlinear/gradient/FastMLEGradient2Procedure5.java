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
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;

/**
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second partial derivatives.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMLEGradient2Procedure5 extends FastMLEGradient2Procedure
{
	/**
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param func
	 *            Gradient function (must produce a strictly positive value, i.e. the mean of a Poisson process)
	 */
	public FastMLEGradient2Procedure5(final double[] x, final Gradient2Function func)
	{
		super(x, func);
		if (n != 5)
			throw new IllegalArgumentException("Function must compute 5 gradients");
	}

	@Override
	protected void reset2()
	{
		d1[0] = 0;
		d1[1] = 0;
		d1[2] = 0;
		d1[3] = 0;
		d1[4] = 0;
		d2[0] = 0;
		d2[1] = 0;
		d2[2] = 0;
		d2[3] = 0;
		d2[4] = 0;
	}

	@Override
	public void execute(double uk, double[] duk_dt, double[] d2uk_dt2)
	{
		u[k] = uk;
		final double xk = x[k++];
		final double xk_uk_minus1 = xk / uk - 1.0;
		final double xk_uk2 = xk / (uk * uk);
		d1[0] += duk_dt[0] * xk_uk_minus1;
		d1[1] += duk_dt[1] * xk_uk_minus1;
		d1[2] += duk_dt[2] * xk_uk_minus1;
		d1[3] += duk_dt[3] * xk_uk_minus1;
		d1[4] += duk_dt[4] * xk_uk_minus1;
		d2[0] += d2uk_dt2[0] * xk_uk_minus1 - duk_dt[0] * duk_dt[0] * xk_uk2;
		d2[1] += d2uk_dt2[1] * xk_uk_minus1 - duk_dt[1] * duk_dt[1] * xk_uk2;
		d2[2] += d2uk_dt2[2] * xk_uk_minus1 - duk_dt[2] * duk_dt[2] * xk_uk2;
		d2[3] += d2uk_dt2[3] * xk_uk_minus1 - duk_dt[3] * duk_dt[3] * xk_uk2;
		d2[4] += d2uk_dt2[4] * xk_uk_minus1 - duk_dt[4] * duk_dt[4] * xk_uk2;
	}

	/**
	 * Reset the first derivative vector
	 */
	@Override
	protected void reset1()
	{
		d1[0] = 0;
		d1[1] = 0;
		d1[2] = 0;
		d1[3] = 0;
		d1[4] = 0;
	}

	@Override
	public void execute(double uk, double[] duk_dt)
	{
		u[k] = uk;
		final double xk = x[k++];
		final double xk_uk_minus1 = xk / uk - 1.0;
		d1[0] += duk_dt[0] * xk_uk_minus1;
		d1[1] += duk_dt[1] * xk_uk_minus1;
		d1[2] += duk_dt[2] * xk_uk_minus1;
		d1[3] += duk_dt[3] * xk_uk_minus1;
		d1[4] += duk_dt[4] * xk_uk_minus1;
	}
}
