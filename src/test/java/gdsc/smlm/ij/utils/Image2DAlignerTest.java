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
package gdsc.smlm.ij.utils;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.function.StandardFloatValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class Image2DAlignerTest
{
	private FloatImage2D createData(int x, int y, double cx, double cy)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, x, y, GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE,
				null);
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = 1;
		a[Gaussian2DFunction.X_POSITION] = cx;
		a[Gaussian2DFunction.Y_POSITION] = cy;
		a[Gaussian2DFunction.X_SD] = 1.2;
		a[Gaussian2DFunction.Y_SD] = 1.1;
		StandardFloatValueProcedure p = new StandardFloatValueProcedure();
		return new FloatImage2D(x, y, p.getValues(f, a));
	}

	@Test
	public void canCorrelateSquareImage()
	{
		canCorrelate(16, 16, false);
	}

	@Test
	public void canCorrelateNonSquareImage()
	{
		canCorrelate(16, 32, false);
	}

	@Test
	public void canCorrelateNonPow2Image()
	{
		canCorrelate(17, 29, false);
	}

	@Test
	public void canCorrelateSquareImageUsingIJImage()
	{
		canCorrelate(16, 16, true);
	}

	@Test
	public void canCorrelateNonSquareImageUsingIJImage()
	{
		canCorrelate(16, 32, true);
	}

	@Test
	public void canCorrelateNonPow2ImageUsingIJImage()
	{
		canCorrelate(17, 29, true);
	}

	private void canCorrelate(int maxx, int maxy, boolean ijMode)
	{
		double cx = (maxx - 1) / 2.0;
		double cy = (maxy - 1) / 2.0;

		double[] shift = new double[] { -1, -0.5, 0, 1.5, 2 };

		Image2D reference = createData(maxx, maxy, cx, cy);
		Image2DAligner a = new Image2DAligner();
		if (ijMode)
			a.setReference(reference.getImageProcessor());
		else
			a.setReference(reference);

		for (double sy : shift)
			for (double sx : shift)
			{
				canCorrelate(a, maxx, maxy, cx, cy, cx + sx, cy + sy, 0.25, 0, 1, ijMode);
				canCorrelate(a, maxx, maxy, cx, cy, cx + sx, cy + sy, 0.25, 5, 0.25, ijMode);
			}
	}

	private void canCorrelate(Image2DAligner a, int maxx, int maxy, double cx1, double cy1, double cx2, double cy2,
			double window, int refinements, double tolerance, boolean ijMode)
	{
		double[] result;
		Image2D c;
		Image2D target = createData(maxx, maxy, cx2, cy2);

		//Utils.display("Ref", reference.getImageProcessor());
		//Utils.display("Tar", target.getImageProcessor());

		a.setEdgeWindow(window);
		//a.setSearchMode(gdsc.smlm.ij.utils.StackAligner.SearchMode.BINARY);
		//a.setRelativeThreshold(1e-6);

		double[] e = new double[] { cx1 - cx2, cy1 - cy2 };
		int cx = maxx / 2 - (int) Math.round(e[0]);
		int cy = maxy / 2 - (int) Math.round(e[1]);
		int index = target.getIndex(cx, cy);

		// Debug the convergence
		//for (int i = 0; i <= refinements; i++)
		//{
		//	result = a.align(target.copy(), i, error);
		//	c = a.getCorrelation();
		//
		//	System.out.printf("e %s %g, o %s\n", java.util.Arrays.toString(e), c.get(index),
		//			java.util.Arrays.toString(result));
		//}

		// Test
		if (ijMode)
			result = a.align(target.getImageProcessor(), refinements);
		else
			result = a.align(target, refinements);
		c = a.getCorrelation();
		System.out.printf("e %s %g, o %s\n", java.util.Arrays.toString(e), c.get(index),
				java.util.Arrays.toString(result));

		for (int i = 0; i < 2; i++)
			Assert.assertEquals(e[i], result[i], tolerance);
	}
}
