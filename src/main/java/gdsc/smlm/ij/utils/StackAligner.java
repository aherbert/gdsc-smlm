package gdsc.smlm.ij.utils;

import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PositionChecker;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.BFGSOptimizer;

import gdsc.core.math.interpolation.CubicSplinePosition;
import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.cspline.CubicSplineCalculator;
import ij.ImageStack;
import pl.edu.icm.jlargearrays.LargeArray;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Perform 3D image alignment using cross-correlation
 */
public class StackAligner implements Cloneable
{
	private double edgeWindow;

	/** The number of slices (max z) of the discrete Hartley transform. */
	private int ns;
	/** The number of rows (max y) of the discrete Hartley transform. */
	private int nr;
	/** The number of columns (max x) of the discrete Hartley transform. */
	private int nc;

	private DHT3D reference;
	private float[] buffer, region;

	// Allow cached window weights
	private double[] wx = null;
	private double[] wy = null;
	private double[] wz = null;

	private CubicSplineCalculator calc;

	/**
	 * Instantiates a new stack aligner with a default edge window of 0.25
	 */
	public StackAligner()
	{
		this(0.25);
	}

	/**
	 * Instantiates a new stack aligner.
	 *
	 * @param edgeWindow
	 *            the alpha value for the Tukey edge window
	 */
	public StackAligner(double edgeWindow)
	{
		this.edgeWindow = edgeWindow;
	}

	/**
	 * Sets the reference stack and assumes the target stack will be the same size.
	 * <p>
	 * The dimension are converted to the next power of 2 for speed. The combined size must fit within the maximum size
	 * of a single array.
	 *
	 * @param stack
	 *            the stack (destructively modified)
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined target dimensions is too large for an array
	 */
	public void setReference(ImageStack stack)
	{
		setReference(stack, stack.getWidth(), stack.getHeight(), stack.getSize());
	}

	/**
	 * Sets the reference stack and the size of the target stack.
	 * <p>
	 * The dimension are converted to the next power of 2 for speed. The combined size must fit within the maximum size
	 * of a single array.
	 *
	 * @param stack
	 *            the stack (destructively modified)
	 * @param w
	 *            the width of the target stack
	 * @param h
	 *            the height of the target stack
	 * @param d
	 *            the depth of the target stack
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined target dimensions is too large for an array
	 */
	public void setReference(ImageStack stack, int w, int h, int d)
	{
		check3D(stack);
		if (w < 2 || h < 2 || d < 2)
			throw new IllegalArgumentException("Require a 3D target stack");
		nc = Maths.nextPow2(Math.max(w, stack.getWidth()));
		nr = Maths.nextPow2(Math.max(h, stack.getHeight()));
		ns = Maths.nextPow2(Math.max(d, stack.getSize()));
		long size = (long) ns * nr * nc;
		// Don't support using large arrays for simplicity
		if (size > LargeArray.getMaxSizeOf32bitArray())
			throw new IllegalArgumentException("3D data too large");
		// Window and pad the reference
		reference = createDHT(stack);
	}

	private void check3D(ImageStack stack)
	{
		if (stack.getWidth() < 2 || stack.getHeight() < 2 || stack.getSize() < 2)
			throw new IllegalArgumentException("Require a 3D stack");
	}

	private DHT3D createDHT(ImageStack stack)
	{
		// Apply window
		int w = stack.getWidth(), h = stack.getHeight(), d = stack.getSize();
		if (edgeWindow > 0)
		{
			double[] wx = createXWindow(w);
			double[] wy = createYWindow(h);
			double[] wz = createZWindow(d);
			for (int z = 0; z < d; z++)
			{
				float[] pixels = (float[]) stack.getPixels(1 + z);
				if (wz[z] == 0)
				{
					// Special case happens with Tukey window at the ends
					Arrays.fill(pixels, 0);
				}
				else
				{
					applyWindow(pixels, w, h, wx, wy, wz[z]);
				}
			}
		}

		DHT3D dht;
		if (w < nc || h < nr || d < ns)
		{
			// Pad into the desired data size
			int size = nc * nr;
			float[] dest = new float[ns * size];

			int ix = getInsert(nc, w);
			int iy = getInsert(nr, h);
			int iz = getInsert(ns, d);
			int insertPos = iy * nc + ix;
			for (int z = 0, slice = 0; z < ns; z++)
			{
				if (z >= iz && slice < d)
				{
					// Stack uses 1-based index
					float[] source = (float[]) stack.getPixels(++slice);
					insert(source, w, dest, nc, z * size + insertPos);
				}
			}
			dht = new DHT3D(nc, nr, ns, dest, false);
		}
		else
		{
			// This will just copy the data
			dht = new DHT3D(stack);
		}

		dht.transform();
		return dht;
	}

	private double[] createXWindow(int n)
	{
		return wx = createWindow(wx, n);
	}

	private double[] createYWindow(int n)
	{
		return wy = createWindow(wy, n);
	}

	private double[] createZWindow(int n)
	{
		return wz = createWindow(wz, n);
	}

	private double[] createWindow(double[] w, int n)
	{
		if (w == null || w.length != n)
			w = ImageWindow.tukey(n, edgeWindow);
		return w;
	}

	private static void applyWindow(float[] image, int maxx, int maxy, double[] wx, double[] wy, double wz)
	{
		for (int y = 0, i = 0; y < maxy; y++)
		{
			double w = wy[y] * wz;
			for (int x = 0; x < maxx; x++, i++)
			{
				image[i] *= wx[x] * w;
			}
		}
	}

	private static int getInsert(int maxN, int n)
	{
		// Note the FHT power spectrum centre is at n/2 of an even sized image.
		// So we must insert the centre at that point. To do this we check for odd/even
		// and offset if necessary. 
		int diff = maxN - n;
		return ((diff & 1) == 1) ? (diff + 1) / 2 : diff / 2;
	}

	private static void insert(float[] source, int sw, float[] dest, int dw, int to)
	{
		for (int from = 0; from < source.length; from += sw, to += dw)
		{
			System.arraycopy(source, from, dest, to, sw);
		}
	}

	/**
	 * Align the stack with the reference. Compute the translation required to move the target stack onto the reference
	 * stack for maximum correlation.
	 *
	 * @param stack
	 *            the stack
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(ImageStack stack)
	{
		return align(stack, 0, 0);
	}

	/**
	 * Align the stack with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * stack onto the reference stack for maximum correlation.
	 *
	 * @param stack
	 *            the stack
	 * @param refinements
	 *            the refinements for sub-pixel accuracy
	 * @param error
	 *            the error for sub-pixel accuracy (i.e. stop when improvements are less than this error)
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(ImageStack stack, int refinements, double error)
	{
		check3D(stack);
		int w = stack.getWidth(), h = stack.getHeight(), d = stack.getSize();
		if (w > nc || h > nr || d > ns)
			throw new IllegalArgumentException("Stack is larger than the initialised reference");

		DHT3D target = createDHT(stack);
		DHT3D correlation = target.conjugateMultiply(reference, buffer);
		buffer = correlation.getData();
		correlation.inverseTransform();
		correlation.swapOctants();
		float[] data = correlation.getData();
		int maxi = SimpleArrayUtils.findMaxIndex(data);
		int[] xyz = correlation.getXYZ(maxi);

		// Report the shift required to move from the centre of the target image to the reference
		// @formatter:off
		double[] result = new double[] {
			nc/2 - xyz[0],
			nr/2 - xyz[1],
			ns/2 - xyz[2],
			data[maxi]
		};
		// @formatter:on

		if (refinements > 0)
		{
			// Perform sub-pixel alignment
			// Create a cubic spline using a small region of pixels around the maximum
			if (calc == null)
				calc = new CubicSplineCalculator();
			// Avoid out-of-bounds errors
			int x = Math.max(0, xyz[0] - 1);
			int y = Math.max(0, xyz[1] - 1);
			int z = Math.max(0, xyz[2] - 1);
			Image3D crop = correlation.crop(x, y, z, 4, 4, 4, region);
			region = crop.getData();
			CustomTricubicFunction f = CustomTricubicFunction.create(calc.compute(region));

			// Find the maximum starting at the current origin
			int ox = xyz[0] - x;
			int oy = xyz[1] - y;
			int oz = xyz[2] - z;

			// Scale to the cubic spline dimensions of 0-1
			double[] origin = new double[] { ox / 3.0, oy / 3.0, oz / 3.0 };

			try
			{
				final SplineFunction sf = new SplineFunction(f, origin);

				// @formatter:off				
				// Scale the error for the position check
				BFGSOptimizer optimiser = new BFGSOptimizer(
						new PositionChecker(-1, error / 3.0, refinements));

				PointValuePair opt = optimiser.optimize(
						GoalType.MAXIMIZE,
						maxEvaluations,
						bounds, 
						gradientTolerance,
						stepLength,
						new InitialGuess(origin),
						new ObjectiveFunction(new MultivariateFunction(){
							public double value(double[] point)
							{
								return sf.value(point);
							}}),
						new ObjectiveFunctionGradient(new MultivariateVectorFunction(){
							public double[] value(double[] point) throws IllegalArgumentException
							{
								// This must be new each time
								double[] df_da = new double[3];
								sf.value(point, df_da);
								return df_da;
							}}));
				// @formatter:on

				// Check it is higher
				if (opt.getValue() > result[3])
				{
					result[3] = opt.getValue();
					// Convert the maximum back with scaling
					double[] optimum = opt.getPointRef();
					for (int i = 0; i < 3; i++)
						result[i] -= (optimum[i] - origin[i]) * 3.0;
					return result;
				}
			}
			catch (Exception e)
			{
				// Ignore this
			}
		}

		return result;
	}

	// For optimisation

	/** Do not have a maximum evaluations as we will converge using the refinements parameter. */
	private static final MaxEval maxEvaluations = new MaxEval(Integer.MAX_VALUE);

	/** The bounds of the spline are always 0-1 */
	private static final SimpleBounds bounds = new SimpleBounds(new double[3], SimpleArrayUtils.newDoubleArray(3, 1));

	/** Set a maximum step length at 1 pixel scaled to the spline dimensions. */
	private static final BFGSOptimizer.StepLength stepLength = new BFGSOptimizer.StepLength(
			SimpleArrayUtils.newDoubleArray(3, 1.0 / 3));
	/**
	 * This is the cut-off for the maximum gradient relative to the function value. When gradients are too small then
	 * the optimisation will end.
	 */
	private static final BFGSOptimizer.GradientTolerance gradientTolerance = new BFGSOptimizer.GradientTolerance(1e-6);

	private static class SplineFunction
	{
		final CustomTricubicFunction f;
		CubicSplinePosition[] s = new CubicSplinePosition[3];
		double[] table;

		SplineFunction(CustomTricubicFunction f, double[] origin)
		{
			this.f = f;
			for (int i = 0; i < 3; i++)
				s[i] = new CubicSplinePosition(origin[i]);
		}

		double value(double[] point)
		{
			initialise(point);
			return f.value(table);
		}

		void value(double[] point, double[] df_da)
		{
			initialise(point);
			f.gradient(table, df_da);
		}

		void initialise(double[] point)
		{
			// Allow caching the the spline positions and the table
			for (int i = 0; i < 3; i++)
			{
				if (s[i].getX() != point[i])
				{
					s[i] = new CubicSplinePosition(point[i]);
					table = null;
				}
			}
			if (table == null)
				table = CustomTricubicFunction.computePowerTable(s[0], s[1], s[2]);
		}
	}

	/**
	 * Copy the aligner. This copies the initialised state for use in alignment on multiple threads concurrently.
	 *
	 * @return the stack aligner
	 */
	public StackAligner copy()
	{
		StackAligner copy;
		try
		{
			copy = (StackAligner) clone();
			// Reset objects that are not thread safe
			copy.calc = null;
			copy.buffer = null;
			copy.region = null;
			return copy;
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}

	/**
	 * Gets the edge window.
	 *
	 * @return the edge window
	 */
	public double getEdgeWindow()
	{
		return edgeWindow;
	}

	/**
	 * Sets the edge window.
	 *
	 * @param edgeWindow
	 *            the new edge window
	 */
	public void setEdgeWindow(double edgeWindow)
	{
		this.edgeWindow = edgeWindow;
	}
}