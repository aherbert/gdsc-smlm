package gdsc.smlm.ij.utils;

import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PositionChecker;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.BFGSOptimizer;

import gdsc.core.ij.Utils;
import gdsc.core.math.interpolation.CubicSplinePosition;
import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.cspline.CubicSplineCalculator;
import ij.ImageStack;

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
 * Perform 3D image alignment using normalised cross-correlation.
 * <p>
 * Uses the following formula:
 * 
 * <pre>
 *  ( Σ xiyi - nx̄ӯ ) / ( (Σ xi^2 - nx̄^2) (Σ yi^2 - nӯ^2) )^0.5
 * </pre>
 * 
 * The summation in the numerator is computed using a conjugate multiplication in the frequency domain. The summation
 * terms are computed using rolling sum tables. Images are converted to the full range of an unsigned 16-bit integer
 * before computation to avoid errors in the rolling sum tables. This should have minimal impact on the
 * correlation value since it is normalised.
 * 
 * @see <a href=
 *      "https://en.wikipedia.org/wiki/Pearson_correlation_coefficient">https://en.wikipedia.org/wiki/Pearson_correlation_coefficient</a>
 * @see <a href="http://scribblethink.org/Work/nvisionInterface/nip.html">Fast Normalized Cross-Correlation by J.P.
 *      Lewis</a>
 */
public class StackAligner implements Cloneable
{
	private static final int X = 0;
	private static final int XX = 1;
	private static final int Y = 0;
	private static final int YY = 1;

	private static class DHTData
	{
		DHT3D dht;
		long[] s_;
		long[] ss;
		// Original dimensions
		int w, h, d;
		// Insert position
		int ix, iy, iz;

		DHTData(DHT3D dht, int w, int h, int d)
		{
			setDHT(dht, w, h, d);
		}

		void setDHT(DHT3D dht, int w, int h, int d)
		{
			this.dht = dht;
			s_ = resize(s_);
			ss = resize(ss);
			this.w = w;
			this.h = h;
			this.d = d;
			ix = getInsert(dht.nc, w);
			iy = getInsert(dht.nr, h);
			iz = getInsert(dht.ns, d);
		}

		private long[] resize(long[] data)
		{
			return (data == null || data.length != dht.getDataLength()) ? new long[dht.getDataLength()] : data;
		}
	}

	/**
	 * The search mode for sub-pixel refinement.
	 */
	public enum SearchMode
	{
		/**
		 * Perform a binary search by condensing the cube vertices around the highest value of a tricubic interpolation
		 */
		BINARY,
		/** Use the local gradient of a tricubic interpolation to find the maximum. */
		GRADIENT
	}

	private double edgeWindow;
	private double relativeThreshold = 1e-3;
	public SearchMode searchMode = SearchMode.GRADIENT;

	/** The number of slices (max z) of the discrete Hartley transform. */
	private int ns;
	/** The number of rows (max y) of the discrete Hartley transform. */
	private int nr;
	/** The number of columns (max x) of the discrete Hartley transform. */
	private int nc;
	/** The number of rows by columns of the discrete Hartley transform. */
	private int nr_by_nc;

	private DHTData reference;

	// Not thread safe as they are used for the target image
	private DHTData target;
	private float[] buffer, region;

	// Allow cached window weights
	private double[] wx = null;
	private double[] wy = null;
	private double[] wz = null;

	private CubicSplineCalculator calc;
	private Image3D ref, tar;

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
	 *            the stack (may be destructively modified)
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
		// Check the stack will fit in an Image3D
		Image3D.checkSize(nc, nr, ns, true);
		// Window and pad the reference
		setReference(createDHT(stack, reference));
	}

	private void setReference(DHTData dhtData)
	{
		reference = dhtData;
		// We could pre-compute the mean and sum-of-squares for each overlap position here
		// But that would require more storage again
	}

	private void check3D(ImageStack stack)
	{
		if (stack.getWidth() < 2 || stack.getHeight() < 2 || stack.getSize() < 2)
			throw new IllegalArgumentException("Require a 3D stack");
	}

	private DHTData createDHT(ImageStack stack, DHTData dhtData)
	{
		if (stack.getBitDepth() != 32)
			return createDHT(new Image3D(stack), dhtData);

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
		//if (w < nc || h < nr || d < ns)
		//{

		// Pad into the desired data size.
		// We always do this so the data is reused
		int size = ns * nr * nc;
		float[] dest;
		if (dhtData == null || dhtData.dht.getDataLength() != size)
		{
			dest = new float[size];
		}
		else
		{
			// Re-use space
			dest = dhtData.dht.getData();
			Arrays.fill(dest, 0f);
		}
		dht = new DHT3D(nc, nr, ns, dest, false);
		int ix = getInsert(nc, w);
		int iy = getInsert(nr, h);
		int iz = getInsert(ns, d);
		dht.insert(ix, iy, iz, stack);

		//}
		//else
		//{
		//	// This will just copy the data
		//	dht = new DHT3D(stack);
		//}

		if (dhtData == null)
			dhtData = new DHTData(dht, w, h, d);
		else
			dhtData.setDHT(dht, w, h, d);

		return prepareDHT(dhtData);
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
		applyWindow(image, 0, maxx, maxy, wx, wy, wz);
	}

	private static void applyWindow(float[] image, int i, int maxx, int maxy, double[] wx, double[] wy, double wz)
	{
		for (int y = 0; y < maxy; y++)
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

	/**
	 * Prepare the DHT.
	 * <p>
	 * Converts the data to the full range of a 16-bit signed integer with a mean of zero.
	 * This may reduce the precision slightly but allows the computation of a rolling sum table with no errors. The
	 * rolling sum and sum-of-squares table is computed and the DHT is transformed to the frequency domain.
	 *
	 * @param dhtData
	 *            the dht data
	 * @return the DHT data
	 */
	private DHTData prepareDHT(DHTData dhtData)
	{
		DHT3D dht = dhtData.dht;
		long[] s_ = dhtData.s_;
		long[] ss = dhtData.ss;

		// XXX Update this
		
		// Convert to a 16-bit signed integer
		// This makes it possible to compute the rolling tables without error
		float[] data = dht.getData();
		float[] limits = Maths.limits(data);
		double min = limits[0];
		double max = limits[1];
		if ((max - min) == 0.0)
		{
			// No data variation so just zero fill!
			Arrays.fill(s_, 0);
			Arrays.fill(ss, 0);
		}
		else
		{
			double scale = LIMIT / (max - min);

			// Compute the rolling sum tables
			int nr_by_nc = dht.nr_by_nc;
			int nc = dht.nc;
			int nr = dht.nr;
			int ns = dht.ns;

			// This has been adapted from Image3D to compute two rolling sum table at once

			// First build a table for each XY slice
			for (int s = 0; s < ns; s++)
			{
				long sum = 0, sum2 = 0;
				int i = s * nr_by_nc;
				// Initialise first row sum
				// sum = rolling sum of (0 - colomn)
				for (int c = 0; c < nc; c++, i++)
				{
					int v = transform(data[i], min, scale);
					data[i] = v;
					sum += v;
					sum2 += v * v;
					s_[i] = sum;
					ss[i] = sum2;
				}
				// Remaining rows
				// sum = rolling sum of (0 - colomn) + sum of same position above
				for (int r = 1, ii = i - nc; r < nr; r++)
				{
					sum = 0;
					sum2 = 0;
					for (int c = 0; c < nc; c++, i++, ii++)
					{
						int v = transform(data[i], min, scale);
						data[i] = v;
						sum += v;
						sum2 += v * v;
						// Add the sum from the previous row
						s_[i] = sum + s_[ii];
						ss[i] = sum2 + ss[ii];
					}
				}
			}

			// Now sum across slices
			// sum = rolling sum of (0,0 to column,row) + sum of same position above
			// => rolling sum of (0,0,0 to column,row,slice)
			for (int s = 1; s < ns; s++)
			{
				int i = s * nr_by_nc;
				int ii = i - nr_by_nc;
				for (int j = 0; j < nr_by_nc; j++, i++, ii++)
				{
					s_[i] += s_[ii];
					ss[i] += ss[ii];
				}
			}
		}
		if (ref == null)
			ref = new Image3D(nc, nr, ns, data.clone());
		else
			tar = new Image3D(nc, nr, ns, data.clone());
		// Transform the data
		dht.transform();
		return dhtData;
	}

	private static double LIMIT = 1023; //65535.0;

	private static int transform(float f, double min, double scale)
	{
		double value = (f - min) * scale;
		if (value < 0.0)
			return 0;
		if (value > LIMIT)
			return (int) LIMIT;
		return (int) (value + 0.5);
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
	public void setReference(Image3D stack)
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
	 *            the stack (may be destructively modified)
	 * @param w
	 *            the width of the target stack
	 * @param h
	 *            the height of the target stack
	 * @param d
	 *            the depth of the target stack
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined target dimensions is too large for an array
	 */
	public void setReference(Image3D stack, int w, int h, int d)
	{
		check3D(stack);
		if (w < 2 || h < 2 || d < 2)
			throw new IllegalArgumentException("Require a 3D target stack");
		nc = Maths.nextPow2(Math.max(w, stack.getWidth()));
		nr = Maths.nextPow2(Math.max(h, stack.getHeight()));
		ns = Maths.nextPow2(Math.max(d, stack.getSize()));
		nr_by_nc = nr * nc;
		// Window and pad the reference
		setReference(createDHT(stack, reference));
	}

	private void check3D(Image3D stack)
	{
		if (stack.getWidth() < 2 || stack.getHeight() < 2 || stack.getSize() < 2)
			throw new IllegalArgumentException("Require a 3D stack");
	}

	private DHTData createDHT(Image3D stack, DHTData dhtData)
	{
		// Apply window
		int w = stack.getWidth(), h = stack.getHeight(), d = stack.getSize();
		if (edgeWindow > 0)
		{
			double[] wx = createXWindow(w);
			double[] wy = createYWindow(h);
			double[] wz = createZWindow(d);
			float[] pixels = stack.getData();
			int inc = stack.nr_by_nc;
			for (int z = 0, i = 0; z < d; z++)
			{
				if (wz[z] == 0)
				{
					// Special case happens with Tukey window at the ends
					for (int j = 0; j < inc; j++)
						pixels[i++] = 0;
				}
				else
				{
					applyWindow(pixels, i, w, h, wx, wy, wz[z]);
					i += inc;
				}
			}
		}

		DHT3D dht;
		if (w < nc || h < nr || d < ns)
		{
			// Pad into the desired data size
			int size = ns * nr * nc;
			float[] dest;
			if (dhtData == null || dhtData.dht.getDataLength() != size)
			{
				dest = new float[size];
			}
			else
			{
				// Re-use space
				dest = dhtData.dht.getData();
				Arrays.fill(dest, 0f);
			}
			dht = new DHT3D(nc, nr, ns, dest, false);
			int ix = getInsert(nc, w);
			int iy = getInsert(nr, h);
			int iz = getInsert(ns, d);
			dht.insert(ix, iy, iz, stack);
		}
		else
		{
			// This will just use the data
			dht = new DHT3D(w, h, d, stack.getData(), false);
		}

		if (dhtData == null)
			dhtData = new DHTData(dht, w, h, d);
		else
			dhtData.setDHT(dht, w, h, d);

		return prepareDHT(dhtData);
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

		target = createDHT(stack, target);
		return align(target, refinements, error);
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
	public double[] align(Image3D stack)
	{
		return align(stack, 0, 0);
	}

	/**
	 * Align the stack with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * stack onto the reference stack for maximum correlation.
	 * <p>
	 * Refinement uses a default sub-pixel accuracy of 1e-2;
	 *
	 * @param stack
	 *            the stack
	 * @param refinements
	 *            the refinements for sub-pixel accuracy
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(Image3D stack, int refinements)
	{
		check3D(stack);
		int w = stack.getWidth(), h = stack.getHeight(), d = stack.getSize();
		if (w > nc || h > nr || d > ns)
			throw new IllegalArgumentException("Stack is larger than the initialised reference");

		target = createDHT(stack, target);
		return align(target, refinements, 1e-2);
	}

	/**
	 * Align the stack with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * stack onto the reference stack for maximum correlation.
	 *
	 * @param stack
	 *            the stack
	 * @param refinements
	 *            the maximum number of refinements for sub-pixel accuracy
	 * @param error
	 *            the error for sub-pixel accuracy (i.e. stop when improvements are less than this error)
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(Image3D stack, int refinements, double error)
	{
		check3D(stack);
		int w = stack.getWidth(), h = stack.getHeight(), d = stack.getSize();
		if (w > nc || h > nr || d > ns)
			throw new IllegalArgumentException("Stack is larger than the initialised reference");

		target = createDHT(stack, target);
		return align(target, refinements, error);
	}

	/**
	 * Align the stack with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * stack onto the reference stack for maximum correlation.
	 *
	 * @param target
	 *            the target
	 * @param refinements
	 *            the maximum number of refinements for sub-pixel accuracy
	 * @param error
	 *            the error for sub-pixel accuracy (i.e. stop when improvements are less than this error)
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	private double[] align(DHTData target, int refinements, double error)
	{
		// Multiply by the reference. This allows the reference to be shared across threads.
		DHT3D correlation = target.dht.conjugateMultiply(reference.dht, buffer);
		buffer = correlation.getData(); // Store for reuse
		correlation.inverseTransform();
		correlation.swapOctants();

		// Normalise:
		//  ( Σ xiyi - nx̄ӯ ) / ( (Σ xi^2 - nx̄^2) (Σ yi^2 - nӯ^2) )^0.5
		// 
		// (sumXY - sumX*sumY/n) / sqrt( (sumXX - sumX^2 / n) * (sumYY - sumY^2 / n) )

		// Only do this over the range where at least half the original images overlap.
		int ix = Math.max(reference.ix, target.ix);
		int iy = Math.max(reference.iy, target.iy);
		int iz = Math.max(reference.iz, target.iz);
		int iw = Math.max(reference.ix + reference.w, target.ix + target.w);
		int ih = Math.max(reference.iy + reference.h, target.iy + target.h);
		int id = Math.max(reference.iz + reference.d, target.iz + target.d);

		float[] data = correlation.getData();

		// Compute sum from rolling sum using:
		// sum(x,y,z,w,h,d) = 
		// + s(x+w-1,y+h-1,z+d-1) 
		// - s(x-1,y+h-1,z+d-1)
		// - s(x+w-1,y-1,z+d-1)
		// + s(x-1,y-1,z+d-1)
		// /* Stack above must be subtracted so reverse sign*/
		// - s(x+w-1,y+h-1,z-1) 
		// + s(x-1,y+h-1,z-1)
		// + s(x+w-1,y-1,z-1)
		// - s(x-1,y-1,z-1)
		// Note: 
		// s(i,j,k) = 0 when either i,j,k < 0
		// i = imax when i>imax 
		// j = jmax when j>jmax 
		// k = kmax when k>kmax

		// Note: The correlation is for the movement of the reference over the target
		int nc_2 = nc / 2;
		int nr_2 = nr / 2;
		int ns_2 = ns / 2;

		// Compute the shift from the centre
		int dx = nc_2 - ix;
		int dy = nr_2 - iy;
		int dz = ns_2 - iz;

		// For the reference (moved -dx,-dy,-dz over the target)
		int rx = -dx;
		int ry = -dy;
		int rz = -dz;

		// For the target (moved dx,dy,dz over the reference)
		int tx = dx;
		int ty = dy;
		int tz = dz;

		// Precompute the x-1,x+w-1,y-1,y+h-1
		int nx = iw - ix;
		int[] rx_1 = new int[nx];
		int[] rx_w_1 = new int[nx];
		int[] tx_1 = new int[nx];
		int[] tx_w_1 = new int[nx];
		int[] w = new int[nx];
		for (int c = ix, i = 0; c < iw; c++, i++)
		{
			rx_1[i] = Math.max(-1, rx - 1);
			rx_w_1[i] = Math.min(nc, rx + nc) - 1;
			rx++;
			tx_1[i] = Math.max(-1, tx - 1);
			tx_w_1[i] = Math.min(nc, tx + nc) - 1;
			tx--;
			w[i] = rx_w_1[i] - rx_1[i];
		}
		int ny = ih - iy;
		int[] ry_1 = new int[ny];
		int[] ry_h_1 = new int[ny];
		int[] ty_1 = new int[ny];
		int[] ty_h_1 = new int[ny];
		int[] h = new int[ny];
		for (int r = iy, j = 0; r < ih; r++, j++)
		{
			ry_1[j] = Math.max(-1, ry - 1);
			ry_h_1[j] = Math.min(nr, ry + nr) - 1;
			ry++;
			ty_1[j] = Math.max(-1, ty - 1);
			ty_h_1[j] = Math.min(nr, ty + nr) - 1;
			ty--;
			h[j] = ry_h_1[j] - ry_1[j];
		}

		long[] rs_ = reference.s_;
		long[] rss = reference.ss;
		long[] ts_ = target.s_;
		long[] tss = target.ss;
		long[] rsum = new long[2];
		long[] tsum = new long[2];

		for (int s = iz; s < id; s++)
		{
			// Compute the z-1,z+d-1
			int rz_1 = Math.max(-1, rz - 1);
			int rz_d_1 = Math.min(ns, rz + ns) - 1;
			rz++;
			int tz_1 = Math.max(-1, tz - 1);
			int tz_d_1 = Math.min(ns, tz + ns) - 1;
			tz--;
			int d = rz_d_1 - rz_1;

			for (int r = iy, j = 0; r < ih; r++, j++)
			{
				int base = s * nr_by_nc + r * nc;
				int hd = h[j] * d;
				for (int c = ix, i = 0; c < iw; c++, i++)
				{
					double sumXY = data[base + c];

					compute(rx_1[i], ry_1[j], rz_1, rx_w_1[i], ry_h_1[j], rz_d_1, w[i], h[j], d, rs_, rss, rsum);
					compute(tx_1[i], ty_1[j], tz_1, tx_w_1[i], ty_h_1[j], tz_d_1, w[i], h[j], d, ts_, tss, tsum);

					// XXX debug this ...
					if (w[i] == nc && h[j] == nr && d == ns)
					{
						double rs = ref.computeSum(rx_1[i] + 1, ry_1[j] + 1, rz_1 + 1, w[i], h[j], d);
						double ts = tar.computeSum(tx_1[i] + 1, ty_1[j] + 1, tz_1 + 1, w[i], h[j], d);
						sumXY = 0;
						long sumXX = 0;
						long sumYY = 0;
						for (int k = ref.getDataLength(); k-- > 0;)
						{
							int a = (int) ref.getData()[k];
							int b = (int) tar.getData()[k];
							sumXY += a * b;
							sumXX += a * a;
							sumYY += b * b;
						}

						// The actual sumXY from the correlation is incorrect!
						// Is it because the data have been converted to int?

						System.out.printf("%g vs %g, %d vs %g, %d vs %g, %d vs %d, %d vs %d\n", data[base + c], sumXY,
								rsum[X], rs, tsum[Y], ts, sumXX, rsum[XX], sumYY, tsum[YY]);
					}

					// Compute the correlation
					double one_over_n = 1.0 / (w[i] * hd);
					// (sumXY - sumX*sumY/n) / sqrt( (sumXX - sumX^2 / n) * (sumYY - sumY^2 / n) )

					double pearsons1 = sumXY - (one_over_n * rsum[X] * tsum[Y]);
					double pearsons2 = rsum[XX] - (one_over_n * rsum[X] * rsum[X]);
					double pearsons3 = tsum[YY] - (one_over_n * tsum[Y] * tsum[Y]);

					double R;
					if (pearsons2 == 0 || pearsons3 == 0)
					{
						// If there is data and all the variances are the same then correlation is perfect
						if (rsum[XX] == tsum[YY] && rsum[XX] == sumXY && rsum[XX] > 0)
						{
							R = 1;
						}
						else
						{
							R = 0;
						}
					}
					else
					{
						R = pearsons1 / (Math.sqrt(pearsons2 * pearsons3));
						// Leave as raw for debugging
						//R = Maths.clip(-1, 1, R);
					}
					data[base + c] = (float) R;
				}
			}
		}

		//Utils.display("corr", correlation.getImageStack());
		int maxi = correlation.findMaxIndex(ix, iy, iz, iw - ix, ih - iy, id - iz);
		int[] xyz = correlation.getXYZ(maxi);

		// Q. What if the correlation surface is flat? This will pick the start and not
		// the centre of the plateau.

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
			int x = Maths.clip(0, correlation.getWidth() - 4, xyz[0] - 1);
			int y = Maths.clip(0, correlation.getHeight() - 4, xyz[1] - 1);
			int z = Maths.clip(0, correlation.getSize() - 4, xyz[2] - 1);
			Image3D crop = correlation.crop(x, y, z, 4, 4, 4, region);
			region = crop.getData();
			CustomTricubicFunction f = CustomTricubicFunction.create(calc.compute(region));

			// Find the maximum starting at the current origin
			int ox = xyz[0] - x;
			int oy = xyz[1] - y;
			int oz = xyz[2] - z;

			// Scale to the cubic spline dimensions of 0-1
			double[] origin = new double[] { ox / 3.0, oy / 3.0, oz / 3.0 };

			// Simple condensing search
			if (searchMode == SearchMode.BINARY)
			{
				// Can this use the current origin as a start point?
				// Currently we evaluate 8-cube vertices. A better search
				// would evaluate 27 points around the optimum, pick the best then condense 
				// the range.
				double[] optimum = f.search(true, refinements, relativeThreshold, -1);
				double value = optimum[3];
				if (value > result[3])
				{
					result[3] = value;
					// Convert the maximum back with scaling
					for (int i = 0; i < 3; i++)
						result[i] -= (optimum[i] - origin[i]) * 3.0;
					return result;
				}
			}
			else
			{
				// Gradient search
				try
				{
					final SplineFunction sf = new SplineFunction(f, origin);

				// @formatter:off				
				BFGSOptimizer optimiser = new BFGSOptimizer(
						// Use a simple check on the relative value change and 
						// set the number of refinements
						new SimpleValueChecker(relativeThreshold, -1, refinements));

				PointValuePair opt = optimiser.optimize(
						maxEvaluations,
						bounds, 
						gradientTolerance,
						stepLength,
						new InitialGuess(origin),
						// Scale the error for the position check
						new PositionChecker(-1, error / 3.0),
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

					// Check it is higher. Invert since we did a minimisation.
					double value = -opt.getValue();
					if (value > result[3])
					{
						result[3] = value;
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
		}

		return result;
	}

	/**
	 * Compute the sum from the rolling sum tables
	 *
	 * @param x_1
	 *            the x value -1
	 * @param y_1
	 *            the y value -1
	 * @param z_1
	 *            the z value -1
	 * @param x_w_1
	 *            the x value +w -1
	 * @param y_h_1
	 *            the y value +h -1
	 * @param z_d_1
	 *            the z value +d -1
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param d
	 *            the depth
	 * @param s_
	 *            the sum table
	 * @param ss
	 *            the sum-of-squares table
	 * @param sum
	 *            the sum (output = [sum, sum-of-squares])
	 */
	private void compute(int x_1, int y_1, int z_1, int x_w_1, int y_h_1, int z_d_1, int w, int h, int d, long[] s_,
			long[] ss, long[] sum)
	{
		// Compute sum from rolling sum using:
		// sum(x,y,z,w,h,d) = 
		// + s(x+w-1,y+h-1,z+d-1) 
		// - s(x-1,y+h-1,z+d-1)
		// - s(x+w-1,y-1,z+d-1)
		// + s(x-1,y-1,z+d-1)
		// /* Stack above must be subtracted so reverse sign*/
		// - s(x+w-1,y+h-1,z-1) 
		// + s(x-1,y+h-1,z-1)
		// + s(x+w-1,y-1,z-1)
		// - s(x-1,y-1,z-1)
		// Note: 
		// s(i,j,k) = 0 when either i,j,k < 0
		// i = imax when i>imax 
		// j = jmax when j>jmax 
		// k = kmax when k>kmax

		// This has been adapted from Image3D to compute the twos sums together

		int xw_yh_zd = reference.dht.getIndex(x_w_1, y_h_1, z_d_1);
		//int xw_yh_zd = z_d_1 * nr_by_nc + y_h_1 * nc + x_w_1;
		sum[0] = 0;
		sum[1] = 0;
		if (z_1 >= 0)
		{
			int xw_yh_z = xw_yh_zd - d * nr_by_nc;
			if (y_1 >= 0)
			{
				int h_ = h * nc;
				if (x_1 >= 0)
				{
					sum[0] = s_[xw_yh_zd - w - h_] - s_[xw_yh_z - w - h_] - s_[xw_yh_zd - w] + s_[xw_yh_z - w];
					sum[1] = ss[xw_yh_zd - w - h_] - ss[xw_yh_z - w - h_] - ss[xw_yh_zd - w] + ss[xw_yh_z - w];
				}
				sum[0] = sum[0] + s_[xw_yh_z - h_] - s_[xw_yh_zd - h_];
				sum[1] = sum[1] + ss[xw_yh_z - h_] - ss[xw_yh_zd - h_];
			}
			else if (x_1 >= 0)
			{
				sum[0] = s_[xw_yh_z - w] - s_[xw_yh_zd - w];
				sum[1] = ss[xw_yh_z - w] - ss[xw_yh_zd - w];
			}
			sum[0] = sum[0] + s_[xw_yh_zd] - s_[xw_yh_z];
			sum[1] = sum[1] + ss[xw_yh_zd] - ss[xw_yh_z];
		}
		else
		{
			if (y_1 >= 0)
			{
				int h_ = h * nc;
				if (x_1 >= 0)
				{
					sum[0] = s_[xw_yh_zd - w - h_] - s_[xw_yh_zd - w];
					sum[1] = ss[xw_yh_zd - w - h_] - ss[xw_yh_zd - w];
				}
				sum[0] -= s_[xw_yh_zd - h_];
				sum[1] -= ss[xw_yh_zd - h_];
			}
			else if (x_1 >= 0)
			{
				sum[0] = -s_[xw_yh_zd - w];
				sum[1] = -ss[xw_yh_zd - w];
			}
			sum[0] = sum[0] + s_[xw_yh_zd];
			sum[1] = sum[1] + ss[xw_yh_zd];
		}
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
			// BFGS algorithm minimimises so invert
			return -f.value(table);
		}

		void value(double[] point, double[] df_da)
		{
			initialise(point);
			f.gradient(table, df_da);
			// BFGS algorithm minimimises so invert
			for (int i = 0; i < 3; i++)
				df_da[i] = -df_da[i];
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
			copy.target = null;
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

	/**
	 * Gets the relative threshold for change in the correlation value for halting refinement. If this is negative it is
	 * disabled.
	 *
	 * @return the relative threshold
	 */
	public double getRelativeThreshold()
	{
		return relativeThreshold;
	}

	/**
	 * Sets the relative threshold for change in the correlation value for halting refinement. Set to negative to
	 * disable. Refinement will then only be halted by the number of refinement steps or the position error.
	 *
	 * @param relativeThreshold
	 *            the new relative threshold
	 */
	public void setRelativeThreshold(double relativeThreshold)
	{
		this.relativeThreshold = relativeThreshold;
	}
}