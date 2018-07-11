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
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.cspline.CubicSplineCalculator;
import ij.ImageStack;
import ij.process.ImageProcessor;

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
public class Image3DAligner implements Cloneable
{
	private static final int X = 0;
	private static final int XX = 1;
	private static final int Y = 0;
	private static final int YY = 1;

	private class DHTData
	{
		DoubleDHT3D dht;
		double[] input;
		double[] s_;
		double[] ss;
		// Original dimensions and 3D size
		int w, h, d, size;
		// Insert position
		int ix, iy, iz;

		DHTData(DoubleDHT3D dht, int w, int h, int d)
		{
			setDHT(dht, w, h, d);
		}

		void setDHT(DoubleDHT3D dht, int w, int h, int d)
		{
			this.dht = dht;
			s_ = resize(s_);
			ss = resize(ss);
			this.w = w;
			this.h = h;
			this.d = d;
			size = w * h * d;
			ix = getInsert(dht.nc, w);
			iy = getInsert(dht.nr, h);
			iz = getInsert(dht.ns, d);
			// Make storage of the original data optional. It is just used for
			// the spatial domain correlation check
			if (isCheckCorrelation())
			{
				input = resize(input);
			}
		}

		private double[] resize(double[] data)
		{
			return (data == null || data.length != dht.getDataLength()) ? new double[dht.getDataLength()] : data;
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
	private double relativeThreshold = 1e-6;
	private SearchMode searchMode = SearchMode.GRADIENT;
	private boolean checkCorrelation = true;
	private double minimumOverlap = 0.5;
	private double minimumDimensionOverlap = 0.75;
	private boolean fastMultiply = true;

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
	private double[] buffer, region;
	private double frequencyDomainCorrelationError;
	private int[] crop;

	// Allow cached window weights
	private double[] wx = null;
	private double[] wy = null;
	private double[] wz = null;

	private CubicSplineCalculator calc;

	/**
	 * Instantiates a new image aligner with a default edge window of 0.25
	 */
	public Image3DAligner()
	{
		this(0.25);
	}

	/**
	 * Instantiates a new image aligner.
	 *
	 * @param edgeWindow
	 *            the alpha value for the Tukey edge window
	 */
	public Image3DAligner(double edgeWindow)
	{
		setEdgeWindow(edgeWindow);
	}

	/**
	 * Sets the reference image and assumes the target image will be the same size.
	 * <p>
	 * The dimension are converted to the next power of 2 for speed. The combined size must fit within the maximum size
	 * of a single array.
	 *
	 * @param image
	 *            the image (destructively modified)
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined target dimensions is too large for an array
	 */
	public void setReference(ImageStack image)
	{
		setReference(image, image.getWidth(), image.getHeight(), image.getSize());
	}

	/**
	 * Sets the reference image and the size of the target image.
	 * <p>
	 * The dimension are converted to the next power of 2 for speed. The combined size must fit within the maximum size
	 * of a single array.
	 *
	 * @param image
	 *            the image (may be destructively modified)
	 * @param w
	 *            the width of the target image
	 * @param h
	 *            the height of the target image
	 * @param d
	 *            the depth of the target image
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined target dimensions is too large for an array
	 */
	public void setReference(ImageStack image, int w, int h, int d)
	{
		check3D(image);
		if (w < 2 || h < 2 || d < 2)
			throw new IllegalArgumentException("Require a 3D target image");
		nc = Maths.nextPow2(Math.max(w, image.getWidth()));
		nr = Maths.nextPow2(Math.max(h, image.getHeight()));
		ns = Maths.nextPow2(Math.max(d, image.getSize()));
		// Check the image will fit in an Image3D
		Image3D.checkSize(nc, nr, ns, true);
		nr_by_nc = nr * nc;
		// Window and pad the reference
		setReference(createDHT(image, reference));
	}

	private void setReference(DHTData dhtData)
	{
		reference = dhtData;
		if (fastMultiply)
			reference.dht.initialiseFastMultiply();
	}

	private void check3D(ImageStack image)
	{
		if (image.getWidth() < 2 || image.getHeight() < 2 || image.getSize() < 2)
			throw new IllegalArgumentException("Require a 3D image");
		// Check for data
		int size = image.getWidth() * image.getHeight();
		for (int s = 1; s <= image.getSize(); s++)
		{
			ImageProcessor ip = image.getProcessor(s);
			for (int i = 0; i < size; i++)
				if (ip.getf(i) != 0)
					return;
		}
		throw new IllegalArgumentException("No data in 3D image");
	}

	private DHTData createDHT(ImageStack image, DHTData dhtData)
	{
		if (image.getBitDepth() != 32)
			return createDHT(new FloatImage3D(image), dhtData);

		// Shift mean to 0 with optional window
		int w = image.getWidth(), h = image.getHeight(), d = image.getSize();
		double[] wx = createXWindow(w);
		double[] wy = createYWindow(h);
		double[] wz = createZWindow(d);

		// We need to compute the weighted centre
		double[] sum = new double[2];

		for (int z = 0; z < d; z++)
		{
			float[] pixels = (float[]) image.getPixels(1 + z);
			if (wz[z] == 0)
			{
				// Special case happens with Tukey window at the ends
			}
			else
			{
				calculateWeightedCentre(pixels, w, h, wx, wy, wz[z], sum);
			}
		}

		double shift = sum[0] / sum[1];

		for (int z = 0; z < d; z++)
		{
			float[] pixels = (float[]) image.getPixels(1 + z);
			if (wz[z] == 0)
			{
				// Special case happens with Tukey window at the ends
				Arrays.fill(pixels, 0);
			}
			else
			{
				applyWindow(pixels, w, h, wx, wy, wz[z], shift);
			}
		}

		DoubleDHT3D dht;

		// Pad into the desired data size.
		// We always do this so the data is reused
		int size = ns * nr * nc;
		double[] dest;
		if (dhtData == null || dhtData.dht.getDataLength() != size)
		{
			dest = new double[size];
		}
		else
		{
			// Re-use space
			dest = dhtData.dht.getData();
			Arrays.fill(dest, 0);
		}
		dht = new DoubleDHT3D(nc, nr, ns, dest, false);
		int ix = getInsert(nc, w);
		int iy = getInsert(nr, h);
		int iz = getInsert(ns, d);
		dht.insert(ix, iy, iz, image);

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

	private static void calculateWeightedCentre(float[] image, int maxx, int maxy, double[] wx, double[] wy, double wz,
			double[] sum)
	{
		calculateWeightedCentre(image, 0, maxx, maxy, wx, wy, wz, sum);
	}

	private static void calculateWeightedCentre(float[] image, int i, int maxx, int maxy, double[] wx, double[] wy,
			double wz, double[] sum)
	{
		for (int y = 0; y < maxy; y++)
		{
			double wyz = wy[y] * wz;
			for (int x = 0; x < maxx; x++, i++)
			{
				double w = wx[x] * wyz;
				sum[0] += image[i] * w;
				sum[1] += w;
			}
		}
	}

	private static void calculateWeightedCentre(Image3D image, int i, int maxx, int maxy, double[] wx, double[] wy,
			double wz, double[] sum)
	{
		for (int y = 0; y < maxy; y++)
		{
			double wyz = wy[y] * wz;
			for (int x = 0; x < maxx; x++, i++)
			{
				double w = wx[x] * wyz;
				sum[0] += image.get(i) * w;
				sum[1] += w;
			}
		}
	}

	private static void applyWindow(float[] image, int maxx, int maxy, double[] wx, double[] wy, double wz,
			double shift)
	{
		applyWindow(image, 0, maxx, maxy, wx, wy, wz, shift);
	}

	private static void applyWindow(float[] image, int i, int maxx, int maxy, double[] wx, double[] wy, double wz,
			double shift)
	{
		for (int y = 0; y < maxy; y++)
		{
			double wyz = wy[y] * wz;
			for (int x = 0; x < maxx; x++, i++)
			{
				image[i] = (float) ((image[i] - shift) * wx[x] * wyz);
			}
		}
	}

	private static void applyWindow(Image3D image, int i, int maxx, int maxy, double[] wx, double[] wy, double wz,
			double shift)
	{
		for (int y = 0; y < maxy; y++)
		{
			double wyz = wy[y] * wz;
			for (int x = 0; x < maxx; x++, i++)
			{
				image.set(i, (image.get(i) - shift) * wx[x] * wyz);
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
	 * Converts the data to quantised data. Any zero value (from padding
	 * or weighting) remains zero.
	 * <p>
	 * This may reduce the precision slightly but allows the computation of a rolling sum table have minimal errors. The
	 * rolling sum and sum-of-squares table is computed and the DHT is transformed to the frequency domain.
	 *
	 * @param dhtData
	 *            the dht data
	 * @return the DHT data
	 */
	private DHTData prepareDHT(DHTData dhtData)
	{
		DoubleDHT3D dht = dhtData.dht;
		double[] s_ = dhtData.s_;
		double[] ss = dhtData.ss;

		// Note previous versions converted to 10-bit integer data. However the 3D DHT creates very large
		// output values and errors occurred when computing the conjugate multiple in the frequency
		// domain verses the spatial domain. A check has been added to compute the spatial domain
		// correlation for the corresponding max correlation in the frequency domain. This allow
		// the code to report when the correlation value is incorrect.

		double[] data = dht.getData();
		double[] limits = Maths.limits(data);
		double min = limits[0];
		double max = limits[1];

		// Note: The image has been shifted to a mean of 0 so that zero padding
		// for frequency domain transform does not add any information.
		// We need to maintain the sign information and ensure that zero is still
		// zero.

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
			double sum_ = 0, sum2 = 0;
			int i = s * nr_by_nc;
			// Initialise first row sum
			// sum = rolling sum of (0 - colomn)
			for (int c = 0; c < nc; c++, i++)
			{
				double v = transform(data[i], scale);
				data[i] = v;
				sum_ += v;
				sum2 += v * v;
				s_[i] = sum_;
				ss[i] = sum2;
			}
			// Remaining rows
			// sum = rolling sum of (0 - colomn) + sum of same position above
			for (int r = 1, ii = i - nc; r < nr; r++)
			{
				sum_ = 0;
				sum2 = 0;
				for (int c = 0; c < nc; c++, i++, ii++)
				{
					double v = transform(data[i], scale);
					data[i] = v;
					sum_ += v;
					sum2 += v * v;
					// Add the sum from the previous row
					s_[i] = sum_ + s_[ii];
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

		// Store after numerical transform
		if (dhtData.input != null && dhtData.input.length == ss.length)
			System.arraycopy(dht.getData(), 0, dhtData.input, 0, ss.length);

		// Transform the data
		dht.transform();
		return dhtData;
	}

	/**
	 * The limit for the range of the data as an integer.
	 * <p>
	 * When this it too high the sumXY from the DHT conjugate multiplication
	 * does not match the sum from correlation in the spatial domain.
	 * <p>
	 * In theory the largest sumXY should be 2^bits * 2^bits * max integer (the size of the largest array).
	 * 10-bit integer: 2^10 * 2^10 * 2^31 = 2^51. This is smaller than the mantissa of a double (2^52)
	 * so should be represented correctly.
	 */
	private static double LIMIT = 1024;

	private static double transform(double f, double scale)
	{
		// Ensure zero is zero
		if (f == 0.0)
			return 0.0;

		// Maintain the sign information
		double value = f * scale;
		return Math.round(value); // / scale;
	}

	/**
	 * Sets the reference image and assumes the target image will be the same size.
	 * <p>
	 * The dimension are converted to the next power of 2 for speed. The combined size must fit within the maximum size
	 * of a single array.
	 *
	 * @param image
	 *            the image (destructively modified)
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined target dimensions is too large for an array
	 */
	public void setReference(Image3D image)
	{
		setReference(image, image.getWidth(), image.getHeight(), image.getSize());
	}

	/**
	 * Sets the reference image and the size of the target image.
	 * <p>
	 * The dimension are converted to the next power of 2 for speed. The combined size must fit within the maximum size
	 * of a single array.
	 *
	 * @param image
	 *            the image (may be destructively modified)
	 * @param w
	 *            the width of the target image
	 * @param h
	 *            the height of the target image
	 * @param d
	 *            the depth of the target image
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined target dimensions is too large for an array
	 */
	public void setReference(Image3D image, int w, int h, int d)
	{
		check3D(image);
		if (w < 2 || h < 2 || d < 2)
			throw new IllegalArgumentException("Require a 3D target image");
		nc = Maths.nextPow2(Math.max(w, image.getWidth()));
		nr = Maths.nextPow2(Math.max(h, image.getHeight()));
		ns = Maths.nextPow2(Math.max(d, image.getSize()));
		nr_by_nc = nr * nc;
		// Window and pad the reference
		setReference(createDHT(image, reference));
	}

	private void check3D(Image3D image)
	{
		if (image.getWidth() < 2 || image.getHeight() < 2 || image.getSize() < 2)
			throw new IllegalArgumentException("Require a 3D image");
		// Check for data
		for (int i = 0, size = image.getDataLength(); i < size; i++)
			if (image.get(i) != 0)
				return;
		throw new IllegalArgumentException("No data in 3D image");
	}

	private DHTData createDHT(Image3D image, DHTData dhtData)
	{
		// Shift mean to 0 with optional window
		int w = image.getWidth(), h = image.getHeight(), d = image.getSize();
		double[] wx = createXWindow(w);
		double[] wy = createYWindow(h);
		double[] wz = createZWindow(d);
		int inc = image.nr_by_nc;

		// We need to compute the weighted centre
		double[] sum = new double[2];

		for (int z = 0, i = 0; z < d; z++)
		{
			if (wz[z] == 0)
			{
				// Special case happens with Tukey window at the ends
			}
			else
			{
				calculateWeightedCentre(image, i, w, h, wx, wy, wz[z], sum);
			}
			i += inc;
		}

		double shift = sum[0] / sum[1];

		for (int z = 0, i = 0; z < d; z++)
		{
			if (wz[z] == 0)
			{
				// Special case happens with Tukey window at the ends
				for (int j = 0; j < inc; j++)
					image.set(i++, 0);
			}
			else
			{
				applyWindow(image, i, w, h, wx, wy, wz[z], shift);
				i += inc;
			}
		}

		//System.out.printf("Sum = %g => %g\n", sum[0], Maths.sum(pixels));

		DoubleDHT3D dht;

		// Pad into the desired data size.
		// We always do this to handle input of float/double Image3D data.
		int size = ns * nr * nc;
		double[] dest;
		if (dhtData == null || dhtData.dht.getDataLength() != size)
		{
			dest = new double[size];
		}
		else
		{
			// Re-use space
			dest = dhtData.dht.getData();
			Arrays.fill(dest, 0);
		}
		dht = new DoubleDHT3D(nc, nr, ns, dest, false);
		int ix = getInsert(nc, w);
		int iy = getInsert(nr, h);
		int iz = getInsert(ns, d);
		dht.insert(ix, iy, iz, image);

		if (dhtData == null)
			dhtData = new DHTData(dht, w, h, d);
		else
			dhtData.setDHT(dht, w, h, d);

		return prepareDHT(dhtData);
	}

	/**
	 * Align the image with the reference. Compute the translation required to move the target image onto the reference
	 * image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(ImageStack image)
	{
		return align(image, 0, 0);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
	 * <p>
	 * Refinement uses a default sub-pixel accuracy of 1e-2;
	 *
	 * @param image
	 *            the image
	 * @param refinements
	 *            the refinements for sub-pixel accuracy
	 * @param error
	 *            the error for sub-pixel accuracy (i.e. stop when improvements are less than this error)
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(ImageStack image, int refinements)
	{
		check3D(image);
		int w = image.getWidth(), h = image.getHeight(), d = image.getSize();
		if (w > nc || h > nr || d > ns)
			throw new IllegalArgumentException("Image is larger than the initialised reference");

		target = createDHT(image, target);
		return align(target, refinements, 1e-2);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @param refinements
	 *            the refinements for sub-pixel accuracy
	 * @param error
	 *            the error for sub-pixel accuracy (i.e. stop when improvements are less than this error)
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(ImageStack image, int refinements, double error)
	{
		check3D(image);
		int w = image.getWidth(), h = image.getHeight(), d = image.getSize();
		if (w > nc || h > nr || d > ns)
			throw new IllegalArgumentException("Image is larger than the initialised reference");

		target = createDHT(image, target);
		return align(target, refinements, error);
	}

	/**
	 * Align the image with the reference. Compute the translation required to move the target image onto the reference
	 * image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(Image3D image)
	{
		return align(image, 0, 0);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
	 * <p>
	 * Refinement uses a default sub-pixel accuracy of 1e-2;
	 *
	 * @param image
	 *            the image
	 * @param refinements
	 *            the refinements for sub-pixel accuracy
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(Image3D image, int refinements)
	{
		check3D(image);
		int w = image.getWidth(), h = image.getHeight(), d = image.getSize();
		if (w > nc || h > nr || d > ns)
			throw new IllegalArgumentException("Image is larger than the initialised reference");

		target = createDHT(image, target);
		return align(target, refinements, 1e-2);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @param refinements
	 *            the maximum number of refinements for sub-pixel accuracy
	 * @param error
	 *            the error for sub-pixel accuracy (i.e. stop when improvements are less than this error)
	 * @return [x,y,z,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(Image3D image, int refinements, double error)
	{
		check3D(image);
		int w = image.getWidth(), h = image.getHeight(), d = image.getSize();
		if (w > nc || h > nr || d > ns)
			throw new IllegalArgumentException("Image is larger than the initialised reference");

		target = createDHT(image, target);
		return align(target, refinements, error);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
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
		DoubleDHT3D correlation = target.dht.conjugateMultiply(reference.dht, buffer);
		buffer = correlation.getData(); // Store for reuse
		correlation.inverseTransform();
		correlation.swapOctants();

		// Normalise:
		//  ( Σ xiyi - nx̄ӯ ) / ( (Σ xi^2 - nx̄^2) (Σ yi^2 - nӯ^2) )^0.5
		//
		// (sumXY - sumX*sumY/n) / sqrt( (sumXX - sumX^2 / n) * (sumYY - sumY^2 / n) )

		// Only do this over the range where at least half the original images overlap,
		// i.e. the insert point of one will be the middle of the other when shifted.
		int ix = Math.min(reference.ix, target.ix);
		int iy = Math.min(reference.iy, target.iy);
		int iz = Math.min(reference.iz, target.iz);
		int ixw = Math.max(reference.ix + reference.w, target.ix + target.w);
		int iyh = Math.max(reference.iy + reference.h, target.iy + target.h);
		int izd = Math.max(reference.iz + reference.d, target.iz + target.d);

		if (minimumDimensionOverlap > 0)
		{
			double f = (1 - minimumDimensionOverlap) / 2;
			int ux = (int) (Math.round(Math.min(reference.w, target.w) * f));
			int uy = (int) (Math.round(Math.min(reference.h, target.h) * f));
			int uz = (int) (Math.round(Math.min(reference.d, target.d) * f));
			ix += ux;
			ixw -= ux;
			iy += uy;
			iyh -= uy;
			iz += uz;
			izd -= uz;
		}

		crop = new int[] { ix, iy, iz, ixw - ix, iyh - iy, izd - iz };

		// The maximum correlation unnormalised. Since this is unnormalised
		// it will be biased towards the centre of the image. This is used
		// to restrict the bounds for finding the maximum of the normalised correlation
		// which should be close to this.
		int maxi = correlation.findMaxIndex(ix, iy, iz, crop[3], crop[4], crop[5]);
		int[] xyz = correlation.getXYZ(maxi);

		// Check in the spatial domain
		checkCorrelation(target, correlation, maxi);

		// Compute sum from rolling sum using:
		// sum(x,y,z,w,h,d) =
		// + s(x+w-1,y+h-1,z+d-1)
		// - s(x-1,y+h-1,z+d-1)
		// - s(x+w-1,y-1,z+d-1)
		// + s(x-1,y-1,z+d-1)
		// /* Image above must be subtracted so reverse sign*/
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
		int[] centre = new int[] { nc_2, nr_2, ns_2 };

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
		int nx = crop[3];
		int[] rx_1 = new int[nx];
		int[] rx_w_1 = new int[nx];
		int[] tx_1 = new int[nx];
		int[] tx_w_1 = new int[nx];
		int[] w = new int[nx];
		for (int c = ix, i = 0; c < ixw; c++, i++)
		{
			rx_1[i] = Math.max(-1, rx - 1);
			rx_w_1[i] = Math.min(nc, rx + nc) - 1;
			rx++;
			tx_1[i] = Math.max(-1, tx - 1);
			tx_w_1[i] = Math.min(nc, tx + nc) - 1;
			tx--;
			w[i] = rx_w_1[i] - rx_1[i];
		}
		int ny = crop[4];
		int[] ry_1 = new int[ny];
		int[] ry_h_1 = new int[ny];
		int[] ty_1 = new int[ny];
		int[] ty_h_1 = new int[ny];
		int[] h = new int[ny];
		for (int r = iy, j = 0; r < iyh; r++, j++)
		{
			ry_1[j] = Math.max(-1, ry - 1);
			ry_h_1[j] = Math.min(nr, ry + nr) - 1;
			ry++;
			ty_1[j] = Math.max(-1, ty - 1);
			ty_h_1[j] = Math.min(nr, ty + nr) - 1;
			ty--;
			h[j] = ry_h_1[j] - ry_1[j];
		}

		double[] rs_ = reference.s_;
		double[] rss = reference.ss;
		double[] ts_ = target.s_;
		double[] tss = target.ss;
		double[] rsum = new double[2];
		double[] tsum = new double[2];

		int size = Math.min(reference.size, target.size);
		int minimumN = (int) (Math.round(size * minimumOverlap));
		int maxj = -1;
		double max = 0;

		for (int s = iz; s < izd; s++)
		{
			// Compute the z-1,z+d-1
			int rz_1 = Math.max(-1, rz - 1);
			int rz_d_1 = Math.min(ns, rz + ns) - 1;
			rz++;
			int tz_1 = Math.max(-1, tz - 1);
			int tz_d_1 = Math.min(ns, tz + ns) - 1;
			tz--;
			int d = rz_d_1 - rz_1;

			for (int r = iy, j = 0; r < iyh; r++, j++)
			{
				int base = s * nr_by_nc + r * nc;
				int hd = h[j] * d;
				for (int c = ix, i = 0; c < ixw; c++, i++)
				{
					double sumXY = buffer[base + c];

					compute(rx_1[i], ry_1[j], rz_1, rx_w_1[i], ry_h_1[j], rz_d_1, w[i], h[j], d, rs_, rss, rsum);
					compute(tx_1[i], ty_1[j], tz_1, tx_w_1[i], ty_h_1[j], tz_d_1, w[i], h[j], d, ts_, tss, tsum);

					// Compute the correlation
					// (sumXY - sumX*sumY/n) / sqrt( (sumXX - sumX^2 / n) * (sumYY - sumY^2 / n) )

					int n = w[i] * hd;
					double numerator = sumXY - (rsum[X] * tsum[Y] / n);
					double denominator1 = rsum[XX] - (rsum[X] * rsum[X] / n);
					double denominator2 = tsum[YY] - (tsum[Y] * tsum[Y] / n);

					double R;
					if (denominator1 == 0 || denominator2 == 0)
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
						R = numerator / Math.sqrt(denominator1 * denominator2);
						// Leave as raw for debugging
						//R = Maths.clip(-1, 1, R);
					}

					buffer[base + c] = R;

					if (n < minimumN)
						continue;

					if (R > 1.0001) // some margin for error
					{
						// Normalisation has failed.
						// This occurs when the correlation sum XY is incorrect.
						// The other terms are exact due to the quantisation to integer data.
						// It is likely to occur at the bounds.

						System.out.printf("Bad normalisation [%d,%d,%d] = %g  (overlap=%g)\n", c, r, s, R,
								(double) n / size);
						continue;
					}

					if (R > max)
					{
						max = R;
						maxj = base + c;
					}
					else if (R == max)
					{
						// Get shift from centre
						int[] xyz1 = correlation.getXYZ(maxj);
						int[] xyz2 = correlation.getXYZ(base + c);
						int d1 = 0, d2 = 0;
						for (int k = 0; k < 3; k++)
						{
							d1 += Maths.pow2(xyz1[k] - centre[k]);
							d2 += Maths.pow2(xyz2[k] - centre[k]);
						}
						if (d2 < d1)
						{
							max = R;
							maxj = base + c;
						}
					}
				}
			}
		}

		// The maximum correlation with normalisation
		maxi = maxj; //correlation.findMaxIndex(ix, iy, iz, iw - ix, ih - iy, id - iz);
		xyz = correlation.getXYZ(maxi);

		// Report the shift required to move from the centre of the target image to the reference
		// @formatter:off
		double[] result = new double[] {
			nc_2 - xyz[0],
			nr_2 - xyz[1],
			ns_2 - xyz[2],
			buffer[maxi]
		};
		// @formatter:on

		if (refinements > 0)
		{
			// Perform sub-pixel alignment
			// Create a cubic spline using a small region of pixels around the maximum
			if (calc == null)
				calc = new CubicSplineCalculator();
			// Avoid out-of-bounds errors. Only use the range that was normalised
			int x = Maths.clip(ix, ixw - 4, xyz[0] - 1);
			int y = Maths.clip(iy, iyh - 4, xyz[1] - 1);
			int z = Maths.clip(iz, izd - 4, xyz[2] - 1);
			DoubleImage3D crop = correlation.crop(x, y, z, 4, 4, 4, region);
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
							@Override
							public double value(double[] point)
							{
								return sf.value(point);
							}}),
						new ObjectiveFunctionGradient(new MultivariateVectorFunction(){
							@Override
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
	 * Check the correlation in the spatial domain verses the maximum correlation in the frequency domain.
	 *
	 * @param target
	 *            the target
	 * @param correlation
	 *            the correlation
	 * @param maxi
	 *            the index of the maximum correlation
	 */
	private void checkCorrelation(DHTData target, DoubleDHT3D correlation, int maxi)
	{
		if (target.input == null || reference.input == null)
			// No check possible
			return;

		// The maximum correlation without normalisation
		int[] xyz = correlation.getXYZ(maxi);

		// Find the range for the target and reference
		int nc_2 = nc / 2;
		int nr_2 = nr / 2;
		int ns_2 = ns / 2;
		int tx = Math.max(0, xyz[0] - nc_2);
		int ty = Math.max(0, xyz[1] - nr_2);
		int tz = Math.max(0, xyz[2] - ns_2);
		int w = Math.min(nc, xyz[0] + nc_2) - tx;
		int h = Math.min(nr, xyz[1] + nr_2) - ty;
		int d = Math.min(ns, xyz[2] + ns_2) - tz;

		// For the reference we express as a shift relative to the centre
		// and subtract the half-width.
		// Formally: (nc_2 - xyz[0]) // shift
		//           + nc_2          // centre
		//           - nc_2          // Half width
		int rx = Math.max(0, -xyz[0] + nc_2);
		int ry = Math.max(0, -xyz[1] + nr_2);
		int rz = Math.max(0, -xyz[2] + ns_2);

		double[] tar = target.input;
		double[] ref = reference.input;
		double o = correlation.get(maxi);
		double e = 0;
		for (int z = 0; z < d; z++)
		{
			for (int y = 0; y < h; y++)
			{
				int i = (tz + z) * nr_by_nc + (ty + y) * nc + tx;
				int j = (rz + z) * nr_by_nc + (ry + y) * nc + rx;
				for (int x = 0; x < w; x++)
				{
					e += tar[i++] * ref[j++];
				}
			}
		}

		//System.out.printf("Raw %d,%d,%d = %g\n", xyz[0], xyz[1], xyz[1], o);

		frequencyDomainCorrelationError = DoubleEquality.relativeError(o, e);
		if (frequencyDomainCorrelationError > 0.05)
		{
			System.err.printf("3D Correlation Error = %s : Spatial = %s, Freq = %s\n",
					Utils.rounded(frequencyDomainCorrelationError), Double.toString(e), Double.toString(o));
		}
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
	private void compute(int x_1, int y_1, int z_1, int x_w_1, int y_h_1, int z_d_1, int w, int h, int d, double[] s_,
			double[] ss, double[] sum)
	{
		// Compute sum from rolling sum using:
		// sum(x,y,z,w,h,d) =
		// + s(x+w-1,y+h-1,z+d-1)
		// - s(x-1,y+h-1,z+d-1)
		// - s(x+w-1,y-1,z+d-1)
		// + s(x-1,y-1,z+d-1)
		// /* Image above must be subtracted so reverse sign*/
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

		//int xw_yh_zd = reference.dht.getIndex(x_w_1, y_h_1, z_d_1);
		int xw_yh_zd = z_d_1 * nr_by_nc + y_h_1 * nc + x_w_1;
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
	 * @return the image aligner
	 */
	public Image3DAligner copy()
	{
		Image3DAligner copy;
		try
		{
			copy = (Image3DAligner) clone();
			// Reset objects that are not thread safe
			copy.calc = null;
			copy.buffer = null;
			copy.region = null;
			copy.target = null;
			copy.crop = null;
			copy.frequencyDomainCorrelationError = 0;
			return copy;
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}

	/**
	 * Gets the correlation image from the last alignment.
	 *
	 * @return the correlation (or null)
	 */
	public Image3D getCorrelation()
	{
		try
		{
			DoubleImage3D image = new DoubleImage3D(nc, nr, ns, buffer);
			image.fillOutside(crop[0], crop[1], crop[2], crop[3], crop[5], crop[5], 0);
			return image;
		}
		catch (IllegalArgumentException e)
		{
			// Thrown when buffer is null or does not match the dimensions.
			return null;
		}
	}

	/**
	 * Gets the frequency domain correlation error from the last correlation.
	 *
	 * @return the frequency domain correlation error
	 */
	public double getFrequencyDomainCorrelationError()
	{
		return frequencyDomainCorrelationError;
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
		this.edgeWindow = Maths.clip(0, 1, edgeWindow);
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

	/**
	 * Gets the search mode.
	 *
	 * @return the search mode
	 */
	public SearchMode getSearchMode()
	{
		return searchMode;
	}

	/**
	 * Sets the search mode.
	 *
	 * @param searchMode
	 *            the new search mode
	 */
	public void setSearchMode(SearchMode searchMode)
	{
		this.searchMode = searchMode;
	}

	/**
	 * Checks if the spatial domain correlation check is enabled.
	 *
	 * @return true, if the spatial domain correlation check is enabled
	 */
	public boolean isCheckCorrelation()
	{
		return checkCorrelation;
	}

	/**
	 * Sets the spatial domain correlation check flag. If true then the original untransformed data will be stored in
	 * memory. The point of the highest correlation in the frequency domain will be recomputed in the spatial domain.
	 * The error between the two can be returned using {@link #getFrequencyDomainCorrelationError()}.
	 *
	 * @param checkCorrelation
	 *            the new check correlation flag
	 */
	public void setCheckCorrelation(boolean checkCorrelation)
	{
		this.checkCorrelation = checkCorrelation;
	}

	/**
	 * Gets the minimum overlap between the smaller image and the other image.
	 *
	 * @return the minimum overlap
	 */
	public double getMinimumOverlap()
	{
		return minimumOverlap;
	}

	/**
	 * Sets the minimum overlap between the smaller image and the other image.
	 *
	 * @param minimumOverlap
	 *            the new minimum overlap
	 */
	public void setMinimumOverlap(double minimumOverlap)
	{
		this.minimumOverlap = Maths.clip(0, 1, minimumOverlap);
	}

	/**
	 * Gets the minimum overlap between the smaller image and the other image in each dimension.
	 *
	 * @return the minimum dimension overlap
	 */
	public double getMinimumDimensionOverlap()
	{
		return minimumDimensionOverlap;
	}

	/**
	 * Sets the minimum overlap between the smaller image and the other image in each dimension.
	 *
	 * @param minimumDimensionOverlap
	 *            the new minimum dimension overlap
	 */
	public void setMinimumDimensionOverlap(double minimumDimensionOverlap)
	{
		this.minimumDimensionOverlap = Maths.clip(0, 1, minimumDimensionOverlap);
	}

	/**
	 * Checks if is fast multiply.
	 *
	 * @return true, if is fast multiply
	 */
	public boolean isFastMultiply()
	{
		return fastMultiply;
	}

	/**
	 * Sets the fast multiply flag. This initialises the DHT for multiplication at the cost of extra memory storage. The
	 * storage requirements are 2 double arrays and 1 integer array of the same length at the FHT data.
	 *
	 * @param fastMultiply
	 *            the new fast multiply flag
	 */
	public void setFastMultiply(boolean fastMultiply)
	{
		this.fastMultiply = fastMultiply;
	}
}
