package gdsc.smlm.ij.utils;

import java.util.Arrays;

import gdsc.core.ij.Utils;
import gdsc.core.math.interpolation.CachedBicubicInterpolator;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

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
 * Perform 2D image alignment using normalised cross-correlation.
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
public class Image2DAligner implements Cloneable
{
	private static final int X = 0;
	private static final int XX = 1;
	private static final int Y = 0;
	private static final int YY = 1;

	private class DHTData
	{
		DoubleDHT2D dht;
		double[] input;
		double[] s_;
		double[] ss;
		// Original dimensions and 2D size
		int w, h, size;
		// Insert position
		int ix, iy;

		DHTData(DoubleDHT2D dht, int w, int h)
		{
			setDHT(dht, w, h);
		}

		void setDHT(DoubleDHT2D dht, int w, int h)
		{
			this.dht = dht;
			s_ = resize(s_);
			ss = resize(ss);
			this.w = w;
			this.h = h;
			size = w * h;
			ix = getInsert(dht.nc, w);
			iy = getInsert(dht.nr, h);
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

	private double edgeWindow;
	private double relativeThreshold = 1e-6;
	private boolean checkCorrelation = true;
	private double minimumOverlap = 0.5;
	private double minimumDimensionOverlap = 0.75;

	/** The number of rows (max y) of the discrete Hartley transform. */
	private int nr;
	/** The number of columns (max x) of the discrete Hartley transform. */
	private int nc;
	/** The number of rows by columns of the discrete Hartley transform. */
	private int nr_by_nc;

	/** The reference. */
	private DHTData reference;

	// Not thread safe as they are used for the target image
	private DHTData target;
	private double[] buffer, region;
	private double frequencyDomainCorrelationError;
	private int[] crop;

	// Allow cached window weights
	private double[] wx = null;
	private double[] wy = null;

	/**
	 * Instantiates a new image aligner with a default edge window of 0.25
	 */
	public Image2DAligner()
	{
		this(0.25);
	}

	/**
	 * Instantiates a new image aligner.
	 *
	 * @param edgeWindow
	 *            the alpha value for the Tukey edge window
	 */
	public Image2DAligner(double edgeWindow)
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
	 *             If any dimension is less than 2
	 */
	public void setReference(ImageProcessor image)
	{
		setReference(image, image.getWidth(), image.getHeight());
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
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2
	 */
	public void setReference(ImageProcessor image, int w, int h)
	{
		check2D(image);
		if (w < 2 || h < 2)
			throw new IllegalArgumentException("Require a 2D target image");
		nc = Maths.nextPow2(Math.max(w, image.getWidth()));
		nr = Maths.nextPow2(Math.max(h, image.getHeight()));
		nr_by_nc = nr * nc;
		// Window and pad the reference
		setReference(createDHT(image, reference));
	}

	/**
	 * Sets the reference.
	 *
	 * @param dhtData
	 *            the new reference
	 */
	private void setReference(DHTData dhtData)
	{
		reference = dhtData;
		dhtData.dht.initialiseFastMultiply();
	}

	/**
	 * Check the image is 2D and has data.
	 *
	 * @param image
	 *            the image
	 */
	private void check2D(ImageProcessor image)
	{
		if (image.getWidth() < 2 || image.getHeight() < 2)
			throw new IllegalArgumentException("Require a 2D image");
		// Check for data
		int size = image.getWidth() * image.getHeight();
		for (int i = 0; i < size; i++)
			if (image.getf(i) != 0)
				return;
		throw new IllegalArgumentException("No data in 2D image");
	}

	/**
	 * Creates the DHT.
	 *
	 * @param image
	 *            the image
	 * @param dhtData
	 *            the dht data
	 * @return the DHT data
	 */
	private DHTData createDHT(ImageProcessor image, DHTData dhtData)
	{
		if (image.getBitDepth() != 32)
			return createDHT(new FloatImage2D(image), dhtData);

		// Shift mean to 0 with optional window		
		int w = image.getWidth(), h = image.getHeight();
		double[] wx = createXWindow(w);
		double[] wy = createYWindow(h);

		// We need to compute the weighted centre
		double[] sum = new double[2];

		float[] pixels = (float[]) image.getPixels();
		calculateWeightedCentre(pixels, w, h, wx, wy, sum);

		double shift = sum[0] / sum[1];

		applyWindow(pixels, w, h, wx, wy, shift);

		DoubleDHT2D dht;

		// Pad into the desired data size.
		// We always do this so the data is reused
		double[] dest;
		if (dhtData == null || dhtData.dht.getDataLength() != nr_by_nc)
		{
			dest = new double[nr_by_nc];
		}
		else
		{
			// Re-use space
			dest = dhtData.dht.getData();
			Arrays.fill(dest, 0);
		}
		dht = new DoubleDHT2D(nc, nr, dest, false);
		int ix = getInsert(nc, w);
		int iy = getInsert(nr, h);
		dht.insert(ix, iy, image);

		if (dhtData == null)
			dhtData = new DHTData(dht, w, h);
		else
			dhtData.setDHT(dht, w, h);

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

	private double[] createWindow(double[] w, int n)
	{
		if (w == null || w.length != n)
			w = ImageWindow.tukey(n, edgeWindow);
		return w;
	}

	private static void calculateWeightedCentre(float[] image, int maxx, int maxy, double[] wx, double[] wy,
			double[] sum)
	{
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double w = wx[x] * wy[y];
				sum[0] += image[i] * w;
				sum[1] += w;
			}
		}
	}

	private static void calculateWeightedCentre(Image2D image, int maxx, int maxy, double[] wx, double[] wy,
			double[] sum)
	{
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double w = wx[x] * wy[y];
				sum[0] += image.get(i) * w;
				sum[1] += w;
			}
		}
	}

	private static void applyWindow(float[] image, int maxx, int maxy, double[] wx, double[] wy, double shift)
	{
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				image[i] = (float) ((image[i] - shift) * wx[x] * wy[y]);
			}
		}
	}

	private static void applyWindow(Image2D image, int maxx, int maxy, double[] wx, double[] wy, double shift)
	{
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				image.set(i, (image.get(i) - shift) * wx[x] * wy[y]);
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
		DoubleDHT2D dht = dhtData.dht;
		double[] s_ = dhtData.s_;
		double[] ss = dhtData.ss;

		// Note previous versions converted to 10-bit integer data. However the 2D DHT creates very large
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
		int nc = dht.nc;
		int nr = dht.nr;

		// This has been adapted from Image2D to compute two rolling sum table at once

		double sum_ = 0, sum2 = 0;
		int i = 0;
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
		for (int r = 1, ii = 0; r < nr; r++)
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
	 *             If any dimension is less than 2
	 */
	public void setReference(Image2D image)
	{
		setReference(image, image.getWidth(), image.getHeight());
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
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2
	 */
	public void setReference(Image2D image, int w, int h)
	{
		check2D(image);
		if (w < 2 || h < 2)
			throw new IllegalArgumentException("Require a 2D target image");
		nc = Maths.nextPow2(Math.max(w, image.getWidth()));
		nr = Maths.nextPow2(Math.max(h, image.getHeight()));
		nr_by_nc = nr * nc;
		// Window and pad the reference
		setReference(createDHT(image, reference));
	}

	/**
	 * Check the image is 2D and has data.
	 *
	 * @param image
	 *            the image
	 */
	private void check2D(Image2D image)
	{
		if (image.getWidth() < 2 || image.getHeight() < 2)
			throw new IllegalArgumentException("Require a 2D image");
		// Check for data
		for (int i = 0, size = image.getDataLength(); i < size; i++)
			if (image.get(i) != 0)
				return;
		throw new IllegalArgumentException("No data in 2D image");
	}

	private DHTData createDHT(Image2D image, DHTData dhtData)
	{
		// Shift mean to 0 with optional window		
		int w = image.getWidth(), h = image.getHeight();
		double[] wx = createXWindow(w);
		double[] wy = createYWindow(h);

		// We need to compute the weighted centre
		double[] sum = new double[2];

		calculateWeightedCentre(image, w, h, wx, wy, sum);

		double shift = sum[0] / sum[1];

		applyWindow(image, w, h, wx, wy, shift);

		//System.out.printf("Sum = %g => %g\n", sum[0], Maths.sum(pixels));

		DoubleDHT2D dht;

		// Pad into the desired data size.
		// We always do this to handle input of float/double Image2D data.
		double[] dest;
		if (dhtData == null || dhtData.dht.getDataLength() != nr_by_nc)
		{
			dest = new double[nr_by_nc];
		}
		else
		{
			// Re-use space
			dest = dhtData.dht.getData();
			Arrays.fill(dest, 0);
		}
		dht = new DoubleDHT2D(nc, nr, dest, false);
		int ix = getInsert(nc, w);
		int iy = getInsert(nr, h);
		dht.insert(ix, iy, image);

		if (dhtData == null)
			dhtData = new DHTData(dht, w, h);
		else
			dhtData.setDHT(dht, w, h);

		return prepareDHT(dhtData);
	}

	/**
	 * Align the image with the reference. Compute the translation required to move the target image onto the reference
	 * image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @return [x,y,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(ImageProcessor image)
	{
		return align(image, 0);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @param refinements
	 *            the refinements for sub-pixel accuracy
	 * @return [x,y,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(ImageProcessor image, int refinements)
	{
		check2D(image);
		int w = image.getWidth(), h = image.getHeight();
		if (w > nc || h > nr)
			throw new IllegalArgumentException("Image is larger than the initialised reference");

		target = createDHT(image, target);
		return align(target, refinements);
	}

	/**
	 * Align the image with the reference. Compute the translation required to move the target image onto the reference
	 * image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @return [x,y,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(Image2D image)
	{
		return align(image, 0);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
	 *
	 * @param image
	 *            the image
	 * @param refinements
	 *            the refinements for sub-pixel accuracy
	 * @return [x,y,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	public double[] align(Image2D image, int refinements)
	{
		check2D(image);
		int w = image.getWidth(), h = image.getHeight();
		if (w > nc || h > nr)
			throw new IllegalArgumentException("Image is larger than the initialised reference");

		target = createDHT(image, target);
		return align(target, refinements);
	}

	/**
	 * Align the image with the reference with sub-pixel accuracy. Compute the translation required to move the target
	 * image onto the reference image for maximum correlation.
	 *
	 * @param target
	 *            the target
	 * @param refinements
	 *            the maximum number of refinements for sub-pixel accuracy
	 * @return [x,y,value]
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if larger than the initialised reference
	 */
	private double[] align(DHTData target, int refinements)
	{
		// Multiply by the reference. This allows the reference to be shared across threads.
		DoubleDHT2D correlation = target.dht.conjugateMultiply(reference.dht, buffer);
		buffer = correlation.getData(); // Store for reuse
		correlation.inverseTransform();
		correlation.swapQuadrants();

		// Normalise:
		//  ( Σ xiyi - nx̄ӯ ) / ( (Σ xi^2 - nx̄^2) (Σ yi^2 - nӯ^2) )^0.5
		// 
		// (sumXY - sumX*sumY/n) / sqrt( (sumXX - sumX^2 / n) * (sumYY - sumY^2 / n) )

		// Only do this over the range where at least half the original images overlap,
		// i.e. the insert point of one will be the middle of the other when shifted.
		int ix = Math.min(reference.ix, target.ix);
		int iy = Math.min(reference.iy, target.iy);
		int ixw = Math.max(reference.ix + reference.w, target.ix + target.w);
		int iyh = Math.max(reference.iy + reference.h, target.iy + target.h);

		if (minimumDimensionOverlap > 0)
		{
			double f = (1 - minimumDimensionOverlap) / 2;
			int ux = (int) (Math.round(Math.min(reference.w, target.w) * f));
			int uy = (int) (Math.round(Math.min(reference.h, target.h) * f));
			ix += ux;
			ixw -= ux;
			iy += uy;
			iyh -= uy;
		}

		crop = new int[] { ix, iy, ixw - ix, iyh - iy };

		// The maximum correlation unnormalised. Since this is unnormalised
		// it will be biased towards the centre of the image. This is used
		// to restrict the bounds for finding the maximum of the normalised correlation
		// which should be close to this.
		int maxi = correlation.findMaxIndex(ix, iy, crop[2], crop[3]);
		int[] xy = correlation.getXY(maxi);

		// Check in the spatial domain
		checkCorrelation(target, correlation, maxi);

		// Compute sum from rolling sum using:
		// sum(x,y,w,h) = 
		// + s(x+w-1,y+h-1) 
		// - s(x-1,y+h-1)
		// - s(x+w-1,y-1)
		// + s(x-1,y-1)
		// Note: 
		// s(i,j) = 0 when either i,j < 0
		// i = imax when i>imax 
		// j = jmax when j>jmax 

		// Note: The correlation is for the movement of the reference over the target
		int nc_2 = nc / 2;
		int nr_2 = nr / 2;
		int[] centre = new int[] { nc_2, nr_2 };

		// Compute the shift from the centre
		int dx = nc_2 - ix;
		int dy = nr_2 - iy;

		// For the reference (moved -dx,-dy over the target)
		int rx = -dx;
		int ry = -dy;

		// For the target (moved dx,dy over the reference)
		int tx = dx;
		int ty = dy;

		// Precompute the x-1,x+w-1
		int nx = crop[2];
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

		for (int r = iy; r < iyh; r++)
		{
			// Compute the y-1,y+h-1
			int ry_1 = Math.max(-1, ry - 1);
			int ry_h_1 = Math.min(nr, ry + nr) - 1;
			ry++;
			int ty_1 = Math.max(-1, ty - 1);
			int ty_h_1 = Math.min(nr, ty + nr) - 1;
			ty--;
			int h = ry_h_1 - ry_1;

			int base = r * nc;
			for (int c = ix, i = 0; c < ixw; c++, i++)
			{
				double sumXY = buffer[base + c];

				compute(rx_1[i], ry_1, rx_w_1[i], ry_h_1, w[i], h, rs_, rss, rsum);
				compute(tx_1[i], ty_1, tx_w_1[i], ty_h_1, w[i], h, ts_, tss, tsum);

				// Compute the correlation
				// (sumXY - sumX*sumY/n) / sqrt( (sumXX - sumX^2 / n) * (sumYY - sumY^2 / n) )

				int n = w[i] * h;
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

					System.out.printf("Bad normalisation [%d,%d] = %g  (overlap=%g)\n", c, r, R, (double) n / size);
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
					int[] xyz1 = correlation.getXY(maxj);
					int[] xyz2 = correlation.getXY(base + c);
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

		// The maximum correlation with normalisation
		maxi = maxj; //correlation.findMaxIndex(ix, iy, iz, iw - ix, ih - iy, id - iz);
		xy = correlation.getXY(maxi);

		// Report the shift required to move from the centre of the target image to the reference
		// @formatter:off
		double[] result = new double[] {
			nc_2 - xy[0],
			nr_2 - xy[1],
			buffer[maxi]
		};
		// @formatter:on

		if (refinements > 0)
		{
			// Perform sub-pixel alignment.
			// Create a cubic spline using a small region of pixels around the maximum.
			// Avoid out-of-bounds errors. Only use the range that was normalised.
			int x = Maths.clip(ix, ixw - 5, xy[0] - 2);
			int y = Maths.clip(iy, iyh - 5, xy[1] - 2);
			DoubleImage2D crop = correlation.crop(x, y, 5, 5, region);
			FloatProcessor fp = new FloatProcessor(5, 5, crop.getData());

			// Find the maximum starting at the current origin
			int ox = xy[0] - x;
			int oy = xy[1] - y;

			double[] optimum = performCubicFit(fp, ox, oy, refinements, getRelativeThreshold());

			// Shift the result
			result[0] -= (optimum[0] - ox);
			result[1] -= (optimum[1] - oy);
			result[2] = optimum[2];
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
	private void checkCorrelation(DHTData target, DoubleDHT2D correlation, int maxi)
	{
		if (target.input == null || reference.input == null)
			// No check possible
			return;

		// The maximum correlation without normalisation 
		int[] xy = correlation.getXY(maxi);

		// Find the range for the target and reference
		int nc_2 = nc / 2;
		int nr_2 = nr / 2;
		int tx = Math.max(0, xy[0] - nc_2);
		int ty = Math.max(0, xy[1] - nr_2);
		int w = Math.min(nc, xy[0] + nc_2) - tx;
		int h = Math.min(nr, xy[1] + nr_2) - ty;

		// For the reference we express as a shift relative to the centre
		// and subtract the half-width.
		// Formally: (nc_2 - xy[0]) // shift 
		//           + nc_2          // centre
		//           - nc_2          // Half width
		int rx = Math.max(0, -xy[0] + nc_2);
		int ry = Math.max(0, -xy[1] + nr_2);

		double[] tar = target.input;
		double[] ref = reference.input;
		double o = correlation.get(maxi);
		double e = 0;
		for (int y = 0; y < h; y++)
		{
			int i = (ty + y) * nc + tx;
			int j = (ry + y) * nc + rx;
			for (int x = 0; x < w; x++)
			{
				e += tar[i++] * ref[j++];
			}
		}

		//System.out.printf("Raw %d,%d = %g\n", xy[0], xy[1], o);

		frequencyDomainCorrelationError = DoubleEquality.relativeError(o, e);
		if (frequencyDomainCorrelationError > 0.05)
		{
			System.err.printf("2D Correlation Error = %s : Spatial = %s, Freq = %s\n",
					Utils.rounded(frequencyDomainCorrelationError), Double.toString(e), Double.toString(o));
		}
	}

	/**
	 * Compute the sum from the rolling sum tables.
	 *
	 * @param x_1
	 *            the x value -1
	 * @param y_1
	 *            the y value -1
	 * @param x_w_1
	 *            the x value +w -1
	 * @param y_h_1
	 *            the y value +h -1
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param s_
	 *            the sum table
	 * @param ss
	 *            the sum-of-squares table
	 * @param sum
	 *            the sum (output = [sum, sum-of-squares])
	 */
	private void compute(int x_1, int y_1, int x_w_1, int y_h_1, int w, int h, double[] s_, double[] ss, double[] sum)
	{
		// Compute sum from rolling sum using:
		// sum(x,y,w,h) = 
		// + s(x+w-1,y+h-1) 
		// - s(x-1,y+h-1)
		// - s(x+w-1,y-1)
		// + s(x-1,y-1)
		// Note: 
		// s(i,j) = 0 when either i,j < 0
		// i = imax when i>imax 
		// j = jmax when j>jmax 

		// This has been adapted from Image2D to compute the twos sums together

		//int xw_yh = reference.dht.getIndex(x_w_1, y_h_1);
		int xw_yh = y_h_1 * nc + x_w_1;
		sum[0] = 0;
		sum[1] = 0;
		if (y_1 >= 0)
		{
			int h_ = h * nc;
			if (x_1 >= 0)
			{
				sum[0] = s_[xw_yh - w - h_] - s_[xw_yh - w];
				sum[1] = ss[xw_yh - w - h_] - ss[xw_yh - w];
			}
			sum[0] -= s_[xw_yh - h_];
			sum[1] -= ss[xw_yh - h_];
		}
		else if (x_1 >= 0)
		{
			sum[0] = -s_[xw_yh - w];
			sum[1] = -ss[xw_yh - w];
		}
		sum[0] = sum[0] + s_[xw_yh];
		sum[1] = sum[1] + ss[xw_yh];
	}

	/**
	 * Iteratively search the cubic spline surface around the given pixel
	 * to maximise the value.
	 * <p>
	 * At each round each of 8 points around the current maximum (+/- range) are evaluated. The optimum is picked and
	 * the range is halved. The initial range is 0.5 so the maximum distance that can be walked in any direction is 1
	 * pixel when the number of refinements is unlimited. With refinements = 3 the distance is 0.5 + 0.25 + 0.125 =
	 * 0.875.
	 *
	 * @param fp
	 *            Float processor containing a peak surface
	 * @param i
	 *            The peak x position
	 * @param j
	 *            The peak y position
	 * @param refinements
	 *            the maximum number of refinements
	 * @param relativeThreshold
	 *            Sets the relative threshold for change in the value for halting refinement. This applies
	 *            only if the position moved during the refinement step.
	 * @return The peak location with sub-pixel accuracy [x,y,value]
	 */
	public static double[] performCubicFit(FloatProcessor fp, int i, int j, int refinements, double relativeThreshold)
	{
		double[] centre = new double[] { i, j, fp.getf(i, j) };
		// Working space
		double[] xrange = new double[3];
		double[] yrange = new double[3];
		// This value will be progressively halved. 
		// Start with a value that allows the number of iterations to fully cover the region +/- 1 pixel
		// 0.5 will result in an minimum range of 0.5 / 2^9 = 0.000976
		double range = 0.5;
		while (refinements-- > 0)
		{
			double previous = centre[2];
			if (performCubicFit(fp, range, centre, xrange, yrange))
			{
				// The centre moved. Check convergence.
				if ((centre[2] - previous) / centre[2] < relativeThreshold)
					break;
			}
			range /= 2;
		}
		return centre;
	}

	/**
	 * Perform a cubic fit refinement.
	 *
	 * @param fp
	 *            Float processor containing a peak surface
	 * @param range
	 *            the range
	 * @param centre
	 *            the centre
	 * @param xrange
	 *            the xrange working space
	 * @param yrange
	 *            the yrange working space
	 * @return true, if the centre moved
	 */
	private static boolean performCubicFit(FloatProcessor fp, double range, double[] centre, double[] xrange,
			double[] yrange)
	{
		boolean moved = false;
		xrange[0] = centre[0] - range;
		xrange[1] = centre[0];
		xrange[2] = centre[0] + range;
		yrange[0] = centre[1] - range;
		yrange[1] = centre[1];
		yrange[2] = centre[1] + range;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (i == 1 && j == 1)
					// Current maximum
					continue;
				double v = fp.getBicubicInterpolatedPixel(xrange[i], yrange[j], fp);
				if (centre[2] < v)
				{
					centre[0] = xrange[i];
					centre[1] = yrange[j];
					centre[2] = v;
					moved = true;
				}
			}
		}
		return moved;
	}

	/**
	 * Iteratively search the cubic spline surface around the given pixel
	 * to maximise the value.
	 * <p>
	 * At each round each of 8 points around the current maximum (+/- range) are evaluated. The optimum is picked and
	 * the range is halved. The initial range is 0.5 so the maximum distance that can be walked in any direction is 1
	 * pixel when the number of refinements is unlimited. With refinements = 3 the distance is 0.5 + 0.25 + 0.125 =
	 * 0.875.
	 *
	 * @param fp
	 *            Float processor containing a peak surface
	 * @param i
	 *            The peak x position
	 * @param j
	 *            The peak y position
	 * @param refinements
	 *            the maximum number of refinements
	 * @param relativeThreshold
	 *            Sets the relative threshold for change in the value for halting refinement. This applies
	 *            only if the position moved during the refinement step.
	 * @return The peak location with sub-pixel accuracy [x,y,value]
	 */
	public static double[] performCubicSearch(FloatProcessor fp, int i, int j, int refinements, double relativeThreshold)
	{
		// TODO - implement this
		// We compute these dynamically as required.
		// Move this functionality to gdsc.core.math.interpolation.BicubicInterpolatingFunction
		
		CachedBicubicInterpolator[][] nodes = new CachedBicubicInterpolator[fp.getWidth()][fp.getHeight()];
		
		double[] centre = new double[] { i, j, fp.getf(i, j) };
		// Working space
		double[] xrange = new double[3];
		double[] yrange = new double[3];
		// This value will be progressively halved. 
		// Start with a value that allows the number of iterations to fully cover the region +/- 1 pixel
		// 0.5 will result in an minimum range of 0.5 / 2^9 = 0.000976
		double range = 0.5;
		while (refinements-- > 0)
		{
			double previous = centre[2];
			if (performCubicFit(fp, range, centre, xrange, yrange))
			{
				// The centre moved. Check convergence.
				if ((centre[2] - previous) / centre[2] < relativeThreshold)
					break;
			}
			range /= 2;
		}
		return centre;
	}	
	
	/**
	 * Copy the aligner. This copies the initialised state for use in alignment on multiple threads concurrently.
	 *
	 * @return the image aligner
	 */
	public Image2DAligner copy()
	{
		Image2DAligner copy;
		try
		{
			copy = (Image2DAligner) clone();
			// Reset objects that are not thread safe
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
	public Image2D getCorrelation()
	{
		try
		{
			DoubleImage2D image = new DoubleImage2D(nc, nr, buffer);
			image.fillOutside(crop[0], crop[1], crop[2], crop[3], 0);
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
}