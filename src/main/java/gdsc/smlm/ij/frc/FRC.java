package gdsc.smlm.ij.frc;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.logging.NullTrackProgress;
import gdsc.core.logging.TrackProgress;
import gdsc.core.utils.FloatEquality;
import gdsc.core.utils.Maths;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import ij.ImageStack;
import ij.process.FHT2;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Compute the Fourier Ring Correlation, a measure of the resolution of a microscopy image.
 * <p>
 * Adapted by Alex Herbert from the FIRE (Fourier Image REsolution) plugin produced as part of the paper:<br>
 * Niewenhuizen, et al (2013). Measuring image resolution in optical nanoscopy. Nature Methods, 10, 557<br>
 * http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2448.html
 * 
 * @author Alex Herbert
 * @author Bernd Rieger, b.rieger@tudelft.nl
 */
public class FRC
{
	public enum ThresholdMethod
	{
		//@formatter:off
		FIXED_1_OVER_7{ public String getName() { return "Fixed 1/7"; }}, 
		HALF_BIT{ public String getName() { return "Half-bit"; }}, 
		ONE_BIT{ public String getName() { return "One-bit"; }}, 
		TWO_BIT{ public String getName() { return "Two-bit"; }}, 
		ONE_SIGMA{ public String getName() { return "One sigma"; }},
		TWO_SIGMA{ public String getName() { return "Two sigma"; }},
		THREE_SIGMA{ public String getName() { return "Three sigma"; }},
		FOUR_SIGMA{ public String getName() { return "Four sigma"; }};
		//@formatter:on

		@Override
		public String toString()
		{
			return getName();
		}

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();
	};

	// Properties controlling the algorithm

	/**
	 * The correlation is computed using intervals around the circle circumference. The number of samples for half the
	 * circle is computed as: Pi * radius * sampling factor
	 */
	public double perimeterSamplingFactor = 1;

	/**
	 * The correlation is computed using intervals around the circle circumference of the Fourier transform. The Fourier
	 * image is 2-fold radially symmetric and so the calculation can use only half the circle for speed. Note: The
	 * results will differ slightly due to the implementation of the Fourier image not being exactly symmetric and the
	 * sample points used on the perimeter not matching between the two semi-circles.
	 */
	public boolean useHalfCircle = true;

	/** Used to track the progess within {@link #calculateFrcCurve(ImageProcessor, ImageProcessor)}. */
	public TrackProgress progress = null;

	private TrackProgress getTrackProgress()
	{
		if (progress == null)
			return new NullTrackProgress();
		return progress;
	}

	private static final double THIRD = 1.0 / 3.0;
	private static final double LAST_THIRD = 1.0 - 2 * THIRD;

	/**
	 * Calculate the Fourier Ring Correlation curve and numerator for two images.
	 * 
	 * @param ip1
	 *            The first image
	 * @param ip2
	 *            The second image
	 * @return An array of septuples representing [][radius,correlation,N,sum1,sum2,sum3,denominator] where:
	 *         radius is the distance from the centre of the Fourier transform;
	 *         correlation is the FRC at the given radius from the centre of the Fourier transform (i.e. 1/spatial
	 *         frequency);
	 *         N is the number of samples used to compute the correlation;
	 *         sum1 is the numerator of the correlation (the conjugate multiple of the two FFT images);
	 *         sum2 and sum3 are the two sums of the absolute FFT transform of each input image;
	 *         and denominator is Math.sqrt(sum2*sum3).
	 *         Note that correlation = sum1 / denominator.
	 */
	public double[][] calculateFrcCurve(ImageProcessor ip1, ImageProcessor ip2)
	{
		// Allow a progress tracker to be input
		TrackProgress progess = getTrackProgress();
		progess.incrementProgress(0);
		progess.status("Calculating complex FFT images...");

		// Pad images to the same size
		final int maxWidth = FastMath.max(ip1.getWidth(), ip2.getWidth());
		final int maxHeight = FastMath.max(ip1.getHeight(), ip2.getHeight());
		ip1 = pad(ip1, maxWidth, maxHeight);
		ip2 = pad(ip2, maxWidth, maxHeight);

		FloatProcessor[] fft1, fft2;

		boolean basic = false;
		if (basic)
		{
			fft1 = getComplexFFT(ip1);
			progess.incrementProgress(THIRD);
			fft2 = getComplexFFT(ip2);
			progess.incrementProgress(THIRD);
		}
		else
		{
			// Speed up by reusing the FHT object which performs pre-computation
			FHT2 fht = new FHT2();
			fht.setShowProgress(false);

			ip1 = getSquareTaperedImage(ip1);
			float[] f1 = (float[]) ip1.getPixels();
			fht.rc2DFHT(f1, false, ip1.getWidth());
			FHT2 fht1 = new FHT2(ip1, true);
			fft1 = getProcessors(fht1.getComplexTransform2());
			progess.incrementProgress(THIRD);

			ip2 = getSquareTaperedImage(ip2);
			float[] f2 = (float[]) ip2.getPixels();
			fht.rc2DFHT(f2, false, ip2.getWidth());
			FHT2 fht2 = new FHT2(ip2, true);
			fft2 = getProcessors(fht2.getComplexTransform2());
			progess.incrementProgress(THIRD);
		}

		progess.status("Preparing FRC curve calculation...");

		final int size = fft1[0].getWidth();
		final int centre = size / 2;

		// In-line for speed 
		float[] dataA1 = (float[]) fft1[0].getPixels();
		float[] dataB1 = (float[]) fft1[1].getPixels();
		float[] dataA2 = (float[]) fft2[0].getPixels();
		float[] dataB2 = (float[]) fft2[1].getPixels();
		float[] numerator = new float[dataA1.length];
		float[] absFFT1 = new float[dataA1.length];
		float[] absFFT2 = new float[dataA1.length];

		if (basic)
		{
			compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
		}
		else
		{
			computeMirroredFast(size, numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
		}

		int radius = 1;
		final int max = centre - 1;

		progess.status("Calculating FRC curve...");

		double[][] frcCurve = new double[(int) max][];

		// Radius zero is always 1. Set number of pixel to 1 when r==0.
		// Avoids divide by zero error.
		frcCurve[0] = new double[] { 0, 1, 1, 0, 0, 0, 0 };

		float[][] images = new float[][] { numerator, absFFT1, absFFT2 };

		final double limit = (useHalfCircle) ? Math.PI : 2 * Math.PI;
		while (radius < max)
		{
			// Inline the calculation for speed
			double sum1 = 0;
			double sum2 = 0;
			double sum3 = 0;

			final double angleStep = 1 / (perimeterSamplingFactor * radius);

			double angle = 0D;
			int numSum = 0;

			while (angle < limit)
			{
				double cosA = FastMath.cos(angle);
				double x = centre + radius * cosA;
				//double sinA = FastMath.sin(angle);
				double y = centre + radius * getSine(angle, cosA);
				double[] values = getInterpolatedValues(x, y, images, size);
				sum1 += values[0];
				sum2 += values[1];
				sum3 += values[2];

				numSum++;
				angle += angleStep;
			}

			double denominator = Math.sqrt(sum2 * sum3);
			double correlation = sum1 / denominator;

			frcCurve[radius] = new double[] { radius, correlation, numSum, sum1, sum2, sum3, denominator };

			radius++;
		}

		progess.incrementProgress(LAST_THIRD);
		progess.status("Finished calculating FRC curve...");

		return frcCurve;
	}

	// Package level to allow JUnit test
	static void compute(float[] numerator, float[] absFFT1, float[] absFFT2, float[] dataA1, float[] dataB1,
			float[] dataA2, float[] dataB2)
	{
		for (int i = dataA1.length; i-- > 0;)
		{
			compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i);
		}
	}

	// Package level to allow JUnit test
	static void computeMirrored(int size, float[] numerator, float[] absFFT1, float[] absFFT2, float[] dataA1,
			float[] dataB1, float[] dataA2, float[] dataB2)
	{
		// Note: Since this is symmetric around the centre we could compute half of it.
		// This is non-trivial since the centre is greater than half of the image, i.e.
		// not (size-1)/2;
		// So we compute up to the centre and copy back to the other half but must not miss
		// the edge pixels.
		final int centre = size / 2;

		// Do the first row, This is not mirrored
		int i = 0;
		while (i < size)
		{
			compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i++);
		}

		// Compute remaining rows up to the centre. These are mirrored
		int j = numerator.length - 1;
		for (int y = 1; y < centre; y++)
		{
			// The first entry in each row is not mirrored so compute and increment i
			compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i++);
			for (int x = 1; x < size; x++, i++, j--)
			{
				compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i);
				// Mirror
				numerator[j] = numerator[i];
				absFFT1[j] = absFFT1[i];
				absFFT2[j] = absFFT2[i];
			}
			// The last entry in each reverse row is not mirrored so compute and decrement j
			compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, j--);
		}

		// Do the centre row. This is mirrored with itself
		compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i++);
		for (int x = 1; x <= centre; x++, i++, j--)
		{
			compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i);
			// Mirror
			numerator[j] = numerator[i];
			absFFT1[j] = absFFT1[i];
			absFFT2[j] = absFFT2[i];
		}
	}

	// Package level to allow JUnit test
	static void computeMirroredFast(int size, float[] numerator, float[] absFFT1, float[] absFFT2, float[] dataA1,
			float[] dataB1, float[] dataA2, float[] dataB2)
	{
		// The same as computeMirrored but ignores the pixels that are not a mirror since
		// these are not used in the FRC calculation. 
		
		// Note: Since this is symmetric around the centre we could compute half of it.
		// This is non-trivial since the centre is greater than half of the image, i.e.
		// not (size-1)/2;
		// So we compute up to the centre and copy back to the other half.
		final int centre = size / 2;

		// Ignore the first row since this is not mirrored
		int i = size;

		// Compute remaining rows up to the centre. These are mirrored
		int j = numerator.length - 1;
		for (int y = 1; y < centre; y++)
		{
			// The first entry in each row is not mirrored so just increment i
			i++;
			for (int x = 1; x < size; x++, i++, j--)
			{
				compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i);
				// Mirror
				numerator[j] = numerator[i];
				absFFT1[j] = absFFT1[i];
				absFFT2[j] = absFFT2[i];
			}
			// The last entry in each reverse row is not mirrored so just decrement j
			j--;
		}

		// Do the centre row. This is mirrored with itself
		i++;
		for (int x = 1; x <= centre; x++, i++, j--)
		{
			compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2, i);
			// Mirror
			numerator[j] = numerator[i];
			absFFT1[j] = absFFT1[i];
			absFFT2[j] = absFFT2[i];
		}
	}

	private static void compute(float[] numerator, float[] absFFT1, float[] absFFT2, float[] dataA1, float[] dataB1,
			float[] dataA2, float[] dataB2, int i)
	{
		final float a1i = dataA1[i];
		final float a2i = dataA2[i];
		final float b1i = dataB1[i];
		final float b2i = dataB2[i];
		numerator[i] = a1i * a2i + b1i * b2i;
		absFFT1[i] = a1i * a1i + b1i * b1i;
		absFFT2[i] = a2i * a2i + b2i * b2i;
	}

	static boolean checkSymmetry(float[] data, int size)
	{
		// Symmetry is around the centre
		int centre = size / 2;

		float maxRelativeError = 1e-10f, maxAbsoluteError = 1e-16f;
		int error = 0;

		for (int y = centre, y2 = centre; y >= 0 && y2 < size; y--, y2++)
		{
			for (int x = centre, x2 = centre, i = size * y + x, j = size * y2 + x2; x >= 0 &&
					x2 < size; x--, x2++, i--, j++)
			{
				if (data[i] != data[j] || !FloatEquality.almostEqualRelativeOrAbsolute(data[i], data[j],
						maxRelativeError, maxAbsoluteError))
				{
					//System.out.printf("[%d] %f != [%d] %f\n", i, data[i], j, data[j]);
					if (--error < 0)
					{
						//gdsc.core.ij.Utils.display("check", new FloatProcessor(size, size, data));
						return false;
					}
				}
			}
		}
		return true;
	}

	/**
	 * Gets the sine of the angle given the cosine value.
	 *
	 * @param angle
	 *            the angle (in radians between 0 and 2*pi)
	 * @param cosA
	 *            the cosine of the angle
	 * @return the sine
	 */
	public static double getSine(double angle, double cosA)
	{
		return Math.sqrt(1 - (cosA * cosA)) * ((angle > Math.PI) ? -1 : 1); // Place in correct domain
	}

	private ImageProcessor pad(ImageProcessor ip, int width, int height)
	{
		if (ip.getWidth() != width || ip.getHeight() != height)
		{
			ImageProcessor ip2 = ip.createProcessor(width, height);
			ip2.insert(ip, 0, 0);
			ip = ip2;
		}
		return ip;
	}

	/**
	 * Convert an image into a Fourier image with real and imaginary parts
	 * 
	 * @param ip
	 * @return the real and imaginary parts
	 */
	public FloatProcessor[] getComplexFFT(ImageProcessor ip)
	{
		FloatProcessor taperedDataImage = getSquareTaperedImage(ip);

		FHT2 fht = new FHT2(taperedDataImage);
		fht.setShowProgress(false);
		fht.transform();

		ImageStack stack1 = fht.getComplexTransform();
		return getProcessors(stack1);
	}

	private FloatProcessor[] getProcessors(ImageStack stack1)
	{
		FloatProcessor[] ret = new FloatProcessor[2];
		ret[0] = ((FloatProcessor) stack1.getProcessor(1));
		ret[1] = ((FloatProcessor) stack1.getProcessor(2));
		return ret;
	}

	/**
	 * Applies a Tukey window function to the image and then pads it to the next square size power of two.
	 * 
	 * @param dataImage
	 * @return The square tapered image
	 */
	public FloatProcessor getSquareTaperedImage(ImageProcessor dataImage)
	{
		// Use a Tukey window function
		float[] taperX = getWindowFunctionX(dataImage.getWidth());
		float[] taperY = getWindowFunctionY(dataImage.getHeight());

		final int size = FastMath.max(dataImage.getWidth(), dataImage.getHeight());

		// Pad to a power of 2
		int newSize = 0;
		for (int i = 4; i < 15; i++)
		{
			newSize = (int) Math.pow(2.0, i);
			if (size <= newSize)
			{
				break;
			}
		}

		if (size > newSize)
			return null; // Error

		dataImage = dataImage.toFloat(0, null);
		float[] data = (float[]) dataImage.getPixels();
		float[] pixels = new float[newSize * newSize];
		// Note that the limits at 0 and size-1 the taper is zero so this can be ignored
		final int maxy_1 = dataImage.getHeight() - 1;
		final int maxx_1 = dataImage.getWidth() - 1;
		final int oldWidth = dataImage.getWidth();
		for (int y = 1; y < maxy_1; y++)
		{
			final float yTmp = taperY[y];
			for (int x = 1, i = y * oldWidth + 1, ii = y * newSize + 1; x < maxx_1; x++, i++, ii++)
			{
				pixels[ii] = data[i] * taperX[x] * yTmp;
			}
		}

		return new FloatProcessor(newSize, newSize, pixels, null);
	}

	// Cache the Tukey window function.
	// Using methods to check the length should make it thread safe since we create an instance reference
	// to an array of the correct length.
	private static float[] taperX = new float[0], taperY = new float[0];

	private static float[] getWindowFunctionX(int size)
	{
		float[] taper = getWindowFunction(taperX, size);
		taperX = taper;
		return taper;
	}

	private static float[] getWindowFunctionY(int size)
	{
		float[] taper = getWindowFunction(taperY, size);
		taperY = taper;
		return taper;
	}

	private static float[] getWindowFunction(float[] taper, int size)
	{
		if (taper.length != size)
		{
			// Re-use cached values
			taper = check(taperX, size);
			if (taper != null)
				return taper;
			taper = check(taperY, size);
			if (taper != null)
				return taper;
			taper = getTukeyWindowFunction(size);
		}
		return taper;
	}

	private static float[] check(float[] taper, int size)
	{
		return (taper.length == size) ? taper : null;
	}

	private static float[] getTukeyWindowFunction(int size)
	{
		float[] taper = new float[size];

		//// Original code. This created a non-symmetric window
		//final int boundary = size / 8;
		//final int upperBoundary = size - boundary;
		//for (int i = 0; i < size; i++)
		//{
		//	if ((i < boundary) || (i >= upperBoundary))
		//	{
		//		final double d = Math.sin(12.566370614359172D * i / size);
		//		taper[i] = (float) (d * d);
		//	}
		//	else
		//	{
		//		taper[i] = 1;
		//	}
		//}

		// New optimised code. This matches ImageWindow.tukey(size, 0.25) 
		final int boundary = size / 8;
		final int middle = size / 2;
		final double FOUR_PI_OVER_SIZE = 12.566370614359172D / (size - 1);
		{
			int i = 1, j = size - 2;
			while (i <= boundary)
			{
				final double d = Math.sin(FOUR_PI_OVER_SIZE * i);
				taper[i++] = taper[j--] = (float) (d * d);
			}
			while (i <= middle)
			{
				taper[i++] = taper[j--] = 1f;
			}
		}

		//// Use the GDSC ImageWindow class:
		//// FRC has Image width / Width of edge region = 8 so use alpha 0.25.
		//double[] w = gdsc.core.utils.ImageWindow.tukey(size, 0.25);
		//for (int i = 0; i < size; i++)
		//{
		//	System.out.printf("[%d] %f  %f\n", i, w[i], taper[i]);
		//	taper[i] = (float) w[i];
		//}

		return taper;
	}

	/**
	 * Adapted from ij.process.ImageProcessor.getInterpolatedValue(int,int).
	 * <p>
	 * Removed bounds checking and compute multiple values at the same time for multiple images.
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	private double[] getInterpolatedValues(final double x, final double y, float[][] images, final int maxx)
	{
		final int xbase = (int) x;
		final int ybase = (int) y;
		double xFraction = x - xbase;
		double yFraction = y - ybase;
		if (xFraction < 0.0)
			xFraction = 0.0;
		if (yFraction < 0.0)
			yFraction = 0.0;

		final int lowerLeftIndex = ybase * maxx + xbase;
		final int lowerRightIndex = lowerLeftIndex + 1;
		final int upperLeftIndex = lowerLeftIndex + maxx;
		final int upperRightIndex = upperLeftIndex + 1;

		final int noImages = 3; //images.length;
		double[] values = new double[noImages];
		for (int i = 0; i < noImages; i++)
		{
			final float[] image = images[i];
			final double lowerLeft = image[lowerLeftIndex];
			final double lowerRight = image[lowerRightIndex];
			final double upperRight = image[upperLeftIndex];
			final double upperLeft = image[upperRightIndex];

			final double upperAverage = upperLeft + xFraction * (upperRight - upperLeft);
			final double lowerAverage = lowerLeft + xFraction * (lowerRight - lowerLeft);
			values[i] = lowerAverage + yFraction * (upperAverage - lowerAverage);
		}
		return values;
	}

	/**
	 * Perform LOESS smoothing on the FRC curve.
	 * <p>
	 * The input curve is copied and then the correlation values are smoothed using a LOESS interpolation with the given
	 * parameters. If smoothing fails the original curve values are returned.
	 * 
	 * @param frcCurve
	 * @param bandwidth
	 * @param robustness
	 * @return A new FRC curve
	 */
	public static double[][] getSmoothedCurve(double[][] frcCurve, double bandwidth, int robustness)
	{
		double[][] sCurve = new double[frcCurve.length][3];

		double[] xVals = new double[frcCurve.length];
		double[] yVals = new double[frcCurve.length];

		for (int i = 0; i < frcCurve.length; i++)
		{
			xVals[i] = frcCurve[i][0];
			yVals[i] = frcCurve[i][1];

			sCurve[i][0] = frcCurve[i][0];
			sCurve[i][1] = frcCurve[i][1];
			sCurve[i][2] = frcCurve[i][2];
		}

		double[] ySmoothed = new double[frcCurve.length];

		try
		{
			LoessInterpolator loess = new LoessInterpolator(bandwidth, robustness);
			ySmoothed = loess.smooth(xVals, yVals);

			for (int i = 0; i < frcCurve.length; i++)
				sCurve[i][1] = ySmoothed[i];
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}

		return sCurve;
	}

	/**
	 * Perform LOESS smoothing on the FRC curve
	 * <p>
	 * The input curve is copied and then the correlation values are smoothed using a LOESS interpolation with a
	 * bandwidth of 0.0707 and robustness of 0. If smoothing fails the original curve values are returned.
	 * 
	 * @param frcCurve
	 * @return A new FRC curve
	 */
	public static double[][] getSmoothedCurve(double[][] frcCurve)
	{
		double bandwidth = 0.0707;
		int robustness = 0;
		return getSmoothedCurve(frcCurve, bandwidth, robustness);
	}

	private final static double TWO_PI = 2.0 * Math.PI;

	/**
	 * Calculate the curve representing the minimum correlation required to distinguish two images for each resolution
	 * in the input FRC curve.
	 * 
	 * @param frcCurve
	 * @param method
	 * @return The threshold curve representing the threshold for each input spatial frequency
	 */
	public static double[] calculateThresholdCurve(double[][] frcCurve, ThresholdMethod method)
	{
		double[] threshold = new double[frcCurve.length];

		// ADH: 
		// Half-Bit and 3 Sigma are explained in Supp section 5.4. However this is based on 
		// Heel, et al (2005) which gives a better explanation so I have updated using their 
		// equations.
		// Equation S.84 has an error compared to equation (13) in Heel. In fact equation (17) 
		// from Heel is what should be implemented for Half-bit.

		// Note: The original code used frcCurve[i][2] which holds the number of samples that 
		// were taken from the circle. This makes the curve dependent on the number of samples 
		// taken (e.g. half-circle/full-circle with different sampling factors).
		// To make the curve sampling independent I assume 2*pi*r samples were taken. 

		// See: Heel, M. v. & Schatz, M. Fourier shell correlation threshold criteria. J. Struct. Bio. 151, 250–262 (2005).

		// Fourier Shell Radius:
		// This is set to 1 for r==0 in Heel
		double nr = 1;

		int sigma = 0; // For the N-sigma methods

		switch (method)
		{
			case FIXED_1_OVER_7:
				Arrays.fill(threshold, 1.0 / 7.0);
				break;

			// Note: The bit curve approach unity when r==0, i.e. nr=1	

			case HALF_BIT:
				// This is actually equation (17) from Heel:
				calculateBitCurve(threshold, 0.5);
				break;

			case ONE_BIT:
				// This is equation (14) from Heel:
				calculateBitCurve(threshold, 1);
				break;

			case TWO_BIT:
				calculateBitCurve(threshold, 2);
				break;

			// We use fall through here to increment sigma to the appropriate level.
			case FOUR_SIGMA:
				sigma++;
			case THREE_SIGMA:
				sigma++;
			case TWO_SIGMA:
				sigma++;
			case ONE_SIGMA:
				sigma++;
				for (int i = 0; i < threshold.length; i++, nr = TWO_PI * i)
				{
					// Note: frcCurve[i][2] holds the number of samples that were taken from the circle.
					//threshold[i] = (3.0 / Math.sqrt(frcCurve[i][2] / 2.0));

					// Heel, Equation (2):
					// We actually want to know the number of pixels contained in the Fourier shell of radius r.
					// We can compute this assuming we sampled the full circle 2*pi*r.
					threshold[i] = sigma / Math.sqrt(nr / 2.0);
				}
				break;

			default:
				Arrays.fill(threshold, 1.0 / 7.0);
		}

		return threshold;
	}

	/**
	 * Compute the threshold curve for the given number of bits
	 * 
	 * @param threshold
	 * @param bits
	 */
	private static void calculateBitCurve(final double[] threshold, double bits)
	{
		// This is adapted from equation (13) from Heel
		// See: Heel, M. v. & Schatz, M. Fourier shell correlation threshold criteria. J. Struct. Bio. 151, 250–262 (2005).

		// Approach unity when r -> 0:
		threshold[0] = 1;

		// Find the SNR in each half of the dataset: 
		// "because the total reconstruction, being the sum of the two half-data 
		// set reconstructions, will have twice the SNR value of each of the half 
		// data sets"
		// Eq. (15) = log2(SNR+1) = n-bits
		final double snr = (Math.pow(2, bits) - 1) / 2;
		final double snr1 = snr + 1;
		final double twoRootSnr = 2 * Math.sqrt(snr);
		final double twoRootSnr1 = twoRootSnr + 1;

		// Sense check: 
		// 1/2-bit is equation (17) from Heel:
		// snr = 0.2071, twoRootSnr = 0.9102
		// 1-bit is equation (14) from Heel:
		// snr = 0.5, twoRootSnr = 1.4142

		for (int i = 1; i < threshold.length; i++)
		{
			// nr = number of samples in Fourier circle = 2*pi*r
			final double sqrtNr = Math.sqrt(TWO_PI * i);
			threshold[i] = ((snr + twoRootSnr1 / sqrtNr) / (snr1 + twoRootSnr / sqrtNr));
		}
	}

	/**
	 * Computes the crossing points of the FRC curve and the threshold curve. The intersections can be used to determine
	 * the image resolution using {@link #getCorrectIntersection(ArrayList, ThresholdMethod)}
	 * 
	 * @param frcCurve
	 * @param thresholdCurve
	 * @param max
	 *            The maximum number of intersections to compute
	 * @return The crossing points
	 */
	public static double[][] getIntersections(double[][] frcCurve, double[] thresholdCurve, int max)
	{
		if (frcCurve.length != thresholdCurve.length)
		{
			System.err.println("Error: Unable to calculate FRC curve intersections due to input length mismatch.");
			return null;
		}

		double[][] intersections = new double[Math.min(max, frcCurve.length - 1)][];
		int count = 0;

		for (int i = 1; i < frcCurve.length; i++)
		{
			// http://en.wikipedia.org/wiki/Line-line_intersection
			//
			//     x1,y1            x4,y4      
			//         **        ++ 
			//           **    ++
			//             **++ P(x,y)
			//            ++ **
			//          ++     **
			//        ++         **
			//    x3,y3            ** 
			//                       x2,y2  

			final double y1 = frcCurve[i - 1][1];
			final double y2 = frcCurve[i][1];
			final double y3 = thresholdCurve[i - 1];
			final double y4 = thresholdCurve[i];

			// Check if they cross
			if (!((y3 >= y1 && y4 < y2) || (y1 >= y3 && y2 < y4)))
			{
				continue;
			}

			final double x1 = frcCurve[i - 1][0];
			final double x2 = frcCurve[i][0];
			final double x3 = x1;
			final double x4 = x2;

			final double x1_x2 = x1 - x2;
			final double x3_x4 = x3 - x4;
			final double y1_y2 = y1 - y2;
			final double y3_y4 = y3 - y4;

			// Check if lines are parallel
			if (x1_x2 * y3_y4 - y1_y2 * x3_x4 == 0)
			{
				if (y1 == y3)
					// The lines are the same
					intersections[count++] = new double[] { x1, y1 };
			}
			else
			{
				// Find intersection
				double px = ((x1 * y2 - y1 * x2) * x3_x4 - x1_x2 * (x3 * y4 - y3 * x4)) /
						(x1_x2 * y3_y4 - y1_y2 * x3_x4);
				//double px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) /
				//		((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
				//double py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) /
				//		((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));

				// Check if the intersection is within the two points
				// Q. Is this necessary given the intersection check above?
				if (px >= x1 && px < x2)
				{
					double py = Maths.interpolateY(x1, y3, x2, y4, px);
					intersections[count++] = new double[] { px, py };
				}
			}
			if (count >= max)
				break;
		}

		return Arrays.copyOf(intersections, count);
	}

	/**
	 * Get the correction intersection representing the image resolution. The intersection chosen depends on the method
	 * used to calculate the threshold curve using {@link #calculateThresholdCurve(double[][], ThresholdMethod)}
	 * <p>
	 * The intersection corresponds the lowest spatial frequency at which there is no significant correlation between
	 * the images.
	 * 
	 * @param intersections
	 * @param method
	 * @return The intersection (or null if no crossings)
	 */
	public static double[] getCorrectIntersection(double[][] intersections, ThresholdMethod method)
	{
		if (intersections == null || intersections.length == 0)
			return null;

		int pos = 0;
		switch (method)
		{
			case FIXED_1_OVER_7:
				// always use the first intersection
				break;

			// The N-sigma curves are above 1 at close to zero spatial frequency.
			// The bit curves are 1 at zero spatial frequency.
			// This means that any FRC curve starting around 1 (due to smoothing)
			// may cross the line twice. 
			// If so the second crossing is the one that is desired.
			default:
				if (intersections.length == 2)
				{
					// Choose the intersection. Just discard an intersection with a correlation above 0.9
					if (intersections[0][1] > 0.9)
						pos++;
				}
		}

		return intersections[pos];
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided images.
	 * 
	 * @param ip1
	 * @param ip2
	 * @param method
	 * @return The FIRE number (in pixels) and the correlation
	 */
	public double calculateFireNumber(ImageProcessor ip1, ImageProcessor ip2, ThresholdMethod method)
	{
		double[][] frcCurve = calculateFrcCurve(ip1, ip2);
		return calculateFireNumber(frcCurve, method);
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided FRC curve data.
	 * 
	 * @param frcCurve
	 * @param method
	 * @return The FIRE number (in pixels)
	 */
	public static double calculateFireNumber(double[][] frcCurve, ThresholdMethod method)
	{
		double[] result = calculateFire(frcCurve, method);
		if (result == null)
			return Double.NaN;
		return result[0];
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided images.
	 * 
	 * @param ip1
	 * @param ip2
	 * @param method
	 * @return The FIRE number (in pixels) and the correlation (null if computation failed)
	 */
	public double[] calculateFire(ImageProcessor ip1, ImageProcessor ip2, ThresholdMethod method)
	{
		double[][] frcCurve = calculateFrcCurve(ip1, ip2);
		return calculateFire(frcCurve, method);
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided FRC curve data.
	 * 
	 * @param frcCurve
	 * @param method
	 * @return The FIRE number (in pixels) and the correlation (null if computation failed)
	 */
	public static double[] calculateFire(double[][] frcCurve, ThresholdMethod method)
	{
		double[] thresholdCurve = calculateThresholdCurve(frcCurve, method);
		double[][] intersections = getIntersections(frcCurve, thresholdCurve, 2);

		if (intersections != null && intersections.length != 0)
		{
			double[] intersection = getCorrectIntersection(intersections, method);
			// Since the Fourier calculation only uses half of the image (from centre to the edge) 
			// we must double the curve length to get the original maximum image width. In addition
			// the computation was up to the edge-1 pixels so add back a pixel to the curve length.
			double fire = 2 * (frcCurve.length + 1) / intersection[0];
			return new double[] { fire, intersection[1] };
		}
		return null;
	}

	/**
	 * Apply spurious correlation correction using the Q-factor. Follows the method described in Niewenhuizen, et al
	 * (2013), on-line methods. This method computes a value that is subtracted from the numerator of the FRC and added
	 * to the denominator of the FRC before computing the correlation. The correlation in the FRC curve will be updated.
	 * The sums of the FRC curve are unchanged.
	 * <p>
	 * The correlation can be reset by calling this method with a Q-value of zero.
	 * <p>
	 * Note: Spurious correlation correction is only useful when computing the resolution of a single set of
	 * localisations split into two images. In this case the same emitter can be present in both images leading to
	 * spurious contribution to the correlation. Correction can be omitted providing the number of emitters is
	 * sufficiently high and the sample has spectral signal content at the computed resolution for the uncorrected
	 * curve (see Niewenhuizen, et al (2013), Nature Methods, 10, 557, Supplementary Material p.22).
	 *
	 * @param frcCurve
	 *            the frc curve
	 * @param nmPerPixel
	 *            the nm per pixel in the images used to compute the FRC curve
	 * @param qValue
	 *            the q value
	 * @param mean
	 *            the mean of the localisation precision
	 * @param sigma
	 *            the standard deviation of the localisation precision
	 */
	public static void applyQCorrection(double[][] frcCurve, double nmPerPixel, double qValue, double mean,
			double sigma)
	{
		if (qValue <= 0)
		{
			// Reset the correlation 
			for (int i = 1; i < frcCurve.length; i++)
				frcCurve[i][1] = frcCurve[i][3] / frcCurve[i][6];
			return;
		}

		double[] q = computeQ(frcCurve, nmPerPixel);

		// H(q) is the factor in the correlation averages related to the localization
		// uncertainties that depends on the mean and width of the
		// distribution of localization uncertainties
		double[] hq = computeHq(q, mean, sigma);

		// Subtract the average residual correlation from the numerator and add to the denominator 
		for (int i = 1; i < frcCurve.length; i++)
		{
			double numerator = frcCurve[i][3];
			double denominator = frcCurve[i][6];
			double residual = qValue * hq[i];
			frcCurve[i][1] = (numerator - residual) / (denominator + residual);
		}
	}

	/**
	 * Compute the size of the field of view for the FRC curve (this is named L)
	 *
	 * @param frcCurve
	 *            the frc curve
	 * @return L
	 */
	public static int computeL(double[][] frcCurve)
	{
		// Since the Fourier calculation only uses half of the image (from centre to the edge) 
		// we must double the curve length to get the original maximum image width. In addition
		// the computation was up to the edge-1 pixels so add back a pixel to the curve length.
		// Note: frcCurveLength == L == Size of field of view.
		return ((int) (frcCurve[(frcCurve.length - 1)][0]) + 1) * 2;
	}

	/**
	 * Compute q. This is defined as 1/L, 2/L, ..., for all the spatial frequencies in the FRC curve where L is the size
	 * of the field of view. This is converted to nm using the pixel size of the input image.
	 *
	 * @param frcCurve
	 *            the frc curve
	 * @param nmPerPixel
	 *            the nm per pixel in the images used to compute the FRC curve
	 * @return the q array
	 */
	public static double[] computeQ(double[][] frcCurve, double nmPerPixel)
	{
		final double L = computeL(frcCurve);

		double[] q = new double[frcCurve.length];
		double conversion = 1.0 / (L * nmPerPixel);
		for (int i = 0; i < q.length; i++)
		{
			final double radius = frcCurve[i][0];
			q[i] = radius * conversion;
		}
		return q;
	}

	/**
	 * Compute the localization PDF factor H(q) for all q. This is the integral of the distribution function of the
	 * localisation uncertainty. It is assumed to Gaussian with the specified mean and width.
	 *
	 * @param q
	 *            the q array
	 * @param mean
	 *            the mean of the Gaussian
	 * @param sigma
	 *            the width of the Gaussian
	 * @return the Hq array
	 */
	public static double[] computeHq(double[] q, double mean, double sigma)
	{
		// H(q) is the factor in the correlation averages related to the localization
		// uncertainties that depends on the mean and width of the
		// distribution of localization uncertainties
		double[] hq = new double[q.length];
		final double four_pi2 = 4 * Math.PI * Math.PI;
		double eight_pi2_s2 = 2 * four_pi2 * sigma * sigma;
		hq[0] = 1; // TODO - what should this be at zero?
		for (int i = 1; i < q.length; i++)
		{
			// Q. Should q be in the correct units?
			double q2 = q[i] * q[i];
			//double q2 = i * i;
			double d = 1 + eight_pi2_s2 * q2;
			hq[i] = FastMath.exp((-four_pi2 * mean * mean * q2) / d) / Math.sqrt(d);
		}
		return hq;
	}
}
