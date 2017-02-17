package gdsc.smlm.ij.frc;

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

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.logging.NullTrackProgress;
import gdsc.core.logging.TrackProgress;
import gdsc.core.utils.ImageWindow;

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
		THREE_SIGMA{ public String getName() { return "Three sigma"; }};
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
	 * @param ip2
	 * @return An array of sextuples representing [][radius,correlation,N,sum1,sum2,sum3] where correlation is the FRC
	 *         at the given radius from the centre of the Fourier transform (i.e. 1/spatial frequency) and N is the
	 *         number of samples used to compute the correlation. sum1 is the numerator of the correlation (the
	 *         conjugate multiple of the two FFT images), sum2 and sum3 are the two sums of the absolute FFT transform
	 *         of each input image. correlation = sum1 / Math.sqrt(sum2*sum3).
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

		// In-line for speed
		float[] numerator = new float[size * size];
		float[] absFFT1 = new float[numerator.length];
		float[] absFFT2 = new float[numerator.length];

		float[] dataA1 = (float[]) fft1[0].getPixels();
		float[] dataB1 = (float[]) fft1[1].getPixels();
		float[] dataA2 = (float[]) fft2[0].getPixels();
		float[] dataB2 = (float[]) fft2[1].getPixels();

		for (int y = size, i = size * size - 1; y-- > 0;)
		{
			for (int x = size; x-- > 0; i--)
			{
				numerator[i] = dataA1[i] * dataA2[i] + dataB1[i] * dataB2[i];
				absFFT1[i] = dataA1[i] * dataA1[i] + dataB1[i] * dataB1[i];
				absFFT2[i] = dataA2[i] * dataA2[i] + dataB2[i] * dataB2[i];
			}
		}

		int radius = 1;
		final int centre = size / 2;
		final int max = centre - 1;

		progess.status("Calculating FRC curve...");

		double[][] frcCurve = new double[(int) max][];

		// Radius zero is always 1
		// Avoid divide by zero errors. Not sure if this is OK
		frcCurve[0] = new double[] { 0, 1, 1, 0, 0, 0 };

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
				double x = centre + radius * Math.cos(angle);
				double y = centre + radius * Math.sin(angle);

				double[] values = getInterpolatedValues(x, y, images, size);
				sum1 += values[0];
				sum2 += values[1];
				sum3 += values[2];

				numSum++;
				angle += angleStep;
			}

			double val = sum1 / Math.sqrt(sum2 * sum3);

			frcCurve[radius] = new double[] { radius, val, numSum, sum1, sum2, sum3 };

			radius++;
		}

		progess.incrementProgress(LAST_THIRD);
		progess.status("Finished calculating FRC curve...");

		return frcCurve;
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

		// Use the GDSC ImageWindow class:
		// FRC has Image width / Width of edge region = 8 so use alpha 0.25.
		double[] w = ImageWindow.tukey(size, 0.25);
		for (int i = 0; i < size; i++)
		{
			//System.out.printf("[%d] %f  %f\n", i, w[i], taper[i]);
			taper[i] = (float) w[i];
		}

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
	public double[][] getSmoothedCurve(double[][] frcCurve, double bandwidth, int robustness)
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
	public double[][] getSmoothedCurve(double[][] frcCurve)
	{
		double bandwidth = 0.0707;
		int robustness = 0;
		return getSmoothedCurve(frcCurve, bandwidth, robustness);
	}

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

		switch (method)
		{
			case HALF_BIT:
				for (int i = 0; i < threshold.length; i++)
					threshold[i] = ((0.2071 * Math.sqrt(frcCurve[i][2]) + 1.9102) /
							(1.2701 * Math.sqrt(frcCurve[i][2]) + 0.9102));
				break;

			case THREE_SIGMA:
				for (int i = 0; i < threshold.length; i++)
					threshold[i] = (3.0 / Math.sqrt(frcCurve[i][2] / 2.0));
				break;

			case FIXED_1_OVER_7:
			default:
				Arrays.fill(threshold, 1.0 / 7.0);
		}

		return threshold;
	}

	/**
	 * Computes the crossing points of the FRC curve and the threshold curve. The intersections can be used to determine
	 * the image resolution using {@link #getCorrectIntersection(ArrayList, ThresholdMethod)}
	 * 
	 * @param frcCurve
	 * @param thresholdCurve
	 * @return The crossing points
	 */
	public double[] getIntersections(double[][] frcCurve, double[] thresholdCurve)
	{
		if (frcCurve.length != thresholdCurve.length)
		{
			getTrackProgress().log("Error: Unable to calculate FRC curve intersections due to input length mismatch.");
			return null;
		}

		double[] intersections = new double[frcCurve.length - 1];
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
					intersections[count++] = x1;
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
					intersections[count++] = px;
				}
			}
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
	 * @return The intersection (or zero if no crossings)
	 */
	public double getCorrectIntersection(double[] intersections, ThresholdMethod method)
	{
		if (intersections == null || intersections.length == 0)
			return 0;

		switch (method)
		{
			// The half-bit and 3-sigma curves are often above 1 at close to zero spatial frequency.
			// This means that any FRC curve starting at 1 may cross the line twice. 
			// If so the second crossing is the one that is desired.
			case HALF_BIT:
			case THREE_SIGMA:
				return (intersections.length > 1) ? intersections[1] : intersections[0];

			case FIXED_1_OVER_7:
			default:
				return intersections[0];
		}
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided images.
	 * 
	 * @param ip1
	 * @param ip2
	 * @param method
	 * @return The FIRE number (in pixels)
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
	public double calculateFireNumber(double[][] frcCurve, ThresholdMethod method)
	{
		double[] thresholdCurve = calculateThresholdCurve(frcCurve, method);
		double[] intersections = getIntersections(frcCurve, thresholdCurve);

		double fire = Double.NaN;
		if (intersections == null || intersections.length != 0)
		{
			double spatialFrequency = getCorrectIntersection(intersections, method);
			// Since the Fourier calculation only uses half of the image (from centre to the edge) 
			// we must double the curve length to get the original maximum image width. In addition
			// the computation was up to the edge-1 pixels so add back a pixel to the curve length.
			fire = 2 * (frcCurve.length + 1) / spatialFrequency;
		}
		return fire;
	}
}
