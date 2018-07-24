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
package uk.ac.sussex.gdsc.smlm.ij.frc;

import java.util.Arrays;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.util.FastMath;
import org.jtransforms.fft.FloatFFT_2D;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import uk.ac.sussex.gdsc.core.ij.process.FHT2;
import uk.ac.sussex.gdsc.core.logging.NullTrackProgress;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.math.RadialStatistics;
import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.Maths;

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
	/** Flag indicating if the JTransforms library is available. */
	private static boolean JTRANSFORMS;
	static
	{
		try
		{
			final int size = 8;
			final FloatFFT_2D fft = new FloatFFT_2D(size, size);
			final float[] data = new float[size * size * 2];
			fft.realForwardFull(data);
			JTRANSFORMS = true;
		}
		catch (final Throwable t)
		{
			JTRANSFORMS = false;
		}
	}

	/**
	 * Specify the method to create the FRC threshold curve. The intersection between the observed curve and the
	 * threshold curve determines the resolution.
	 */
	public enum ThresholdMethod
	{
		/** The fixed level threshold using 1 over 7. */
		//@formatter:off
		FIXED_1_OVER_7{ @Override
		public String getName() { return "Fixed 1/7"; }},

		/** The half bit threshold. */
		HALF_BIT{ @Override
		public String getName() { return "Half-bit"; }},

		/** The one bit threshold. */
		ONE_BIT{ @Override
		public String getName() { return "One-bit"; }},

		/** The two bit threshold. */
		TWO_BIT{ @Override
		public String getName() { return "Two-bit"; }},

		/** The one sigma threshold. */
		ONE_SIGMA{ @Override
		public String getName() { return "One sigma"; }},

		/** The two sigma threshold. */
		TWO_SIGMA{ @Override
		public String getName() { return "Two sigma"; }},

		/** The three sigma threshold. */
		THREE_SIGMA{ @Override
		public String getName() { return "Three sigma"; }},

		/** The four sigma threshold. */
		FOUR_SIGMA{ @Override
		public String getName() { return "Four sigma"; }};
		//@formatter:on

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Enum#toString()
		 */
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
	}

	/**
	 * Specify the sampling method to compute the Fourier Ring Correlation.
	 */
	public enum SamplingMethod
	{
		//@formatter:off
		/**
		 * Compute using all pixels from radius n (inclusive) to n+1 (exclusive) assigned to ring n.
		 * This does not require interpolation.
		 */
		RADIAL_SUM{ @Override
		public String getName() { return "Radial Sum"; }},
		/**
		 * Compute using sampled points on a circle circumference of radius n for each ring.
		 * Values are computed using interpolation of the surrounding pixels. The number of points
		 * on the circumference can be controlled using the perimeter sampling factor.
		 */
		INTERPOLATED_CIRCLE{ @Override
		public String getName() { return "Interpolated Circle"; }};
		//@formatter:on

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Enum#toString()
		 */
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
	}

	/**
	 * Specify the method used to compute the Fourier transform.
	 */
	public enum FourierMethod
	{
		//@formatter:off
		/** Use the JTransforms Java library. */
		JTRANSFORMS{ @Override
		public String getName() { return "JTransforms"; }},

		/** Use ImageJ's Fast Hartley Transform. */
		FHT{ @Override
		public String getName() { return "FHT"; }};
		//@formatter:on

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Enum#toString()
		 */
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
	}

	/**
	 * Contains the result of a single ring from the Fourier Ring Correlation (FRC) between two images.
	 */
	public static class FRCCurveResult implements Cloneable
	{
		/** The radius of the ring. */
		private final int radius;

		/** The number of samples on the ring. */
		private final int nSamples;

		/** The denominator. */
		private final double numerator, sum1, sum2, denominator;

		/** The correlation. */
		private double correlation;

		/**
		 * Instantiates a new FRC curve result.
		 *
		 * @param radius
		 *            the radius
		 * @param nSamples
		 *            the n samples
		 * @param sum0
		 *            the sum 0
		 * @param sum1
		 *            the sum 1
		 * @param sum2
		 *            the sum 2
		 */
		private FRCCurveResult(int radius, int nSamples, double sum0, double sum1, double sum2)
		{
			this.radius = radius;
			this.nSamples = nSamples;
			this.numerator = sum0;
			this.sum1 = sum1;
			this.sum2 = sum2;
			denominator = Math.sqrt(sum1 * sum2);
			setCorrelation(numerator / denominator);
		}

		/**
		 * Gets the radius of the ring.
		 *
		 * @return the radius
		 */
		public int getRadius()
		{
			return radius;
		}

		/**
		 * Gets the number of samples taken from the ring used to compute the correlation.
		 *
		 * @return the number of samples
		 */
		public int getNumberOfSamples()
		{
			return nSamples;
		}

		/**
		 * Gets the correlation (the normalised conjugate multiple of the two FFT images). This is the numerator divided
		 * by the denominator.
		 *
		 * @return the correlation
		 */
		public double getCorrelation()
		{
			return correlation;
		}

		/**
		 * Sets the correlation.
		 *
		 * @param correlation
		 *            the new correlation
		 */
		private void setCorrelation(double correlation)
		{
			this.correlation = Maths.clip(-1, 1, correlation);
		}

		/**
		 * Gets the numerator of the correlation. This is the sum of the conjugate multiple of the FFT of input image 1
		 * and 2.
		 *
		 * @return the numerator
		 */
		public double getNumerator()
		{
			return numerator;
		}

		/**
		 * Gets the denominator of the correlation.
		 * <p>
		 * Note: the denominator is Math.sqrt(sum1*sum2)
		 *
		 * @return the denominator
		 */
		public double getDenominator()
		{
			return denominator;
		}

		/**
		 * Gets the sum of the absolute FFT of input image 1.
		 *
		 * @return the sum
		 */
		public double getSum1()
		{
			return sum1;
		}

		/**
		 * Gets the sum of the absolute FFT of input image 2.
		 *
		 * @return the sum
		 */
		public double getSum2()
		{
			return sum2;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Object#clone()
		 */
		@Override
		public FRCCurveResult clone()
		{
			try
			{
				return (FRCCurveResult) super.clone();
			}
			catch (final CloneNotSupportedException e)
			{
				return null; // This should not happen
			}
		}
	}

	/**
	 * Contains the result of computing all the rings of the Fourier Ring Correlation (FRC) between two images.
	 */
	public static class FRCCurve implements Cloneable
	{
		/** The nm per pixel for the super-resolution images. */
		final public double nmPerPixel;

		/** The size of the field of view in the Fourier image (named L). */
		final public int fieldOfView;

		/** The mean of the first input image after application of the Tukey window taper. */
		final public double mean1;
		/** The mean of the second input image after application of the Tukey window taper. */
		final public double mean2;

		/** The results. */
		private FRCCurveResult[] results;

		/**
		 * Instantiates a new FRC curve.
		 *
		 * @param nmPerPixel
		 *            the nm per pixel
		 * @param fieldOfView
		 *            the field of view
		 * @param mean1
		 *            the mean 1
		 * @param mean2
		 *            the mean 2
		 * @param results
		 *            the results
		 */
		private FRCCurve(double nmPerPixel, int fieldOfView, double mean1, double mean2, FRCCurveResult[] results)
		{
			this.nmPerPixel = nmPerPixel;
			this.fieldOfView = fieldOfView;
			this.mean1 = mean1;
			this.mean2 = mean2;
			this.results = results;
		}

		/**
		 * Gets the FRC curve result.
		 *
		 * @param index
		 *            the index
		 * @return the FRC curve result
		 */
		public FRCCurveResult get(int index)
		{
			return results[index];
		}

		/**
		 * Gets the number of results.
		 *
		 * @return the size
		 */
		public int getSize()
		{
			return results.length;
		}

		/**
		 * Return a deep copy of this curve.
		 *
		 * @return the copy
		 */
		public FRCCurve copy()
		{
			// Let the Java framework copy the primitives
			final FRCCurve curve = this.clone();
			// Clone the curve entries
			curve.results = new FRCCurveResult[results.length];
			for (int i = results.length; i-- > 0;)
				curve.results[i] = results[i].clone();
			return curve;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Object#clone()
		 */
		@Override
		protected FRCCurve clone()
		{
			try
			{
				return (FRCCurve) super.clone();
			}
			catch (final CloneNotSupportedException e)
			{
				return null; // This should not happen
			}
		}

		/**
		 * Gets the radius values for each ring.
		 *
		 * @return the radius values
		 */
		public double[] getRadiusValues()
		{
			final double[] values = new double[getSize()];
			for (int i = values.length; i-- > 0;)
				values[i] = get(i).getRadius();
			return values;
		}

		/**
		 * Gets the correlation values for each ring.
		 *
		 * @return the correlation values
		 */
		public double[] getCorrelationValues()
		{
			final double[] values = new double[getSize()];
			for (int i = values.length; i-- > 0;)
				values[i] = get(i).getCorrelation();
			return values;
		}
	}

	/**
	 * Contains the Fourier Image Resolution (FIRE) result computed from the intersection of the FRC curve and a
	 * threshold curve.
	 */
	public static class FIREResult
	{
		/** The fire number (in nm). */
		final public double fireNumber;

		/** The correlation. */
		final public double correlation;

		/**
		 * Instantiates a new FIRE result.
		 *
		 * @param fireNumber
		 *            the fire number
		 * @param correlation
		 *            the correlation
		 */
		private FIREResult(double fireNumber, double correlation)
		{
			this.fireNumber = fireNumber;
			this.correlation = correlation;
		}
	}

	// Properties controlling the algorithm

	/**
	 * Depending on the sampling method, the correlation is computed using interpolated values from intervals around the
	 * circle circumference. The number of samples for half the circle is computed as: Pi * radius * sampling factor.
	 */
	private double perimeterSamplingFactor = 1;

	/**
	 * Gets the perimeter sampling factor.
	 *
	 * @return the perimeter sampling factor
	 */
	public double getPerimeterSamplingFactor()
	{
		return perimeterSamplingFactor;
	}

	/**
	 * Sets the perimeter sampling factor.
	 * <p>
	 * Depending on the sampling method, the correlation is computed using interpolated values from intervals around the
	 * circle circumference. The number of samples for half the circle is computed as: Pi * radius * sampling factor.
	 *
	 * @param perimeterSamplingFactor
	 *            the new perimeter sampling factor
	 */
	public void setPerimeterSamplingFactor(double perimeterSamplingFactor)
	{
		if (perimeterSamplingFactor <= 0)
			perimeterSamplingFactor = 1;
		this.perimeterSamplingFactor = perimeterSamplingFactor;
	}

	/**
	 * Control the method for generating the Fourier circle.
	 * <p>
	 * The correlation is computed using intervals around the circle circumference of the Fourier transform.
	 * <p>
	 * Note in the case of using interpolated pixels on the perimeter the Fourier image is 2-fold radially symmetric and
	 * so the calculation can use only half the circle for speed.
	 */
	private SamplingMethod samplingMethod = SamplingMethod.RADIAL_SUM;

	/**
	 * Gets the sampling method.
	 *
	 * @return the sampling method
	 */
	public SamplingMethod getSamplingMethod()
	{
		return samplingMethod;
	}

	/**
	 * Control the method for generating the Fourier circle.
	 * <p>
	 * The correlation is computed using intervals around the circle circumference of the Fourier transform. The radial
	 * sum does not use interpolation but instead assigns all pixels at radius r to the interval {@code n<=r<n+1} for all n
	 * from 0 to the max radius.
	 * <p>
	 * Note in the case of using interpolated pixels on the perimeter the Fourier image is 2-fold radially symmetric and
	 * so the calculation can use only half the circle for speed.
	 *
	 * @param
	 * 		   samplingMethod
	 *            the new sampling method
	 */
	public void setSamplingMethod(SamplingMethod samplingMethod)
	{
		if (samplingMethod == null)
			samplingMethod = SamplingMethod.RADIAL_SUM;
		this.samplingMethod = samplingMethod;
	}

	/** The fourier method. */
	private FourierMethod fourierMethod = FourierMethod.JTRANSFORMS;

	/**
	 * Gets the Fourier method.
	 *
	 * @return the Fourier method
	 */
	public FourierMethod getFourierMethod()
	{
		return fourierMethod;
	}

	/**
	 * Sets the Fourier method.
	 *
	 * @param fourierMethod
	 *            the new Fourier method
	 */
	public void setFourierMethod(FourierMethod fourierMethod)
	{
		this.fourierMethod = fourierMethod;
	}

	/** Used to track the progress within {@link #calculateFrcCurve(ImageProcessor, ImageProcessor, double)}. */
	private TrackProgress progress = null;

	/**
	 * Sets the track progress.
	 * <p>
	 * Used to track the progress within {@link #calculateFrcCurve(ImageProcessor, ImageProcessor, double)}.
	 *
	 * @param progress
	 *            the new track progress
	 */
	public void setTrackProgress(TrackProgress progress)
	{
		this.progress = progress;
	}

	/**
	 * Gets the track progress.
	 *
	 * @return the track progress
	 */
	private TrackProgress getTrackProgress()
	{
		return progress = NullTrackProgress.createIfNull(progress);
	}

	/** Constant containing the value of a third. */
	private static final double THIRD = 1.0 / 3.0;

	/**
	 * Constant containing the value of a third so that value of 2 * {@link FRC#THIRD} + {@link FRC#LAST_THIRD} == 1.
	 */
	private static final double LAST_THIRD = 1.0 - 2 * THIRD;

	/**
	 * Calculate the Fourier Ring Correlation curve for two images.
	 *
	 * @param ip1
	 *            The first image
	 * @param ip2
	 *            The second image
	 * @param nmPerPixel
	 *            the nm per pixel for the super-resolution images
	 * @return the FRC curve
	 */
	public FRCCurve calculateFrcCurve(ImageProcessor ip1, ImageProcessor ip2, double nmPerPixel)
	{
		// Allow a progress tracker to be input
		final TrackProgress progess = getTrackProgress();
		progess.incrementProgress(0);
		progess.status("Calculating complex FFT images...");

		// Pad images to the same size if different
		final int maxWidth = FastMath.max(ip1.getWidth(), ip2.getWidth());
		final int maxHeight = FastMath.max(ip1.getHeight(), ip2.getHeight());
		if (FastMath.max(maxWidth, maxHeight) > MAX_SIZE)
		{
			progess.status("Error calculating FRC curve...");
			progess.incrementProgress(1);
			return null;
		}
		final int fieldOfView = FastMath.max(maxWidth, maxHeight);
		ip1 = pad(ip1, maxWidth, maxHeight);
		ip2 = pad(ip2, maxWidth, maxHeight);

		// The mean of each image after applying the taper
		double mean1, mean2;

		// Real and imaginary components
		float[] re1 = null, im1, re2, im2;

		// Do the first image
		ip1 = getSquareTaperedImage(ip1);
		mean1 = taperedImageMean;
		final int size = ip1.getWidth();

		if (JTRANSFORMS && fourierMethod == FourierMethod.JTRANSFORMS)
		{
			// Speed up by reusing the FFT object which performs pre-computation
			final float[] data = new float[size * size * 2];
			final FloatFFT_2D fft = new FloatFFT_2D(size, size);

			float[] pixels = (float[]) ip1.getPixels();
			System.arraycopy(pixels, 0, data, 0, pixels.length);
			fft.realForwardFull(data);

			// Get the data
			re1 = pixels;
			im1 = new float[pixels.length];
			for (int i = 0, j = 0; i < data.length; j++)
			{
				re1[j] = data[i++];
				im1[j] = data[i++];
			}
			FHT2.swapQuadrants(new FloatProcessor(size, size, re1));
			FHT2.swapQuadrants(new FloatProcessor(size, size, im1));
			progess.incrementProgress(THIRD);

			ip2 = getSquareTaperedImage(ip2);
			mean2 = taperedImageMean;

			pixels = (float[]) ip2.getPixels();
			System.arraycopy(pixels, 0, data, 0, pixels.length);
			for (int i = pixels.length; i < data.length; i++)
				data[i] = 0;
			fft.realForwardFull(data);

			// Get the data
			re2 = pixels;
			im2 = new float[pixels.length];
			for (int i = 0, j = 0; i < data.length; j++)
			{
				re2[j] = data[i++];
				im2[j] = data[i++];
			}
			FHT2.swapQuadrants(new FloatProcessor(size, size, re2));
			FHT2.swapQuadrants(new FloatProcessor(size, size, im2));
			progess.incrementProgress(THIRD);
		}
		else
		{
			// Simple implementation. This is left for testing.
			//FloatProcessor[] fft = getComplexFFT(ip1);
			//mean1 = taperedImageMean;
			//re1 = (float[]) fft[0].getPixels();
			//im1 = (float[]) fft[1].getPixels();
			//progess.incrementProgress(THIRD);
			//
			//fft = getComplexFFT(ip2);
			//mean2 = taperedImageMean;
			//re2 = (float[]) fft[0].getPixels();
			//im2 = (float[]) fft[1].getPixels();
			//progess.incrementProgress(THIRD);

			// Speed up by reusing the FHT object which performs pre-computation

			final float[] f1 = (float[]) ip1.getPixels();
			final FHT2 fht1 = new FHT2(f1, ip1.getWidth(), false);
			fht1.transform();
			FloatProcessor[] fft = fht1.getComplexTransformProcessors();
			re1 = (float[]) fft[0].getPixels();
			im1 = (float[]) fft[1].getPixels();
			progess.incrementProgress(THIRD);

			ip2 = getSquareTaperedImage(ip2);
			mean2 = taperedImageMean;

			final float[] f2 = (float[]) ip2.getPixels();
			final FHT2 fht2 = new FHT2(f2, ip2.getWidth(), false);
			fht2.copyTables(fht1);
			fft = fht2.getComplexTransformProcessors();
			re2 = (float[]) fft[0].getPixels();
			im2 = (float[]) fft[1].getPixels();
			progess.incrementProgress(THIRD);
		}

		progess.status("Preparing FRC curve calculation...");

		final int centre = size / 2;

		// In-line for speed
		final float[] conjMult = new float[re1.length];
		final float[] absFFT1 = new float[re1.length];
		final float[] absFFT2 = new float[re1.length];

		// Normalise the FFT to the field of view, i.e. normalise by 1/sqrt(N) for each dimension
		final double norm = 1.0 / fieldOfView;
		for (int i = 0; i < re1.length; i++)
		{
			re1[i] *= norm;
			im1[i] *= norm;
			re2[i] *= norm;
			im2[i] *= norm;
		}

		final boolean basic = false;
		if (basic)
			compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2);
		else
			computeMirroredFast(size, conjMult, absFFT1, absFFT2, re1, im1, re2, im2);

		progess.status("Calculating FRC curve...");

		final int max = centre - 1;
		final FRCCurveResult[] results = new FRCCurveResult[max];

		if (samplingMethod == SamplingMethod.INTERPOLATED_CIRCLE)
		{
			// Set the results for the centre pixel
			final int cx = size * centre + centre;
			results[0] = new FRCCurveResult(0, 1, conjMult[cx], absFFT1[cx], absFFT2[cx]);

			final float[][] images = new float[][] { conjMult, absFFT1, absFFT2 };
			for (int radius = 1; radius < max; radius++)
			{
				// Inline the calculation for speed
				double sum0 = 0;
				double sum1 = 0;
				double sum2 = 0;

				// Note: The image has 2-fold radial symmetry. So we only need to sample
				// angles from 0-pi. To sample the perimeter at pixel intervals we need
				// pi*r samples. So the angle step is max_angle / samples == pi / (pi*r) == 1 / r.
				// The number of samples is increased using the sampling factor.

				final double angleStep = 1 / (perimeterSamplingFactor * radius);

				double angle = 0;
				int numSum = 0;

				while (angle < Math.PI)
				{
					final double cosA = FastMath.cos(angle);
					final double x = centre + radius * cosA;
					//double sinA = FastMath.sin(angle);
					final double sinA = getSine(angle, cosA);
					final double y = centre + radius * sinA;
					final double[] values = getInterpolatedValues(x, y, images, size);
					sum0 += values[0];
					sum1 += values[1];
					sum2 += values[2];

					numSum++;
					angle += angleStep;
				}

				results[radius] = new FRCCurveResult(radius, numSum, sum0, sum1, sum2);
			}
		}
		else
		{
			// Compute the radial sum as per the DIP image Matlab toolbox
			final double[][] sum = RadialStatistics.radialSumAndCount(size, conjMult, absFFT1, absFFT2);
			for (int radius = 0; radius < max; radius++)
				results[radius] = new FRCCurveResult(radius, (int) sum[3][radius], sum[0][radius], sum[1][radius],
						sum[2][radius]);
		}

		progess.incrementProgress(LAST_THIRD);
		progess.status("Finished calculating FRC curve...");

		return new FRCCurve(nmPerPixel, fieldOfView, mean1, mean2, results);
	}

	/**
	 * Compute the conjugate multiple of two FFT images.
	 *
	 * @param conjMult
	 *            the conjugate multiplication of FFT 1 and FFT 2
	 * @param absFFT1
	 *            the absolute magnitude of FFT 1
	 * @param absFFT2
	 *            the absolute magnitude of FFT 2
	 * @param re1
	 *            the real part of FFT 1
	 * @param im1
	 *            the imaginary part of FFT 1
	 * @param re2
	 *            the real part of FFT 2
	 * @param im2
	 *            the imaginary part of FFT 2
	 */
	// Package level to allow JUnit test
	static void compute(float[] conjMult, float[] absFFT1, float[] absFFT2, float[] re1, float[] im1, float[] re2,
			float[] im2)
	{
		for (int i = re1.length; i-- > 0;)
			compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i);
	}

	/**
	 * Compute the conjugate multiple of two FFT images.
	 *
	 * @param size
	 *            the size of the FFT image
	 * @param conjMult
	 *            the conjugate multiplication of FFT 1 and FFT 2
	 * @param absFFT1
	 *            the absolute magnitude of FFT 1
	 * @param absFFT2
	 *            the absolute magnitude of FFT 2
	 * @param re1
	 *            the real part of FFT 1
	 * @param im1
	 *            the imaginary part of FFT 1
	 * @param re2
	 *            the real part of FFT 2
	 * @param im2
	 *            the imaginary part of FFT 2
	 */
	// Package level to allow JUnit test
	static void computeMirrored(int size, float[] conjMult, float[] absFFT1, float[] absFFT2, float[] re1, float[] im1,
			float[] re2, float[] im2)
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
			compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i++);

		// Compute remaining rows up to the centre. These are mirrored
		int j = conjMult.length - 1;
		for (int y = 1; y < centre; y++)
		{
			// The first entry in each row is not mirrored so compute and increment i
			compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i++);
			for (int x = 1; x < size; x++, i++, j--)
			{
				compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i);
				// Mirror
				conjMult[j] = conjMult[i];
				absFFT1[j] = absFFT1[i];
				absFFT2[j] = absFFT2[i];
			}
			// The last entry in each reverse row is not mirrored so compute and decrement j
			compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, j--);
		}

		// Do the centre row. This is mirrored with itself
		compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i++);
		for (int x = 1; x <= centre; x++, i++, j--)
		{
			compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i);
			// Mirror
			conjMult[j] = conjMult[i];
			absFFT1[j] = absFFT1[i];
			absFFT2[j] = absFFT2[i];
		}
	}

	/**
	 * Compute the conjugate multiple of two FFT images.
	 *
	 * @param size
	 *            the size
	 * @param conjMult
	 *            the conjugate multiplication of FFT 1 and FFT 2
	 * @param absFFT1
	 *            the absolute magnitude of FFT 1
	 * @param absFFT2
	 *            the absolute magnitude of FFT 2
	 * @param re1
	 *            the real part of FFT 1
	 * @param im1
	 *            the imaginary part of FFT 1
	 * @param re2
	 *            the real part of FFT 2
	 * @param im2
	 *            the imaginary part of FFT 2
	 */
	// Package level to allow JUnit test
	static void computeMirroredFast(int size, float[] conjMult, float[] absFFT1, float[] absFFT2, float[] re1,
			float[] im1, float[] re2, float[] im2)
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
		int j = conjMult.length - 1;
		for (int y = 1; y < centre; y++)
		{
			// The first entry in each row is not mirrored so just increment i
			i++;
			for (int x = 1; x < size; x++, i++, j--)
			{
				compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i);
				// Mirror
				conjMult[j] = conjMult[i];
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
			compute(conjMult, absFFT1, absFFT2, re1, im1, re2, im2, i);
			// Mirror
			conjMult[j] = conjMult[i];
			absFFT1[j] = absFFT1[i];
			absFFT2[j] = absFFT2[i];
		}
	}

	/**
	 * Compute the conjugate multiple of two FFT images at the given index.
	 *
	 * @param conjMult
	 *            the conjugate multiplication of FFT 1 and FFT 2
	 * @param absFFT1
	 *            the absolute magnitude of FFT 1
	 * @param absFFT2
	 *            the absolute magnitude of FFT 2
	 * @param re1
	 *            the real part of FFT 1
	 * @param im1
	 *            the imaginary part of FFT 1
	 * @param re2
	 *            the real part of FFT 2
	 * @param im2
	 *            the imaginary part of FFT 2
	 * @param i
	 *            the index
	 */
	private static void compute(float[] conjMult, float[] absFFT1, float[] absFFT2, float[] re1, float[] im1,
			float[] re2, float[] im2, int i)
	{
		final float re1i = re1[i];
		final float im1i = im1[i];
		final float re2i = re2[i];
		final float im2i = im2[i];
		conjMult[i] = re1i * re2i + im1i * im2i;
		absFFT1[i] = re1i * re1i + im1i * im1i;
		absFFT2[i] = re2i * re2i + im2i * im2i;
	}

	/**
	 * Check symmetry.
	 *
	 * @param data
	 *            the data
	 * @param size
	 *            the size
	 * @return true, if successful
	 */
	static boolean checkSymmetry(float[] data, int size)
	{
		// Symmetry is around the centre
		final int centre = size / 2;

		final float maxRelativeError = 1e-10f, maxAbsoluteError = 1e-16f;
		int error = 0;

		for (int y = centre, y2 = centre; y >= 0 && y2 < size; y--, y2++)
			for (int x = centre, x2 = centre, i = size * y + x, j = size * y2 + x2; x >= 0 &&
					x2 < size; x--, x2++, i--, j++)
				if (data[i] != data[j] || !FloatEquality.almostEqualRelativeOrAbsolute(data[i], data[j],
						maxRelativeError, maxAbsoluteError))
					//System.out.printf("[%d] %f != [%d] %f\n", i, data[i], j, data[j]);
					if (--error < 0)
						//uk.ac.sussex.gdsc.core.ij.Utils.display("check", new FloatProcessor(size, size, data));
						return false;
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
		final double sine = Math.sqrt(1 - (cosA * cosA));
		return ((angle > Math.PI) ? -sine : sine); // Place in correct domain
	}

	/**
	 * Pad.
	 *
	 * @param ip
	 *            the image
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @return the image processor
	 */
	private static ImageProcessor pad(ImageProcessor ip, int width, int height)
	{
		if (ip.getWidth() != width || ip.getHeight() != height)
		{
			final ImageProcessor ip2 = ip.createProcessor(width, height);
			ip2.insert(ip, 0, 0);
			return ip2;
		}
		return ip;
	}

	/**
	 * Convert an image into a Fourier image with real and imaginary parts.
	 *
	 * @param ip
	 *            The image
	 * @return the real and imaginary parts
	 */
	public FloatProcessor[] getComplexFFT(ImageProcessor ip)
	{
		final FloatProcessor taperedDataImage = getSquareTaperedImage(ip);

		final FHT2 fht = new FHT2(taperedDataImage);
		fht.transform();

		return fht.getComplexTransformProcessors();
	}

	/** Max size for Fourier images. Must be a power of 2 that is smaller that Math.sqrt(Integer.MAX_VALUE) */
	private static final int MAX_SIZE = 32768;

	/**
	 * Returns the closest power-of-two number greater than or equal to x.
	 * <p>
	 * Copied from the JTransforms library class edu.emory.mathcs.utils.ConcurrencyUtils.
	 *
	 * @param x
	 *            the x
	 * @return the closest power-of-two number greater than or equal to x
	 */
	public static int nextPow2(int x)
	{
		if (x < 1)
			throw new IllegalArgumentException("x must be greater or equal 1");
		if ((x & (x - 1)) == 0)
			return x; // x is already a power-of-two number
		x |= (x >>> 1);
		x |= (x >>> 2);
		x |= (x >>> 4);
		x |= (x >>> 8);
		x |= (x >>> 16);
		x |= (x >>> 32);
		return x + 1;
	}

	/** The tapered image mean. */
	private double taperedImageMean;

	/**
	 * Applies a Tukey window function to the image and then pads it to the next square size power of two.
	 *
	 * @param dataImage
	 *            the data image
	 * @return The square tapered image
	 */
	public FloatProcessor getSquareTaperedImage(ImageProcessor dataImage)
	{
		taperedImageMean = 0;

		final int size = FastMath.max(dataImage.getWidth(), dataImage.getHeight());
		if (size > MAX_SIZE)
			return null; // Too large so error

		// Use a Tukey window function
		final float[] taperX = getWindowFunctionX(dataImage.getWidth());
		final float[] taperY = getWindowFunctionY(dataImage.getHeight());

		// Pad to a power of 2
		final int newSize = nextPow2(size);

		dataImage = dataImage.toFloat(0, null);
		final float[] data = (float[]) dataImage.getPixels();
		final float[] pixels = new float[newSize * newSize];
		// Note that the limits at 0 and size-1 the taper is zero so this can be ignored
		final int maxy_1 = dataImage.getHeight() - 1;
		final int maxx_1 = dataImage.getWidth() - 1;
		final int oldWidth = dataImage.getWidth();
		for (int y = 1; y < maxy_1; y++)
		{
			final float yTmp = taperY[y];
			for (int x = 1, i = y * oldWidth + 1, ii = y * newSize + 1; x < maxx_1; x++, i++, ii++)
			{
				final float v = data[i] * taperX[x] * yTmp;
				taperedImageMean += v;
				pixels[ii] = v;
			}
		}
		// Take the mean over the non-padded image size
		taperedImageMean /= (dataImage.getPixelCount());

		return new FloatProcessor(newSize, newSize, pixels, null);
	}

	// Cache the Tukey window function.
	// Using methods to check the length should make it thread safe since we create an instance reference
	// to an array of the correct length.
	private static float[] taperX = new float[0], taperY = new float[0];

	/**
	 * Gets the window function X.
	 *
	 * @param size
	 *            the size
	 * @return the window function X
	 */
	private static float[] getWindowFunctionX(int size)
	{
		final float[] taper = getWindowFunction(taperX, size);
		taperX = taper;
		return taper;
	}

	/**
	 * Gets the window function Y.
	 *
	 * @param size
	 *            the size
	 * @return the window function Y
	 */
	private static float[] getWindowFunctionY(int size)
	{
		final float[] taper = getWindowFunction(taperY, size);
		taperY = taper;
		return taper;
	}

	/**
	 * Gets the window function.
	 *
	 * @param taper
	 *            the taper
	 * @param size
	 *            the size
	 * @return the window function
	 */
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

	/**
	 * Check.
	 *
	 * @param taper
	 *            the taper
	 * @param size
	 *            the size
	 * @return the float[]
	 */
	private static float[] check(float[] taper, int size)
	{
		return (taper.length == size) ? taper : null;
	}

	/**
	 * Gets the tukey window function.
	 *
	 * @param size
	 *            the size
	 * @return the tukey window function
	 */
	private static float[] getTukeyWindowFunction(int size)
	{
		final float[] taper = new float[size];

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

		// New optimised code. This matches uk.ac.sussex.gdsc.core.utils.ImageWindow.tukey(size, 0.25)
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
				taper[i++] = taper[j--] = 1f;
		}

		//// Use the GDSC ImageWindow class:
		//// FRC has Image width / Width of edge region = 8 so use alpha 0.25.
		//double[] w = uk.ac.sussex.gdsc.core.utils.ImageWindow.tukey(size, 0.25);
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
	 *            the x
	 * @param y
	 *            the y
	 * @param images
	 *            the images
	 * @param maxx
	 *            the maxx
	 * @return the interpolated values
	 */
	private static double[] getInterpolatedValues(final double x, final double y, float[][] images, final int maxx)
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
		final double[] values = new double[noImages];
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
	 * The correlation values are smoothed using a LOESS interpolation with the given parameters. If smoothing fails the
	 * original curve values are returned. If successful then the input curve is optionally copied and updated with the
	 * smoothed correlation.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param bandwidth
	 *            the bandwidth
	 * @param robustness
	 *            the robustness
	 * @param inPlace
	 *            Set to true to modify the correlation in place (do not create a copy)
	 * @return A new FRC curve
	 */
	private static FRCCurve getSmoothedCurve(FRCCurve frcCurve, double bandwidth, int robustness, boolean inPlace)
	{
		final double[] xVals = new double[frcCurve.getSize()];
		final double[] yVals = new double[frcCurve.getSize()];

		for (int i = 0; i < frcCurve.getSize(); i++)
		{
			xVals[i] = frcCurve.get(i).getRadius();
			yVals[i] = frcCurve.get(i).getCorrelation();
		}

		double[] ySmoothed = new double[frcCurve.getSize()];

		try
		{
			final LoessInterpolator loess = new LoessInterpolator(bandwidth, robustness);
			ySmoothed = loess.smooth(xVals, yVals);

			if (!inPlace)
				frcCurve = frcCurve.copy();

			for (int i = 0; i < frcCurve.getSize(); i++)
				frcCurve.get(i).setCorrelation(ySmoothed[i]);
		}
		catch (final Exception e)
		{
			e.printStackTrace();
		}

		return frcCurve;
	}

	/**
	 * Perform LOESS smoothing on the FRC curve
	 * <p>
	 * | * The correlation values are smoothed using a LOESS interpolation with bandwidth of 0.0707 and robustness of 0.
	 * If smoothing fails the original curve values are returned. If successful then the input curve is optionally
	 * copied and updated with the smoothed correlation.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param inPlace
	 *            Set to true to modify the correlation in place (do not create a copy)
	 * @return A new FRC curve
	 */
	public static FRCCurve getSmoothedCurve(FRCCurve frcCurve, boolean inPlace)
	{
		final double bandwidth = 0.0707;
		final int robustness = 0;
		return getSmoothedCurve(frcCurve, bandwidth, robustness, inPlace);
	}

	/** The constant 2 * Math.PI */
	private final static double TWO_PI = 2.0 * Math.PI;

	/**
	 * Calculate the curve representing the minimum correlation required to distinguish two images for each
	 * resolution
	 * in the input FRC curve.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param thresholdMethod
	 *            the threshold method
	 * @return The threshold curve representing the threshold for each input spatial frequency
	 */
	public static double[] calculateThresholdCurve(FRCCurve frcCurve, ThresholdMethod thresholdMethod)
	{
		final double[] threshold = new double[frcCurve.getSize()];

		// ADH:
		// Half-Bit and 3 Sigma are explained in Supp section 5.4. However this is based on
		// Heel, et al (2005) which gives a better explanation so I have updated using their
		// equations.
		// Equation S.84 has an error compared to equation (13) in Heel. In fact equation (17)
		// from Heel is what should be implemented for Half-bit.

		// Note: The original code used frcCurve.get(i)[2] which holds the number of samples that
		// were taken from the circle. This makes the curve dependent on the number of samples
		// taken (e.g. half-circle/full-circle with different sampling factors).
		// To make the curve sampling independent I assume 2*pi*r samples were taken.

		// See: Heel, M. v. & Schatz, M. Fourier shell correlation threshold criteria. J. Struct. Bio. 151, 250–262 (2005).

		// Fourier Shell Radius:
		// This is set to 1 for r==0 in Heel
		double nr = 1;

		int sigma = 0; // For the N-sigma methods

		switch (thresholdMethod)
		{
			// This is first as it is the most commonly used
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
					// Heel, Equation (2):
					// We actually want to know the number of pixels contained in the Fourier shell of radius r.
					// We can compute this assuming we sampled the full circle 2*pi*r.
					threshold[i] = sigma / Math.sqrt(nr / 2.0);
				break;

			default:
				Arrays.fill(threshold, 1.0 / 7.0);
		}

		return threshold;
	}

	/**
	 * Compute the threshold curve for the given number of bits.
	 *
	 * @param threshold
	 *            the threshold
	 * @param bits
	 *            the bits
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
		final double snr = (FastMath.pow(2, bits) - 1) / 2;
		final double snr_p_1 = snr + 1;
		final double twoRootSnr = 2 * Math.sqrt(snr);
		final double twoRootSnr_p_1 = twoRootSnr + 1;

		// Sense check:
		// 1/2-bit is equation (17) from Heel:
		// snr = 0.2071, twoRootSnr = 0.9102
		// 1-bit is equation (14) from Heel:
		// snr = 0.5, twoRootSnr = 1.4142

		for (int i = 1; i < threshold.length; i++)
		{
			// nr = number of samples in Fourier circle = 2*pi*r
			final double sqrtNr = Math.sqrt(TWO_PI * i);
			threshold[i] = ((snr + twoRootSnr_p_1 / sqrtNr) / (snr_p_1 + twoRootSnr / sqrtNr));
		}
	}

	/**
	 * Computes the crossing points of the FRC curve and the threshold curve. The intersections can be used to
	 * determine the image resolution using {@link #getCorrectIntersection(double[][], ThresholdMethod)}
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param thresholdCurve
	 *            the threshold curve
	 * @param max
	 *            The maximum number of intersections to compute
	 * @return The crossing points
	 */
	public static double[][] getIntersections(FRCCurve frcCurve, double[] thresholdCurve, int max)
	{
		if (frcCurve.getSize() != thresholdCurve.length)
		{
			System.err.println("Error: Unable to calculate FRC curve intersections due to input length mismatch.");
			return null;
		}

		final double[][] intersections = new double[Math.min(max, frcCurve.getSize() - 1)][];
		int count = 0;

		for (int i = 1; i < frcCurve.getSize(); i++)
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

			final double y1 = frcCurve.get(i - 1).getCorrelation();
			final double y2 = frcCurve.get(i).getCorrelation();
			final double y3 = thresholdCurve[i - 1];
			final double y4 = thresholdCurve[i];

			// Check if they cross
			if (!((y3 >= y1 && y4 < y2) || (y1 >= y3 && y2 < y4)))
				continue;

			final double x1 = frcCurve.get(i - 1).getRadius();
			final double x2 = frcCurve.get(i).getRadius();
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
				final double px = ((x1 * y2 - y1 * x2) * x3_x4 - x1_x2 * (x3 * y4 - y3 * x4)) /
						(x1_x2 * y3_y4 - y1_y2 * x3_x4);
				//double px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) /
				//		((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
				//double py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) /
				//		((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));

				// Check if the intersection is within the two points
				// Q. Is this necessary given the intersection check above?
				if (px >= x1 && px < x2)
				{
					final double py = Maths.interpolateY(x1, y3, x2, y4, px);
					intersections[count++] = new double[] { px, py };
				}
			}
			if (count >= max)
				break;
		}

		return Arrays.copyOf(intersections, count);
	}

	/**
	 * Get the correction intersection representing the image resolution. The intersection chosen depends on the
	 * method used to calculate the threshold curve using {@link #calculateThresholdCurve(FRCCurve, ThresholdMethod)}
	 * <p>
	 * The intersection corresponds the lowest spatial frequency at which there is no significant correlation
	 * between the images.
	 *
	 * @param intersections
	 *            the intersections
	 * @param thresholdMethod
	 *            the threshold method
	 * @return The intersection (or null if no crossings)
	 */
	public static double[] getCorrectIntersection(double[][] intersections, ThresholdMethod thresholdMethod)
	{
		if (intersections == null || intersections.length == 0)
			return null;

		int pos = 0;
		switch (thresholdMethod)
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
					// Choose the intersection. Just discard an intersection with a correlation above 0.9
					if (intersections[0][1] > 0.9)
						pos++;
		}

		return intersections[pos];
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided images.
	 *
	 * @param ip1
	 *            the first image
	 * @param ip2
	 *            the second image
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @param thresholdMethod
	 *            the threshold method
	 * @return The FIRE number (in pixels) and the correlation
	 */
	public double calculateFireNumber(ImageProcessor ip1, ImageProcessor ip2, double nmPerPixel,
			ThresholdMethod thresholdMethod)
	{
		final FRCCurve frcCurve = calculateFrcCurve(ip1, ip2, nmPerPixel);
		if (frcCurve == null)
			return Double.NaN;
		return calculateFireNumber(frcCurve, thresholdMethod);
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided FRC curve
	 * data.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param thresholdMethod
	 *            the threshold method
	 * @return The FIRE number (in pixels)
	 */
	public static double calculateFireNumber(FRCCurve frcCurve, ThresholdMethod thresholdMethod)
	{
		final FIREResult result = calculateFire(frcCurve, thresholdMethod);
		if (result == null)
			return Double.NaN;
		return result.fireNumber;
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided images.
	 *
	 * @param ip1
	 *            the first image
	 * @param ip2
	 *            the second image
	 * @param nmPerPixel
	 *            the nm per pixel
	 * @param thresholdMethod
	 *            the method
	 * @return The FIRE result (null if computation failed)
	 */
	public FIREResult calculateFire(ImageProcessor ip1, ImageProcessor ip2, double nmPerPixel,
			ThresholdMethod thresholdMethod)
	{
		FRCCurve frcCurve = calculateFrcCurve(ip1, ip2, nmPerPixel);
		if (frcCurve == null)
			return null;
		frcCurve = getSmoothedCurve(frcCurve, false);
		return calculateFire(frcCurve, thresholdMethod);
	}

	/**
	 * Utility function that calculates the Fourier Image Resolution (FIRE) number using the provided FRC curve
	 * data.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param thresholdMethod
	 *            the threshold method
	 * @return The FIRE result (null if computation failed)
	 */
	public static FIREResult calculateFire(FRCCurve frcCurve, ThresholdMethod thresholdMethod)
	{
		final double[] thresholdCurve = calculateThresholdCurve(frcCurve, thresholdMethod);
		final double[][] intersections = getIntersections(frcCurve, thresholdCurve, 2);

		if (intersections != null && intersections.length != 0)
		{
			final double[] intersection = getCorrectIntersection(intersections, thresholdMethod);
			final double fireNumber = frcCurve.fieldOfView / intersection[0];
			return new FIREResult(frcCurve.nmPerPixel * fireNumber, intersection[1]);
		}
		// Edge case where the entire curve has a correlation of 1.
		// This happens when the two split images are the same, e.g. comparing a dataset to itself.
		if (perfect(frcCurve))
			return new FIREResult(0, 1);
		return null;
	}

	/**
	 * Perfect.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @return true, if successful
	 */
	private static boolean perfect(FRCCurve frcCurve)
	{
		for (int i = 0; i < frcCurve.getSize(); i++)
			if (frcCurve.get(i).getCorrelation() != 1)
				return false;
		return true;
	}

	/**
	 * Apply spurious correlation correction using the Q-factor. Follows the method described in Niewenhuizen, et al
	 * (2013), on-line methods. This method computes a value that is subtracted from the numerator of the FRC and
	 * added to the denominator of the FRC before computing the correlation. The correlation in the FRC curve will be
	 * updated. The sums of the FRC curve are unchanged.
	 * <p>
	 * The localisation precision is input in units of nm. This is converted using the nmPerPixel stored in the FRC
	 * curve object into units of super-resolution pixels.
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
	 *            the FRC curve
	 * @param qValue
	 *            the q value
	 * @param mean
	 *            the mean of the Gaussian localisation precision (in units of nm)
	 * @param sigma
	 *            the width of the Gaussian localisation precision (in units of nm)
	 */
	public static void applyQCorrection(FRCCurve frcCurve, double qValue, double mean, double sigma)
	{
		if (qValue <= 0)
		{
			// Reset the correlation
			for (int i = 1; i < frcCurve.getSize(); i++)
				frcCurve.get(i).setCorrelation(frcCurve.get(i).getNumerator() / frcCurve.get(i).getDenominator());
			return;
		}

		// q must be in pixel units
		final double[] q = computeQ(frcCurve, false);

		// H(q) is the factor in the correlation averages related to the localization
		// uncertainties that depends on the mean and width of the
		// distribution of localization uncertainties
		final double[] hq = computeHq(q, mean / frcCurve.nmPerPixel, sigma / frcCurve.nmPerPixel);

		// Precompute Q normalisation
		final double qNorm = (1 / frcCurve.mean1 + 1 / frcCurve.mean2);
		qValue /= qNorm;

		// Subtract the average residual correlation from the numerator and add to the denominator
		for (int i = 1; i < frcCurve.getSize(); i++)
		{
			// The numerator and denominator are computed using the radial sum.
			// Convert this to the radial mean by dividing by the number of samples.
			final int n = frcCurve.get(i).getNumberOfSamples();
			final double numerator = frcCurve.get(i).getNumerator() / n;
			final double denominator = frcCurve.get(i).getDenominator() / n;
			// Matlab code provided by Bernd Reiger computes the residual as:
			// Q/Qnorm*sinc(q).^2.*exp(-4*pi^2*sigma_mean^2*q.^2./(1+8*pi^2*sigma_std^2*q.^2))./sqrt(1+8*pi^2*sigma_std^2*q.^2);
			// Qnorm = (1/mean1+1/mean2);
			// Matlab sinc is sin(pi*x) / pi*x
			final double sinc_q = sinc(Math.PI * q[i]);
			final double residual = qValue * sinc_q * sinc_q * hq[i];
			frcCurve.get(i).setCorrelation((numerator - residual) / (denominator + residual));
		}
	}

	/**
	 * Sinc.
	 *
	 * @param x
	 *            the x
	 * @return the double
	 */
	private static double sinc(double x)
	{
		return Math.sin(x) / x;
	}

	/**
	 * Compute q. This is defined as 1/L, 2/L, ..., for all the spatial frequencies in the FRC curve where L is the
	 * size of the field of view. This is converted to nm using the pixel size of the input image.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param nmPerPixel
	 *            the nm per pixel in the images used to compute the FRC curve
	 * @return the q array (in nm^-1)
	 */
	private static double[] computeQ(FRCCurve frcCurve, double nmPerPixel)
	{
		final double L = frcCurve.fieldOfView;

		final double[] q = new double[frcCurve.getSize()];
		final double conversion = 1.0 / (L * nmPerPixel);
		for (int i = 0; i < q.length; i++)
			q[i] = frcCurve.get(i).radius * conversion;
		return q;
	}

	/**
	 * Compute q. This is defined as 1/L, 2/L, ..., for all the spatial frequencies in the FRC curve where L is the
	 * size of the field of view. This is optionally converted to nm using the pixel size units in the FRC curve.
	 *
	 * @param frcCurve
	 *            the FRC curve
	 * @param useUnits
	 *            Set to true to convert the pixel units to nm
	 * @return the q array
	 */
	public static double[] computeQ(FRCCurve frcCurve, boolean useUnits)
	{
		return computeQ(frcCurve, (useUnits) ? frcCurve.nmPerPixel : 1);
	}

	/**
	 * Compute the localization PDF factor H(q) for all q. This is the integral of the distribution function of the
	 * localisation uncertainty. It is assumed to Gaussian with the specified mean and width.
	 *
	 * @param q
	 *            the q array
	 * @param mean
	 *            the mean of the Gaussian (in units of super-resolution pixels)
	 * @param sigma
	 *            the width of the Gaussian (in units of super-resolution pixels)
	 * @return the Hq array
	 */
	public static double[] computeHq(double[] q, double mean, double sigma)
	{
		// H(q) is the factor in the correlation averages related to the localization
		// uncertainties that depends on the mean and width of the
		// distribution of localization uncertainties
		final double[] hq = new double[q.length];
		final double four_pi2 = 4 * Math.PI * Math.PI;
		final double eight_pi2_s2 = 2 * four_pi2 * sigma * sigma;
		hq[0] = 1;
		for (int i = 1; i < q.length; i++)
		{
			final double q2 = q[i] * q[i];
			final double d = 1 + eight_pi2_s2 * q2;
			hq[i] = FastMath.exp((-four_pi2 * mean * mean * q2) / d) / Math.sqrt(d);
		}
		return hq;
	}
}
