package gdsc.smlm.ij.utils;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.util.Tools;

import java.awt.Rectangle;

/**
 * Aligns an image stack to a reference image using XY translation to maximise the correlation. Takes in:
 * <ul>
 * <li>The reference image
 * <li>The image/stack to align.
 * <li>Optional Max/Min values for the X and Y translation
 * <li>Window function to reduce edge artifacts in frequency space
 * </ul>
 * <p>
 * The alignment is calculated using the maximum correlation between the images. The correlation is computed using the
 * frequency domain (note that conjugate multiplication in the frequency domain is equivalent to correlation in the
 * space domain).
 * <p>
 * Output new stack with the best alignment with optional sub-pixel accuracy.
 * <p>
 * By default restricts translation so that at least half of the smaller image width/height is within the larger image
 * (half-max translation). This can be altered by providing a translation bounds. Note that when using normalised
 * correlation all scores are set to zero outside the half-max translation due to potential floating-point summation
 * error during normalisation.
 */
public class AlignImagesFFT
{
	public enum WindowMethod
	{
		NONE("None"), HANNING("Hanning"), COSINE("Cosine"), TUKEY("Tukey");

		private String name;

		private WindowMethod(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	public enum SubPixelMethod
	{
		NONE("None"), CUBIC("Cubic"), GAUSSIAN("Gaussian");

		private String name;

		private SubPixelMethod(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	private double lastXOffset = 0;
	private double lastYOffset = 0;
	private boolean doTranslation = true;

	// This is used for debugging the normalisation
	private FloatProcessor normalisedRefIp;
	// The location where the reference/target was inserted into the normalised FFT image
	private Rectangle refImageBounds = new Rectangle();
	private Rectangle targetImageBounds = new Rectangle();

	/**
	 * Aligns all images in the target stack to the current processor in the reference.
	 * <p>
	 * If no target is provided then all slices are aligned to the current processor in the reference.
	 * 
	 * @param refImp
	 * @param targetImp
	 * @param windowMethod
	 * @param bounds
	 * @param subPixelMethod
	 * @param interpolationMethod
	 *            see {@link ij.process.ImageProcessor#getInterpolationMethods() }
	 * @param normalised
	 * @param showCorrelationImage
	 * @param showNormalisedImage
	 * @param clipOutput
	 * @return
	 */
	public ImagePlus align(ImagePlus refImp, ImagePlus targetImp, WindowMethod windowMethod, Rectangle bounds,
			SubPixelMethod subPixelMethod, int interpolationMethod, boolean normalised, boolean showCorrelationImage,
			boolean showNormalisedImage, boolean clipOutput)
	{
		ImageProcessor refIp = refImp.getProcessor();
		if (targetImp == null)
			targetImp = refImp;

		// Check same size
		if (!isValid(refIp, targetImp))
			return null;

		// Fourier transforms use the largest power-two dimension that covers both images
		int maxN = Math.max(refIp.getWidth(), refIp.getHeight());
		int maxM = Math.max(targetImp.getWidth(), targetImp.getHeight());
		maxN = Math.max(maxN, maxM);

		this.normalisedRefIp = padAndZero(refIp, maxN, windowMethod, refImageBounds);
		if (showNormalisedImage)
			new ImagePlus(refImp.getTitle() + " Normalised Ref", normalisedRefIp).show();
		maxN = normalisedRefIp.getWidth(); // Update with the power-two result

		// Set up the output stack
		ImageStack outStack = new ImageStack(targetImp.getWidth(), targetImp.getHeight());
		ImageStack correlationStack = null;
		ImageStack normalisedStack = null;
		FloatProcessor fpCorrelation = null;
		FloatProcessor fpNormalised = null;
		if (showCorrelationImage)
		{
			correlationStack = new ImageStack(maxN, maxN);
			fpCorrelation = new FloatProcessor(maxN, maxN);
		}
		if (showNormalisedImage)
		{
			normalisedStack = new ImageStack(maxN, maxN);
			fpNormalised = new FloatProcessor(maxN, maxN);
		}

		// Subtract mean to normalise the numerator of the cross-correlation.
		// ---
		// The effectively normalises the numerator of the correlation but does not address the denominator.
		// The denominator should be calculated using rolling sums for each offset position.
		// See: http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html
		// Following the computation of the correlation each offset (u,v) position should then be divided
		// by the energy of the reference image under the target image. This equates to:
		//   Sum(x,y) [ f(x,y) - f_(u,v) ]^2
		// where f_(u,v) is the mean of the region under the target feature
		// ---

		// Calculate rolling sum of squares
		double[] s = null;
		double[] ss = null;
		if (normalised)
		{
			s = new double[normalisedRefIp.getPixelCount()];
			ss = new double[s.length];
			calculateRollingSums(normalisedRefIp, s, ss);
		}

		FHT refFHT = fft(normalisedRefIp, maxN);

		if (bounds == null)
		{
			bounds = createHalfMaxBounds(refImp.getWidth(), refImp.getHeight(), targetImp.getWidth(),
					targetImp.getHeight());
		}

		// Process each image in the target stack
		ImageStack stack = targetImp.getStack();
		for (int slice = 1; slice <= stack.getSize(); slice++)
		{
			ImageProcessor targetIp = stack.getProcessor(slice);
			outStack.addSlice(
					null,
					alignImages(refFHT, s, ss, targetIp, slice, windowMethod, bounds, fpCorrelation, fpNormalised,
							subPixelMethod, interpolationMethod, clipOutput));
			if (showCorrelationImage)
				correlationStack.addSlice(null, fpCorrelation.duplicate());
			if (showNormalisedImage)
				normalisedStack.addSlice(null, fpNormalised.duplicate());
			if (Utils.isInterrupted())
				return null;
		}

		if (showCorrelationImage)
			new ImagePlus(targetImp.getTitle() + " Correlation", correlationStack).show();
		if (showNormalisedImage)
			new ImagePlus(targetImp.getTitle() + " Normalised Target", normalisedStack).show();

		return new ImagePlus(targetImp.getTitle() + " Aligned", outStack);
	}

	private ImageProcessor refIp = null;
	private double[] s = null;
	private double[] ss = null;
	private FHT refFHT = null;

	/**
	 * Initialises the reference image for batch alignment. All target images should be equal or smaller than the
	 * reference.
	 * 
	 * @param refImp
	 * @param windowMethod
	 * @param normalised
	 *            True if the correlation should be normalised (score of -1 to 1)
	 */
	public void init(ImagePlus refImp, WindowMethod windowMethod, boolean normalised)
	{
		refIp = null;
		s = null;
		ss = null;
		refFHT = null;

		if (refImp == null)
			return;

		init(refImp.getProcessor(), windowMethod, normalised);
	}

	/**
	 * Initialises the reference image for batch alignment. All target images should be equal or smaller than the
	 * reference.
	 * 
	 * @param refIp
	 * @param windowMethod
	 * @param normalised
	 *            True if the correlation should be normalised (score of -1 to 1)
	 */
	public void init(ImageProcessor refIp, WindowMethod windowMethod, boolean normalised)
	{
		s = null;
		ss = null;
		refFHT = null;

		if (refIp == null || noValue(refIp))
			return;

		// Fourier transforms use the largest power-two dimension that covers both images
		int maxN = Math.max(refIp.getWidth(), refIp.getHeight());

		this.normalisedRefIp = padAndZero(refIp, maxN, windowMethod, refImageBounds);

		// Subtract mean to normalise the numerator of the cross-correlation.
		// ---
		// The effectively normalises the numerator of the correlation but does not address the denominator.
		// The denominator should be calculated using rolling sums for each offset position.
		// See: http://www.idiom.com/~zilla/Papers/nvisionInterface/nip.html
		// Following the computation of the correlation each offset (u,v) position should then be divided
		// by the energy of the reference image under the target image. This equates to:
		//   Sum(x,y) [ f(x,y) - f_(u,v) ]^2
		// where f_(u,v) is the mean of the region under the target feature
		// ---

		// Calculate rolling sum of squares
		s = null;
		ss = null;
		if (normalised)
		{
			s = new double[normalisedRefIp.getPixelCount()];
			ss = new double[s.length];
			calculateRollingSums(normalisedRefIp, s, ss);
		}

		refFHT = fft(normalisedRefIp, maxN);
	}

	/**
	 * Aligns all images in the target stack to the pre-initialised reference.
	 * 
	 * @param targetImp
	 * @param windowMethod
	 * @param bounds
	 * @param subPixelMethod
	 * @param interpolationMethod
	 *            see {@link ij.process.ImageProcessor#getInterpolationMethods() }
	 * @param clipOutput
	 *            Set to true to ensure the output image has the same max as the input. Applies to bicubic
	 *            interpolation
	 * @return
	 */
	public ImagePlus align(ImagePlus targetImp, WindowMethod windowMethod, Rectangle bounds,
			SubPixelMethod subPixelMethod, int interpolationMethod, boolean clipOutput)
	{
		if (refFHT == null || targetImp == null)
			return null;

		int maxN = refFHT.getWidth();

		// Check correct size
		if (targetImp.getWidth() > maxN || targetImp.getHeight() > maxN)
			return null;

		// Set up the output stack
		ImageStack outStack = new ImageStack(targetImp.getWidth(), targetImp.getHeight());

		if (bounds == null)
		{
			bounds = createHalfMaxBounds(refIp.getWidth(), refIp.getHeight(), targetImp.getWidth(),
					targetImp.getHeight());
		}

		// Process each image in the target stack
		ImageStack stack = targetImp.getStack();
		for (int slice = 1; slice <= stack.getSize(); slice++)
		{
			ImageProcessor targetIp = stack.getProcessor(slice);
			outStack.addSlice(
					null,
					alignImages(refFHT, s, ss, targetIp, slice, windowMethod, bounds, null, null, subPixelMethod,
							interpolationMethod, clipOutput));
			if (Utils.isInterrupted())
				return null;
		}

		return new ImagePlus(targetImp.getTitle() + " Aligned", outStack);
	}

	/**
	 * Aligns the target image to the pre-initialised reference and return the shift and score for the alignment
	 * 
	 * @param targetImp
	 * @param windowMethod
	 * @param bounds
	 * @param subPixelMethod
	 * @return [ x_shift, y_shift, score ]
	 */
	public double[] align(ImageProcessor targetIp, WindowMethod windowMethod, Rectangle bounds,
			SubPixelMethod subPixelMethod)
	{
		if (refFHT == null || targetIp == null)
			return null;

		int maxN = refFHT.getWidth();

		// Check correct size
		if (targetIp.getWidth() > maxN || targetIp.getHeight() > maxN)
			return null;

		if (bounds == null)
		{
			bounds = createHalfMaxBounds(refIp.getWidth(), refIp.getHeight(), targetIp.getWidth(), targetIp.getHeight());
		}

		return alignImages(refFHT, s, ss, targetIp, windowMethod, bounds, subPixelMethod);
	}

	private void calculateRollingSums(FloatProcessor ip, double[] s_, double[] ss)
	{
		// Compute the rolling sum and sum of squares
		// s(u,v) = f(u,v) + s(u-1,v) + s(u,v-1) - s(u-1,v-1) 
		// ss(u,v) = f(u,v) * f(u,v) + ss(u-1,v) + ss(u,v-1) - ss(u-1,v-1)
		// where s(u,v) = ss(u,v) = 0 when either u,v < 0

		int maxx = ip.getWidth();
		int maxy = ip.getHeight();
		float[] originalData = (float[]) ip.getPixels();
		double[] data = Tools.toDouble(originalData);

		// First row
		double cs_ = 0; // Column sum
		double css = 0; // Column sum-squares
		for (int i = 0; i < maxx; i++)
		{
			cs_ += data[i];
			css += data[i] * data[i];
			s_[i] = cs_;
			ss[i] = css;
		}

		// Remaining rows:
		// sum = rolling sum of row + sum of row above
		for (int y = 1; y < maxy; y++)
		{
			int i = y * maxx;
			cs_ = 0;
			css = 0;

			// Remaining columns
			for (int x = 0; x < maxx; x++, i++)
			{
				cs_ += data[i];
				css += data[i] * data[i];

				s_[i] = s_[i - maxx] + cs_;
				ss[i] = ss[i - maxx] + css;
			}
		}
	}

	/**
	 * Normalise the correlation matrix using the standard deviation of the region from the reference that is covered by
	 * the target
	 * 
	 * @param subCorrMat
	 * @param s
	 * @param ss
	 * @param targetIp
	 */
	private void normalise(FloatProcessor subCorrMat, double[] s, double[] ss, ImageProcessor targetIp)
	{
		int maxx = subCorrMat.getWidth();
		int maxy = subCorrMat.getHeight();
		Rectangle imageBounds = new Rectangle(0, 0, maxx, maxy); //refImageBounds;

		int NU = targetIp.getWidth();
		int NV = targetIp.getHeight();

		// Locate where the target image was inserted when padding
		int x = targetImageBounds.x; // (maxx - NU) / 2;
		int y = targetImageBounds.y; // (maxy - NV) / 2;

		//IJ.log(String.format("maxx=%d,  maxy=%d, NU=%d, NV=%d, x=%d, y=%d", maxx, maxy, NU, NV, x, y));

		// Calculate overlap:
		// Assume a full size target image relative to the reference and then compensate with the insert location
		int halfNU = maxx / 2 - x;
		int halfNV = maxy / 2 - y;

		// Normalise within the bounds of the largest image (i.e. only allow translation 
		// up to half of the longest edge from the reference or target).
		// The further the translation from the half-max translation the more likely there 
		// can be errors in the normalisation score due to floating point summation errors. 
		// This is observed mainly at the very last pixel overlap between images.
		// To see this set: 
		// union = imageBounds;
		// TODO - More analysis to determine under what conditions this occurs.
		Rectangle union = refImageBounds.union(targetImageBounds);

		// Normalise using the denominator
		float[] data = (float[]) subCorrMat.getPixels();
		float[] newData = new float[data.length];
		for (int yyy = union.y; yyy < union.y + union.height; yyy++)
		{
			int i = yyy * maxx + union.x;
			for (int xxx = union.x; xxx < union.x + union.width; xxx++, i++)
			{
				double sum = 0;
				double sumSquares = 0;

				int minU = xxx - halfNU - 1;
				int maxU = Math.min(minU + NU, maxx - 1);
				int minV = yyy - halfNV - 1;
				int maxV = Math.min(minV + NV, maxy - 1);

				// Compute sum from rolling sum using:
				// sum(u,v) = 
				// + s(u+N-1,v+N-1) 
				// - s(u-1,v+N-1)
				// - s(u+N-1,v-1)
				// + s(u-1,v-1)
				// Note: 
				// s(u,v) = 0 when either u,v < 0
				// s(u,v) = s(umax,v) when u>umax
				// s(u,v) = s(u,vmax) when v>vmax
				// s(u,v) = s(umax,vmax) when u>umax,v>vmax
				// Likewise for ss

				// + s(u+N-1,v+N-1) 
				int index = maxV * maxx + maxU;
				sum += s[index];
				sumSquares += ss[index];

				if (minU >= 0)
				{
					// - s(u-1,v+N-1)
					index = maxV * maxx + minU;
					sum -= s[index];
					sumSquares -= ss[index];
				}
				if (minV >= 0)
				{
					// - s(u+N-1,v-1)
					index = minV * maxx + maxU;
					sum -= s[index];
					sumSquares -= ss[index];

					if (minU >= 0)
					{
						// + s(u-1,v-1)
						index = minV * maxx + minU;
						sum += s[index];
						sumSquares += ss[index];
					}
				}

				// Reset to bounds to calculate the number of pixels
				if (minU < 0)
					minU = 0;
				if (minV < 0)
					minV = 0;

				Rectangle regionBounds = new Rectangle(xxx - halfNU, yyy - halfNV, NU, NV);
				Rectangle r = imageBounds.intersection(regionBounds);

				//int n = (maxU - minU + 1) * (maxV - minV + 1);
				int n = r.width * r.height;

				if (n < 1)
					continue;

				// Get the sum of squared differences
				double residuals = sumSquares - sum * sum / n;

				//				// Check using the original data
				//				double sx = 0;
				//				double ssx = 0;
				//				int nn = 0;
				//				for (int yy = yyy - halfNV; yy < yyy - halfNV + NV; yy++)
				//					for (int xx = xxx - halfNU; xx < xxx - halfNU + NU; xx++)
				//					{
				//						if (xx >= 0 && xx < maxx && yy >= 0 && yy < maxy)
				//						{
				//							float value = normalisedRefIp.getf(xx, yy);
				//							sx += value;
				//							ssx += value * value;
				//							nn++;
				//						}
				//					}
				//				gdsc.fitting.utils.DoubleEquality eq = new gdsc.fitting.utils.DoubleEquality(8, 1e-16);
				//				if (n != nn)
				//				{
				//					System.out.printf("Wrong @ %d,%d %d <> %d\n", xxx, yyy, n, nn);
				//					residuals = ssx - sx * sx / nn;
				//				}
				//				else if (!eq.almostEqualComplement(sx, sum) || !eq.almostEqualComplement(ssx, sumSquares))
				//				{
				//					System.out.printf("Wrong @ %d,%d %g <> %g : %g <> %g\n", xxx, yyy, sx, sum, ssx, sumSquares);
				//					residuals = ssx - sx * sx / nn;
				//				}

				double normalisation = (residuals > 0) ? Math.sqrt(residuals) : 0;

				if (normalisation > 0)
				{
					newData[i] = (float) (data[i] / normalisation);
					// Watch out for normalisation errors which cause problems when displaying the image data.
					if (newData[i] < -1.1f)
						newData[i] = -1.1f;
					if (newData[i] > 1.1f)
						newData[i] = 1.1f;
				}
			}
		}
		subCorrMat.setPixels(newData);
	}

	public static Rectangle createHalfMaxBounds(int width1, int height1, int width2, int height2)
	{
		// Restrict translation so that at least half of the smaller image width/height 
		// is within the larger image (half-max translation)
		int maxx = Math.max(width1, width2);
		int maxy = Math.max(height1, height2);
		maxx /= 2;
		maxy /= 2;
		return new Rectangle(-maxx, -maxy, 2 * maxx, 2 * maxy);
	}

	public static Rectangle createBounds(int minXShift, int maxXShift, int minYShift, int maxYShift)
	{
		int w = maxXShift - minXShift;
		int h = maxYShift - minYShift;
		Rectangle bounds = new Rectangle(minXShift, minYShift, w, h);
		return bounds;
	}

	private boolean isValid(ImageProcessor refIp, ImagePlus targetImp)
	{
		if (refIp == null || targetImp == null)
			return false;

		// Check images have values. No correlation is possible with 
		if (noValue(refIp))
			return false;

		return true;
	}

	/**
	 * @param ip
	 * @return true if the image has not pixels with a value
	 */
	private boolean noValue(ImageProcessor ip)
	{
		for (int i = 0; i < ip.getPixelCount(); i++)
			if (ip.getf(i) != 0)
				return false;
		return true;
	}

	private ImageProcessor alignImages(FHT refFHT, double[] s, double[] ss, ImageProcessor targetIp, int slice,
			WindowMethod windowMethod, Rectangle bounds, FloatProcessor fpCorrelation, FloatProcessor fpNormalised,
			SubPixelMethod subPixelMethod, int interpolationMethod, boolean clipOutput)
	{
		lastXOffset = lastYOffset = 0;

		if (noValue(targetIp))
		{
			// Zero correlation with empty image
			IJ.log(String.format("Best Slice %d  x %g  y %g = %g", slice, 0, 0, 0));
			if (fpCorrelation != null)
				fpCorrelation.setPixels(new float[refFHT.getPixelCount()]);
			if (fpNormalised != null)
				fpNormalised.setPixels(new float[refFHT.getPixelCount()]);
			return targetIp.duplicate();
		}

		// Perform correlation analysis in Fourier space (A and B transform to F and G) 
		// using the complex conjugate of G multiplied by F:
		//   C(u,v) = F(u,v) G*(u,v)		

		int maxN = refFHT.getWidth();

		ImageProcessor paddedTargetIp = padAndZero(targetIp, maxN, windowMethod, targetImageBounds);
		FloatProcessor normalisedTargetIp = normalise(paddedTargetIp);
		FHT targetFHT = fft(normalisedTargetIp, maxN);
		FloatProcessor subCorrMat = correlate(refFHT, targetFHT);

		//new ImagePlus("Unnormalised correlation", subCorrMat.duplicate()).show();

		int originX = (maxN / 2);
		int originY = (maxN / 2);

		// Normalise using the denominator
		if (s != null)
			normalise(subCorrMat, s, ss, targetIp);

		// Copy back result images
		if (fpCorrelation != null)
			fpCorrelation.setPixels(subCorrMat.getPixels());
		if (fpNormalised != null)
			fpNormalised.setPixels(normalisedTargetIp.getPixels());

		Rectangle intersect = new Rectangle(0, 0, subCorrMat.getWidth(), subCorrMat.getHeight());

		// Restrict the translation
		if (bounds != null)
		{
			// Restrict bounds to image limits
			intersect = intersect.intersection(new Rectangle(originX + bounds.x, originY + bounds.y, bounds.width,
					bounds.height));
		}

		int[] iCoord = getPeak(subCorrMat, intersect.x, intersect.y, intersect.width, intersect.height);
		float scoreMax = subCorrMat.getf(iCoord[0], iCoord[1]);
		double[] dCoord = new double[] { iCoord[0], iCoord[1] };

		String estimatedScore = "";
		if (subPixelMethod != SubPixelMethod.NONE)
		{
			double[] centre = null;
			if (subPixelMethod == SubPixelMethod.CUBIC)
			{
				centre = performCubicFit(subCorrMat, iCoord[0], iCoord[1]);
			}
			else
			{
				// Perform sub-peak analysis using the method taken from Jpiv
				centre = performGaussianFit(subCorrMat, iCoord[0], iCoord[1]);
				// Check the centre has not moved too far
				if (!(Math.abs(dCoord[0] - iCoord[0]) < intersect.width / 2 && Math.abs(dCoord[1] - iCoord[1]) < intersect.height / 2))
				{
					centre = null;
				}
			}

			if (centre != null)
			{
				dCoord[0] = centre[0];
				dCoord[1] = centre[1];

				double score = subCorrMat.getBicubicInterpolatedPixel(centre[0], centre[1], subCorrMat);
				//				if (score < -1)
				//					score = -1;
				//				if (score > 1)
				//					score = 1;
				estimatedScore = String.format(" (interpolated score %g)", score);
			}
		}
		else if (IJ.debugMode)
		{
			// Used for debugging - Check if interpolation rounds to a different integer 
			double[] centre = performCubicFit(subCorrMat, iCoord[0], iCoord[1]);
			if (centre != null)
			{
				centre[0] = Math.round(centre[0]);
				centre[1] = Math.round(centre[1]);

				if (centre[0] != iCoord[0] || centre[1] != iCoord[1])
					IJ.log(String.format("Cubic rounded to different integer: %d,%d => %d,%d", iCoord[0], iCoord[1],
							(int) centre[0], (int) centre[1]));
			}

			centre = performGaussianFit(subCorrMat, iCoord[0], iCoord[1]);
			if (centre != null)
			{
				centre[0] = Math.round(centre[0]);
				centre[1] = Math.round(centre[1]);

				if (centre[0] != iCoord[0] || centre[1] != iCoord[1])
					IJ.log(String.format("Gaussian rounded to different integer: %d,%d => %d,%d", iCoord[0], iCoord[1],
							(int) centre[0], (int) centre[1]));
			}
		}

		// The correlation image is the size of the reference.
		// Offset from centre of reference
		lastXOffset = dCoord[0] - originX;
		lastYOffset = dCoord[1] - originY;

		IJ.log(String.format("Best Slice %d  x %g  y %g = %g%s", slice, lastXOffset, lastYOffset, scoreMax,
				estimatedScore));

		// Translate the result and crop to the original size
		if (!doTranslation)
			return targetIp;

		ImageProcessor resultIp = translate(interpolationMethod, targetIp, lastXOffset, lastYOffset, clipOutput);
		return resultIp;
	}

	private double[] alignImages(FHT refFHT, double[] s, double[] ss, ImageProcessor targetIp,
			WindowMethod windowMethod, Rectangle bounds, SubPixelMethod subPixelMethod)
	{
		lastXOffset = lastYOffset = 0;

		if (noValue(targetIp))
		{
			// Zero correlation with empty image
			return new double[] { 0, 0, 0 };
		}

		// Perform correlation analysis in Fourier space (A and B transform to F and G) 
		// using the complex conjugate of G multiplied by F:
		//   C(u,v) = F(u,v) G*(u,v)		

		int maxN = refFHT.getWidth();

		// Allow the input target to be a FHT
		FHT targetFHT;
		if (targetIp instanceof FHT && targetIp.getWidth() == maxN)
		{
			targetFHT = (FHT) targetIp;
		}
		else
		{
			targetFHT = transformTarget(targetIp, windowMethod);
		}
		FloatProcessor subCorrMat = correlate(refFHT, targetFHT);

		int originX = (maxN / 2);
		int originY = (maxN / 2);

		// Normalise using the denominator
		if (s != null)
			normalise(subCorrMat, s, ss, targetIp);

		Rectangle intersect = new Rectangle(0, 0, subCorrMat.getWidth(), subCorrMat.getHeight());

		// Restrict the translation
		if (bounds != null)
		{
			// Restrict bounds to image limits
			intersect = intersect.intersection(new Rectangle(originX + bounds.x, originY + bounds.y, bounds.width,
					bounds.height));
		}

		int[] iCoord = getPeak(subCorrMat, intersect.x, intersect.y, intersect.width, intersect.height);
		double scoreMax = subCorrMat.getf(iCoord[0], iCoord[1]);
		double[] dCoord = new double[] { iCoord[0], iCoord[1] };

		double[] centre = null;
		switch (subPixelMethod)
		{
			case CUBIC:
				centre = performCubicFit(subCorrMat, iCoord[0], iCoord[1]);
				break;
			case GAUSSIAN:
				// Perform sub-peak analysis using the method taken from Jpiv
				centre = performGaussianFit(subCorrMat, iCoord[0], iCoord[1]);
				// Check the centre has not moved too far
				if (!(Math.abs(dCoord[0] - iCoord[0]) < intersect.width / 2 && Math.abs(dCoord[1] - iCoord[1]) < intersect.height / 2))
				{
					centre = null;
				}
				break;
			default:
				break;
		}

		if (centre != null)
		{
			dCoord[0] = centre[0];
			dCoord[1] = centre[1];

			double score = subCorrMat.getBicubicInterpolatedPixel(centre[0], centre[1], subCorrMat);
			if (score < -1)
				score = -1;
			if (score > 1)
				score = 1;
			scoreMax = score;
		}

		// The correlation image is the size of the reference.
		// Offset from centre of reference
		lastXOffset = dCoord[0] - originX;
		lastYOffset = dCoord[1] - originY;

		return new double[] { lastXOffset, lastYOffset, scoreMax };
	}

	/**
	 * Transforms a target image processor for alignment with the initialised reference. The FHT can be passed to the
	 * {@link #align(ImageProcessor, WindowMethod, Rectangle, SubPixelMethod)} method
	 * <p>
	 * If the {@link #init(ImageProcessor, WindowMethod, boolean)} method has not been called this returns null.
	 * 
	 * @param targetIp
	 * @param windowMethod
	 * @return The FHT
	 */
	public FHT transformTarget(ImageProcessor targetIp, WindowMethod windowMethod)
	{
		if (refFHT == null || targetIp == null)
			return null;
		int maxN = refFHT.getWidth();
		FHT targetFHT;
		ImageProcessor paddedTargetIp = padAndZero(targetIp, maxN, windowMethod, targetImageBounds);
		FloatProcessor normalisedTargetIp = normalise(paddedTargetIp);
		targetFHT = fft(normalisedTargetIp, maxN);
		return targetFHT;
	}

	/**
	 * Convert to unit length, return a float processor
	 * 
	 * @param ip
	 * @return
	 */
	public static FloatProcessor normalise(ImageProcessor ip)
	{
		float[] pixels = new float[ip.getPixelCount()];

		// Normalise to unit length and subtract mean
		double sum = 0;
		for (int i = 0; i < pixels.length; i++)
		{
			sum += ip.getf(i) * ip.getf(i);
		}
		if (sum > 0)
		{
			double factor = 1.0 / Math.sqrt(sum);
			for (int i = 0; i < pixels.length; i++)
			{
				pixels[i] = (float) (ip.getf(i) * factor);
			}
		}

		return new FloatProcessor(ip.getWidth(), ip.getHeight(), pixels, null);
	}

	/**
	 * Duplicate and translate the image processor
	 * 
	 * @param interpolationMethod
	 * @param ip
	 * @param xOffset
	 * @param yOffset
	 * @param clipOutput
	 *            Set to true to ensure the output image has the same max as the input. Applies to bicubic
	 *            interpolation
	 * @return New translated processor
	 */
	public static ImageProcessor translate(int interpolationMethod, ImageProcessor ip, double xOffset, double yOffset,
			boolean clipOutput)
	{
		ImageProcessor newIp = ip.duplicate();
		translateProcessor(interpolationMethod, newIp, xOffset, yOffset, clipOutput);
		return newIp;
	}

	/**
	 * Translate the image processor in place
	 * 
	 * @param interpolationMethod
	 * @param ip
	 * @param xOffset
	 * @param yOffset
	 * @param clipOutput
	 *            Set to true to ensure the output image has the same max as the input. Applies to bicubic
	 *            interpolation
	 */
	public static void translateProcessor(int interpolationMethod, ImageProcessor ip, double xOffset, double yOffset,
			boolean clipOutput)
	{
		// Check if interpolation is needed
		if (xOffset == (int) xOffset && yOffset == (int) yOffset)
		{
			interpolationMethod = ImageProcessor.NONE;
		}

		// Bicubic interpolation can generate values outside the input range. 
		// Optionally clip these. This is not applicable for ColorProcessors.
		ImageStatistics stats = null;
		if (interpolationMethod == ImageProcessor.BICUBIC && clipOutput && !(ip instanceof ColorProcessor))
			stats = ImageStatistics.getStatistics(ip, ImageStatistics.MIN_MAX, null);

		ip.setInterpolationMethod(interpolationMethod);
		ip.translate(xOffset, yOffset);

		if (interpolationMethod == ImageProcessor.BICUBIC && clipOutput && !(ip instanceof ColorProcessor))
		{
			float max = (float) stats.max;
			for (int i = ip.getPixelCount(); i-- > 0;)
			{
				if (ip.getf(i) > max)
					ip.setf(i, max);
			}
		}
	}

	/**
	 * Iteratively search the cubic spline surface around the given pixel
	 * to maximise the value.
	 * 
	 * @param fp
	 *            Float processor containing a peak surface
	 * @param i
	 *            The peak x position
	 * @param j
	 *            The peak y position
	 * @return The peak location with sub-pixel accuracy
	 */
	public static double[] performCubicFit(FloatProcessor fp, int i, int j)
	{
		double[] centre = new double[] { i, j };
		// This value will be progressively halved. 
		// Start with a value that allows the number of iterations to fully cover the region +/- 1 pixel
		// TODO - Test if 0.67 is better as this can cover +/- 1 pixel in 2 iterations
		double range = 0.5;
		for (int c = 10; c-- > 0;)
		{
			centre = performCubicFit(fp, centre[0], centre[1], range);
			range /= 2;
		}
		return centre;
	}

	private static double[] performCubicFit(FloatProcessor fp, double x, double y, double range)
	{
		double[] centre = new double[2];
		double peakValue = Double.NEGATIVE_INFINITY;
		for (double x0 : new double[] { x - range, x, x + range })
		{
			for (double y0 : new double[] { y - range, y, y + range })
			{
				double v = fp.getBicubicInterpolatedPixel(x0, y0, fp);
				if (peakValue < v)
				{
					peakValue = v;
					centre[0] = x0;
					centre[1] = y0;
				}
			}
		}
		return centre;
	}

	/**
	 * Perform an interpolated Gaussian fit.
	 * <p>
	 * The following functions for peak finding using Gaussian fitting have been extracted from Jpiv:
	 * http://www.jpiv.vennemann-online.de/
	 * 
	 * @param fp
	 *            Float processor containing a peak surface
	 * @param i
	 *            The peak x position
	 * @param j
	 *            The peak y position
	 * @return The peak location with sub-pixel accuracy
	 */
	public static double[] performGaussianFit(FloatProcessor fp, int i, int j)
	{
		// Extract Pixel block
		float[][] pixelBlock = new float[fp.getWidth()][fp.getHeight()];
		for (int x = pixelBlock.length; x-- > 0;)
		{
			for (int y = pixelBlock[0].length; y-- > 0;)
			{
				if (Float.isNaN(fp.getf(x, y)))
				{
					pixelBlock[x][y] = -1;
				}
				else
				{
					pixelBlock[x][y] = fp.getf(x, y);
				}
			}
		}

		// Extracted as per the code in Jpiv2.PivUtils:
		int x = 0, y = 0, w = fp.getWidth(), h = fp.getHeight();
		int[] iCoord = new int[2];
		double[] dCoord = new double[2];
		// This will weight the function more towards the centre of the correlation pixels.
		// I am not sure if this is necessary.
		//pixelBlock = divideByWeightingFunction(pixelBlock, x, y, w, h);
		iCoord = getPeak(pixelBlock);
		dCoord = gaussianPeakFit(pixelBlock, iCoord[0], iCoord[1]);
		double[] ret = null;
		// more or less acceptable peak fit
		if (Math.abs(dCoord[0] - iCoord[0]) < w / 2 && Math.abs(dCoord[1] - iCoord[1]) < h / 2)
		{
			dCoord[0] += x;
			dCoord[1] += y;
			// Jpiv block is in [Y,X] format (not [X,Y])
			ret = new double[] { dCoord[1], dCoord[0] };

			//    		IJ.log(String.format("Fitted x %d -> %g   y %d -> %g",  
			//    				i, ret[0],
			//    				j, ret[1]));
		}
		return (ret);
	}

	/**
	 * Divides the correlation matrix by a pyramid weighting function.
	 * 
	 * @param subCorrMat
	 *            The biased correlation function
	 * @param xOffset
	 *            If this matrix is merely a search area within a larger
	 *            correlation matrix, this is the offset of the search area.
	 * @param yOffset
	 *            If this matrix is merely a search area within a larger
	 *            correlation matrix, this is the offset of the search area.
	 * @param w
	 *            Width of the original correlation matrix.
	 * @param h
	 *            Height of the original correlation matrix.
	 * @return The corrected correlation function
	 */
	@SuppressWarnings("unused")
	private static float[][] divideByWeightingFunction(float[][] subCorrMat, int xOffset, int yOffset, int w, int h)
	{
		for (int i = 0; i < subCorrMat.length; i++)
		{
			for (int j = 0; j < subCorrMat[0].length; j++)
			{
				subCorrMat[i][j] = subCorrMat[i][j] *
						(Math.abs(j + xOffset - w / 2) / w * 2 + Math.abs(i + yOffset - h / 2) / h * 2 + 1);
			}
		}
		return subCorrMat;
	}

	/**
	 * Finds the highest value in a correlation matrix.
	 * 
	 * @param subCorrMat
	 *            A single correlation matrix.
	 * @return The indices of the highest value {i,j} or {y,x}.
	 */
	private static int[] getPeak(float[][] subCorrMat)
	{
		int[] coord = new int[2];
		float peakValue = 0;
		for (int i = 0; i < subCorrMat.length; ++i)
		{
			for (int j = 0; j < subCorrMat[0].length; ++j)
			{
				if (subCorrMat[i][j] > peakValue)
				{
					peakValue = subCorrMat[i][j];
					coord[0] = j;
					coord[1] = i;
				}
			}
		}
		return (coord);
	}

	/**
	 * Gaussian peak fit.
	 * See Raffel, Willert, Kompenhans;
	 * Particle Image Velocimetry;
	 * 3rd printing;
	 * S. 131
	 * for details
	 * 
	 * @param subCorrMat
	 *            some two dimensional data containing a correlation peak
	 * @param x
	 *            the horizontal peak position
	 * @param y
	 *            the vertical peak position
	 * @return a double array containing the peak position
	 *         with sub pixel accuracy
	 */
	private static double[] gaussianPeakFit(float[][] subCorrMat, int x, int y)
	{
		double[] coord = new double[2];
		// border values
		if (x == 0 || x == subCorrMat[0].length - 1 || y == 0 || y == subCorrMat.length - 1)
		{
			coord[0] = x;
			coord[1] = y;
		}
		else
		{
			coord[0] = x +
					(Math.log(subCorrMat[y][x - 1]) - Math.log(subCorrMat[y][x + 1])) /
					(2 * Math.log(subCorrMat[y][x - 1]) - 4 * Math.log(subCorrMat[y][x]) + 2 * Math
							.log(subCorrMat[y][x + 1]));
			coord[1] = y +
					(Math.log(subCorrMat[y - 1][x]) - Math.log(subCorrMat[y + 1][x])) /
					(2 * Math.log(subCorrMat[y - 1][x]) - 4 * Math.log(subCorrMat[y][x]) + 2 * Math
							.log(subCorrMat[y + 1][x]));
		}
		return (coord);
	}

	private int[] getPeak(FloatProcessor subCorrMat, int minX, int minY, int w, int h)
	{
		int width = subCorrMat.getWidth();
		float max = Float.NEGATIVE_INFINITY;
		int maxi = 0;
		float[] data = (float[]) subCorrMat.getPixels();
		for (int y = minY; y < minY + h; y++)
			for (int x = 0, i = y * width + minX; x < w; x++, i++)
			{
				if (max < data[i])
				{
					max = data[i];
					maxi = i;
				}
			}
		return new int[] { maxi % width, maxi / width };
	}

	private FloatProcessor correlate(FHT refComplex, FHT targetComplex)
	{
		FHT fht = refComplex.conjugateMultiply(targetComplex);
		fht.setShowProgress(false);
		fht.inverseTransform();
		fht.swapQuadrants();
		fht.resetMinAndMax();
		ImageProcessor ip = fht;
		return ip.toFloat(0, null);
	}

	// The following Fast Fourier Transform routines have been extracted from the ij.plugins.FFT class
	FHT fft(ImageProcessor ip, int maxN)
	{
		FHT fht = new FHT(ip);
		fht.setShowProgress(false);
		fht.transform();
		return fht;
	}

	/**
	 * Centre image on zero, padding if necessary to next square power-two above the given max dimension.
	 * <p>
	 * Optionally apply a window function so the image blends smoothly to zero background.
	 * 
	 * @param ip
	 * @param maxN
	 * @return
	 */
	FloatProcessor padAndZero(ImageProcessor ip, int maxN, WindowMethod windowMethod, Rectangle padBounds)
	{
		boolean pad = true;
		int i = 2;
		while (i < maxN)
			i *= 2;
		if (i == maxN && ip.getWidth() == maxN && ip.getHeight() == maxN)
		{
			pad = false;
		}
		maxN = i;

		// This should shift the image so it smoothly blends with a zero background
		// Ideally this would window the image so that the result has an average of zero with smooth edges transitions.
		// However this involves shifting the image and windowing. The average depends on both
		// and so would have to be solved iteratively.		

		if (windowMethod != WindowMethod.NONE)
		{
			// Use separable for speed. 
			//ip = applyWindow(ip, windowMethod);
			ip = applyWindowSeparable(ip, windowMethod);
		}

		// Get average
		double sum = 0;
		for (int ii = 0; ii < ip.getPixelCount(); ii++)
			sum += ip.getf(ii);
		double av = sum / ip.getPixelCount();

		// Create the result image
		FloatProcessor ip2 = new FloatProcessor(maxN, maxN);
		float[] data = (float[]) ip2.getPixels();

		padBounds.width = ip.getWidth();
		padBounds.height = ip.getHeight();
		if (pad)
		{
			// Place into middle of image => Correlation is centre-to-centre alignment
			int x = (maxN - ip.getWidth()) / 2;
			int y = (maxN - ip.getHeight()) / 2;

			padBounds.x = x;
			padBounds.y = y;

			for (int yy = 0, index = 0; yy < ip.getHeight(); yy++)
			{
				int ii = (yy + y) * maxN + x;
				for (int xx = 0; xx < ip.getWidth(); xx++, index++, ii++)
					data[ii] = (float) (ip.getf(index) - av);
			}
		}
		else
		{
			padBounds.x = 0;
			padBounds.y = 0;

			// Copy pixels
			for (int ii = 0; ii < ip.getPixelCount(); ii++)
				data[ii] = (float) (ip.getf(ii) - av);
		}

		return ip2;
	}

	/**
	 * @return the lastXOffset
	 */
	public double getLastXOffset()
	{
		return lastXOffset;
	}

	/**
	 * @return the lastYOffset
	 */
	public double getLastYOffset()
	{
		return lastYOffset;
	}

	/**
	 * Apply a window function to reduce edge artifacts.
	 * <p>
	 * Applied as two 1-dimensional window functions. Faster than the nonseparable form but has direction dependent
	 * corners.
	 * 
	 * @param ip
	 * @param windowMethod
	 * @return
	 */
	public static FloatProcessor applyWindowSeparable(ImageProcessor ip, WindowMethod windowMethod)
	{
		int maxx = ip.getWidth();
		int maxy = ip.getHeight();
		double[] wx = null;
		double[] wy = null;

		switch (windowMethod)
		{
			case HANNING:
				wx = hanning(maxx);
				wy = hanning(maxy);
				break;
			case COSINE:
				wx = cosine(maxx);
				wy = cosine(maxy);
				break;
			case TUKEY:
				wx = tukey(maxx, ALPHA);
				wy = tukey(maxy, ALPHA);
				break;
			default:
				return ip.toFloat(0, null);
		}

		float[] data = new float[ip.getPixelCount()];

		// Calculate total signal of window function applied to image (t1).
		// Calculate total signal of window function applied to a flat signal of intensity 1 (t2).
		// Divide t1/t2 => Result is the mean shift for image so that the average is zero.

		double sumWindowFunction = 0;
		double sumImage = 0;
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double w = wx[x] * wy[y];
				sumWindowFunction += w;
				sumImage += ip.getf(i) * w;
			}
		}

		// Shift to zero
		double shift = sumImage / sumWindowFunction;
		//double sum = 0;
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double value = (ip.getf(i) - shift) * wx[x] * wy[y];
				//sum += value;
				data[i] = (float) value;
			}
		}

		return new FloatProcessor(maxx, maxy, data, null);
	}

	/**
	 * Apply a window function to reduce edge artifacts
	 * <p>
	 * Applied as a nonseparable form.
	 * 
	 * @param ip
	 * @param windowMethod
	 * @return
	 */
	public static FloatProcessor applyWindow(ImageProcessor ip, WindowMethod windowMethod)
	{
		//if (true)
		//	return applyWindowSeparable(ip, windowMethod, duplicate);

		int maxx = ip.getWidth();
		int maxy = ip.getHeight();

		WindowFunction wf = null;
		switch (windowMethod)
		{
			case HANNING: //
				wf = instance.new Hanning();
				break;
			case COSINE:
				wf = instance.new Cosine();
				break;
			case TUKEY:
				wf = instance.new Tukey(ALPHA);
				break;
			default:
				return ip.toFloat(0, null);
		}

		float[] data = new float[ip.getPixelCount()];

		double cx = maxx * 0.5;
		double cy = maxy * 0.5;
		double maxDistance = Math.sqrt(maxx * maxx + maxy * maxy);

		// Calculate total signal of window function applied to image (t1).
		// Calculate total signal of window function applied to a flat signal of intensity 1 (t2).
		// Divide t1/t2 => Result is the mean shift for image so that the average is zero.

		double sumWindowFunction = 0;
		double sumImage = 0;
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double distance = Math.sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
				double w = wf.weight(0.5 - (distance / maxDistance));
				sumWindowFunction += w;
				sumImage += ip.getf(i) * w;
			}
		}

		// Shift to zero
		double shift = sumImage / sumWindowFunction;
		//double sum = 0;
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double distance = Math.sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
				double w = wf.weight(0.5 - (distance / maxDistance));
				double value = (ip.getf(i) - shift) * w;
				//sum += value;
				data[i] = (float) value;
			}
		}

		return new FloatProcessor(maxx, maxy, data, null);
	}

	private static AlignImagesFFT instance = new AlignImagesFFT();
	private static double ALPHA = 0.5;

	interface WindowFunction
	{
		/**
		 * Return the weight for the window at a fraction of the distance from the edge of the window.
		 * 
		 * @param fractionDistance
		 *            (range 0-1)
		 * @return
		 */
		double weight(double fractionDistance);
	}

	class Hanning implements WindowFunction
	{
		public double weight(double fractionDistance)
		{
			return 0.5 * (1 - Math.cos(Math.PI * 2 * fractionDistance));
		}
	}

	class Cosine implements WindowFunction
	{
		public double weight(double fractionDistance)
		{
			return Math.sin(Math.PI * fractionDistance);
		}
	}

	class Tukey implements WindowFunction
	{
		final double alpha;

		public Tukey(double alpha)
		{
			this.alpha = alpha;
		}

		public double weight(double fractionDistance)
		{
			if (fractionDistance < alpha / 2)
				return 0.5 * (1 + Math.cos(Math.PI * (2 * fractionDistance / alpha - 1)));
			if (fractionDistance > 1 - alpha / 2)
				return 0.5 * (1 + Math.cos(Math.PI * (2 * fractionDistance / alpha - 2 / alpha + 1)));
			return 1;
		}
	}

	// Should these be replaced with periodic functions as per use in spectral analysis:
	// http://en.wikipedia.org/wiki/Window_function

	private static double[] window(WindowFunction wf, int N)
	{
		double N_1 = N - 1;
		double[] w = new double[N];
		for (int i = 0; i < N; i++)
			w[i] = wf.weight(i / N_1);
		return w;
	}

	private static double[] hanning(int N)
	{
		return window(instance.new Hanning(), N);
	}

	private static double[] cosine(int N)
	{
		return window(instance.new Cosine(), N);
	}

	private static double[] tukey(int N, double alpha)
	{
		return window(instance.new Tukey(alpha), N);
	}

	/**
	 * @return if false the image will not be translated
	 */
	public boolean isDoTranslation()
	{
		return doTranslation;
	}

	/**
	 * Set to false to prevent the image processor from being translated. The translation can be retrieved using the
	 * lastOffset properties.
	 * 
	 * @param doTranslation
	 *            if false the image will not be translated
	 */
	public void setDoTranslation(boolean doTranslation)
	{
		this.doTranslation = doTranslation;
	}
}
