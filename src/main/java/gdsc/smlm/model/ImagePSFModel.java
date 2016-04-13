package gdsc.smlm.model;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.Sort;

/**
 * Generates a Point Spread Function using an image constructed from diffraction limited spots imaged axially through
 * the plane of focus.
 * <p>
 * The input image must be square. The X/Y centre is the middle of the square. The stack can have any number of slices.
 * The z-centre must be identified. Any pixels below zero will be set to zero.
 * <p>
 * The model can be used to draw a PSF for a point on an image. If the z coordinate is positive then the PSF image from
 * a positive index (above the z-centre) is used. If the z-coordinate is negative then the PSF image from a negative
 * index (below the z-centre) is used. I.e. the input stack is assumed to be imaged axially with increasing z-stage
 * position moving the stage closer to the objective. Thus higher z coordinates correspond to further into the stack.
 */
public class ImagePSFModel extends PSFModel
{
	public static final double DEFAULT_NOISE_FRACTION = 5e-2;
	private static final boolean COM_CHECK = false;

	//private double[][] inputImage;
	private double[][] sumImage;
	private double[][] cumulativeImage;
	private int psfWidth, zCentre;
	private double[][] xyCentre;
	private double unitsPerPixel;
	private double unitsPerSlice;
	private double fwhm;

	/**
	 * Construct the ImagePSF.
	 * <p>
	 * The noise fraction parameter can specify how to remove noise. All pixels below the fraction are set to zero. The
	 * remaining pixels are normalised to 1 to create a PDF for the image.
	 * 
	 * @param image
	 *            The image consisting of a stack of square pixel buffers. The buffers are stored in YX order.
	 * @param zCentre
	 *            The centre of the PSF image
	 * @param unitsPerPixel
	 *            The distance between adjacent X/Y pixels
	 * @param unitsPerSlice
	 *            The distance between adjacent Z pixels
	 * @param fwhm
	 *            The full-width at half-maximum for the z-centre
	 * @param noiseFraction
	 *            The noise fraction
	 */
	public ImagePSFModel(float[][] image, int zCentre, double unitsPerPixel, double unitsPerSlice, double fwhm,
			double noiseFraction)
	{
		super();
		init(image, zCentre, unitsPerPixel, unitsPerSlice, fwhm, noiseFraction);
	}

	/**
	 * Construct the ImagePSF.
	 * 
	 * @param image
	 *            The image consisting of a stack of square pixel buffers. The buffers are stored in YX order.
	 * @param zCentre
	 *            The centre of the PSF image
	 * @param unitsPerPixel
	 *            The distance between adjacent X/Y pixels
	 * @param unitsPerSlice
	 *            The distance between adjacent Z pixels
	 * @param fwhm
	 *            The full-width at half-maximum for the z-centre
	 */
	public ImagePSFModel(float[][] image, int zCentre, double unitsPerPixel, double unitsPerSlice, double fwhm)
	{
		this(image, zCentre, unitsPerPixel, unitsPerSlice, fwhm, DEFAULT_NOISE_FRACTION);
	}

	/**
	 * Private constructor used in the {@link #copy()} method
	 */
	private ImagePSFModel()
	{
		super();
	}

	private void init(float[][] image, int zCentre, double unitsPerPixel, double unitsPerSlice, double fwhm,
			double noiseFraction)
	{
		if (image == null || image.length == 0)
			throw new IllegalArgumentException("Image cannot be null/empty");
		for (int i = 0; i < image.length; i++)
			if (image[i] == null)
				throw new IllegalArgumentException("Image contains null plane");
		if (zCentre < 0 || zCentre >= image.length)
			throw new IllegalArgumentException("z-centre is not within the bounds of the image stack");
		final int size = image[0].length;
		double edge = Math.sqrt(size);
		if (edge != (int) edge)
			throw new IllegalArgumentException("Image planes are not square");
		psfWidth = (int) edge;
		//if (psfWidth % 2 != 1)
		//	throw new IllegalArgumentException("Image edge length is not an odd number");
		xyCentre = new double[image.length][];
		Arrays.fill(xyCentre, new double[] { psfWidth * 0.5, psfWidth * 0.5 });
		for (int i = 1; i < image.length; i++)
			if (image[i].length != size)
				throw new IllegalArgumentException("Image planes are not the same size");
		this.zCentre = zCentre;
		if (unitsPerPixel <= 0 || unitsPerPixel > 1)
			throw new IllegalArgumentException("Units per pixel must be between 0 and 1: " + unitsPerPixel);
		if (image.length > 1 && (unitsPerSlice <= 0 || unitsPerSlice > 1))
			throw new IllegalArgumentException("Units per slice must be between 0 and 1: " + unitsPerSlice);
		this.unitsPerPixel = unitsPerPixel;
		this.unitsPerSlice = unitsPerSlice;

		// Duplicate and convert to double
		this.sumImage = duplicate(image);

		if (noiseFraction > 0)
		{
			//double[] scratch = new double[size];
			for (int i = 0; i < sumImage.length; i++)
			{
				//subtractNoise(sumImage[i], scratch, noiseFraction);
				subtractNoise(sumImage[i], noiseFraction);
			}
		}

		// Debugging display
		//Utils.display("floor", sumImage, psfWidth, psfWidth);

		// Normalise so that the highest intensity frame sums to 1.
		normalise(this.sumImage);

		// Debugging display
		//Utils.display("norm", sumImage, psfWidth, psfWidth);

		// Used for debugging
		if (COM_CHECK)
		{
			// Find X/Y centre using centre of mass
			double sx = 0;
			double sy = 0;
			double s = 0;
			double[] data = sumImage[zCentre];
			double cx = xyCentre[zCentre][0];
			double cy = xyCentre[zCentre][1];
			for (int y = 0; y < psfWidth; y++)
			{
				if (Math.abs(y + 0.5 - cy) > fwhm)
					continue;
				for (int x = 0, j = y * psfWidth; x < psfWidth; x++, j++)
				{
					if (Math.abs(x + 0.5 - cx) > fwhm)
						continue;
					// Centre in middle of pixel
					sx += (x + 0.5) * data[j];
					sy += (y + 0.5) * data[j];
					s += data[j];
				}
			}
			sx = sx / s - cx;
			sy = sy / s - cy;
			System.out.printf("%dx%d centre [ %f %f ] ( %f %f )\n", psfWidth, psfWidth, sx, sy, sx / unitsPerPixel, sy /
					unitsPerPixel);
		}

		// Create a cumulative sum image
		cumulativeImage = new double[sumImage.length][];
		for (int i = 0; i < sumImage.length; i++)
		{
			cumulativeImage[i] = calculateCumulativeImage(sumImage[i]);
		}

		// Debugging display
		//Utils.display("cum", cumulativeImage, psfWidth, psfWidth);

		// Then create a rolling sum table
		for (int i = 0; i < sumImage.length; i++)
		{
			calculateRollingSums(sumImage[i]);
		}
	}

	/**
	 * The noise fraction parameter can specify how to remove noise. Pixels are sorted in descending order and
	 * cumulatively summed. All pixels below the fraction of the total sum are set to zero. The remaining pixels are
	 * normalised to 1 to create a PDF for the image.
	 * 
	 * @param image
	 * @param scratch
	 * @param noiseFraction
	 */
	@SuppressWarnings("unused")
	private void subtractNoise(double[] image, double[] scratch, double noiseFraction)
	{
		// Sort ascending and store the original indices
		double sum = 0;
		for (int i = 0; i < image.length; i++)
		{
			sum += image[i];
			scratch[i] = image[i];
		}
		Arrays.sort(scratch);
		Sort.reverse(scratch);
		// Find the cut-off using the sum and the noise fraction
		final double cutoff = sum - noiseFraction * sum;
		sum = 0;
		int cut = -1;
		for (int i = 0; i < image.length; i++)
		{
			sum += scratch[i];
			if (sum > cutoff)
			{
				// We will zero all pixels after the cut-off that have a lower pixel value
				cut = i;
				while (cut < image.length && scratch[i] == scratch[cut])
					cut++;
				break;
			}
		}
		System.out.printf("Cut = %d, cutoff = %f\n", cut, cutoff);
		if (cut == -1)
			return;
		// All pixels to be included must subtract the noise floor.
		// All pixels below the noise floor are zeroed.
		final double floor = scratch[cut];
		for (int i = 0; i < image.length; i++)
		{
			final double newValue = image[i] - floor;
			image[i] = (newValue > 0) ? newValue : 0;
		}
	}

	/**
	 * The noise fraction parameter can specify how to remove noise. All pixels below the fraction of the maximum are
	 * set to zero. The remaining pixels are adjusted by the height of the first pixel below the cutoff.
	 * 
	 * @param image
	 * @param noiseFraction
	 */
	private void subtractNoise(double[] image, double noiseFraction)
	{
		double max = 0, sum = 0;
		for (final double v : image)
		{
			if (max < v)
				max = v;
			sum += v;
		}
		if (max <= 0)
			return;
		final double cutoff = max * noiseFraction;
		// Find highest value below cutoff
		double floor = 0;
		for (final double v : image)
		{
			if (v < cutoff && floor < v)
				floor = v;
		}
		// All pixels to be included must subtract the noise floor.
		// All pixels below the noise floor are zeroed.
		double sum2 = 0;
		for (int i = 0; i < image.length; i++)
		{
			final double newValue = image[i] - floor;
			image[i] = (newValue > 0) ? newValue : 0;
			sum2 += image[i];
		}
		// Re-normalise to the same intensity
		double scale = sum / sum2;
		for (int i = 0; i < image.length; i++)
			image[i] *= scale;
	}

	private double[][] duplicate(float[][] image)
	{
		final int size = image[0].length;
		double[][] duplicate = new double[image.length][size];
		for (int i = 0; i < image.length; i++)
			for (int j = 0; j < size; j++)
			{
				final float f = image[i][j];
				// Ignore negative pixels
				if (f > 0)
					duplicate[i][j] = f;
			}
		return duplicate;
	}

	/**
	 * Normalise the image so that the brightest frame has a sum of 1
	 * 
	 * @param image
	 */
	private void normalise(double[][] image)
	{
		if (image == null || image.length == 0)
			return;

		double max = 0;
		for (int i = 0; i < image.length; i++)
			max = FastMath.max(max, sum(image[i]));

		if (max <= 0)
			return;

		for (int i = 0; i < image.length; i++)
		{
			final double[] data = image[i];
			for (int j = 0; j < data.length; j++)
			{
				data[j] /= max;
			}
		}
	}

	private double sum(double[] data)
	{
		double sum = 0;
		for (double f : data)
			sum += f;
		return sum;
	}

	private double[] calculateCumulativeImage(double[] s)
	{
		boolean normalised = true;
		double[] c = new double[s.length + 1];
		if (normalised)
		{
			// Normalised image as input
			double sum = 0;
			for (int i = 0; i < s.length; i++)
			{
				sum += s[i];
				c[i + 1] = sum;
			}
		}
		else
		{
			// Normalise as we go
			int n = s.length;
			double mean = 0, sum = 0;

			for (int i = 0; i < n; i++)
			{
				mean += (s[i] - mean) / ((double) (i + 1));
			}

			c[0] = 0;

			for (int i = 0; i < n; i++)
			{
				sum += (s[i] / mean) / n;
				c[i + 1] = sum;
			}
		}
		return c;
	}

	private void calculateRollingSums(double[] s)
	{
		// Compute the rolling sum
		// s(u,v) = f(u,v) + s(u-1,v) + s(u,v-1) - s(u-1,v-1) 
		// where s(u,v) = 0 when either u,v < 0

		final int maxx = psfWidth;
		final int maxy = psfWidth;

		// First row
		double cs = 0; // Column sum
		for (int i = 0; i < maxx; i++)
		{
			cs += s[i];
			s[i] = cs;
		}

		// Remaining rows:
		// sum = rolling sum of row + sum of row above
		for (int y = 1, i = maxx; y < maxy; y++)
		{
			cs = 0;

			// Remaining columns
			for (int x = 0; x < maxx; x++, i++)
			{
				cs += s[i];
				s[i] = (s[i - maxx] + cs);
			}
		}
	}

	/**
	 * @param randomGenerator
	 * @param image
	 *            The image consisting of a stack of square pixel buffers. The buffers are stored in YX order.
	 * @param zCentre
	 *            The centre of the PSF image
	 * @param unitsPerPixel
	 *            The distance between adjacent X/Y pixels
	 * @param unitsPerSlice
	 *            The distance between adjacent Z pixels
	 * @param fwhm
	 *            The full-width at half-maximum for the z-centre
	 */
	public ImagePSFModel(RandomGenerator randomGenerator, float[][] image, int zCentre, double unitsPerPixel,
			double unitsPerSlice, double fwhm)
	{
		super(randomGenerator);
		init(image, zCentre, unitsPerPixel, unitsPerSlice, fwhm, DEFAULT_NOISE_FRACTION);
	}

	/**
	 * @param randomDataGenerator
	 * @param image
	 *            The image consisting of a stack of square pixel buffers. The buffers are stored in YX order.
	 * @param zCentre
	 *            The centre of the PSF image
	 * @param unitsPerPixel
	 *            The distance between adjacent X/Y pixels
	 * @param unitsPerSlice
	 *            The distance between adjacent Z pixels
	 * @param fwhm
	 *            The full-width at half-maximum for the z-centre
	 */
	public ImagePSFModel(RandomDataGenerator randomDataGenerator, float[][] image, int zCentre, double unitsPerPixel,
			double unitsPerSlice, double fwhm)
	{
		super(randomDataGenerator);
		init(image, zCentre, unitsPerPixel, unitsPerSlice, fwhm, DEFAULT_NOISE_FRACTION);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#create3D(float[], int, int, double, double, double, double, boolean)
	 */
	public double create3D(float[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		try
		{
			return drawPSF(data, width, height, sum, x0, x1, x2, poissonNoise);
		}
		catch (IllegalArgumentException e)
		{
			return 0;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#create3D(double[], int, int, double, double, double, double, boolean)
	 */
	public double create3D(double[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		try
		{
			return drawPSF(data, width, height, sum, x0, x1, x2, poissonNoise);
		}
		catch (IllegalArgumentException e)
		{
			return 0;
		}
	}

	public double drawPSF(float[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		final int slice = getSlice(x2);
		if (slice < 0 || slice >= xyCentre.length)
		{
			return insert(data, 0, 0, 0, 0, 0, null, false);
		}

		// Parameter check
		if (width < 1)
			throw new IllegalArgumentException("Width cannot be less than 1");
		if (height < 1)
			throw new IllegalArgumentException("Height cannot be less than 1");
		if (data == null)
			data = new float[width * height];
		else if (data.length < width * height)
			throw new IllegalArgumentException("Data length cannot be smaller than width * height");

		// Evaluate the PSF over the full range about the PSF centre
		final double cx = xyCentre[slice][0] * unitsPerPixel;
		final double cy = xyCentre[slice][1] * unitsPerPixel;
		final int x0min = clip((int) (x0 - cx), width);
		final int x1min = clip((int) (x1 - cy), height);
		final int x0max = clip((int) Math.ceil(x0 - cx + psfWidth * unitsPerPixel), width);
		final int x1max = clip((int) Math.ceil(x1 - cy + psfWidth * unitsPerPixel), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Shift centre to origin and draw the PSF
		double[] psf = drawPSF(x0range, x1range, sum, x0 - x0min, x1 - x1min, x2, true);

		return insert(data, x0min, x1min, x0max, x1max, width, psf, poissonNoise);
	}

	public double drawPSF(double[] data, final int width, final int height, final double sum, double x0, double x1,
			double x2, boolean poissonNoise)
	{
		final int slice = getSlice(x2);
		if (slice < 0 || slice >= xyCentre.length)
		{
			return insert(data, 0, 0, 0, 0, 0, null, false);
		}

		// Parameter check
		if (width < 1)
			throw new IllegalArgumentException("Width cannot be less than 1");
		if (height < 1)
			throw new IllegalArgumentException("Height cannot be less than 1");
		if (data == null)
			data = new double[width * height];
		else if (data.length < width * height)
			throw new IllegalArgumentException("Data length cannot be smaller than width * height");

		// Evaluate the PSF over the full range about the PSF centre
		final double cx = xyCentre[slice][0] * unitsPerPixel;
		final double cy = xyCentre[slice][1] * unitsPerPixel;
		final int x0min = clip((int) (x0 - cx), width);
		final int x1min = clip((int) (x1 - cy), height);
		final int x0max = clip((int) Math.ceil(x0 - cx + psfWidth * unitsPerPixel), width);
		final int x1max = clip((int) Math.ceil(x1 - cy + psfWidth * unitsPerPixel), height);

		final int x0range = x0max - x0min;
		final int x1range = x1max - x1min;

		// min should always be less than max
		if (x0range < 1)
			throw new IllegalArgumentException("Dimension 0 range not within data bounds");
		if (x1range < 1)
			throw new IllegalArgumentException("Dimension 1 range not within data bounds");

		// Shift centre to origin and draw the PSF
		double[] psf = drawPSF(x0range, x1range, sum, x0 - x0min, x1 - x1min, x2, true);

		return insert(data, x0min, x1min, x0max, x1max, width, psf, poissonNoise);
	}

	/**
	 * Construct a PSF function based at the origin using the specified range in each dimension.
	 * 
	 * @param x0range
	 *            The maximum range in dimension 0 (width)
	 * @param x1range
	 *            The maximum range in dimension 1 (height)
	 * @param sum
	 *            The integral
	 * @param x0
	 *            The centre in dimension 0
	 * @param x1
	 *            The centre in dimension 1
	 * @param x2
	 *            The centre in dimension 2
	 * @param interpolate
	 *            Compute bilinear interpolation. Not using this option can result in artifacts if the PSF and image
	 *            pixels are similar sizes (i.e. unitsPerPixel is close to 1)
	 * @return The data (packed in yx order, length = x0range * x1range)
	 * @return
	 */
	public double[] drawPSF(int x0range, int x1range, double sum, double x0, double x1, double x2, boolean interpolate)
	{
		double[] data = new double[x0range * x1range];

		final int slice = getSlice(x2);
		if (slice < 0 || slice >= sumImage.length)
			return data;
		final double[] sumPsf = sumImage[slice];

		// Determine PSF blocks. 
		// We need to map vertices of each pixel in the PSF onto the output image.
		// Use (psfX+1) to describe upper bounds of pixel => mapping back describes lower bounds of PSF pixel 
		// outX = (psfX + 1 - psfX0) * unitsPerPixel + origin
		// =>
		// psfX = (outX - origin) / unitsPerPixel + psfX0 - 1

		// Note that if the PSF pixels do not exactly fit into the image pixels, then the lookup index 
		// will be rounded. This will cause the inserted PSF to be incorrect. The correct method is
		// to do linear interpolation between the pixel sums.

		double[] u = createInterpolationLookup(x0range, x0, xyCentre[slice][0]);
		double[] v = createInterpolationLookup(x1range, x1, xyCentre[slice][1]);

		int[] lu = new int[u.length];
		int[] lv = new int[v.length];

		if (interpolate)
		{
			for (int i = 0; i < u.length; i++)
				lu[i] = (int) u[i];
			for (int i = 0; i < v.length; i++)
				lv[i] = (int) v[i];
		}
		else
		{
			// If we are not interpolating then we use the centre of the pixel so convert 
			// the vertices to the pixel centre using a 0.5 pixel offset.
			for (int i = 0; i < u.length; i++)
				lu[i] = (int) (u[i] + 0.5);
			for (int i = 0; i < v.length; i++)
				lv[i] = (int) (v[i] + 0.5);
		}

		// Data will be the sum of the input pixels from (u,v) to (u+1,v+1)

		// Check data can be inserted from the PSF
		if (lu[0] > psfWidth - 1 || lv[0] > psfWidth - 1 || lu[x0range] < 0 || lv[x1range] < 0)
			return data;

		for (int y = 0; y < x1range; y++)
		{
			if (lv[y] > psfWidth - 1)
				break;
			if (lv[y + 1] < 0)
				continue;
			final int lowerV = lv[y];
			final int upperV = FastMath.min(lv[y + 1], psfWidth - 1);
			for (int x = 0, i = y * x0range; x < x0range; x++, i++)
			{
				if (lu[x] > psfWidth - 1)
					break;
				if (lu[x + 1] < 0)
					continue;

				if (interpolate)
				{
					// Do linear interpolation
					// sum = 
					// + s(upperU,upperV) 
					// - s(lowerU,upperV)
					// - s(upperU,lowerV)
					// + s(lowerU,lowerV)
					double uUuV = interpolate(sumPsf, lu[x + 1], lv[y + 1], u[x + 1], v[y + 1]);
					double lUuV = interpolate(sumPsf, lu[x], lv[y + 1], u[x], v[y + 1]);
					double uUlV = interpolate(sumPsf, lu[x + 1], lv[y], u[x + 1], v[y]);
					double lUlV = interpolate(sumPsf, lu[x], lv[y], u[x], v[y]);
					data[i] = uUuV - lUuV - uUlV + lUlV;
				}
				else
				{
					// No interpolation
					data[i] = sum(sumPsf, lu[x], lowerV, lu[x + 1], upperV);

					//// Check using the original input image that the sum is correct
					//double s = 0;
					//for (int vv = lowerV + 1; vv <= upperV; vv++)
					//{
					//	if (vv > psfWidth-1)
					//		break;
					//	if (vv < 0)
					//		continue;
					//	for (int uu = u[x] + 1, index = vv * psfWidth + u[x] + 1; uu <= u[x + 1]; uu++, index++)
					//	{
					//		if (uu > psfWidth-1)
					//			break;
					//		if (uu < 0)
					//			continue;
					//		s += psf[index];
					//	}
					//}
				}
			}
		}

		// The PSF is normalised so the brightest plane is 1. Just multiply by the sum to create the integral.
		for (int i = 0; i < data.length; i++)
			data[i] *= sum;

		return data;
	}

	private double[] createInterpolationLookup(int range, double origin, double xyCentre)
	{
		double[] pixel = new double[range + 1];
		for (int i = 0; i < pixel.length; i++)
		{
			pixel[i] = ((i - origin) / unitsPerPixel + xyCentre - 1);
		}
		return pixel;
	}

	private double sum(double[] s, int lowerU, int lowerV, int upperU, int upperV)
	{
		// Compute sum from rolling sum using:
		// sum = 
		// + s(upperU,upperV) 
		// - s(lowerU,upperV)
		// - s(upperU,lowerV)
		// + s(lowerU,lowerV)
		// Note: 
		// s(u,v) = 0 when either u,v < 0
		// s(u,v) = s(umax,v) when u>umax
		// s(u,v) = s(u,vmax) when v>vmax
		// s(u,v) = s(umax,vmax) when u>umax,v>vmax

		//		if (lowerU > psfWidth - 1 || lowerV > psfWidth - 1)
		//			return 0;
		//		if (upperU < 0 || upperV < 0)
		//			return 0;

		upperU = FastMath.min(upperU, psfWidth - 1);
		//upperV = FastMath.min(upperV, psfWidth - 1);

		int index = upperV * psfWidth + upperU;
		double sum = s[index];

		if (lowerU >= 0)
		{
			index = upperV * psfWidth + lowerU;
			sum -= s[index];
		}
		if (lowerV >= 0)
		{
			index = lowerV * psfWidth + upperU;
			sum -= s[index];

			if (lowerU >= 0)
			{
				// + s(u-1,v-1)
				index = lowerV * psfWidth + lowerU;
				sum += s[index];
			}
		}
		return sum;
	}

	private double safeSum(double[] s, int lowerU, int lowerV, int upperU, int upperV)
	{
		// Compute sum from rolling sum using:
		// sum = 
		// + s(upperU,upperV) 
		// - s(lowerU,upperV)
		// - s(upperU,lowerV)
		// + s(lowerU,lowerV)
		// Note: 
		// s(u,v) = 0 when either u,v < 0
		// s(u,v) = s(umax,v) when u>umax
		// s(u,v) = s(u,vmax) when v>vmax
		// s(u,v) = s(umax,vmax) when u>umax,v>vmax

		if (lowerU > psfWidth - 1 || lowerV > psfWidth - 1)
			return 0;
		if (upperU < 0 || upperV < 0)
			return 0;

		upperU = FastMath.min(upperU, psfWidth - 1);
		upperV = FastMath.min(upperV, psfWidth - 1);

		int index = upperV * psfWidth + upperU;
		double sum = s[index];

		if (lowerU >= 0)
		{
			index = upperV * psfWidth + lowerU;
			sum -= s[index];
		}
		if (lowerV >= 0)
		{
			index = lowerV * psfWidth + upperU;
			sum -= s[index];

			if (lowerU >= 0)
			{
				// + s(u-1,v-1)
				index = lowerV * psfWidth + lowerU;
				sum += s[index];
			}
		}
		return sum;
	}

	private double safeSum(double[] s, int upperU, int upperV)
	{
		if (upperU < 0 || upperV < 0)
			return 0;
		upperU = FastMath.min(upperU, psfWidth - 1);
		upperV = FastMath.min(upperV, psfWidth - 1);

		int index = upperV * psfWidth + upperU;
		return s[index];
	}

	/**
	 * Find the sum to the point x,y. x,y is assumed to be within the range x0,y0 to x0+1,y0+1.
	 * 
	 * @param sum
	 * @param x0
	 * @param y0
	 * @param x
	 * @param y
	 * @return The sum
	 */
	private double interpolate(double[] sum, int x0, int y0, double x, double y)
	{
		double sum_00_x0y0 = safeSum(sum, x0, y0);
		double sum_00_x1y0 = safeSum(sum, x0 + 1, y0);
		double sum_00_x0y1 = safeSum(sum, x0, y0 + 1);
		double sum_x0y0_x1y1 = safeSum(sum, x0, y0, x0 + 1, y0 + 1);

		x -= x0;
		y -= y0;
		return x * (sum_00_x1y0 - sum_00_x0y0) + y * (sum_00_x0y1 - sum_00_x0y0) + sum_00_x0y0 + x * y * sum_x0y0_x1y1;
	}

	private int clip(int x, int max)
	{
		if (x < 0)
			x = 0;
		if (x > max)
			x = max;
		return x;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.model.PSFModel#getFwhm()
	 */
	public double getFwhm()
	{
		return fwhm;
	}

	/**
	 * Produce a shallow copy of this object. This shares the pre-computed PSF image data but will allow
	 * the copy to store its own version of the most recently created PSF. The copy has a new random data generator.
	 * 
	 * @return A shallow copy of this object
	 */
	public ImagePSFModel copy()
	{
		ImagePSFModel model = new ImagePSFModel();
		model.sumImage = sumImage;
		model.cumulativeImage = cumulativeImage;
		//model.inputImage = inputImage;
		model.psfWidth = psfWidth;
		model.xyCentre = xyCentre;
		model.zCentre = zCentre;
		model.unitsPerPixel = unitsPerPixel;
		model.unitsPerSlice = unitsPerSlice;
		model.fwhm = fwhm;
		return model;
	}

	@Override
	public int sample3D(float[] data, int width, int height, int n, double x0, double x1, double x2)
	{
		if (n <= 0)
			return insertSample(data, width, height, null, null);
		double[][] sample = sample(n, x0, x1, x2);
		return insertSample(data, width, height, sample[0], sample[1]);
	}

	@Override
	public int sample3D(double[] data, int width, int height, int n, double x0, double x1, double x2)
	{
		if (n <= 0)
			return insertSample(data, width, height, null, null);
		double[][] sample = sample(n, x0, x1, x2);
		return insertSample(data, width, height, sample[0], sample[1]);
	}

	private double[][] sample(final int n, double x0, double x1, double x2)
	{
		final int slice = getSlice(x2);
		if (slice < 0 || slice >= sumImage.length)
			return new double[][] { null, null };

		final double[] sumPsf = cumulativeImage[slice];

		final RandomGenerator randomX, randomY;

		// Use the input generator
		randomX = rand.getRandomGenerator();
		randomY = rand.getRandomGenerator();

		//// Debugging - Use a uniform distribution to sample x
		//randomX = new AbstractRandomGenerator()
		//{
		//	int pos = 0;
		//
		//	@Override
		//	public double nextDouble()
		//	{
		//		double p = (double) pos / n;
		//		if (pos++ >= n)
		//			pos = 0;
		//		return p;
		//	}
		//
		//	@Override
		//	public void setSeed(long seed)
		//	{
		//		pos = Math.abs((int) seed) % n;
		//	}
		//};
		//// Debugging - Use a fixed distribution to sample y
		//randomY = new AbstractRandomGenerator()
		//{
		//	public double nextDouble()
		//	{
		//		return 0.5;
		//	}
		//
		//	@Override
		//	public void setSeed(long seed)
		//	{
		//	}
		//};

		// Ensure the generated index is adjusted to the correct position
		// The index will be generated at 0,0 of a pixel in the PSF image.
		// We must subtract the PSF centre so that the middle coords are zero.

		x0 -= xyCentre[slice][0] * unitsPerPixel;
		x1 -= xyCentre[slice][1] * unitsPerPixel;
		//x0 -= 0.5 * psfWidth * unitsPerPixel;
		//x1 -= 0.5 * psfWidth * unitsPerPixel;

		final double max = sumPsf[sumPsf.length - 1];
		double[] x = new double[n];
		double[] y = new double[n];
		int count = 0;
		double sx = 0, sy = 0, s = 0;
		for (int i = 0; i < n; i++)
		{
			final double p = randomX.nextDouble();
			// If outside the observed PSF then skip 
			if (p > max)
				continue;
			final int index = findIndex(sumPsf, p);

			// Interpolate xi using the fraction of the pixel
			double xi = index % psfWidth;
			xi += (p - sumPsf[index]) / (sumPsf[index + 1] - sumPsf[index]);
			// Add random dither within pixel for y
			final double yi = randomY.nextDouble() + (index / psfWidth);

			if (COM_CHECK)
			{
				final double v = 1;
				sx += xi * v;
				sy += yi * v;
				s += v;
			}

			x[count] = x0 + (xi * this.unitsPerPixel);
			y[count] = x1 + (yi * this.unitsPerPixel);
			count++;
		}

		if (COM_CHECK)
		{
			sx = sx / s - xyCentre[slice][0];
			sy = sy / s - xyCentre[slice][1];
			System.out.printf("%dx%d sample centre [ %f %f ] ( %f %f )\n", psfWidth, psfWidth, sx, sy, sx /
					unitsPerPixel, sy / unitsPerPixel);
		}

		x = Arrays.copyOf(x, count);
		y = Arrays.copyOf(y, count);

		return new double[][] { x, y };
	}

	private int getSlice(double x2)
	{
		// We assume the PSF was imaged axially with increasing z-stage position (moving the stage 
		// closer to the objective). Thus higher z-coordinate are for higher slice numbers.
		final int slice = (int) Math.round(x2 / unitsPerSlice) + zCentre;
		return slice;
	}

	/**
	 * Find the index such that sum[index] <= p < sum[index+1]
	 * 
	 * @param sum
	 * @param p
	 * @return the index (or -1)
	 */
	private int findIndex(double[] sum, double p)
	{
		/* perform binary search */
		int upper = sum.length - 1;
		int lower = 0;

		while (upper - lower > 1)
		{
			//final int mid = (upper + lower) / 2;
			final int mid = upper + lower >>> 1;

			if (p >= sum[mid])
			{
				lower = mid;
			}
			else
			{
				upper = mid;
			}
		}

		///* sanity check the result */
		// We do not need this:
		// The lowest value for p is always 0.
		// The check against the max has already been performed.
		//if (p < sum[lower] || p >= sum[lower + 1])
		//{
		//	return -1;
		//}

		return lower;
	}

	/**
	 * Set the PSF centre for the given slice. The centre must be with the width/size of the PSF.
	 * 
	 * @param slice
	 * @param x0
	 * @param x1
	 * @return True if set
	 */
	public boolean setCentre(int slice, double x0, double x1)
	{
		if (slice < 0 || slice >= xyCentre.length)
			return false;
		if (x0 < 0 || x0 >= psfWidth)
			return false;
		if (x1 < 0 || x1 >= psfWidth)
			return false;
		xyCentre[slice] = new double[] { x0, x1 };
		return true;
	}

	/**
	 * Set the PSF centre for the given slice. The centre is relative to 0,0.
	 * 
	 * @param slice
	 * @param x0
	 * @param x1
	 * @return True if set
	 */
	public boolean setRelativeCentre(int slice, double x0, double x1)
	{
		x0 += 0.5 * psfWidth;
		x1 += 0.5 * psfWidth;
		return setCentre(slice, x0, x1);
	}
}
