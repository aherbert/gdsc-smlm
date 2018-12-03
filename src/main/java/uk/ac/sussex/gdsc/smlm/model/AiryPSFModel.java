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
package uk.ac.sussex.gdsc.smlm.model;

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;
import org.apache.commons.math3.util.FastMath;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Contains methods for generating models of a Point Spread Function using a Airy pattern.
 * <p>
 * Out-of-focus regions are computed using a width spreading of the Airy pattern. A true diffraction model for
 * out-of-focus regions is not implemented.
 */
public class AiryPSFModel extends PSFModel
{
    private double zeroW0, zeroW1;
    private double w0;
    private double w1;
    private double zDepth = 0;
    private int ring = 2;
    private boolean singlePixelApproximation = false;
    private int minSamplesPerDimension = 2;
    private int maxSamplesPerDimension = 50;

    // Used for the random sampling of the Airy function
    private static int SAMPLE_RINGS = 4;
    private static double[] r;
    private static double[] sum;
    private static PolynomialSplineFunction spline;

    /**
     * The zeros of J1(x) corresponding to the rings of the Airy pattern.
     */
    public static double[] RINGS = { 0, 3.8317, 7.0156, 10.1735, 13.3237, 16.4706 };

    /**
     * The Airy power corresponding to the rings of the Airy pattern.
     */
    public static double[] POWER;
    static
    {
        POWER = new double[RINGS.length];
        for (int i = 1; i < POWER.length; i++)
            POWER[i] = AiryPattern.power(RINGS[i]);
    }

    /**
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     */
    public AiryPSFModel(double w0, double w1)
    {
        super();
        this.zeroW0 = w0;
        this.zeroW1 = w1;
    }

    /**
     * Instantiates a new airy PSF model.
     *
     * @param randomGenerator
     *            the random generator
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     */
    public AiryPSFModel(RandomGenerator randomGenerator, double w0, double w1)
    {
        super(randomGenerator);
        this.zeroW0 = w0;
        this.zeroW1 = w1;
    }

    /**
     * Instantiates a new airy PSF model.
     *
     * @param randomGenerator
     *            the random generator
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     * @param zDepth
     *            the Z-depth where the 3D PSF is sqrt(2) the width (1.41 x FWHM)
     */
    public AiryPSFModel(RandomGenerator randomGenerator, double w0, double w1, double zDepth)
    {
        super(randomGenerator);
        this.zeroW0 = w0;
        this.zeroW1 = w1;
        setzDepth(zDepth);
    }

    /**
     * Instantiates a new airy PSF model.
     *
     * @param randomDataGenerator
     *            the random data generator
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     */
    public AiryPSFModel(RandomDataGenerator randomDataGenerator, double w0, double w1)
    {
        super(randomDataGenerator);
        this.zeroW0 = w0;
        this.zeroW1 = w1;
    }

    /**
     * Instantiates a new airy PSF model.
     *
     * @param randomDataGenerator
     *            the random data generator
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     * @param zDepth
     *            the Z-depth where the 3D PSF is sqrt(2) the width (1.41 x FWHM)
     */
    public AiryPSFModel(RandomDataGenerator randomDataGenerator, double w0, double w1, double zDepth)
    {
        super(randomDataGenerator);
        this.zeroW0 = w0;
        this.zeroW1 = w1;
        setzDepth(zDepth);
    }

    /**
     * Private constructor used in the {@link #copy()} method
     */
    private AiryPSFModel()
    {
        super();
    }

    /** {@inheritDoc} */
    @Override
    public double create3D(float[] data, final int width, final int height, final double sum, double x0, double x1,
            double x2, boolean poissonNoise)
    {
        if (sum == 0)
            return 0;
        final double scale = createWidthScale(x2);
        try
        {
            final double d = airy2D(data, width, height, sum, x0, x1, scale * zeroW0, scale * zeroW1, poissonNoise);
            return d;
        }
        catch (final IllegalArgumentException e)
        {
            //System.out.println(e.getMessage());
            return 0;
        }
    }

    /** {@inheritDoc} */
    @Override
    public double create3D(double[] data, final int width, final int height, final double sum, double x0, double x1,
            double x2, boolean poissonNoise)
    {
        if (sum == 0)
            return 0;
        final double scale = createWidthScale(x2);
        try
        {
            return airy2D(data, width, height, sum, x0, x1, scale * zeroW0, scale * zeroW1, poissonNoise);
        }
        catch (final IllegalArgumentException e)
        {
            //System.out.println(e.getMessage());
            return 0;
        }
    }

    /**
     * Generate a scale so that at the configured zDepth the scale is sqrt(2).
     *
     * @param z
     *            the z
     * @return The scale
     */
    private double createWidthScale(double z)
    {
        if (zDepth == 0) // Not 3D data
            return 1;

        // Holtzer z-model:
        // Ref: Holtzer, L., Meckel, T. & Schmidt, T. Nanometric three-dimensional tracking of individual quantum dots in cells.
        // Applied Physics Letters 90, 1–3 (2007).
        //   width = sqrt(1 + z^2 / d^2)

        z /= zDepth; // Scale so z=1 at the configured z-depth
        return Math.sqrt(1.0 + z * z);
    }

    /**
     * Construct a Airy pattern on the provided data. Only evaluates the function up to the configured dark ring.
     *
     * @param data
     *            The data (can be null)
     * @param width
     *            The data width
     * @param height
     *            The data height
     * @param sum
     *            The integral
     * @param x0
     *            The centre in dimension 0
     * @param x1
     *            The centre in dimension 1
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     * @param poissonNoise
     *            Add Poisson noise
     * @return The total sum added to the image (useful when poissonNoise is added)
     */
    public double airy2D(float[] data, final int width, final int height, final double sum, double x0, double x1,
            double w0, double w1, boolean poissonNoise)
    {
        if (sum == 0)
            return 0;
        // Parameter check
        checkSize(width, height);
        if (data == null)
            data = new float[width * height];
        else if (data.length < width * height)
            throw new IllegalArgumentException("Data length cannot be smaller than width * height");

        w0 = Math.abs(w0);
        w1 = Math.abs(w1);

        // The second zero (dark ring of an Airy pattern is at 7.0156 of the width
        final int x0min = clip((int) (x0 - RINGS[ring] * w0), width);
        final int x1min = clip((int) (x1 - RINGS[ring] * w1), height);
        final int x0max = clip((int) Math.ceil(x0 + RINGS[ring] * w0), width);
        final int x1max = clip((int) Math.ceil(x1 + RINGS[ring] * w1), height);

        final int x0range = x0max - x0min;
        final int x1range = x1max - x1min;

        // min should always be less than max
        if (x0range < 1)
            throw new IllegalArgumentException("Dimension 0 range not within data bounds");
        if (x1range < 1)
            throw new IllegalArgumentException("Dimension 1 range not within data bounds");

        // Shift centre to origin and compute gaussian
        final double[] gauss = airy2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, w0, w1);

        return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
    }

    /**
     * Construct a Airy pattern on the provided data. Only evaluates the function up to the configured dark ring.
     *
     * @param data
     *            The data (can be null)
     * @param width
     *            The data width
     * @param height
     *            The data height
     * @param sum
     *            The integral
     * @param x0
     *            The centre in dimension 0
     * @param x1
     *            The centre in dimension 1
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     * @param poissonNoise
     *            Add Poisson noise
     * @return The total sum added to the image (useful when poissonNoise is added)
     */
    public double airy2D(double[] data, final int width, final int height, final double sum, double x0, double x1,
            double w0, double w1, boolean poissonNoise)
    {
        if (sum == 0)
            return 0;
        // Parameter check
        checkSize(width, height);
        if (data == null)
            data = new double[width * height];
        else if (data.length < width * height)
            throw new IllegalArgumentException("Data length cannot be smaller than width * height");

        w0 = Math.abs(w0);
        w1 = Math.abs(w1);

        // The second zero (dark ring of an Airy pattern is at 7.0156 of the width
        final int x0min = clip((int) (x0 - RINGS[ring] * w0), width);
        final int x1min = clip((int) (x1 - RINGS[ring] * w1), height);
        final int x0max = clip((int) Math.ceil(x0 + RINGS[ring] * w0), width);
        final int x1max = clip((int) Math.ceil(x1 + RINGS[ring] * w1), height);

        final int x0range = x0max - x0min;
        final int x1range = x1max - x1min;

        // min should always be less than max
        if (x0range < 1)
            throw new IllegalArgumentException("Dimension 0 range not within data bounds");
        if (x1range < 1)
            throw new IllegalArgumentException("Dimension 1 range not within data bounds");

        // Shift centre to origin and compute gaussian
        final double[] gauss = airy2D(x0range, x1range, sum, x0 - x0min, x1 - x1min, w0, w1);

        return insert(data, x0min, x1min, x0max, x1max, width, gauss, poissonNoise);
    }

    /**
     * Construct a Airy pattern on the provided data. Only evaluates the function up to the configured dark ring.
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
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     * @return The data (packed in yx order, length = x0range * x1range)
     */
    public double[] airy2D(int x0range, int x1range, double sum, double x0, double x1, double w0, double w1)
    {
        w0 = Math.abs(w0);
        w1 = Math.abs(w1);

        this.w0 = w0;
        this.w1 = w1;

        // Limit to nth dark ring
        final double limit = RINGS[ring] * RINGS[ring];
        final double[] data = new double[x0range * x1range];

        // Store if the Airy pattern has been clipped
        final boolean clipped = (x0 - RINGS[ring] * w0 < 0) || (x1 - RINGS[ring] * w1 < 0) ||
                (x0 + RINGS[ring] * w0 > x0range) || (x1 + RINGS[ring] * w0 > x1range);

        // Offset by pixel centres by 0.5
        x0 -= 0.5;
        x1 -= 0.5;

        // Pre-compute the Airy intensity used for interpolation.
        // Find the maximum distance from the centre to the edge of the image (normalised using the widths)
        final double max = MathUtils.max(x0 / w0, x1 / w1, (x0range - x0) / w0, (x1range - x1) / w1);
        // Find the maximum distance needed to evaluate the Airy pattern
        final double maxD = FastMath.min(RINGS[ring], Math.sqrt(2 * max * max));
        // Limit the total samples used for interpolation but always sample at least every pixel:
        final double samplesPerPixel = FastMath.max(200 / maxD, 1);
        final int maxR = (int) Math.ceil(maxD * samplesPerPixel);
        final double[] radius = new double[maxR + 1];
        final double[] intensity = new double[maxR + 1];
        for (int r = 0; r <= maxR; r++)
        {
            // TODO - To simulate out of focus planes the intensity function can be pre-computed using
            // a different equation, e.g. Born-Wolf model.
            // Note that the pixel width for evaluation (e.g. the dark rings) would need to be calculated
            // and a different normalisation factor would have to be calculated for clipped data. This may be
            // achieved by pre-calculation of widths and normalisation factors for different z-depths.
            intensity[r] = AiryPattern.intensity(r / samplesPerPixel);
            radius[r] = r / samplesPerPixel;
        }

        double integral = 0;

        // Pre-calculate x offset
        final double[] d0 = new double[x0range];
        final double[] d02 = new double[x0range];
        for (int x = 0; x < x0range; x++)
        {
            d0[x] = (x - x0) / w0;
            d02[x] = d0[x] * d0[x];
        }

        if (singlePixelApproximation)
            // Single point approximation
            for (int y = 0, i = 0; y < x1range; y++)
            {
                double d1 = (y - x1) / w1;
                d1 *= d1;

                for (int x = 0; x < x0range; x++, i++)
                {
                    final double distance2 = d02[x] + d1;
                    if (distance2 < limit)
                    {
                        final double a = intensity(d02[x], d1, limit, samplesPerPixel, intensity, radius);
                        data[i] = a;
                        integral += a;
                    }
                }
            }
        else
        {
            // Integration using Simpson's composite interval

            // Set the number of subintervals adaptively, i.e. for small widths use more samples per pixel.
            final double nPixels = Math.PI * maxD * maxD;
            final double number = Math.sqrt(1000 / nPixels); // Approx 1000 (or so) samples across the image
            final int nSubintervals = FastMath.max(minSamplesPerDimension,
                    FastMath.min(maxSamplesPerDimension, (int) Math.ceil(number * 0.5) * 2));

            final double range0 = 0.5 / w0;
            final double range1 = 0.5 / w1;
            // Allow any point of the square pixel to be within the limit
            final double pixelLimit = limit + Math.sqrt(0.5);
            final double rescale = w0 * w1;

            for (int y = 0, i = 0; y < x1range; y++)
            {
                final double d1 = (y - x1) / w1;
                final double d12 = d1 * d1;

                for (int x = 0; x < x0range; x++, i++)
                {
                    final double distance2 = d02[x] + d12;
                    if (distance2 < pixelLimit)
                    {
                        final double a = integral(d0[x] - range0, d0[x] + range0, d1 - range1, d1 + range1, limit,
                                samplesPerPixel, intensity, radius, nSubintervals) * rescale;
                        data[i] = a;
                        integral += a;
                    }
                }
            }
        }

        //System.out.printf("Integral = %g (nPixels=%g) w0=%g, power = %g (norm = %g) = %g (clipped=%b)\n", integral,
        //		Math.PI * RINGS[ring] * w0 * RINGS[ring] * w1, w0, POWER[ring], integral / POWER[ring], (4 * Math.PI *
        //				w0 * w1), clipped);

        // We must normalise the integral we calculated to the correct power of the Airy pattern,
        // i.e. make the function we calculated a probability density that sums to 1.
        if (clipped)
            // Analysis has shown on unclipped data that the integral up to the nth ring is:
            // integral ~ POWER[ring] * (Math.PI * 4 * w0 * w1)
            // i.e. the full power of the Airy pattern is (Math.PI * 4 * w0 * w1)
            sum *= 1.0 / (4 * Math.PI * w0 * w1);
        else
            // The integral we calculated corresponds to the power at the nth ring
            sum *= POWER[ring] / integral;

        for (int i = 0; i < data.length; i++)
            data[i] *= sum;

        return data;
    }

    /**
     * Calculate the intensity of the Airy pattern at the given distances by interpolation using the lookup table.
     *
     * @param d0
     *            squared distance in dimension 0
     * @param d1
     *            squared distance in dimension 1
     * @param limit
     *            The squared distance limit of the Airy pattern
     * @param samplesPerPixel
     *            The number of samples per pixel of the pattern
     * @param intensity
     *            The Airy intensity at the provided radii
     * @param radius
     *            The radii
     * @return The intensity
     */
    private static double intensity(final double d0, final double d1, final double limit, final double samplesPerPixel,
            final double[] intensity, final double[] radius)
    {
        final double distance2 = d0 + d1;
        if (distance2 < limit)
        {
            final double r = Math.sqrt(distance2);
            //return AiryPattern.intensity(r);

            // Interpolate the intensity at this pixel
            final int index = (int) (r * samplesPerPixel);
            return intensity[index] + (intensity[index + 1] - intensity[index]) * (r - radius[index]) * samplesPerPixel;
        }
        return 0;
    }

    /**
     * Calculate the intensity of the Airy pattern between the specified ranges using the composite Simpson's rule
     *
     * @param ax
     *            Lower limit of x
     * @param bx
     *            Upper limit of x
     * @param ay
     *            Lower limit of y
     * @param by
     *            Upper limit of y
     * @param limit
     *            The squared distance limit of the Airy pattern
     * @param samplesPerPixel
     *            The number of samples per pixel of the pattern
     * @param intensity
     *            The Airy intensity at the provided radii
     * @param radius
     *            The radii
     * @param N
     *            The number of subintervals
     * @return the integral
     */
    private static double integral(final double ax, final double bx, final double ay, final double by,
            final double limit, final double samplesPerPixel, final double[] intensity, final double[] radius,
            final int N)
    {
        final double h = (bx - ax) / N;
        // TODO - The upper and lower bounds can be pre-computed since they are used for each pixel boundary
        double s = integral(ax * ax, ay, by, limit, samplesPerPixel, intensity, radius, N) +
                integral(bx * bx, ay, by, limit, samplesPerPixel, intensity, radius, N);
        for (int n = 1; n < N; n += 2)
        {
            final double x = ax + n * h;
            s += 4 * integral(x * x, ay, by, limit, samplesPerPixel, intensity, radius, N);
        }
        for (int n = 2; n < N; n += 2)
        {
            final double x = ax + n * h;
            s += 2 * integral(x * x, ay, by, limit, samplesPerPixel, intensity, radius, N);
        }
        return s * h / 3;
    }

    /**
     * Calculate the intensity of the Airy pattern between the specified ranges using the composite Simpson's rule
     *
     * @param x2
     *            The squared x distance
     * @param ay
     *            Lower limit of y
     * @param by
     *            Upper limit of y
     * @param limit
     *            The squared distance limit of the Airy pattern
     * @param samplesPerPixel
     *            The number of samples per pixel of the pattern
     * @param intensity
     *            The Airy intensity at the provided radii
     * @param radius
     *            The radii
     * @param N
     *            The number of subintervals
     * @return the integral
     */
    private static double integral(final double x2, final double ay, final double by, final double limit,
            final double samplesPerPixel, final double[] intensity, final double[] radius, final int N)
    {
        final double h = (by - ay) / N;
        // TODO - The upper and lower bounds can be pre-computed since they are used for each pixel boundary
        double s = intensity(x2, ay * ay, limit, samplesPerPixel, intensity, radius) +
                intensity(x2, by * by, limit, samplesPerPixel, intensity, radius);
        for (int n = 1; n < N; n += 2)
        {
            final double y = ay + n * h;
            s += 4 * intensity(x2, y * y, limit, samplesPerPixel, intensity, radius);
        }
        for (int n = 2; n < N; n += 2)
        {
            final double y = ay + n * h;
            s += 2 * intensity(x2, y * y, limit, samplesPerPixel, intensity, radius);
        }
        return s * h / 3;
    }

    private static int clip(int x, int max)
    {
        if (x < 0)
            x = 0;
        if (x > max)
            x = max;
        return x;
    }

    /**
     * @return the Z-depth where the 3D PSF is sqrt(2) the width (1.41 x FWHM)
     */
    public double getzDepth()
    {
        return zDepth;
    }

    /**
     * @param zDepth
     *            the Z-depth where the 3D PSF is sqrt(2) the width (1.41 x FWHM)
     */
    public void setzDepth(double zDepth)
    {
        this.zDepth = Math.abs(zDepth);
    }

    /**
     * @return The width in dimension 0 for the last drawn Airy pattern.
     */
    public double getW0()
    {
        return w0;
    }

    /**
     * @return The width in dimension 1 for the last drawn Airy pattern.
     */
    public double getW1()
    {
        return w1;
    }

    /** {@inheritDoc} */
    @Override
    public AiryPSFModel copy()
    {
        final AiryPSFModel model = new AiryPSFModel();
        model.zeroW0 = zeroW0;
        model.zeroW1 = zeroW1;
        model.zDepth = zDepth;
        model.ring = ring;
        //model.inputImage = inputImage;
        model.singlePixelApproximation = singlePixelApproximation;
        model.minSamplesPerDimension = minSamplesPerDimension;
        model.maxSamplesPerDimension = maxSamplesPerDimension;
        return model;
    }

    /**
     * Gets the ring.
     *
     * @return the ring limit for the calculated Airy pattern
     */
    public int getRing()
    {
        return ring;
    }

    /**
     * Set the limit of the Airy pattern, defined by the dark rings where the pattern is zero. Allowed values are 1-5.
     *
     * @param ring
     *            the ring limit for the calculated Airy pattern
     */
    public void setRing(int ring)
    {
        if (ring < RINGS.length && ring > 1)
            this.ring = ring;
    }

    /**
     * Checks if is single pixel approximation.
     *
     * @return True if the Airy pattern is evaluated once per pixel, otherwise use Simpson's integration
     */
    public boolean isSinglePixelApproximation()
    {
        return singlePixelApproximation;
    }

    /**
     * Sets the single pixel approximation.
     *
     * @param singlePixelApproximation
     *            True if the Airy pattern is evaluated once per pixel, otherwise use Simpson's integration
     */
    public void setSinglePixelApproximation(boolean singlePixelApproximation)
    {
        this.singlePixelApproximation = singlePixelApproximation;
    }

    /**
     * Gets the min samples per dimension.
     *
     * @return The minimum number of samples per dimension for Simpson's integration over each pixel
     */
    public int getMinSamplesPerDimension()
    {
        return minSamplesPerDimension;
    }

    /**
     * Set the minimum number of samples per dimension for Simpson's integration over each pixel. Must be above 0 and is
     * set to the next even number.
     *
     * @param n
     *            The minimum number of samples per dimension for Simpson's integration over each pixel
     */
    public void setMinSamplesPerDimension(int n)
    {
        if (n >= 2)
            this.minSamplesPerDimension = (n % 2 == 0) ? n : n + 1;
    }

    /**
     * Gets the max samples per dimension.
     *
     * @return The maximum number of samples per dimension for Simpson's integration over each pixel
     */
    public int getMaxSamplesPerDimension()
    {
        return maxSamplesPerDimension;
    }

    /**
     * Set the maximum number of samples per dimension for Simpson's integration over each pixel. Must be above 0 and is
     * set to the next even number.
     *
     * @param n
     *            The maximum number of samples per dimension for Simpson's integration over each pixel
     */
    public void setMaxSamplesPerDimension(int n)
    {
        if (n >= 2)
            this.maxSamplesPerDimension = (n % 2 == 0) ? n : n + 1;
    }

    @Override
    public int sample3D(float[] data, int width, int height, int n, double x0, double x1, double x2)
    {
        if (n <= 0)
            return insertSample(data, width, height, null, null);
        final double scale = createWidthScale(x2);
        final double[][] sample = sample(n, x0, x1, scale * zeroW0, scale * zeroW1);
        return insertSample(data, width, height, sample[0], sample[1]);
    }

    @Override
    public int sample3D(double[] data, int width, int height, int n, double x0, double x1, double x2)
    {
        if (n <= 0)
            return insertSample(data, width, height, null, null);
        final double scale = createWidthScale(x2);
        final double[][] sample = sample(n, x0, x1, scale * zeroW0, scale * zeroW1);
        return insertSample(data, width, height, sample[0], sample[1]);
    }

    /**
     * Sample from an Airy distribution.
     *
     * @param n
     *            The number of samples
     * @param x0
     *            The centre in dimension 0
     * @param x1
     *            The centre in dimension 1
     * @param w0
     *            The Airy width for dimension 0
     * @param w1
     *            The Airy width for dimension 1
     * @return The sample x and y values
     */
    public double[][] sample(final int n, final double x0, final double x1, final double w0, final double w1)
    {
        this.w0 = w0;
        this.w1 = w1;
        if (spline == null)
            createAiryDistribution();
        double[] x = new double[n];
        double[] y = new double[n];

        final RandomGenerator random = rand.getRandomGenerator();
        final UnitSphereRandomVectorGenerator vg = new UnitSphereRandomVectorGenerator(2, random);
        int c = 0;
        for (int i = 0; i < n; i++)
        {
            final double p = random.nextDouble();
            if (p > POWER[SAMPLE_RINGS])
                // TODO - We could add a simple interpolation here using a spline from AiryPattern.power()
                continue;
            final double r = spline.value(p);

            // Convert to xy using a random vector generator
            final double[] v = vg.nextVector();
            x[c] = v[0] * r * w0 + x0;
            y[c] = v[1] * r * w1 + x1;
            c++;
        }

        if (c < n)
        {
            x = Arrays.copyOf(x, c);
            y = Arrays.copyOf(y, c);
        }
        return new double[][] { x, y };
    }

    private static synchronized void createAiryDistribution()
    {
        if (spline != null)
            return;

        final double relativeAccuracy = 1e-4;
        final double absoluteAccuracy = 1e-8;
        final int minimalIterationCount = 3;
        final int maximalIterationCount = 32;

        final UnivariateIntegrator integrator = new SimpsonIntegrator(relativeAccuracy, absoluteAccuracy,
                minimalIterationCount, maximalIterationCount);
        final UnivariateFunction f = new UnivariateFunction()
        {
            @Override
            public double value(double x)
            {
                // The pattern profile is in one dimension.
                // Multiply by the perimeter of a circle to convert to 2D volume then normalise by 4 pi
                //return AiryPattern.intensity(x) * 2 * Math.PI * x / (4 * Math.PI);
                return AiryPattern.intensity(x) * 0.5 * x;
            }
        };

        // Integrate up to a set number of dark rings
        final int samples = 1000;
        final double step = RINGS[SAMPLE_RINGS] / samples;
        double to = 0, from = 0;
        r = new double[samples + 1];
        sum = new double[samples + 1];
        for (int i = 1; i < sum.length; i++)
        {
            from = to;
            r[i] = to = step * i;
            sum[i] = integrator.integrate(2000, f, from, to) + sum[i - 1];
        }

        if (DoubleEquality.relativeError(sum[samples], POWER[SAMPLE_RINGS]) > 1e-3)
            throw new RuntimeException("Failed to create the Airy distribution");

        final SplineInterpolator si = new SplineInterpolator();
        spline = si.interpolate(sum, r);
    }

    @Override
    protected boolean computeValueAndGradient(int width, int height, double x0, double x1, double x2, double[] value,
            double[][] jacobian)
    {
        final double[] dx = new double[] { 1e-4, 1e-4, 1e-4 };
        return computeValueAndGradient(width, height, x0, x1, x2, value, jacobian, dx);
    }
}
