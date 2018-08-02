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
package uk.ac.sussex.gdsc.smlm.ij.utils;

import java.util.Arrays;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.math.interpolation.CachedBicubicInterpolator;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.ImageWindow;
import uk.ac.sussex.gdsc.core.utils.Maths;

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
                input = resize(input);
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
    private boolean fastMultiply = true;

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
        if (fastMultiply)
            dhtData.dht.initialiseFastMultiply();
    }

    /**
     * Check the image is 2D and has data.
     *
     * @param image
     *            the image
     */
    private static void check2D(ImageProcessor image)
    {
        if (image.getWidth() < 2 || image.getHeight() < 2)
            throw new IllegalArgumentException("Require a 2D image");
        // Check for data
        final int size = image.getWidth() * image.getHeight();
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
        final int w = image.getWidth(), h = image.getHeight();
        final double[] wx = createXWindow(w);
        final double[] wy = createYWindow(h);

        // We need to compute the weighted centre
        final double[] sum = new double[2];

        final float[] pixels = (float[]) image.getPixels();
        calculateWeightedCentre(pixels, w, h, wx, wy, sum);

        final double shift = sum[0] / sum[1];

        applyWindow(pixels, w, h, wx, wy, shift);

        DoubleDHT2D dht;

        // Pad into the desired data size.
        // We always do this so the data is reused
        double[] dest;
        if (dhtData == null || dhtData.dht.getDataLength() != nr_by_nc)
            dest = new double[nr_by_nc];
        else
        {
            // Re-use space
            dest = dhtData.dht.getData();
            Arrays.fill(dest, 0);
        }
        dht = new DoubleDHT2D(nc, nr, dest, false);
        final int ix = getInsert(nc, w);
        final int iy = getInsert(nr, h);
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
            for (int x = 0; x < maxx; x++, i++)
            {
                final double w = wx[x] * wy[y];
                sum[0] += image[i] * w;
                sum[1] += w;
            }
    }

    private static void calculateWeightedCentre(Image2D image, int maxx, int maxy, double[] wx, double[] wy,
            double[] sum)
    {
        for (int y = 0, i = 0; y < maxy; y++)
            for (int x = 0; x < maxx; x++, i++)
            {
                final double w = wx[x] * wy[y];
                sum[0] += image.get(i) * w;
                sum[1] += w;
            }
    }

    private static void applyWindow(float[] image, int maxx, int maxy, double[] wx, double[] wy, double shift)
    {
        for (int y = 0, i = 0; y < maxy; y++)
            for (int x = 0; x < maxx; x++, i++)
                image[i] = (float) ((image[i] - shift) * wx[x] * wy[y]);
    }

    private static void applyWindow(Image2D image, int maxx, int maxy, double[] wx, double[] wy, double shift)
    {
        for (int y = 0, i = 0; y < maxy; y++)
            for (int x = 0; x < maxx; x++, i++)
                image.set(i, (image.get(i) - shift) * wx[x] * wy[y]);
    }

    private static int getInsert(int maxN, int n)
    {
        // Note the FHT power spectrum centre is at n/2 of an even sized image.
        // So we must insert the centre at that point. To do this we check for odd/even
        // and offset if necessary.
        final int diff = maxN - n;
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
    private static DHTData prepareDHT(DHTData dhtData)
    {
        final DoubleDHT2D dht = dhtData.dht;
        final double[] s_ = dhtData.s_;
        final double[] ss = dhtData.ss;

        // Note previous versions converted to 10-bit integer data. However the 2D DHT creates very large
        // output values and errors occurred when computing the conjugate multiple in the frequency
        // domain verses the spatial domain. A check has been added to compute the spatial domain
        // correlation for the corresponding max correlation in the frequency domain. This allow
        // the code to report when the correlation value is incorrect.

        final double[] data = dht.getData();
        final double[] limits = Maths.limits(data);
        final double min = limits[0];
        final double max = limits[1];

        // Note: The image has been shifted to a mean of 0 so that zero padding
        // for frequency domain transform does not add any information.
        // We need to maintain the sign information and ensure that zero is still
        // zero.

        final double scale = LIMIT / (max - min);

        // Compute the rolling sum tables
        final int nc = dht.nc;
        final int nr = dht.nr;

        // This has been adapted from Image2D to compute two rolling sum table at once

        double sum_ = 0, sum2 = 0;
        int i = 0;
        // Initialise first row sum
        // sum = rolling sum of (0 - colomn)
        for (int c = 0; c < nc; c++, i++)
        {
            final double v = transform(data[i], scale);
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
                final double v = transform(data[i], scale);
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
        final double value = f * scale;
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
    private static void check2D(Image2D image)
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
        final int w = image.getWidth(), h = image.getHeight();
        final double[] wx = createXWindow(w);
        final double[] wy = createYWindow(h);

        // We need to compute the weighted centre
        final double[] sum = new double[2];

        calculateWeightedCentre(image, w, h, wx, wy, sum);

        final double shift = sum[0] / sum[1];

        applyWindow(image, w, h, wx, wy, shift);

        //System.out.printf("Sum = %g => %g\n", sum[0], Maths.sum(pixels));

        DoubleDHT2D dht;

        // Pad into the desired data size.
        // We always do this to handle input of float/double Image2D data.
        double[] dest;
        if (dhtData == null || dhtData.dht.getDataLength() != nr_by_nc)
            dest = new double[nr_by_nc];
        else
        {
            // Re-use space
            dest = dhtData.dht.getData();
            Arrays.fill(dest, 0);
        }
        dht = new DoubleDHT2D(nc, nr, dest, false);
        final int ix = getInsert(nc, w);
        final int iy = getInsert(nr, h);
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
        final int w = image.getWidth(), h = image.getHeight();
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
        final int w = image.getWidth(), h = image.getHeight();
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
        final DoubleDHT2D correlation = target.dht.conjugateMultiply(reference.dht, buffer);
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
            final double f = (1 - minimumDimensionOverlap) / 2;
            final int ux = (int) (Math.round(Math.min(reference.w, target.w) * f));
            final int uy = (int) (Math.round(Math.min(reference.h, target.h) * f));
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
        final int nc_2 = nc / 2;
        final int nr_2 = nr / 2;
        final int[] centre = new int[] { nc_2, nr_2 };

        // Compute the shift from the centre
        final int dx = nc_2 - ix;
        final int dy = nr_2 - iy;

        // For the reference (moved -dx,-dy over the target)
        int rx = -dx;
        int ry = -dy;

        // For the target (moved dx,dy over the reference)
        int tx = dx;
        int ty = dy;

        // Precompute the x-1,x+w-1
        final int nx = crop[2];
        final int[] rx_1 = new int[nx];
        final int[] rx_w_1 = new int[nx];
        final int[] tx_1 = new int[nx];
        final int[] tx_w_1 = new int[nx];
        final int[] w = new int[nx];
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

        final double[] rs_ = reference.s_;
        final double[] rss = reference.ss;
        final double[] ts_ = target.s_;
        final double[] tss = target.ss;
        final double[] rsum = new double[2];
        final double[] tsum = new double[2];

        final int size = Math.min(reference.size, target.size);
        final int minimumN = (int) (Math.round(size * minimumOverlap));
        int maxj = -1;
        double max = 0;

        for (int r = iy; r < iyh; r++)
        {
            // Compute the y-1,y+h-1
            final int ry_1 = Math.max(-1, ry - 1);
            final int ry_h_1 = Math.min(nr, ry + nr) - 1;
            ry++;
            final int ty_1 = Math.max(-1, ty - 1);
            final int ty_h_1 = Math.min(nr, ty + nr) - 1;
            ty--;
            final int h = ry_h_1 - ry_1;

            final int base = r * nc;
            for (int c = ix, i = 0; c < ixw; c++, i++)
            {
                final double sumXY = buffer[base + c];

                compute(rx_1[i], ry_1, rx_w_1[i], ry_h_1, w[i], h, rs_, rss, rsum);
                compute(tx_1[i], ty_1, tx_w_1[i], ty_h_1, w[i], h, ts_, tss, tsum);

                // Compute the correlation
                // (sumXY - sumX*sumY/n) / sqrt( (sumXX - sumX^2 / n) * (sumYY - sumY^2 / n) )

                final int n = w[i] * h;
                final double numerator = sumXY - (rsum[X] * tsum[Y] / n);
                final double denominator1 = rsum[XX] - (rsum[X] * rsum[X] / n);
                final double denominator2 = tsum[YY] - (tsum[Y] * tsum[Y] / n);

                double R;
                if (denominator1 == 0 || denominator2 == 0)
                {
                    // If there is data and all the variances are the same then correlation is perfect
                    if (rsum[XX] == tsum[YY] && rsum[XX] == sumXY && rsum[XX] > 0)
                        R = 1;
                    else
                        R = 0;
                }
                else
                    R = numerator / Math.sqrt(denominator1 * denominator2);
                // Leave as raw for debugging
                //R = Maths.clip(-1, 1, R);

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
                    final int[] xyz1 = correlation.getXY(maxj);
                    final int[] xyz2 = correlation.getXY(base + c);
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
		final double[] result = new double[] {
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
            final int x = Maths.clip(ix, ixw - 5, xy[0] - 2);
            final int y = Maths.clip(iy, iyh - 5, xy[1] - 2);
            final DoubleImage2D crop = correlation.crop(x, y, 5, 5, region);

            // Find the maximum starting at the current origin
            final int ox = xy[0] - x;
            final int oy = xy[1] - y;

            double[] optimum;
            if (ox == 2 && oy == 2 && crop.getWidth() == 5 && crop.getHeight() == 5)
                optimum = performCubicSearch(crop, refinements, getRelativeThreshold());
            else
            {
                final FloatProcessor fp = new FloatProcessor(5, 5, crop.getData());
                optimum = performCubicFit(fp, ox, oy, refinements, getRelativeThreshold());
            }

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
        final int[] xy = correlation.getXY(maxi);

        // Find the range for the target and reference
        final int nc_2 = nc / 2;
        final int nr_2 = nr / 2;
        final int tx = Math.max(0, xy[0] - nc_2);
        final int ty = Math.max(0, xy[1] - nr_2);
        final int w = Math.min(nc, xy[0] + nc_2) - tx;
        final int h = Math.min(nr, xy[1] + nr_2) - ty;

        // For the reference we express as a shift relative to the centre
        // and subtract the half-width.
        // Formally: (nc_2 - xy[0]) // shift
        //           + nc_2          // centre
        //           - nc_2          // Half width
        final int rx = Math.max(0, -xy[0] + nc_2);
        final int ry = Math.max(0, -xy[1] + nr_2);

        final double[] tar = target.input;
        final double[] ref = reference.input;
        final double o = correlation.get(maxi);
        double e = 0;
        for (int y = 0; y < h; y++)
        {
            int i = (ty + y) * nc + tx;
            int j = (ry + y) * nc + rx;
            for (int x = 0; x < w; x++)
                e += tar[i++] * ref[j++];
        }

        //System.out.printf("Raw %d,%d = %g\n", xy[0], xy[1], o);

        frequencyDomainCorrelationError = DoubleEquality.relativeError(o, e);
        if (frequencyDomainCorrelationError > 0.05)
            System.err.printf("2D Correlation Error = %s : Spatial = %s, Freq = %s\n",
                    Utils.rounded(frequencyDomainCorrelationError), Double.toString(e), Double.toString(o));
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
        final int xw_yh = y_h_1 * nc + x_w_1;
        sum[0] = 0;
        sum[1] = 0;
        if (y_1 >= 0)
        {
            final int h_ = h * nc;
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
        final double[] centre = new double[] { i, j, fp.getf(i, j) };
        // This value will be progressively halved.
        // Start with a value that allows the number of iterations to fully cover the region +/- 1 pixel
        // 0.5 will result in an minimum range of 0.5 / 2^9 = 0.000976
        double range = 0.5;
        while (refinements-- > 0)
        {
            final double previous = centre[2];
            if (performCubicFit(fp, range, centre))
                // The centre moved. Check convergence.
                if ((centre[2] - previous) / centre[2] < relativeThreshold)
                    break;
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
     * @return true, if the centre moved
     */
    private static boolean performCubicFit(FloatProcessor fp, double range, double[] centre)
    {
        boolean moved = false;
        for (int i = -1; i <= 1; i++)
        {
            final double x = centre[0] + i * range;
            for (int j = -1; j <= 1; j++)
            {
                if (i == 0 && j == 0)
                    // Current maximum
                    continue;
                final double y = centre[1] + j * range;
                final double v = fp.getBicubicInterpolatedPixel(x, y, fp);
                if (centre[2] < v)
                {
                    centre[0] = x;
                    centre[1] = y;
                    centre[2] = v;
                    moved = true;
                }
            }
        }
        return moved;
    }

    /**
     * Iteratively search the cubic spline surface around the centre to maximise the value.
     * <p>
     * At each round each of 8 points around the current maximum (+/- range) are evaluated. The optimum is picked and
     * the range is halved. The initial range is 0.5 so the maximum distance that can be walked in any direction is 1
     * pixel when the number of refinements is unlimited. With refinements = 3 the distance is 0.5 + 0.25 + 0.125 =
     * 0.875.
     *
     * @param surface
     *            A peak surface (must be 5x5)
     * @param refinements
     *            the maximum number of refinements
     * @param relativeThreshold
     *            Sets the relative threshold for change in the value for halting refinement. This applies
     *            only if the position moved during the refinement step.
     * @return The peak location with sub-pixel accuracy [x,y,value]
     */
    public static double[] performCubicSearch(Image2D surface, int refinements, double relativeThreshold)
    {
        if (surface.getWidth() != 5 || surface.getHeight() != 5)
            throw new IllegalArgumentException("Require a 5x5 input surface");

        final CachedBicubicInterpolator[][] nodes = new CachedBicubicInterpolator[2][2];
        final double[] data = new double[16];
        for (int x = 0; x < 2; x++)
            for (int y = 0; y < 2; y++)
            {
                int offset = y * 5 + x;
                for (int k = 0, index = 0; k < 4; k++)
                {
                    for (int l = 0; l < 4; l++)
                        data[index++] = surface.get(offset + l);
                    offset += 5;
                }
                nodes[x][y] = new CachedBicubicInterpolator();
                nodes[x][y].updateCoefficients(data);
            }

        // Offset centre by 1 so it is exactly in the middle of the 2x2 grid of bicubic interpolators
        final double[] centre = new double[] { 1, 1, surface.get(12) };
        final double[] y = new double[9];
        final int[] iy = new int[3];
        // This value will be progressively halved.
        // Start with a value that allows the number of iterations to fully cover the region +/- 1 pixel
        // 0.5 will result in an minimum range of 0.5 / 2^9 = 0.000976
        double range = 0.5;
        while (refinements-- > 0)
        {
            final double previous = centre[2];
            if (performCubicSearch(surface, nodes, range, centre, y, iy))
                // The centre moved. Check convergence.
                if ((centre[2] - previous) / centre[2] < relativeThreshold)
                    break;
            range /= 2;
        }
        // Add back the pixel offset
        centre[0] += 1;
        centre[1] += 1;
        return centre;
    }

    /**
     * Perform a cubic search refinement.
     *
     * @param surface
     *            the surface (to allow debugging)
     * @param nodes
     *            the nodes for interpolation
     * @param range
     *            the range
     * @param centre
     *            the centre
     * @param y
     *            working space for y
     * @param iy
     *            working space for iy
     * @return true, if the centre moved
     */
    private static boolean performCubicSearch(Image2D surface, CachedBicubicInterpolator[][] nodes, double range,
            double[] centre, double[] y, int[] iy)
    {
        // for debugging
        //FloatProcessor fp = (FloatProcessor) surface.getImageProcessor();

        // Pre-compute the node position and the fraction between 0-1 for y values
        for (int j = -1, k = 0, l = 0; j <= 1; j++, k++)
        {
            double yy = centre[1] + j * range;
            iy[k] = (int) yy;
            yy = yy - iy[k];
            y[l++] = yy;
            y[l++] = yy * yy;
            y[l++] = yy * yy * yy;
        }

        boolean moved = false;
        for (int i = -1; i <= 1; i++)
        {
            // Compute the node position and the fraction between 0-1 for y values
            double x = centre[0] + i * range;
            final int ix = (int) x;
            x = x - ix;
            final double x2 = x * x;
            final double x3 = x * x2;
            for (int k = 0; k < 3; k++)
            {
                if (i == 0 && k == 1)
                    // Current maximum
                    continue;
                final double v = nodes[ix][iy[k]].getValue(x, x2, x3, y[k * 3], y[k * 3 + 1], y[k * 3 + 2]);
                //double v2 = fp.getBicubicInterpolatedPixel(ix + 1 + x, iy[k] + 1 + y[k * 3], fp);
                //System.out.printf("%g vs %g @ %g,%g\n", v, v2, x + ix, iy[k] + y[k * 3]);
                if (centre[2] < v)
                {
                    // Add back the node index to get the correct x/y-values
                    centre[0] = x + ix;
                    centre[1] = y[k * 3] + iy[k];
                    centre[2] = v;
                    moved = true;
                }
            }
        }
        //System.out.printf("Centre %g,%g = %g\n", centre[0], centre[1], centre[2]);
        return moved;
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
        catch (final CloneNotSupportedException e)
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
            final DoubleImage2D image = new DoubleImage2D(nc, nr, buffer);
            image.fillOutside(crop[0], crop[1], crop[2], crop[3], 0);
            return image;
        }
        catch (final IllegalArgumentException e)
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
