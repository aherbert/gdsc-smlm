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
package uk.ac.sussex.gdsc.smlm.utils;

import java.awt.Rectangle;

/**
 * Contains methods for converting an image to float data.
 * <p>
 * Handles unsigned byte, unsigned short, float, double and RGB int array data.
 */
public class ImageConverter
{
    private double rWeight = 1d / 3d, gWeight = 1d / 3d, bWeight = 1d / 3d;

    /**
     * Sets the weighting factors used to do colour conversions. The default values are
     * 1/3, 1/3 and 1/3. E.g. for weighted RGB Conversions use 0.299, 0.587 and 0.114.
     *
     * @param rFactor
     *            the r factor
     * @param gFactor
     *            the g factor
     * @param bFactor
     *            the b factor
     * @throws IllegalArgumentException
     *             if the factors do not sum to a positive value
     */
    public void setWeightingFactors(double rFactor, double gFactor, double bFactor) throws IllegalArgumentException
    {
        if (!(rFactor >= 0 && gFactor >= 0 && bFactor >= 0))
            throw new IllegalArgumentException("Weights must sum to a positive value");

        final double sum = rFactor + gFactor + bFactor;
        if (!(sum > 0 && sum != Double.POSITIVE_INFINITY))
            throw new IllegalArgumentException("Weights must sum to a positive finite value");

        rWeight = rFactor / sum;
        gWeight = gFactor / sum;
        bWeight = bFactor / sum;
    }

    /**
     * Returns the three weighting factors used to do colour conversions.
     *
     * @return the weighting factors
     */
    public double[] getWeightingFactors()
    {
        final double[] weights = new double[3];
        weights[0] = rWeight;
        weights[1] = gWeight;
        weights[2] = bWeight;
        return weights;
    }

    /**
     * Converts the specified RGB pixel to greyscale using the formula g=(r+g+b)/3 and returns it as a float.
     * Call {@link #setWeightingFactors(double, double, double)} to specify different conversion factors.
     *
     * @param c
     *            the RGB value
     * @return the greyscale value
     */
    public float rgbToGreyscale(int c)
    {
        final int r = (c & 0xff0000) >> 16;
        final int g = (c & 0xff00) >> 8;
        final int b = c & 0xff;
        return (float) (r * rWeight + g * gWeight + b * bWeight);
    }

    /**
     * Get the data from the image as a float array (include cropping to the ROI). Data is duplicated if the
     * input is already a float array.
     * <p>
     * Allows reuse of an existing buffer if provided. This will not be truncated if it is larger than the
     * ImageProcessor ROI bounds. If smaller then a new buffer will be created.
     * <p>
     * If the object pixels array is incorrect size (it should be width*height) then null will be returned.
     *
     * @param oPixels
     *            the pixels object
     * @param width
     *            the width
     * @param height
     *            the height
     * @param bounds
     *            the bounds
     * @param buffer
     *            the buffer
     * @return The float array data
     */
    public float[] getData(final Object oPixels, final int width, final int height, final Rectangle bounds,
            float[] buffer)
    {
        if (oPixels == null)
            return null;

        // Ignore the bounds if they specify the entire image size

        if (oPixels instanceof float[])
        {
            final float[] pixels = (float[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = pixels[offset2++];
                return pixels2;
            }
            final float[] pixels2 = allocate(buffer, pixels.length);
            System.arraycopy(pixels, 0, pixels2, 0, pixels.length);
            return pixels2;
        }
        else if (oPixels instanceof short[])
        {
            final short[] pixels = (short[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = pixels[offset2++] & 0xffff;
                return pixels2;
            }
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = pixels[i] & 0xffff;
            return pixels2;
        }
        else if (oPixels instanceof byte[])
        {
            final byte[] pixels = (byte[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = pixels[offset2++] & 0xff;
                return pixels2;
            }
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = pixels[i] & 0xff;
            return pixels2;
        }
        else if (oPixels instanceof double[])
        {
            final double[] pixels = (double[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = (float) pixels[offset2++];
                return pixels2;
            }
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = (float) pixels[i];
            return pixels2;
        }
        else if (oPixels instanceof int[])
        {
            // The default processing assumes RGB
            final int[] pixels = (int[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = rgbToGreyscale(pixels[offset2++]);
                return pixels2;
            }
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = rgbToGreyscale(pixels[i]);
            return pixels2;
        }
        return null;
    }

    private static boolean incorrectSize(int length, int width, int height)
    {
        return length != width * height;
    }

    private static float[] allocate(float[] buffer, int size)
    {
        if (buffer == null || buffer.length < size)
            buffer = new float[size];
        return buffer;
    }

    /**
     * Get the data from the image as a double array (include cropping to the ROI). Data is duplicated if the
     * input is already a double array.
     * <p>
     * Allows reuse of an existing buffer if provided. This will not be truncated if it is larger than the
     * ImageProcessor ROI bounds. If smaller then a new buffer will be created.
     * <p>
     * If the object pixels array is incorrect size (it should be width*height) then null will be returned.
     *
     * @param oPixels
     *            the pixels object
     * @param width
     *            the width
     * @param height
     *            the height
     * @param bounds
     *            the bounds
     * @param buffer
     *            the buffer
     * @return The double array data
     */
    public double[] getDoubleData(final Object oPixels, final int width, final int height, final Rectangle bounds,
            double[] buffer)
    {
        if (oPixels == null)
            return null;

        // Ignore the bounds if they specify the entire image size

        if (oPixels instanceof float[])
        {
            final float[] pixels = (float[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = pixels[offset2++];
                return pixels2;
            }
            final double[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = pixels[i];
            return pixels2;
        }
        else if (oPixels instanceof short[])
        {
            final short[] pixels = (short[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = pixels[offset2++] & 0xffff;
                return pixels2;
            }
            final double[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = pixels[i] & 0xffff;
            return pixels2;
        }
        else if (oPixels instanceof byte[])
        {
            final byte[] pixels = (byte[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = pixels[offset2++] & 0xff;
                return pixels2;
            }
            final double[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = pixels[i] & 0xff;
            return pixels2;
        }
        if (oPixels instanceof double[])
        {
            final double[] pixels = (double[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = pixels[offset2++];
                return pixels2;
            }
            final double[] pixels2 = allocate(buffer, pixels.length);
            System.arraycopy(pixels, 0, pixels2, 0, pixels.length);
            return pixels2;
        }
        else if (oPixels instanceof int[])
        {
            // The default processing assumes RGB
            final int[] pixels = (int[]) oPixels;
            if (incorrectSize(pixels.length, width, height))
                return null;
            if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
            {
                final double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
                for (int ys = 0, offset1 = 0; ys < bounds.height; ys++)
                    for (int xs = 0, offset2 = (ys + bounds.y) * width + bounds.x; xs < bounds.width; xs++)
                        pixels2[offset1++] = rgbToGreyscale(pixels[offset2++]);
                return pixels2;
            }
            final double[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = rgbToGreyscale(pixels[i]);
            return pixels2;
        }
        return null;
    }

    private static double[] allocate(double[] buffer, int size)
    {
        if (buffer == null || buffer.length < size)
            buffer = new double[size];
        return buffer;
    }

    /**
     * Get the data from the image pixels as a float array. Data is not duplicated if the
     * input is already a float array unless a buffer is provided.
     * <p>
     * Allows reuse of an existing buffer if provided. This will not be truncated if it is larger than the
     * pixels array. If smaller then a new buffer will be created.
     *
     * @param oPixels
     *            the pixels
     * @param buffer
     *            the buffer
     * @return The float array data
     */
    public float[] getData(final Object oPixels, float[] buffer)
    {
        if (oPixels == null)
            return null;

        if (oPixels instanceof float[])
        {
            final float[] pixels = (float[]) oPixels;
            if (buffer == null)
                return pixels;
            final float[] pixels2 = allocate(buffer, pixels.length);
            System.arraycopy(pixels, 0, pixels2, 0, pixels.length);
            return pixels2;
        }
        else if (oPixels instanceof short[])
        {
            final short[] pixels = (short[]) oPixels;
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = pixels[i] & 0xffff;
            return pixels2;
        }
        else if (oPixels instanceof byte[])
        {
            final byte[] pixels = (byte[]) oPixels;
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = pixels[i] & 0xff;
            return pixels2;
        }
        else if (oPixels instanceof double[])
        {
            final double[] pixels = (double[]) oPixels;
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = (float) pixels[i];
            return pixels2;
        }
        else if (oPixels instanceof int[])
        {
            // The default processing
            final int[] pixels = (int[]) oPixels;
            final float[] pixels2 = allocate(buffer, pixels.length);
            for (int i = 0; i < pixels.length; i++)
                pixels2[i] = rgbToGreyscale(pixels[i]);
            return pixels2;
        }
        return null;
    }

    /**
     * Checks if the primitive array type is supported.
     *
     * @param oPixels
     *            the pixels
     * @return true, if is supported
     */
    public boolean isSupported(Object oPixels)
    {
        if (oPixels instanceof float[])
            return true;
        if (oPixels instanceof short[])
            return true;
        if (oPixels instanceof byte[])
            return true;
        if (oPixels instanceof double[])
            return true;
        if (oPixels instanceof int[])
            return true;
        return false;
    }
}
