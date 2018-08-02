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

import ij.process.ImageProcessor;

/**
 * Store a 2D image in a single double array. Forms a base for 2D DHT transform using the JTransforms library.
 */
public class DoubleImage2D extends Image2D
{
    /**
     * The data packed in yx order: y * nc + x .
     * 
     * @see #index(int, int)
     */
    protected double[] data;

    /**
     * Instantiates a new 2D image.
     *
     * @param nc
     *            the number of columns
     * @param nr
     *            the number of rows
     * @throws IllegalArgumentException
     *             If the combined dimensions is too large for an array
     */
    public DoubleImage2D(int nc, int nr) throws IllegalArgumentException
    {
        super(nc, nr);
    }

    /**
     * Instantiates a new 2D image.
     *
     * @param image
     *            the image
     * @throws IllegalArgumentException
     *             If the combined dimensions is too large for an array
     */
    public DoubleImage2D(ImageProcessor image) throws IllegalArgumentException
    {
        super(image);
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.ij.utils.Image2D#createData(int)
     */
    @Override
    protected void createData(int size)
    {
        data = new double[size];
    }

    /**
     * Instantiates a new 2D image.
     *
     * @param nc
     *            the number of columns
     * @param nr
     *            the number of rows
     * @param data
     *            the data
     * @throws IllegalArgumentException
     *             If the data is not the correct length
     */
    public DoubleImage2D(int nc, int nr, double[] data) throws IllegalArgumentException
    {
        // Avoid constructor that calls createData(int)
        super(nc, nr, false);
        if (data == null || data.length != checkSize(nc, nr, true))
            throw new IllegalArgumentException("Data is not correct length");
        this.data = data;
    }

    /**
     * Instantiates a new 2D image.
     *
     * @param nc
     *            the number of columns
     * @param nr
     *            the number of rows
     * @param data
     *            the data
     * @param dummy
     *            the dummy flag
     */
    protected DoubleImage2D(int nc, int nr, double[] data, boolean dummy)
    {
        // No checks as this is used internally
        super(nc, nr, false);
        this.data = data;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.ij.utils.Image2D#copy()
     */
    @Override
    public DoubleImage2D copy()
    {
        return new DoubleImage2D(nc, nr, data.clone(), false);
    }

    /**
     * Gets the data.
     *
     * @return the data
     */
    public double[] getData()
    {
        return data;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.ij.utils.Image2D#getDataLength()
     */
    @Override
    public int getDataLength()
    {
        return data.length;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.ij.utils.Image2D#crop(int, int, int, int)
     */
    @Override
    public DoubleImage2D crop(int x, int y, int w, int h) throws IllegalArgumentException
    {
        return crop(x, y, w, h, null);
    }

    /**
     * Crop a sub-region of the data. The target dimensions must be positive.
     *
     * @param x
     *            the x index
     * @param y
     *            the y index
     * @param w
     *            the width
     * @param h
     *            the height
     * @param region
     *            the cropped data (will be reused if the correct size)
     * @return the cropped data
     * @throws IllegalArgumentException
     *             if the region is not within the data
     */
    public DoubleImage2D crop(int x, int y, int w, int h, double[] region) throws IllegalArgumentException
    {
        // Check the region range
        if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr)
            throw new IllegalArgumentException("Region not within the data");
        final int size = h * w;
        if (region == null || region.length != size)
            region = new double[size];
        int base = y * nc + x;
        for (int r = 0, i = 0; r < h; r++)
        {
            System.arraycopy(data, base, region, i, w);
            base += nc;
            i += w;
        }
        return new DoubleImage2D(w, h, region, false);
    }

    /**
     * Crop a sub-region of the data. The target dimensions must be positive.
     *
     * @param image
     *            the image
     * @param x
     *            the x index
     * @param y
     *            the y index
     * @param w
     *            the width
     * @param h
     *            the height
     * @param region
     *            the cropped data (will be reused if the correct size)
     * @return the cropped data
     * @throws IllegalArgumentException
     *             if the region is not within the data
     */
    public static DoubleImage2D crop(ImageProcessor image, int x, int y, int w, int h, double[] region)
            throws IllegalArgumentException
    {
        final int nc = image.getWidth();
        final int nr = image.getHeight();

        // Check the region range
        if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr)
            throw new IllegalArgumentException("Region not within the data");
        final int size = checkSize(w, h, true);
        if (region == null || region.length != size)
            region = new double[size];
        int base = y * nc + x;
        for (int r = 0, i = 0; r < h; r++)
        {
            for (int c = 0; c < w; c++)
                region[i++] = image.getf(base + c);
            base += nc;
        }
        return new DoubleImage2D(w, h, region, false);
    }

    @Override
    public void insert(int x, int y, Image2D image) throws IllegalArgumentException
    {
        if (image instanceof DoubleImage2D)
            insert(x, y, (DoubleImage2D) image);
        else
            super.insert(x, y, image);
    }

    /**
     * Insert a sub-region.
     *
     * @param x
     *            the x position
     * @param y
     *            the y position
     * @param image
     *            the image
     * @throws IllegalArgumentException
     *             if the region is not within the data
     */
    public void insert(int x, int y, DoubleImage2D image) throws IllegalArgumentException
    {
        // Check the region range
        final int w = image.getWidth();
        final int h = image.getHeight();
        if (w < 1 || h < 1)
            return;
        if (x < 0 || (long) x + w > nc || y < 0 || (long) y + h > nr)
            throw new IllegalArgumentException("Region not within the data");
        final double[] region = image.data;
        int base = y * nc + x;
        for (int r = 0, i = 0; r < h; r++)
        {
            System.arraycopy(region, i, data, base, w);
            base += nc;
            i += w;
        }
    }

    @Override
    protected void copyTo(int i, float[] buffer, int j, int size)
    {
        while (size-- > 0)
            buffer[j++] = (float) data[i++];
    }

    @Override
    protected void copyFrom(float[] buffer, int j, int size, int i)
    {
        while (size-- > 0)
            data[i++] = buffer[j++];
    }

    @Override
    public double get(int i)
    {
        return data[i];
    }

    @Override
    public void set(int i, double value)
    {
        data[i] = value;
    }

    @Override
    public float getf(int i)
    {
        return (float) data[i];
    }

    @Override
    public void setf(int i, float value)
    {
        data[i] = value;
    }

    @Override
    protected void fill(int i, int size, double value)
    {
        while (size-- > 0)
            data[i++] = value;
    }
}
