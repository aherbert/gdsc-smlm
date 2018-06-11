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
package gdsc.smlm.ij.utils;

import ij.ImageStack;
import ij.process.ImageProcessor;


/**
 * Store a 3D image in a single float array. Forms a base for 3D DHT transform using the JTransforms library.
 */
public class FloatImage3D extends Image3D
{
	protected float[] data;

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public FloatImage3D(int nc, int nr, int ns) throws IllegalArgumentException
	{
		super(nc, nr, ns);
	}

	/**
	 * Instantiates a new 3D image
	 *
	 * @param stack
	 *            the stack
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public FloatImage3D(ImageStack stack) throws IllegalArgumentException
	{
		super(stack);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#createData(int)
	 */
	@Override
	protected void createData(int size)
	{
		data = new float[size];
	}

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param stack
	 *            the stack
	 * @param data
	 *            the data
	 */
	private FloatImage3D(ImageStack stack, float[] data)
	{
		super(stack.getWidth(), stack.getHeight(), stack.getSize(), stack.getWidth() * stack.getHeight());

		// This is used internally so the data is the correct length
		this.data = data;

		if (stack.getBitDepth() == 32)
		{
			for (int s = 0; s < ns; s++)
			{
				System.arraycopy(stack.getPixels(s + 1), 0, data, s * nr_by_nc, nr_by_nc);
			}
		}
		else
		{
			for (int s = 1, i = 0; s <= ns; s++)
			{
				ImageProcessor ip = stack.getProcessor(s);
				for (int j = 0; i < nr_by_nc; j++)
					data[i++] = ip.getf(j);
			}
		}
	}

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @param data
	 *            the data
	 * @throws IllegalArgumentException
	 *             If the data is not the correct length
	 */
	public FloatImage3D(int nc, int nr, int ns, float[] data) throws IllegalArgumentException
	{
		// Avoid constructor that calls createData(int)
		super(nc, nr, ns, nr * nc);
		if (data == null || data.length != checkSize(nc, nr, ns, true))
			throw new IllegalArgumentException("Data is not correct length");
		this.data = data;
	}

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @param nr_by_nc
	 *            the number of rows multiplied by the number of columns
	 * @param data
	 *            the data
	 */
	protected FloatImage3D(int nc, int nr, int ns, int nr_by_nc, float[] data)
	{
		// No checks as this is used internally		
		super(nc, nr, ns, nr_by_nc);
		this.data = data;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#copy()
	 */
	public FloatImage3D copy()
	{
		return new FloatImage3D(nc, nr, ns, nr_by_nc, data.clone());
	}

	/**
	 * Gets the data.
	 *
	 * @return the data
	 */
	public float[] getData()
	{
		return data;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#getDataLength()
	 */
	@Override
	public int getDataLength()
	{
		return data.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#crop(int, int, int, int, int, int)
	 */
	@Override
	public FloatImage3D crop(int x, int y, int z, int w, int h, int d) throws IllegalArgumentException
	{
		return crop(x, y, z, w, h, d, null);
	}

	/**
	 * Crop a sub-region of the data. The target dimensions must be positive.
	 *
	 * @param x
	 *            the x index
	 * @param y
	 *            the y index
	 * @param z
	 *            the z index
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param d
	 *            the depth
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public FloatImage3D crop(int x, int y, int z, int w, int h, int d, float[] region) throws IllegalArgumentException
	{
		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr || z < 0 || d < 1 ||
				(long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		int size = d * h * w;
		if (region == null || region.length != size)
			region = new float[size];
		for (int s = 0, i = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				System.arraycopy(data, base, region, i, w);
				base += nc;
				i += w;
			}
		}
		return new FloatImage3D(w, h, d, w * h, region);
	}

	/**
	 * Crop a sub-region of the data. The target dimensions must be positive.
	 *
	 * @param stack
	 *            the stack
	 * @param x
	 *            the x index
	 * @param y
	 *            the y index
	 * @param z
	 *            the z index
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param d
	 *            the depth
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public static FloatImage3D crop(ImageStack stack, int x, int y, int z, int w, int h, int d, float[] region)
			throws IllegalArgumentException
	{
		int nc = stack.getWidth();
		int nr = stack.getHeight();
		int ns = stack.getSize();

		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr || z < 0 || d < 1 ||
				(long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		int size = checkSize(w, h, d, true);
		if (region == null || region.length != size)
			region = new float[size];

		if (stack.getBitDepth() != 32)
		{
			// Handle non-float data
			for (int s = 0, i = 0; s < d; s++, z++)
			{
				ImageProcessor ip = stack.getProcessor(1 + z);
				int base = y * nc + x;
				for (int r = 0; r < h; r++)
				{
					for (int c = 0; c < w; c++)
					{
						region[i++] = ip.getf(base + c);
					}
					base += nc;
				}
			}
		}
		else
		{
			for (int s = 0, i = 0; s < d; s++, z++)
			{
				float[] data = (float[]) stack.getPixels(1 + z);
				int base = y * nc + x;
				for (int r = 0; r < h; r++)
				{
					System.arraycopy(data, base, region, i, w);
					base += nc;
					i += w;
				}
			}
		}
		return new FloatImage3D(w, h, d, w * h, region);
	}

	@Override
	public void insert(int x, int y, int z, Image3D image) throws IllegalArgumentException
	{
		if (image instanceof FloatImage3D)
		{
			insert(x, y, z, (FloatImage3D) image);
		}
		else
		{
			super.insert(x, y, z, image);
		}
	}

	/**
	 * Insert a sub-region.
	 *
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param z
	 *            the z position
	 * @param image
	 *            the image
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public void insert(int x, int y, int z, FloatImage3D image) throws IllegalArgumentException
	{
		// Check the region range
		int w = image.getWidth();
		int h = image.getHeight();
		int d = image.getSize();
		if (w < 1 || h < 1 || d < 1)
			return;
		if (x < 0 || (long) x + w > nc || y < 0 || (long) y + h > nr || z < 0 || (long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		float[] region = image.data;
		for (int s = 0, i = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				System.arraycopy(region, i, data, base, w);
				base += nc;
				i += w;
			}
		}
	}

	@Override
	protected void copyTo(int i, float[] buffer, int j, int size)
	{
		System.arraycopy(data, i, buffer, j, size);
	}

	@Override
	protected void copyFrom(float[] buffer, int j, int size, int i)
	{
		System.arraycopy(buffer, j, data, i, size);
	}

	@Override
	public double get(int i)
	{
		return data[i];
	}

	@Override
	public void set(int i, double value)
	{
		data[i] = (float) value;
	}

	@Override
	public float getf(int i)
	{
		return data[i];
	}

	@Override
	public void setf(int i, float value)
	{
		data[i] = value;
	}

	@Override
	protected void fill(int i, int size, double value)
	{
		final float v = (float) value;
		while (size-- > 0)
			data[i++] = v;
	}
}
