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

import gdsc.core.utils.Maths;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Store a 3D image in a single array
 */
public abstract class Image3D
{
	/**
	 * The largest array size for which a regular 1D Java array is used to store the
	 * data (2^30)
	 */
	public static final int MAX_SIZE_OF_32_BIT_ARRAY = 1073741824;

	/**
	 * Check the size can fit in a 1D array
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @param raiseException
	 *            Set to true to raise an exception if too large
	 * @return the size, or -1 if too large
	 * @throws IllegalArgumentException
	 *             if too large (optional)
	 */
	public static int checkSize(int nc, int nr, int ns, boolean raiseException) throws IllegalArgumentException
	{
		if (nc < 0 || nr < 0 || ns < 0)
		{
			if (raiseException)
				throw new IllegalArgumentException("Negative dimensions");
			return -1;
		}
		long size = (long) ns * nr * nc;
		if (size > MAX_SIZE_OF_32_BIT_ARRAY)
		{
			if (raiseException)
				throw new IllegalArgumentException("3D data too large");
			return -1;
		}
		return (int) size;
	}

	/** The number of slices (max z). */
	public final int ns;
	/** The number of rows (max y). */
	public final int nr;
	/** The number of columns (max x). */
	public final int nc;

	/** The number of rows multiplied by the number of columns */
	public final int nr_by_nc;

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
	public Image3D(int nc, int nr, int ns) throws IllegalArgumentException
	{
		createData(checkSize(nc, nr, ns, true));
		this.nc = nc;
		this.nr = nr;
		this.ns = ns;
		nr_by_nc = nr * nc;
	}

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param stack
	 *            the stack
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public Image3D(ImageStack stack) throws IllegalArgumentException
	{
		nc = stack.getWidth();
		nr = stack.getHeight();
		ns = stack.getSize();
		createData(checkSize(nc, nr, ns, true));
		nr_by_nc = nr * nc;
		if (stack.getBitDepth() == 32)
		{
			for (int s = 0; s < ns; s++)
			{
				copyFrom((float[]) stack.getPixels(s + 1), 0, nr_by_nc, s * nr_by_nc);
			}
		}
		else
		{
			for (int s = 1, i = 0; s <= ns; s++)
			{
				ImageProcessor ip = stack.getProcessor(s);
				for (int j = 0; i < nr_by_nc; j++)
					setf(i++, ip.getf(j));
			}
		}
	}

	/**
	 * Instantiates a new 3D image. It is assumed that the sub-class will correctly create the data storage.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @param nr_by_nc
	 *            the number of rows multiplied by the number of columns
	 */
	protected Image3D(int nc, int nr, int ns, int nr_by_nc)
	{
		// No checks as this is used internally
		this.nc = nc;
		this.nr = nr;
		this.ns = ns;
		this.nr_by_nc = nr_by_nc;
	}

	/**
	 * Creates the array to store the data.
	 *
	 * @param size
	 *            the size of the array
	 */
	protected abstract void createData(int size);

	/**
	 * Copy the data from the given index to the buffer.
	 * <p>
	 * Utility method to handle conversion with ImageJ ImageProcessor objects.
	 *
	 * @param i
	 *            the index
	 * @param buffer
	 *            the buffer
	 * @param j
	 *            the buffer index
	 * @param size
	 *            the size
	 */
	protected abstract void copyTo(int i, float[] buffer, int j, int size);

	/**
	 * Copy the data from the given buffer to the given index.
	 * <p>
	 * Utility method to handle conversion with ImageJ ImageProcessor objects.
	 *
	 * @param i
	 *            the index
	 * @param buffer
	 *            the buffer
	 * @param j
	 *            the buffer index
	 * @param size
	 *            the size
	 */
	protected abstract void copyFrom(float[] buffer, int j, int size, int i);

	/**
	 * Gets the value at the given index.
	 *
	 * @param i
	 *            the index
	 * @return the value
	 */
	public abstract double get(int i);

	/**
	 * Sets the value at the given index.
	 *
	 * @param i
	 *            the index
	 * @param value
	 *            the value
	 */
	public abstract void set(int i, double value);

	/**
	 * Gets the value at the given index.
	 *
	 * @param i
	 *            the index
	 * @return the value
	 */
	public abstract float getf(int i);

	/**
	 * Sets the value at the given index.
	 *
	 * @param i
	 *            the index
	 * @param value
	 *            the value
	 */
	public abstract void setf(int i, float value);

	/**
	 * Return a copy of the 3D image.
	 *
	 * @return the copy
	 */
	public abstract Image3D copy();

	/**
	 * Gets the width (the number of columns).
	 *
	 * @return the width
	 */
	public int getWidth()
	{
		return nc;
	}

	/**
	 * Gets the height (the number of rows).
	 *
	 * @return the height
	 */
	public int getHeight()
	{
		return nr;
	}

	/**
	 * Gets the size (the number of slices)
	 *
	 * @return the size
	 */
	public int getSize()
	{
		return ns;
	}

	/**
	 * Gets the data length.
	 *
	 * @return the data length
	 */
	public int getDataLength()
	{
		return ns * nr_by_nc;
	}

	/**
	 * Convert to an image stack.
	 *
	 * @return the image stack
	 */
	public ImageStack getImageStack()
	{
		ImageStack stack = new ImageStack(nc, nr);
		for (int s = 0; s < ns; s++)
		{
			float[] pixels = new float[nr_by_nc];
			copyTo(s * nr_by_nc, pixels, 0, nr_by_nc);
			stack.addSlice(null, pixels);
		}
		return stack;
	}

	/**
	 * Gets the xyz components of the index.
	 *
	 * @param i
	 *            the index
	 * @return the xyz components
	 * @throws IllegalArgumentException
	 *             if the index is not within the data
	 */
	public int[] getXYZ(int i) throws IllegalArgumentException
	{
		if (i < 0 || i >= getDataLength())
			throw new IllegalArgumentException("Index in not in the correct range: 0 <= i < " + getDataLength());
		int[] xyz = new int[3];
		xyz[2] = i / nr_by_nc;
		int j = i % nr_by_nc;
		xyz[1] = j / nc;
		xyz[0] = j % nc;
		return xyz;
	}

	/**
	 * Gets the xyz components of the index.
	 *
	 * @param i
	 *            the index
	 * @param xyz
	 *            the xyz components (must be an array of at least length 3)
	 * @throws IllegalArgumentException
	 *             if the index is not within the data
	 */
	public void getXYZ(int i, int[] xyz) throws IllegalArgumentException
	{
		if (i < 0 || i >= getDataLength())
			throw new IllegalArgumentException("Index in not in the correct range: 0 <= i < " + getDataLength());
		xyz[2] = i / nr_by_nc;
		int j = i % nr_by_nc;
		xyz[1] = j / nc;
		xyz[0] = j % nc;
	}

	/**
	 * Gets the index using the xyz components
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @return the index
	 * @throws IllegalArgumentException
	 *             if the index is not within the data
	 */
	public int getIndex(int x, int y, int z) throws IllegalArgumentException
	{
		if (x < 0 || x >= nc || y < 0 || y >= nr || z < 0 || z >= ns)
			throw new IllegalArgumentException("Index in not inside the image");
		return index(x, y, z);
	}

	/**
	 * Gets the index using the xyz components
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @return the index
	 */
	protected int index(int x, int y, int z)
	{
		return z * nr_by_nc + y * nc + x;
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
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public abstract Image3D crop(int x, int y, int z, int w, int h, int d) throws IllegalArgumentException;

	/**
	 * Crop a sub-region of the data into the given image. The target dimensions must be positive.
	 *
	 * @param x
	 *            the x index
	 * @param y
	 *            the y index
	 * @param z
	 *            the z index
	 * @param image
	 *            the image
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public Image3D crop(int x, int y, int z, Image3D image) throws IllegalArgumentException
	{
		// Check the region range
		int w = image.getWidth();
		int h = image.getHeight();
		int d = image.getSize();
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr || z < 0 || d < 1 ||
				(long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		for (int s = 0, i = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int c = 0; c < w; c++)
					image.set(i++, get(base + c));
				base += nc;
			}
		}
		return image;
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
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public ImageStack cropToStack(int x, int y, int z, int w, int h, int d) throws IllegalArgumentException
	{
		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr || z < 0 || d < 1 ||
				(long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		int size = w * h;
		ImageStack stack = new ImageStack(w, h, d);
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			float[] region = new float[size];
			for (int r = 0, i = 0; r < h; r++)
			{
				copyTo(base, region, i, w);
				base += nc;
				i += w;
			}
			stack.setPixels(region, 1 + s);
		}
		return stack;
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
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public static ImageStack cropToStack(ImageStack stack, int x, int y, int z, int w, int h, int d)
			throws IllegalArgumentException
	{
		int nc = stack.getWidth();
		int nr = stack.getHeight();
		int ns = stack.getSize();

		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr || z < 0 || d < 1 ||
				(long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		ImageStack stack2 = new ImageStack(w, h, d);
		for (int s = 0; s < d; s++, z++)
		{
			ImageProcessor ip = stack.getProcessor(1 + z);
			ip.setRoi(x, y, w, h);
			stack2.setPixels(ip.crop().getPixels(), 1 + s);
		}
		return stack2;
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
	public void insert(int x, int y, int z, Image3D image) throws IllegalArgumentException
	{
		// Check the region range
		int w = image.getWidth();
		int h = image.getHeight();
		int d = image.getSize();
		if (w < 1 || h < 1 || d < 1)
			return;
		if (x < 0 || (long) x + w > nc || y < 0 || (long) y + h > nr || z < 0 || (long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		for (int s = 0, i = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int c = 0; c < w; c++)
					set(base + c, image.get(i++));
				base += nc;
			}
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
	 * @param stack
	 *            the image stack
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public void insert(int x, int y, int z, ImageStack stack) throws IllegalArgumentException
	{
		// Check the region range
		int w = stack.getWidth();
		int h = stack.getHeight();
		int d = stack.getSize();
		if (w < 1 || h < 1 || d < 1)
			return;
		if (x < 0 || (long) x + w > nc || y < 0 || (long) y + h > nr || z < 0 || (long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		boolean isFloat = stack.getBitDepth() == 32;
		FloatProcessor fp = (isFloat) ? new FloatProcessor(w, h) : null;
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			float[] region = (float[]) ((isFloat) ? stack.getPixels(1 + s)
					: stack.getProcessor(1 + s).toFloat(0, fp).getPixels());
			for (int r = 0, i = 0; r < h; r++)
			{
				copyFrom(region, i, w, base);
				base += nc;
				i += w;
			}
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
	public void insert(int x, int y, int z, ImageProcessor image) throws IllegalArgumentException
	{
		// Check the region range
		int w = image.getWidth();
		int h = image.getHeight();
		if (w < 1 || h < 1)
			return;
		if (x < 0 || (long) x + w > nc || y < 0 || (long) y + h > nr || z < 0 || z >= ns)
			throw new IllegalArgumentException("Region not within the data");
		boolean isFloat = image.getBitDepth() == 32;
		int base = z * nr_by_nc + y * nc + x;
		float[] region = (float[]) ((isFloat) ? image.getPixels() : image.toFloat(0, null).getPixels());
		for (int r = 0, i = 0; r < h; r++)
		{
			copyFrom(region, i, w, base);
			base += nc;
			i += w;
		}
	}

	/**
	 * Compute 3D intersect with this object.
	 * <p>
	 * If any of w,h,d are negative then the corresponding x,y,z is updated and the w,h,d is inverted.
	 * The maximum bounds of the given dimensions are then computed by adding the
	 * w,h,d to the x,y,z. The bounds are then clipped to the image dimensions and the intersect returned.
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
	 * @return [x,y,z,w,h,d]
	 */
	public int[] computeIntersect(int x, int y, int z, int w, int h, int d)
	{
		if (w < 0)
		{
			w = -w;
			x = subtract(x, w);
		}
		if (h < 0)
		{
			h = -h;
			y = subtract(y, h);
		}
		if (d < 0)
		{
			d = -d;
			z = subtract(z, d);
		}
		// Compute 3D intersect with this object
		int x2 = clip(nc, x, w);
		int y2 = clip(nr, y, h);
		int z2 = clip(ns, z, d);
		x = clip(nc, x);
		y = clip(nr, y);
		z = clip(ns, z);
		return new int[] { x, y, z, x2 - x, y2 - y, z2 - z };
	}

	/**
	 * Compute 3D intersect with this object or throw an exception if the intersect has no volume.
	 * <p>
	 * If any of w,h,d are negative then the corresponding x,y,z is updated and the w,h,d is inverted.
	 * The maximum bounds of the given dimensions are then computed by adding the
	 * w,h,d to the x,y,z. The bounds are then clipped to the image dimensions and the intersect returned.
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
	 * @return [x,y,z,w,h,d]
	 * @throws IllegalArgumentException
	 *             if the intersect has no volume
	 */
	public int[] computeIntersectOrThrow(int x, int y, int z, int w, int h, int d) throws IllegalArgumentException
	{
		if (w < 0)
		{
			w = -w;
			x = subtract(x, w);
		}
		if (h < 0)
		{
			h = -h;
			y = subtract(y, h);
		}
		if (d < 0)
		{
			d = -d;
			z = subtract(z, d);
		}
		// Compute 3D intersect with this object
		int x2 = clip(nc, x, w);
		int y2 = clip(nr, y, h);
		int z2 = clip(ns, z, d);
		x = clip(nc, x);
		y = clip(nr, y);
		z = clip(ns, z);
		w = checkSize(x2 - x);
		h = checkSize(y2 - y);
		d = checkSize(z2 - z);
		return new int[] { x, y, z, w, h, d };
	}

	/**
	 * Check size is above zero or throw an exception.
	 *
	 * @param size
	 *            the size
	 * @return the size
	 * @throws IllegalArgumentException
	 *             If the size if zero
	 */
	private static int checkSize(int size) throws IllegalArgumentException
	{
		if (size == 0)
			throw new IllegalArgumentException("No intersect");
		return size;
	}

	/**
	 * Subtract the value avoiding underflow.
	 *
	 * @param value
	 *            the value
	 * @param subtraction
	 *            the subtraction
	 * @return the value
	 */
	private static int subtract(int value, int subtraction)
	{
		// Avoid underflow
		long v = (long) value - subtraction;
		return (v < Integer.MIN_VALUE) ? Integer.MIN_VALUE : (int) v;
	}

	/**
	 * Return value clipped to within the given bounds (0 - upper).
	 *
	 * @param upper
	 *            the upper limit
	 * @param value
	 *            the value
	 * @param addition
	 *            the addition
	 * @return the clipped value
	 */
	private static int clip(int upper, int value, int addition)
	{
		// Avoid overflow
		long v = (long) value + addition;
		if (v < 0)
			return 0;
		if (v > upper)
			return upper;
		return (int) v;
	}

	/**
	 * Return value clipped to within the given bounds (0 - upper).
	 *
	 * @param upper
	 *            the upper limit
	 * @param value
	 *            the value
	 * @return the clipped value
	 */
	private static int clip(int upper, int value)
	{
		return (value < 0) ? 0 : (value > upper) ? upper : value;
	}

	/**
	 * Find the index of the minimum value in the region.
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
	 * @return the min index
	 * @throws IllegalArgumentException
	 *             if there is no intersect
	 * @throws IllegalStateException
	 *             if the region is just NaN values
	 */
	public int findMinIndex(int x, int y, int z, int w, int h, int d)
			throws IllegalArgumentException, IllegalStateException
	{
		int[] intersect = computeIntersectOrThrow(x, y, z, w, h, d);
		x = intersect[0];
		y = intersect[1];
		z = intersect[2];
		w = intersect[3];
		h = intersect[4];
		d = intersect[5];
		int index = findValueIndex(x, y, z, w, h, d);
		double min = get(index);
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int j = 0; j < w; j++)
					if (get(base + j) < min)
					{
						index = base + j;
						min = get(index);
					}
				base += nc;
			}
		}
		return index;
	}

	/**
	 * Find the index of the minimum value in the region.
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
	 * @return the min index
	 * @throws IllegalArgumentException
	 *             if there is no intersect
	 * @throws IllegalStateException
	 *             if the region is just NaN values
	 */
	public int findMaxIndex(int x, int y, int z, int w, int h, int d)
			throws IllegalArgumentException, IllegalStateException
	{
		int[] intersect = computeIntersectOrThrow(x, y, z, w, h, d);
		x = intersect[0];
		y = intersect[1];
		z = intersect[2];
		w = intersect[3];
		h = intersect[4];
		d = intersect[5];
		int index = findValueIndex(x, y, z, w, h, d);
		double max = get(index);
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int j = 0; j < w; j++)
				{
					if (get(base + j) > max)
					{
						index = base + j;
						max = get(index);
					}
				}
				base += nc;
			}
		}
		return index;
	}

	/**
	 * Find the index of the first non-NaN value in the region.
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
	 * @return the index
	 * @throws IllegalStateException
	 *             if the region is just NaN values
	 */
	private int findValueIndex(int x, int y, int z, int w, int h, int d) throws IllegalStateException
	{
		// Quick check without loops
		if (!Double.isNaN(get(z * nr_by_nc + y * nc + x)))
			return z * nr_by_nc + y * nc + x;

		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int j = 0; j < w; j++)
					if (!Double.isNaN(get(base + j)))
						return base + j;
				base += nc;
			}
		}
		throw new IllegalStateException("Region is NaN");
	}

	/**
	 * Compute the sum of the region.
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
	 * @return the sum
	 */
	public double computeSum(int x, int y, int z, int w, int h, int d)
	{
		int[] intersect = computeIntersect(x, y, z, w, h, d);
		w = intersect[3];
		h = intersect[4];
		d = intersect[5];
		// Recheck bounds
		if (w == 0 || h == 0 || d == 0)
			return 0;
		x = intersect[0];
		y = intersect[1];
		z = intersect[2];
		double sum = 0;
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int j = 0; j < w; j++)
					sum += get(base + j);
				base += nc;
			}
		}
		return sum;
	}

	/**
	 * Compute the rolling sum table for use in {@link #computeSum(double[], int, int, int, int, int, int)}.
	 * <p>
	 * This is a table of the sum of the volume from 0,0,0 to x,y,z inclusive.
	 *
	 * @param table
	 *            the table (will be reused if the correct size)
	 * @return the rolling sum table
	 */
	public double[] computeRollingSumTable(double[] table)
	{
		if (table == null || table.length != getDataLength())
			table = new double[getDataLength()];

		// First build a table for each XY slice
		for (int s = 0; s < ns; s++)
		{
			double sum = 0;
			int i = s * nr_by_nc;
			// Initialise first row sum
			// sum = rolling sum of (0 - colomn)
			for (int c = 0; c < nc; c++, i++)
			{
				sum += get(i);
				table[i] = sum;
			}
			// Remaining rows
			// sum = rolling sum of (0 - colomn) + sum of same position above
			for (int r = 1, ii = i - nc; r < nr; r++)
			{
				sum = 0;
				for (int c = 0; c < nc; c++, i++)
				{
					sum += get(i);
					// Add the sum from the previous row
					table[i] = sum + table[ii++];
				}
			}
		}

		// Now sum across slices
		// sum = rolling sum of (0,0 to column,row) + sum of same position above
		// => rolling sum of (0,0,0 to column,row,slice)
		for (int s = 1; s < ns; s++)
		{
			int i = s * nr_by_nc;
			int ii = i - nr_by_nc;
			for (int j = 0; j < nr_by_nc; j++)
				table[i++] += table[ii++];
		}

		return table;
	}

	/**
	 * Compute the sum of the region using the precomputed rolling sum table.
	 *
	 * @param table
	 *            the rolling sum table
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
	 * @return the sum
	 */
	public double computeSum(double[] table, int x, int y, int z, int w, int h, int d)
	{
		int[] intersect = computeIntersect(x, y, z, w, h, d);
		w = intersect[3];
		h = intersect[4];
		d = intersect[5];
		// Recheck bounds
		if (w == 0 || h == 0 || d == 0)
			return 0;
		//x = intersect[0];
		//y = intersect[1];
		//z = intersect[2];

		// Compute sum from rolling sum using:
		// sum(x,y,z,w,h,d) =
		// + s(x+w-1,y+h-1,z+d-1)
		// - s(x-1,y+h-1,z+d-1)
		// - s(x+w-1,y-1,z+d-1)
		// + s(x-1,y-1,z+d-1)
		// /* Stack above must be subtracted so reverse sign*/
		// - s(x+w-1,y+h-1,z-1)
		// + s(x-1,y+h-1,z-1)
		// + s(x+w-1,y-1,z-1)
		// - s(x-1,y-1,z-1)
		// Note:
		// s(i,j,k) = 0 when either i,j,k < 0
		// i = imax when i>imax
		// j = jmax when j>jmax
		// k = kmax when k>kmax

		int x_1 = intersect[0] - 1;
		int y_1 = intersect[1] - 1;
		int z_1 = intersect[2] - 1;
		// The intersect has already checked the bounds
		//int x_w_1 = Math.min(x_1 + w, nc);
		//int y_h_1 = Math.min(y_1 + h, nr);
		//int z_d_1 = Math.min(z_1 + d, ns);
		int x_w_1 = x_1 + w;
		int y_h_1 = y_1 + h;
		int z_d_1 = z_1 + d;

		//double sum = table[index(x_w_1, y_h_1, z_d_1)];
		//if (y_1 >= 0)
		//{
		//	sum -= table[index(x_w_1, y_1, z_d_1)];
		//	if (x_1 >= 0)
		//		sum = sum + table[index(x_1, y_1, z_d_1)] - table[index(x_1, y_h_1, z_d_1)];
		//}
		//else if (x_1 >= 0)
		//{
		//	sum -= table[index(x_1, y_h_1, z_d_1)];
		//}
		//if (z_1 >= 0)
		//{
		//	sum -= table[index(x_w_1, y_h_1, z_1)];
		//	if (y_1 >= 0)
		//	{
		//		sum += table[index(x_w_1, y_1, z_1)];
		//		if (x_1 >= 0)
		//			sum = sum - table[index(x_1, y_1, z_1)] + table[index(x_1, y_h_1, z_1)];
		//	}
		//	else if (x_1 >= 0)
		//	{
		//		sum += table[index(x_1, y_h_1, z_1)];
		//	}
		//}
		//return sum;

		// This has been ordered to use the smallest sums first (i.e. closer to x,y,z than x+w,y+h,z+d)
		int xw_yh_zd = index(x_w_1, y_h_1, z_d_1);
		if (z_1 >= 0)
		{
			int xw_yh_z = xw_yh_zd - d * nr_by_nc;
			double sum = 0;
			if (y_1 >= 0)
			{
				int h_ = h * nc;
				if (x_1 >= 0)
					sum = table[xw_yh_zd - w - h_] - table[xw_yh_z - w - h_] - table[xw_yh_zd - w] + table[xw_yh_z - w];
				sum = sum + table[xw_yh_z - h_] - table[xw_yh_zd - h_];
			}
			else if (x_1 >= 0)
				sum = table[xw_yh_z - w] - table[xw_yh_zd - w];
			return sum + table[xw_yh_zd] - table[xw_yh_z];
		}
		double sum = 0;
		if (y_1 >= 0)
		{
			int h_ = h * nc;
			if (x_1 >= 0)
				sum = table[xw_yh_zd - w - h_] - table[xw_yh_zd - w];
			sum -= table[xw_yh_zd - h_];
		}
		else if (x_1 >= 0)
			sum = -table[xw_yh_zd - w];
		return sum + table[xw_yh_zd];
	}

	/**
	 * Compute the sum of the region using the precomputed rolling sum table. Assumes x+w,y+h,z+d will not overflow!
	 *
	 * @param table
	 *            the rolling sum table
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
	 * @return the sum
	 */
	public double computeSumFast(double[] table, int x, int y, int z, int w, int h, int d)
	{
		if (w <= 0 || h <= 0 || d <= 0 || x >= nc || y >= nr || z >= ns)
			return 0;

		// Compute sum from rolling sum using:
		// sum(x,y,z,w,h,d) =
		// + s(x+w-1,y+h-1,z+d-1)
		// - s(x-1,y+h-1,z+d-1)
		// - s(x+w-1,y-1,z+d-1)
		// + s(x-1,y-1,z+d-1)
		// /* Stack above must be subtracted so reverse sign*/
		// - s(x+w-1,y+h-1,z-1)
		// + s(x-1,y+h-1,z-1)
		// + s(x+w-1,y-1,z-1)
		// - s(x-1,y-1,z-1)
		// Note:
		// s(i,j,k) = 0 when either i,j,k < 0
		// i = imax when i>imax
		// j = jmax when j>jmax
		// k = kmax when k>kmax

		// Compute bounds assuming w,h,d is small and positive.
		int x_1, y_1, z_1, x_w_1, y_h_1, z_d_1;
		if (x < 0)
		{
			x_1 = 0;
			x_w_1 = Maths.clip(0, nc, x + w);
		}
		else
		{
			x_1 = x;
			x_w_1 = Math.min(nc, x + w);
		}
		w = x_w_1 - x_1;
		if (w == 0)
			return 0;
		if (y < 0)
		{
			y_1 = 0;
			y_h_1 = Maths.clip(0, nr, y + h);
		}
		else
		{
			y_1 = y;
			y_h_1 = Math.min(nr, y + h);
		}
		h = y_h_1 - y_1;
		if (h == 0)
			return 0;
		if (z < 0)
		{
			z_1 = 0;
			z_d_1 = Maths.clip(0, ns, z + d);
		}
		else
		{
			z_1 = z;
			z_d_1 = Math.min(ns, z + d);
		}
		d = z_d_1 - z_1;
		if (d == 0)
			return 0;

		// Adjust for the -1
		x_1--;
		y_1--;
		z_1--;
		x_w_1--;
		y_h_1--;
		z_d_1--;

		// This has been ordered to use the smallest sums first (i.e. closer to x,y,z than x+w,y+h,z+d)
		int xw_yh_zd = index(x_w_1, y_h_1, z_d_1);
		if (z_1 >= 0)
		{
			int xw_yh_z = xw_yh_zd - d * nr_by_nc;
			double sum = 0;
			if (y_1 >= 0)
			{
				int h_ = h * nc;
				if (x_1 >= 0)
					sum = table[xw_yh_zd - w - h_] - table[xw_yh_z - w - h_] - table[xw_yh_zd - w] + table[xw_yh_z - w];
				sum = sum + table[xw_yh_z - h_] - table[xw_yh_zd - h_];
			}
			else if (x_1 >= 0)
				sum = table[xw_yh_z - w] - table[xw_yh_zd - w];
			return sum + table[xw_yh_zd] - table[xw_yh_z];
		}
		double sum = 0;
		if (y_1 >= 0)
		{
			int h_ = h * nc;
			if (x_1 >= 0)
				sum = table[xw_yh_zd - w - h_] - table[xw_yh_zd - w];
			sum -= table[xw_yh_zd - h_];
		}
		else if (x_1 >= 0)
			sum = -table[xw_yh_zd - w];
		return sum + table[xw_yh_zd];
	}

	/**
	 * Fill the image.
	 *
	 * @param value
	 *            the value
	 */
	public void fill(double value)
	{
		fill(0, getDataLength(), value);
	}

	/**
	 * Fill the region.
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
	 * @param value
	 *            the value
	 */
	public void fill(int x, int y, int z, int w, int h, int d, double value)
	{
		int[] intersect = computeIntersect(x, y, z, w, h, d);
		w = intersect[3];
		h = intersect[4];
		d = intersect[5];
		// Recheck bounds
		if (w == 0 || h == 0 || d == 0)
			return;
		x = intersect[0];
		y = intersect[1];
		z = intersect[2];
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				fill(base, w, value);
				base += nc;
			}
		}
	}

	/**
	 * Fill with the given value from the given index.
	 *
	 * @param i
	 *            the index
	 * @param size
	 *            the size to fill
	 * @param value
	 *            the value
	 */
	protected abstract void fill(int i, int size, double value);

	/**
	 * Fill outside the region. If the region is not within the image then the entire image is filled.
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
	 * @param value
	 *            the value
	 */
	public void fillOutside(int x, int y, int z, int w, int h, int d, double value)
	{
		int[] intersect = computeIntersect(x, y, z, w, h, d);
		w = intersect[3];
		h = intersect[4];
		d = intersect[5];
		// Recheck bounds
		if (w == 0 || h == 0 || d == 0)
		{
			fill(value);
			return;
		}
		x = intersect[0];
		y = intersect[1];
		z = intersect[2];

		// Before
		if (z > 0)
			fill(0, z * nr_by_nc, value);
		// After
		if (z + d < ns)
			fill((z + d) * nr_by_nc, (ns - z - d) * nr_by_nc, value);

		int y_p_h = y + h;
		int fillYBefore = y * nc;
		int yAfter = y_p_h * nc;
		int fillYAfter = (nr - y_p_h) * nc;

		int x_p_w = x + w;
		int fillXBefore = x;
		int fillXAfter = (nc - x_p_w);

		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc;

			if (fillYBefore != 0)
				fill(base, fillYBefore, value);
			if (fillYAfter != 0)
				fill(base + yAfter, fillYAfter, value);

			base += fillYBefore;

			for (int r = 0; r < h; r++)
			{
				if (fillXBefore != 0)
					fill(base, fillXBefore, value);
				if (fillXAfter != 0)
					fill(base + x_p_w, fillXAfter, value);
				base += nc;
			}
		}
	}
}
