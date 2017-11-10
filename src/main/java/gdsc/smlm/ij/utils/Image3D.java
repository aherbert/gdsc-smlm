package gdsc.smlm.ij.utils;

import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Store a 3D image in a single float array
 */
public class Image3D
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

	protected final float[] data;

	/**
	 * Instantiates a new 3D image
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

		data = new float[checkSize(nc, nr, ns, true)];

		nr_by_nc = nr * nc;
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
	 * @param stack
	 *            the stack
	 * @param data
	 *            the data
	 */
	private Image3D(ImageStack stack, float[] data)
	{
		nc = stack.getWidth();
		nr = stack.getHeight();
		ns = stack.getSize();

		this.data = data;

		nr_by_nc = nr * nc;
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
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public Image3D(int nc, int nr, int ns) throws IllegalArgumentException
	{
		data = new float[checkSize(nc, nr, ns, true)];
		this.nc = nc;
		this.nr = nr;
		this.ns = ns;
		nr_by_nc = nr * nc;
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
	public Image3D(int nc, int nr, int ns, float[] data) throws IllegalArgumentException
	{
		if (data == null || data.length != checkSize(nc, nr, ns, true))
			throw new IllegalArgumentException("Data is not correct length");
		this.nc = nc;
		this.nr = nr;
		this.ns = ns;
		nr_by_nc = nr * nc;
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
	protected Image3D(int nc, int nr, int ns, int nr_by_nc, float[] data)
	{
		// No checks as this is used internally		
		this.nc = nc;
		this.nr = nr;
		this.ns = ns;
		this.nr_by_nc = nr_by_nc;
		this.data = data;
	}

	/**
	 * Return a copy of the 3D image.
	 *
	 * @return the copy
	 */
	public Image3D copy()
	{
		return new Image3D(nc, nr, ns, nr_by_nc, data.clone());
	}

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
	 * Gets the data.
	 *
	 * @return the data
	 */
	public float[] getData()
	{
		return data;
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
			System.arraycopy(data, s * nr_by_nc, pixels, 0, nr_by_nc);
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
		if (i < 0 || i >= data.length)
			throw new IllegalArgumentException("Index in not in the correct range: 0 <= i < " + data.length);
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
		if (i < 0 || i >= data.length)
			throw new IllegalArgumentException("Index in not in the correct range: 0 <= i < " + data.length);
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
	private int index(int x, int y, int z)
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
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public Image3D crop(int x, int y, int z, int w, int h, int d, float[] region) throws IllegalArgumentException
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
		return new Image3D(w, h, d, w * h, region);
	}

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
		return crop(x, y, z, image.getWidth(), image.getHeight(), image.getSize(), image.data);
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
				System.arraycopy(data, base, region, i, w);
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
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public static Image3D crop(ImageStack stack, int x, int y, int z, int w, int h, int d, float[] region)
			throws IllegalArgumentException
	{
		if (stack.getBitDepth() != 32)
		{
			// Handle non-float data
			checkSize(w, h, d, true);
			stack = cropToStack(stack, x, y, z, w, h, d);
			int size = d * h * w;
			if (region == null || region.length != size)
				region = new float[size];
			return new Image3D(stack, region);
		}

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
		return new Image3D(w, h, d, w * h, region);
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
				System.arraycopy(region, i, data, base, w);
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
			System.arraycopy(region, i, data, base, w);
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
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int j = 0; j < w; j++)
					if (data[base + j] < data[index])
						index = base + j;
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
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int j = 0; j < w; j++)
					if (data[base + j] > data[index])
						index = base + j;
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
		if (!Float.isNaN(data[z * nr_by_nc + y * nc + x]))
			return z * nr_by_nc + y * nc + x;

		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int j = 0; j < w; j++)
					if (!Float.isNaN(data[base + j]))
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
					sum += data[base + j];
				base += nc;
			}
		}
		return sum;
	}

	/**
	 * Compute the rolling sum table for use in computeSum.
	 *
	 * @return the rolling sum table
	 */
	public double[] computeRollingSumTable()
	{
		return null;
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
		return 0;
	}
}