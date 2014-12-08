package gdsc.smlm.ij.utils;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Rectangle;

/**
 * Contains methods for converting an image to float data
 */
public class ImageConverter
{
	/**
	 * Get the data from the image processor as a float array (include cropping to the ROI)
	 * 
	 * @param ip
	 * @return The float array data
	 */
	public static float[] getData(ImageProcessor ip)
	{
		return getData(ip, null);
	}

	/**
	 * Get the data from the image processor as a float array (include cropping to the ROI). Data is duplicated if the
	 * InputImage is a FloatProcessor.
	 * <p>
	 * Allows reuse of an existing buffer if provided. This will not be truncated if it is larger than the
	 * ImageProcessor ROI bounds. If smaller then a new buffer will be created.
	 * 
	 * @param ip
	 * @param buffer
	 * @return The float array data
	 */
	public static float[] getData(ImageProcessor ip, float[] buffer)
	{
		if (ip == null)
			return null;

		return getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), ip.getRoi(), buffer);
	}

	/**
	 * Get the data from the image as a float array (include cropping to the ROI). Data is duplicated if the
	 * input already a float array.
	 * <p>
	 * Allows reuse of an existing buffer if provided. This will not be truncated if it is larger than the
	 * ImageProcessor ROI bounds. If smaller then a new buffer will be created.
	 * 
	 * @param oPixels
	 * @param width
	 * @param height
	 * @param bounds
	 * @param buffer
	 * @return The float array data
	 */
	public static float[] getData(final Object oPixels, final int width, final int height, final Rectangle bounds,
			float[] buffer)
	{
		if (oPixels == null)
			return null;

		// Ignore the bounds if they specify the entire image size

		if (oPixels instanceof float[])
		{
			float[] pixels = (float[]) oPixels;
			if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
			{
				float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
				for (int ys = bounds.y; ys < bounds.y + bounds.height; ys++)
				{
					int offset1 = (ys - bounds.y) * bounds.width;
					int offset2 = ys * width + bounds.x;
					for (int xs = 0; xs < bounds.width; xs++)
						pixels2[offset1++] = pixels[offset2++];
				}
				return pixels2;
			}
			else
			{
				float[] pixels2 = allocate(buffer, pixels.length);
				System.arraycopy(pixels, 0, pixels2, 0, pixels.length);
				return pixels2;
			}
		}
		else if (oPixels instanceof short[])
		{
			short[] pixels = (short[]) oPixels;
			if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
			{
				float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
				for (int ys = bounds.y; ys < bounds.y + bounds.height; ys++)
				{
					int offset1 = (ys - bounds.y) * bounds.width;
					int offset2 = ys * width + bounds.x;
					for (int xs = 0; xs < bounds.width; xs++)
						pixels2[offset1++] = pixels[offset2++] & 0xffff;
				}
				return pixels2;
			}
			else
			{
				float[] pixels2 = allocate(buffer, pixels.length);
				for (int i = 0; i < pixels.length; i++)
					pixels2[i] = pixels[i] & 0xffff;
				return pixels2;
			}
		}
		else if (oPixels instanceof byte[])
		{
			byte[] pixels = (byte[]) oPixels;
			if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
			{
				float[] pixels2 = allocate(buffer, bounds.width * bounds.height);
				for (int ys = bounds.y; ys < bounds.y + bounds.height; ys++)
				{
					int offset1 = (ys - bounds.y) * bounds.width;
					int offset2 = ys * width + bounds.x;
					for (int xs = 0; xs < bounds.width; xs++)
						pixels2[offset1++] = pixels[offset2++] & 0xff;
				}
				return pixels2;
			}
			else
			{
				float[] pixels2 = allocate(buffer, pixels.length);
				for (int i = 0; i < pixels.length; i++)
					pixels2[i] = pixels[i] & 0xff;
				return pixels2;
			}
		}
		else if (oPixels instanceof int[])
		{
			// The default processing
			ImageProcessor ip = new ColorProcessor(width, height, (int[]) oPixels);
			ip.setRoi(bounds);
			FloatProcessor fp = ip.crop().toFloat(0, null);
			return (float[]) fp.getPixels();
		}
		return null;
	}

	private static float[] allocate(float[] buffer, int size)
	{
		if (buffer == null || buffer.length < size)
			buffer = new float[size];
		return buffer;
	}
	

	/**
	 * Get the data from the image processor as a double array (include cropping to the ROI). Data is duplicated if the
	 * InputImage is a FloatProcessor.
	 * <p>
	 * Allows reuse of an existing buffer if provided. This will not be truncated if it is larger than the
	 * ImageProcessor ROI bounds. If smaller then a new buffer will be created.
	 * 
	 * @param ip
	 * @param buffer
	 * @return The double array data
	 */
	public static double[] getDoubleData(ImageProcessor ip, double[] buffer)
	{
		if (ip == null)
			return null;

		return getDoubleData(ip.getPixels(), ip.getWidth(), ip.getHeight(), ip.getRoi(), buffer);
	}

	/**
	 * Get the data from the image as a double array (include cropping to the ROI). Data is duplicated if the
	 * input already a double array.
	 * <p>
	 * Allows reuse of an existing buffer if provided. This will not be truncated if it is larger than the
	 * ImageProcessor ROI bounds. If smaller then a new buffer will be created.
	 * 
	 * @param oPixels
	 * @param width
	 * @param height
	 * @param bounds
	 * @param buffer
	 * @return The double array data
	 */
	public static double[] getDoubleData(final Object oPixels, final int width, final int height, final Rectangle bounds,
			double[] buffer)
	{
		if (oPixels == null)
			return null;

		// Ignore the bounds if they specify the entire image size

		if (oPixels instanceof float[])
		{
			float[] pixels = (float[]) oPixels;
			if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
			{
				double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
				for (int ys = bounds.y; ys < bounds.y + bounds.height; ys++)
				{
					int offset1 = (ys - bounds.y) * bounds.width;
					int offset2 = ys * width + bounds.x;
					for (int xs = 0; xs < bounds.width; xs++)
						pixels2[offset1++] = pixels[offset2++];
				}
				return pixels2;
			}
			else
			{
				double[] pixels2 = allocate(buffer, pixels.length);
				System.arraycopy(pixels, 0, pixels2, 0, pixels.length);
				return pixels2;
			}
		}
		else if (oPixels instanceof short[])
		{
			short[] pixels = (short[]) oPixels;
			if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
			{
				double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
				for (int ys = bounds.y; ys < bounds.y + bounds.height; ys++)
				{
					int offset1 = (ys - bounds.y) * bounds.width;
					int offset2 = ys * width + bounds.x;
					for (int xs = 0; xs < bounds.width; xs++)
						pixels2[offset1++] = pixels[offset2++] & 0xffff;
				}
				return pixels2;
			}
			else
			{
				double[] pixels2 = allocate(buffer, pixels.length);
				for (int i = 0; i < pixels.length; i++)
					pixels2[i] = pixels[i] & 0xffff;
				return pixels2;
			}
		}
		else if (oPixels instanceof byte[])
		{
			byte[] pixels = (byte[]) oPixels;
			if (bounds != null && (bounds.x != 0 || bounds.y != 0 || bounds.width != width || bounds.height != height))
			{
				double[] pixels2 = allocate(buffer, bounds.width * bounds.height);
				for (int ys = bounds.y; ys < bounds.y + bounds.height; ys++)
				{
					int offset1 = (ys - bounds.y) * bounds.width;
					int offset2 = ys * width + bounds.x;
					for (int xs = 0; xs < bounds.width; xs++)
						pixels2[offset1++] = pixels[offset2++] & 0xff;
				}
				return pixels2;
			}
			else
			{
				double[] pixels2 = allocate(buffer, pixels.length);
				for (int i = 0; i < pixels.length; i++)
					pixels2[i] = pixels[i] & 0xff;
				return pixels2;
			}
		}
		else if (oPixels instanceof int[])
		{
			// The default processing
			ImageProcessor ip = new ColorProcessor(width, height, (int[]) oPixels);
			ip.setRoi(bounds);
			FloatProcessor fp = ip.crop().toFloat(0, null);
			return (double[]) fp.getPixels();
		}
		return null;
	}

	private static double[] allocate(double[] buffer, int size)
	{
		if (buffer == null || buffer.length < size)
			buffer = new double[size];
		return buffer;
	}	
}
