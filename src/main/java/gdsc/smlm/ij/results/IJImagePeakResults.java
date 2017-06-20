package gdsc.smlm.ij.results;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Collection;

import gdsc.core.ij.Utils;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.utils.XmlUtils;
import ij.ImagePlus;
import ij.ImageStack;
import ij.MappedImageStack;
import ij.WindowManager;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.LutLoader;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.MappedFloatProcessor;
import ij.process.ShortProcessor;

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

/**
 * Saves the fit results to an ImageJ image
 */
public class IJImagePeakResults extends IJAbstractPeakResults
{
	public static final String IMAGE_SUFFIX = "SuperRes";
	/**
	 * Display the signal of the peak in the image. The default is a count of 1.
	 */
	public static final int DISPLAY_SIGNAL = 1;
	/**
	 * Interpolate the value over multiple pixels. Depending on the location in the containing pixel this is usually the
	 * 3 closest 8-connected neighbours and the containing pixel. It may be less if the containing pixel is at the image
	 * bounds.
	 */
	public static final int DISPLAY_WEIGHTED = 2;
	/**
	 * Equalise the histogram of the output image. Allows showing a high dynamic range by limiting bright pixels.
	 */
	public static final int DISPLAY_EQUALIZED = 4;
	/**
	 * Display the peak number in the image. The default is a count of 1.
	 */
	public static final int DISPLAY_PEAK = 8;
	/**
	 * Display the peak error in the image. The default is a count of 1.
	 */
	public static final int DISPLAY_ERROR = 16;

	/**
	 * Replace the pixels with the new value. This should not be used with {@link #DISPLAY_WEIGHTED} to avoid the value
	 * being interpolated over multiple pixels.
	 */
	public static final int DISPLAY_REPLACE = 32;
	/**
	 * Use the maximum value. This should not be used with {@link #DISPLAY_WEIGHTED} to avoid the value being
	 * interpolated over multiple pixels.
	 */
	public static final int DISPLAY_MAX = 64;
	/**
	 * Use this to support negative values
	 */
	public static final int DISPLAY_NEGATIVES = 128;
	/**
	 * Mapped all non-zero values to 1-255 in the 8-bit displayed image. Zero and below are mapped to 0 in the LUT.
	 * <p>
	 * This cannot be used with {@link #DISPLAY_EQUALIZED} or {@link #DISPLAY_NEGATIVES}.
	 */
	public static final int DISPLAY_MAPPED = 256;
	/**
	 * Mapped even zero to 1-255 in the 8-bit displayed image. -0.0f and below is mapped to 0 in the LUT. This can be
	 * used for example to display the result of a probability calculation where 0 is a valid display value but must be
	 * distinguished from pixels that have no value computed.
	 * <p>
	 * Must be used with {@link #DISPLAY_MAPPED}.
	 */
	public static final int DISPLAY_MAP_ZERO = 512;

	/** The empty value. */
	private double EMPTY = 0.0;

	protected final String title;
	protected final int imageWidth;
	protected final int imageHeight;
	protected final float scale;
	protected int size = 0;
	protected double[] data;
	protected final float xlimit;
	protected final float ylimit;
	protected ImagePlus imp = null;
	protected boolean imageActive = false;
	protected int displayFlags = 0;
	private int rollingWindowSize = 0;
	private boolean displayImage = true;
	private boolean liveImage = true;
	private boolean uncalibrated = false;

	// Used to draw the image
	private int lastPaintSize = 0;
	private int nextRepaintSize = 0;
	private long nextPaintTime = 0;
	private Object pixels;
	private boolean imageLock = false;
	private double repaintInterval = 0.1;
	private long repaintDelay = 1000;
	private int currentFrame;

	private String lutName = "fire";

	/**
	 * @param title
	 *            Title of the image (appended with a suffix)
	 * @param bounds
	 *            Define the bounding rectangle of the image coordinates. Any results outside this will not be
	 *            displayed.
	 * @param scale
	 *            The image scale. Must be strictly positive.
	 */
	public IJImagePeakResults(String title, Rectangle bounds, float scale)
	{
		if (scale <= 0 || Float.isNaN(scale))
			throw new IllegalArgumentException("Invalid scale: " + scale);

		this.title = title + " " + IMAGE_SUFFIX;

		this.bounds = (Rectangle) bounds.clone();
		if (bounds.width < 0)
			bounds.width = 0;
		if (bounds.height < 0)
			bounds.height = 0;
		this.scale = scale;

		imageWidth = ceil(bounds.width * scale);
		imageHeight = ceil(bounds.height * scale);

		// Set the limits used to check if a coordinate has 4 neighbour cells
		xlimit = imageWidth - 1;
		ylimit = imageHeight - 1;
	}

	private int ceil(float f)
	{
		return (int) Math.ceil(f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#begin()
	 */
	public void begin()
	{
		imageActive = false;

		preBegin();

		// Handle invalid bounds with an empty single pixel image
		boolean validBounds = imageWidth > 0 && imageHeight > 0 &&
				(double) imageWidth * (double) imageHeight < Integer.MAX_VALUE;
		int w, h;
		if (validBounds)
		{
			w = imageWidth;
			h = imageHeight;
		}
		else
		{
			Utils.log("ERROR: Unable to create image results '%s' due to invalid dimensions: width=%d, height=%d",
					title, imageWidth, imageHeight);
			w = h = 1;
		}

		// Q. Should this be changed to handle the data in non-pixel distances.
		// At the moment we hope that the results IO can work out the units and convert them during load. 
		boolean validCalibration = isUncalibrated() ||
				(calibration.hasDistanceUnit() && calibration.getDistanceUnit() == DistanceUnit.PIXEL);

		size = 0;
		lastPaintSize = 0;
		nextRepaintSize = 20; // Let some results appear before drawing
		nextPaintTime = System.currentTimeMillis() + repaintDelay;
		imageLock = false;
		data = new double[w * h];

		// Use negative zero so that we know when positive zero has been written to the array.
		if ((displayFlags & (DISPLAY_MAPPED | DISPLAY_MAP_ZERO)) == (DISPLAY_MAPPED | DISPLAY_MAP_ZERO))
			EMPTY = -0.0f;

		resetData();
		imp = WindowManager.getImage(title);
		currentFrame = 1;

		ImageProcessor ip = createNewProcessor(w, h);

		if (imp == null)
		{
			imp = new ImagePlus(title, ip);
			// Apply the fire lookup table
			WindowManager.setTempCurrentImage(imp);
			LutLoader lut = new LutLoader();
			lut.run(lutName);
			WindowManager.setTempCurrentImage(null);

			if (displayImage)
				imp.show();
		}
		else
		{
			// Copy the lookup table
			ip.setColorModel(imp.getProcessor().getColorModel());
			ImageStack stack = createNewImageStack(w, h);
			stack.addSlice(null, ip);
			// If resizing then remove adornments
			if (stack.getWidth() != imp.getWidth() || stack.getHeight() != imp.getHeight())
			{
				imp.setOverlay(null);
				imp.setRoi((Roi) null);
			}
			imp.setStack(stack);

			if (displayImage)
				imp.show();
			else
				imp.hide();
		}

		imp.setProperty("Info", createInfo());

		if (calibration != null)
		{
			Calibration cal = new Calibration();
			
			// This assumes the input data is in pixels
			
			String unit = "nm";
			double unitPerPixel = calibration.getNmPerPixel() / scale;
			if (unitPerPixel > 100)
			{
				unit = "um";
				unitPerPixel /= 1000.0;
			}
			cal.setUnit(unit);
			cal.pixelHeight = cal.pixelWidth = unitPerPixel;
			imp.setCalibration(cal);
		}

		// We cannot draw anything with no bounds or not in pixels
		imageActive = validBounds && validCalibration;
	}

	private String createInfo()
	{
		StringBuilder sb = new StringBuilder();
		if (source != null)
			sb.append("Source: ").append(source.toXML()).append("\n");
		if (bounds != null)
			sb.append("Bounds: ").append(getBoundsString()).append("\n");
		if (calibration != null)
			sb.append("Calibration:\n").append(XmlUtils.toXML(calibration)).append("\n");
		if (configuration != null)
			sb.append("Configuration:\n").append(configuration).append("\n");
		return (sb.length() > 0) ? sb.toString() : null;
	}

	/**
	 * Check the display flags when {@link #begin()} is called to ensure the image settings are OK. Update the flags if
	 * necessary.
	 * <p>
	 * Use to perform any other processing before begin().
	 */
	protected void preBegin()
	{
		// Signal is OK to be equalised
		if ((displayFlags & DISPLAY_SIGNAL) != 0)
		{
		}
		// Peak and localisation should not use equalisation
		else
		{
			displayFlags &= ~DISPLAY_EQUALIZED;
		}

		// Display peaks cannot use weighting and should show the exact frame number so use replace
		if ((displayFlags & DISPLAY_PEAK) != 0)
		{
			displayFlags &= ~DISPLAY_WEIGHTED;
			displayFlags |= DISPLAY_REPLACE;
		}

		// Mapped values (above zero) cannot use equalisation or be negative
		if ((displayFlags & DISPLAY_MAPPED) != 0)
		{
			displayFlags &= ~DISPLAY_EQUALIZED;
			displayFlags &= ~DISPLAY_NEGATIVES;
		}
	}

	private ImageProcessor createNewProcessor(int imageWidth, int imageHeight)
	{
		// Equalised display requires a 16-bit image to allow fast processing of the histogram 
		if ((displayFlags & DISPLAY_EQUALIZED) != 0)
		{
			pixels = new short[data.length];
			return new ShortProcessor(imageWidth, imageHeight, (short[]) pixels, null);
		}
		else
		{
			pixels = new float[data.length];

			// Special float processor that maps all values to 1-255 in the LUT.
			// Zero is mapped to 0 in the LUT.
			if ((displayFlags & DISPLAY_MAPPED) != 0)
			{
				MappedFloatProcessor fp = new MappedFloatProcessor(imageWidth, imageHeight, (float[]) pixels, null);
				fp.setMapZero((displayFlags & DISPLAY_MAP_ZERO) != 0);
				return fp;
			}

			return new FloatProcessor(imageWidth, imageHeight, (float[]) pixels, null);
		}
	}

	private ImageStack createNewImageStack(int w, int h)
	{
		if ((displayFlags & DISPLAY_MAPPED) != 0)
		{
			MappedImageStack stack = new MappedImageStack(w, h);
			stack.setMapZero((displayFlags & DISPLAY_MAP_ZERO) != 0);
			return stack;
		}
		return new ImageStack(w, h);
	}

	/**
	 * Create the image from a clone of the current data. Should only be called by one thread which has the lock so can
	 * use class variables and the actual pixel buffer.
	 * 
	 * @return The size when the image data was cloned
	 */
	private void createImage()
	{
		double[] data;

		synchronized (this.data)
		{
			data = this.data.clone();
			lastPaintSize = this.size;
			setNextRepaintSize(lastPaintSize);
			if (repaintDelay != 0)
				nextPaintTime = System.currentTimeMillis() + repaintDelay;
		}

		if ((displayFlags & DISPLAY_EQUALIZED) != 0)
		{
			// 16-bit image

			// Get the current maximum
			double max = data[0];
			for (int i = 1; i < data.length; i++)
			{
				if (max < data[i])
					max = data[i];
			}

			// Compress into 16-bit image if necessary
			int K = 65535;
			double norm = K / max;

			short[] pixels = (short[]) this.pixels;

			for (int i = 0; i < pixels.length; i++)
			{
				int index = (int) (norm * data[i]);
				if (index > K)
					index = K;
				pixels[i] = (short) index;
			}

			// Get the histogram
			int[] H = new int[K + 1];
			for (int i = 0; i < pixels.length; i++)
				H[pixels[i] & 0xffff]++;

			// Skip empty data
			int start = 1;
			while (H[start] == 0 && start < K)
				start++;
			//System.out.printf("Start = %d\n", start);

			// Perform weighted histogram equalisation
			// See: ij.plugin.ContrastEnhancer
			double[] sqrt = new double[H.length];
			sqrt[start] = Math.sqrt(H[start]);
			double sum = sqrt[start];
			for (int i = start + 1; i < K; i++)
			{
				sqrt[i] = Math.sqrt(H[i]);
				sum += 2 * sqrt[i];
			}
			sum += Math.sqrt(H[K]);

			double scale = K / sum;

			int[] lut = new int[K + 1];

			lut[0] = 0;
			sum = sqrt[start];
			for (int i = start + 1; i < K; i++)
			{
				double delta = sqrt[i];
				sum += delta;
				lut[i] = (int) (sum * scale + 0.5);
				sum += delta;
			}
			lut[K] = K;

			for (int i = 0; i < pixels.length; i++)
				pixels[i] = (short) lut[pixels[i] & 0xffff];

			imp.setDisplayRange(0, K);
		}
		else
		{
			// 32-bit image. Just copy the data but find the maximum
			float[] pixels = (float[]) this.pixels;
			double max = data[0];
			double min = 0;
			for (int i = 0; i < data.length; i++)
			{
				if (max < data[i])
					max = data[i];
				pixels[i] = (float) data[i];
			}

			if ((displayFlags & DISPLAY_NEGATIVES) != 0)
			{
				for (float f : pixels)
					if (min > f)
						min = f;
			}

			imp.setDisplayRange(min, max);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.AbstractPeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsDev)
	{
		if (!imageActive)
			return;

		final float x = mapX(params[PeakResult.X]);
		final float y = mapY(params[PeakResult.Y]);

		// Check bounds
		if (x < 0 || x >= imageWidth || y < 0 || y >= imageHeight)
			return;

		checkAndUpdateToFrame(peak);

		int[] indices = new int[5];
		float[] values = new float[4];

		getValue(peak, params, error, x, y, indices, values);

		addData(1, indices[4], indices, values);

		updateImage();
	}

	private void addData(int nPoints, int nValues, int[] indices, float[] values)
	{
		// Add the values to the configured indices
		synchronized (data)
		{
			size += nPoints;

			if ((displayFlags & DISPLAY_REPLACE) != 0)
			{
				// Replace the data
				for (int i = nValues; i-- > 0;)
					data[indices[i]] = values[i];
			}
			else if ((displayFlags & DISPLAY_MAX) != 0)
			{
				// Use the highest value
				for (int i = nValues; i-- > 0;)
					data[indices[i]] = max(data[indices[i]], values[i]);
			}
			else
			{
				// Add the data
				for (int i = nValues; i-- > 0;)
					data[indices[i]] += values[i];
			}
		}
	}

	private static double max(final double a, final double b)
	{
		// Ignore possible NaNs or infinity

		return (a > b) ? a : b;
	}

	/**
	 * Map x to the location on the output image.
	 *
	 * @param x
	 *            the x
	 * @return the output x
	 */
	public float mapX(float x)
	{
		return (x - bounds.x) * scale;
	}

	/**
	 * Map y to the location on the output image.
	 *
	 * @param y
	 *            the y
	 * @return the output y
	 */
	public float mapY(float y)
	{
		return (y - bounds.y) * scale;
	}

	/**
	 * Gets the value to add to the image data.
	 * <p>
	 * We construct the indices based on the current settings. 1, 2, or 4 indices will be returned with their values.
	 * The number of indices is stored in the 5th position of the indices array.
	 *
	 * @param peak
	 *            the peak
	 * @param params
	 *            the peak params
	 * @param error
	 *            the peak error
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param indices
	 *            the indices
	 * @param value
	 *            the values for the indices
	 * @return the value
	 */
	private void getValue(int peak, float[] params, double error, float x, float y, int[] indices, float[] value)
	{
		final float v;

		// Use the signal for the count
		if ((displayFlags & DISPLAY_SIGNAL) != 0)
		{
			v = params[PeakResult.INTENSITY];
		}
		// Use the peak number for the count
		else if ((displayFlags & DISPLAY_PEAK) != 0)
		{
			v = peak;
		}
		// Use the peak number for the count
		else if ((displayFlags & DISPLAY_ERROR) != 0)
		{
			v = (float) error;
		}
		else
		{
			v = 1;
		}

		getValue(v, x, y, indices, value);
	}

	/**
	 * Gets the value to add to the image data.
	 * <p>
	 * We construct the indices based on the current settings. 1, 2, or 4 indices will be returned with their values.
	 * The number of indices is stored in the 5th position of the indices array.
	 *
	 * @param v
	 *            the value
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param indices
	 *            the indices
	 * @param values
	 *            the values for the indices
	 */
	private void getValue(float v, float x, float y, int[] indices, float[] values)
	{
		final int x1 = (int) x;
		final int y1 = (int) y;
		final int index = y1 * imageWidth + x1;

		if ((displayFlags & DISPLAY_WEIGHTED) == 0)
		{
			// No interpolation. Just put the value on the containing pixel
			indices[0] = index;
			values[0] = v;
			indices[4] = 1;
			return;
		}

		// Note: It is very unlikely that dx and dy will be 0.5f so we ignore this case for speed.
		// It could be added later to test if speed is impacted since we return the number of indices.
		// If a user wants to add data only to one pixel then they can remove the weighted option. 

		indices[4] = 4;

		// Use bilinear weighting

		final float dx = x - x1;
		final float dy = y - y1;

		final float wx; // X weight for the location pixel
		final float wy; // Y weight for the location pixel

		// Get the 4 neighbours and avoid overrun. In this case the edge pixel will get the entire value.
		final int xDelta, yDelta;

		// Note: The image width/height could be zero making the deltas invalid. However in this case the 
		// getValue(...) method will never be called.

		if (dx < 0.5f)
		{
			// Interpolate to the lower x pixel
			wx = 0.5f + dx;
			if (x1 == 0)
			{
				xDelta = 0;
			}
			else
			{
				xDelta = -1;
			}
		}
		else
		{
			// Interpolate to the upper x pixel
			wx = 1.5f - dx;
			if (x1 == xlimit)
			{
				xDelta = 0;
			}
			else
			{
				xDelta = 1;
			}
		}

		if (dy < 0.5f)
		{
			// Interpolate to the lower y pixel
			wy = 0.5f + dy;
			if (y1 == 0)
			{
				yDelta = 0;
			}
			else
			{
				yDelta = -imageWidth;
			}
		}
		else
		{
			// Interpolate to the upper y pixel
			wy = 1.5f - dy;
			if (y1 == ylimit)
			{
				yDelta = 0;
			}
			else
			{
				yDelta = imageWidth;
			}
		}

		indices[0] = index;
		indices[1] = index + xDelta;
		indices[2] = index + yDelta;
		indices[3] = index + xDelta + yDelta;

		final float wxDelta = 1f - wx;
		final float wyDelta = 1f - wy;

		values[0] = v * wx * wy;
		values[1] = v * wxDelta * wy;
		values[2] = v * wx * wyDelta;
		values[3] = v * wxDelta * wyDelta;
	}

	/**
	 * Simplified method to allow the image to be reconstructed using just T,X,Y coordinates and a value
	 * 
	 * @param peak
	 *            The peak frame
	 * @param x
	 *            The X coordinate
	 * @param y
	 *            The Y coordinate
	 * @param v
	 *            The value
	 */
	public void add(int peak, float x, float y, float v)
	{
		if (!imageActive)
			return;

		x = mapX(x);
		y = mapY(y);

		// Check bounds
		if (x < 0 || x >= imageWidth || y < 0 || y >= imageHeight)
			return;

		checkAndUpdateToFrame(peak);

		int[] indices = new int[5];
		float[] values = new float[4];

		getValue(v, x, y, indices, values);

		addData(1, indices[4], indices, values);

		updateImage();
	}

	/**
	 * Simplified method to allow the image to be reconstructed using just X,Y coordinates and a value
	 * 
	 * @param x
	 *            The X coordinate
	 * @param y
	 *            The Y coordinate
	 * @param v
	 *            The value
	 */
	public void add(float x, float y, float v)
	{
		if (!imageActive)
			return;

		x = mapX(x);
		y = mapY(y);

		// Check bounds
		if (x < 0 || x >= imageWidth || y < 0 || y >= imageHeight)
			return;

		int[] indices = new int[5];
		float[] values = new float[4];

		getValue(v, x, y, indices, values);

		addData(1, indices[4], indices, values);

		updateImage();
	}

	/**
	 * Simplified method to allow the image to be reconstructed using just T,X,Y coordinates and a value
	 * 
	 * @param allpeak
	 *            The peak frames
	 * @param allx
	 *            The X coordinates
	 * @param ally
	 *            The Y coordinates
	 * @param allv
	 *            The values
	 */
	public void add(int[] allpeak, float[] allx, float[] ally, float[] allv)
	{
		if (!imageActive)
			return;

		int[] indices = new int[5];
		float[] values = new float[4];

		int nPoints = 0;
		int nValues = 0;

		// Buffer output in batches
		int[] allIndices = new int[100];
		float[] allValues = new float[allIndices.length];

		boolean replace = ((displayFlags & DISPLAY_REPLACE) != 0);

		// We add at most 4 indices for each peak
		int limit = allIndices.length - 4;

		for (int j = 0; j < allx.length; j++)
		{
			float x = mapX(allx[j]);
			float y = mapY(ally[j]);
			int peak = allpeak[j];

			// Check bounds
			if (x < 0 || x >= imageWidth || y < 0 || y >= imageHeight)
				continue;

			if (shouldUpdate(peak))
			{
				addData(nPoints, nValues, allIndices, allValues);
				nPoints = 0;
				nValues = 0;
				updateToFrame(peak);
			}

			getValue(allv[j], x, y, indices, values);

			for (int i = indices[4]; i-- > 0;)
			{
				allIndices[nValues] = indices[i];
				allValues[nValues] = values[i];
				nValues++;
			}

			nPoints++;

			if (nValues > limit || replace)
			{
				addData(nPoints, nValues, allIndices, allValues);
				nPoints = 0;
				nValues = 0;
				updateImage();
				if (!imageActive)
					return;
			}
		}

		// Now add the values to the configured indices
		addData(nPoints, nValues, allIndices, allValues);

		updateImage();
	}

	/**
	 * Simplified method to allow the image to be reconstructed using just T,X,Y coordinates and a value
	 * 
	 * @param allx
	 *            The X coordinates
	 * @param ally
	 *            The Y coordinates
	 * @param allv
	 *            The values
	 */
	public void add(float[] allx, float[] ally, float[] allv)
	{
		if (!imageActive)
			return;

		int[] indices = new int[5];
		float[] values = new float[4];

		int nPoints = 0;
		int nValues = 0;

		// Buffer output in batches
		int[] allIndices = new int[100];
		float[] allValues = new float[allIndices.length];

		boolean replace = ((displayFlags & DISPLAY_REPLACE) != 0);

		// We add at most 4 indices for each peak
		int limit = allIndices.length - 4;

		for (int j = 0; j < allx.length; j++)
		{
			float x = mapX(allx[j]);
			float y = mapY(ally[j]);

			// Check bounds
			if (x < 0 || x >= imageWidth || y < 0 || y >= imageHeight)
				continue;

			getValue(allv[j], x, y, indices, values);

			for (int i = indices[4]; i-- > 0;)
			{
				allIndices[nValues] = indices[i];
				allValues[nValues] = values[i];
				nValues++;
			}

			nPoints++;

			if (nValues > limit || replace)
			{
				addData(nPoints, nValues, allIndices, allValues);
				nPoints = 0;
				nValues = 0;
				updateImage();
				if (!imageActive)
					return;
			}
		}

		// Now add the values to the configured indices
		addData(nPoints, nValues, allIndices, allValues);

		updateImage();
	}

	/**
	 * Check if the stack should be updated to move the rolling window to the given peak.
	 * 
	 * @param peak
	 * @return True if update is required
	 */
	protected boolean shouldUpdate(int peak)
	{
		return (rollingWindowSize > 0) && (peak >= currentFrame + rollingWindowSize);
	}

	/**
	 * Add frames to the current stack to ensure the rolling window size is enforced (i.e. batches of N are drawn as a
	 * frame)
	 * 
	 * @param peak
	 */
	protected void checkAndUpdateToFrame(int peak)
	{
		if (shouldUpdate(peak))
			updateToFrame(peak);
	}

	/**
	 * Add frames to the current stack to ensure the rolling window size is enforced (i.e. batches of N are drawn as a
	 * frame)
	 * 
	 * @param peak
	 */
	protected void updateToFrame(int peak)
	{
		synchronized (data) // Stop other threads adding more data
		{
			int i = 0;
			ImageStack stack = imp.getStack();
			peak -= rollingWindowSize;
			while (peak >= currentFrame)
			{
				//System.out.printf("%d => %d (%d)\n", currentFrame, currentFrame + rollingWindowSize, peak+rollingWindowSize);
				if (i++ == 0)
				{
					// Draw all current data for first time we move forward.
					// Force repaint 
					forceUpdateImage();
				}

				ImageProcessor ip = createNewProcessor(stack.getWidth(), stack.getHeight());
				stack.addSlice(null, ip);
				currentFrame += rollingWindowSize;
			}

			// Check if any frames were added
			if (i > 0)
			{
				imp.setStack(stack);
				imp.setSlice(stack.getSize());

				resetData();
			}
		}
	}

	private void resetData()
	{
		Arrays.fill(data, EMPTY);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResults#add(gdsc.smlm.results.PeakResult)
	 */
	public void add(PeakResult result)
	{
		add(result.getFrame(), result.origX, result.origY, result.origValue, result.error, result.noise,
				result.getParameters(), null);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.AbstractPeakResults#addAll(java.util.Collection)
	 */
	public void addAll(Collection<PeakResult> results)
	{
		addAll(results.toArray(new PeakResult[results.size()]));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResults#addAll(gdsc.smlm.results.PeakResult[])
	 */
	public void addAll(PeakResult[] results)
	{
		if (!imageActive)
			return;

		int[] indices = new int[5];
		float[] values = new float[4];

		int nPoints = 0;
		int nValues = 0;

		// Buffer output in batches
		int[] allIndices = new int[100];
		float[] allValues = new float[allIndices.length];

		boolean replace = ((displayFlags & DISPLAY_REPLACE) != 0);

		// We add at most 4 indices for each peak
		int limit = allIndices.length - 4;

		for (PeakResult result : results)
		{
			float x = mapX(result.getXPosition());
			float y = mapY(result.getYPosition());

			// Check bounds
			if (x < 0 || x >= imageWidth || y < 0 || y >= imageHeight)
				continue;

			if (shouldUpdate(result.getFrame()))
			{
				addData(nPoints, nValues, allIndices, allValues);
				nPoints = 0;
				nValues = 0;
				updateToFrame(result.getFrame());
			}

			getValue(result.getFrame(), result.getParameters(), result.error, x, y, indices, values);

			for (int i = indices[4]; i-- > 0;)
			{
				allIndices[nValues] = indices[i];
				allValues[nValues] = values[i];
				nValues++;
			}

			nPoints++;

			if (nValues > limit || replace)
			{
				addData(nPoints, nValues, allIndices, allValues);
				nPoints = 0;
				nValues = 0;
				updateImage();
				if (!imageActive)
					return;
			}
		}

		// Now add the values to the configured indices
		addData(nPoints, nValues, allIndices, allValues);

		updateImage();
	}

	protected void updateImage()
	{
		if (size < nextRepaintSize || !liveImage || !displayImage)
			return;

		if (!imp.isVisible())
		{
			//System.out.println("Image has been closed");
			imageActive = false;
			return;
		}

		if (repaintDelay != 0)
		{
			long time = System.currentTimeMillis();

			if (time < nextPaintTime)
			{
				// Get the amount of time it took to acquire the data
				int n = size - lastPaintSize;
				long elapsed = time - (nextPaintTime - repaintDelay);

				if (elapsed > 0)
				{
					// Set the next repaint size using linear scaling
					long remaining = nextPaintTime - time;
					int extra = (int) (n * ((double) remaining / elapsed));
					//System.out.printf("Updating next paint size: %d : %d -> %d\n", lastPaintSize, nextRepaintSize,
					//		nextRepaintSize + extra);
					nextRepaintSize += extra;
				}
				else
				{
					setNextRepaintSize(size);
				}
				return;
			}
		}

		drawImage();
	}

	private void setNextRepaintSize(int size)
	{
		nextRepaintSize = (int) Math.ceil(size + size * repaintInterval);
	}

	/**
	 * This forces all the current data to be written to the image. It is used when a rolling window is being drawn.
	 */
	private void forceUpdateImage()
	{
		if (!imp.isVisible())
		{
			imageActive = false;
			return;
		}

		drawImage();
	}

	private void drawImage()
	{
		if (aquireLock())
		{
			try
			{
				createImage();

				// We direct manipulate the pixel buffer so this is not necessary
				//ImageProcessor ip = imp.getProcessor();
				//ip.setPixels(newPixels);
				//imp.setProcessor(ip);

				imp.updateAndDraw();
			}
			finally
			{
				releaseLock();
			}
		}
	}

	private synchronized boolean aquireLock()
	{
		if (imageLock)
			return false;

		imageLock = true;
		return true;
	}

	private synchronized void releaseLock()
	{
		imageLock = false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#size()
	 */
	public int size()
	{
		return size;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#end()
	 */
	public void end()
	{
		// Wait for previous image to finish rendering
		while (!aquireLock())
		{
			try
			{
				System.out.printf("Waiting for final image\n");
				Thread.sleep(50);
			}
			catch (InterruptedException e)
			{
				// Ignore
			}
		}

		releaseLock();

		drawImage();

		if (rollingWindowSize > 0)
		{
			imp.setDimensions(1, 1, imp.getStackSize());
		}

		imageActive = false;
	}

	/**
	 * Image will be repainted when the size is increased by a fraction of the last size painted.
	 * 
	 * @param repaintInterval
	 *            the repaintInterval to set
	 */
	public void setRepaintInterval(double repaintInterval)
	{
		if (repaintInterval < 0)
			repaintInterval = 0;
		this.repaintInterval = repaintInterval;
	}

	/**
	 * @return the repaintInterval
	 */
	public double getRepaintInterval()
	{
		return repaintInterval;
	}

	/**
	 * Sets the repaint delay (time in milliseconds that must elapse before a repaint).
	 *
	 * @param repaintDelay
	 *            the new repaint delay
	 */
	public void setRepaintDelay(long repaintDelay)
	{
		if (repaintDelay < 0)
			repaintDelay = 0;
		this.repaintDelay = repaintDelay;
	}

	/**
	 * Gets the repaint delay.
	 *
	 * @return the repaint delay
	 */
	public long getRepaintDelay()
	{
		return repaintDelay;
	}

	/**
	 * @param displayFlags
	 *            the displayFlags to set
	 */
	public void setDisplayFlags(int displayFlags)
	{
		this.displayFlags = displayFlags;
	}

	/**
	 * @return the displayFlags
	 */
	public int getDisplayFlags()
	{
		return displayFlags;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		return imageActive;
	}

	/**
	 * @return the rollingWindowSize
	 */
	public int getRollingWindowSize()
	{
		return rollingWindowSize;
	}

	/**
	 * Produce a final output image as a stack. Specify the number of peak frames to combine into each stack frame, e.g.
	 * a window size of 10 will combine 10 consecutive fitting frames into 1 plane.
	 * <p>
	 * This setting only applies before the {@link #begin()} method.
	 * 
	 * @param rollingWindowSize
	 *            the rollingWindowSize to set
	 */
	public void setRollingWindowSize(int rollingWindowSize)
	{
		this.rollingWindowSize = rollingWindowSize;
	}

	/**
	 * Over-ridden to ignore any passed in bounds. The bounds must be set when the image is created.
	 * 
	 * @see gdsc.smlm.results.AbstractPeakResults#setBounds(java.awt.Rectangle)
	 */
	public void setBounds(Rectangle bounds)
	{
		// Ignore. Bounds are only valid when the image is created
	}

	/**
	 * @return True if the image should be displayed
	 */
	public boolean isDisplayImage()
	{
		return displayImage;
	}

	/**
	 * Set to true if the image should be displayed. Should be called before {@link #begin()}
	 * 
	 * @param displayImage
	 *            Set to true if the image should be displayed.
	 */
	public void setDisplayImage(boolean displayImage)
	{
		this.displayImage = displayImage;
	}

	/**
	 * @return The IJ image
	 */
	public ImagePlus getImagePlus()
	{
		return imp;
	}

	/**
	 * @return the lutName
	 */
	public String getLutName()
	{
		return lutName;
	}

	/**
	 * @param lutName
	 *            the lutName to set
	 */
	public void setLutName(String lutName)
	{
		this.lutName = lutName;
	}

	/**
	 * Checks if is live image. If true the image will be updated as data is added.
	 *
	 * @return true, if is live image
	 */
	public boolean isLiveImage()
	{
		return liveImage;
	}

	/**
	 * Sets the live image flag. Set to true and the image will be updated as data is added.
	 *
	 * @param liveImage
	 *            the new live image flag
	 */
	public void setLiveImage(boolean liveImage)
	{
		this.liveImage = liveImage;
	}

	/**
	 * Checks if is uncalibrated flag. An uncalibrated image does not require the calibration to be in pixel units.
	 *
	 * @return true, if is uncalibrated
	 */
	public boolean isUncalibrated()
	{
		return uncalibrated;
	}

	/**
	 * Sets the uncalibrated flag. An uncalibrated image does not require the calibration to be in pixel units.
	 *
	 * @param uncalibrated
	 *            the new uncalibrated flag.
	 */
	public void setUncalibrated(boolean uncalibrated)
	{
		this.uncalibrated = uncalibrated;
	}
}
