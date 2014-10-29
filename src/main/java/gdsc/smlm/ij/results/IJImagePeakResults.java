package gdsc.smlm.ij.results;

import gdsc.smlm.results.PeakResult;
import gdsc.smlm.utils.XmlUtils;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.measure.Calibration;
import ij.plugin.LutLoader;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Collection;

import org.apache.commons.math3.util.FastMath;

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
	public static final int DISPLAY_SIGNAL = 1;
	public static final int DISPLAY_WEIGHTED = 2;
	public static final int DISPLAY_EQUALIZED = 4;
	public static final int DISPLAY_PEAK = 8;
	public static final int DISPLAY_ERROR = 16;

	public static final int DISPLAY_REPLACE = 32;
	public static final int DISPLAY_MAX = 64;

	public static final int DISPLAY_NEGATIVES = 128;

	protected String title;
	protected int imageWidth;
	protected int imageHeight;
	protected float scale;
	protected int size = 0;
	protected float[] data;
	protected float xlimit;
	protected float ylimit;
	protected ImagePlus imp = null;
	protected boolean imageActive = false;
	protected int displayFlags = 0;
	private int rollingWindowSize = 0;
	private boolean displayImage = true;

	// Used to draw the image
	private int nextRepaintSize = 0;
	private Object pixels;
	private boolean imageLock = false;
	private double repaintInterval = 0.1;
	private int currentFrame;

	private String lutName = "fire";

	/**
	 * @param title
	 *            Title of the image (appended with a suffix)
	 * @param bounds
	 *            Define the bounding rectangle of the image coordinates. Any results outside this will not be
	 *            displayed.
	 * @param scale
	 */
	public IJImagePeakResults(String title, Rectangle bounds, float scale)
	{
		this.title = title + " " + IMAGE_SUFFIX;
		this.bounds = bounds;
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

		size = 0;
		nextRepaintSize = 20; // Let some results appear before drawing
		imageLock = false;
		data = new float[imageWidth * imageHeight];
		imp = WindowManager.getImage(title);
		currentFrame = 1;

		ImageProcessor ip = createNewProcessor();

		if (imp == null)
		{
			imp = new ImagePlus(title, ip);
			// Apply the fire lookup table
			WindowManager.setTempCurrentImage(imp);
			LutLoader lut = new LutLoader();
			lut.run(lutName);
			WindowManager.setTempCurrentImage(null);
		}
		else
		{
			// Copy the lookup table
			ip.setColorModel(imp.getProcessor().getColorModel());
			ImageStack stack = new ImageStack(imageWidth, imageHeight);
			stack.addSlice(null, ip);
			imp.setStack(stack);
		}

		imp.setProperty("Info", createInfo());

		if (displayImage)
			imp.show();
		else
			imp.hide();

		if (calibration != null)
		{
			Calibration cal = new Calibration();
			String unit = "nm";
			double unitPerPixel = calibration.nmPerPixel / scale;
			if (unitPerPixel > 100)
			{
				unit = "um";
				unitPerPixel /= 1000.0;
			}
			cal.setUnit(unit);
			cal.pixelHeight = cal.pixelWidth = unitPerPixel;
			imp.setCalibration(cal);
		}

		imageActive = true;
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
	}

	private ImageProcessor createNewProcessor()
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
			return new FloatProcessor(imageWidth, imageHeight, (float[]) pixels, null);
		}
	}

	/**
	 * Create the image from the current data. Should only be called by one thread which has the lock so can use class
	 * variables and the actual pixel buffer.
	 * 
	 * @return
	 */
	private void createImage()
	{
		if ((displayFlags & DISPLAY_EQUALIZED) != 0)
		{
			// 16-bit image

			// Get the current maximum
			float max = 0;
			for (int i = 0; i < data.length; i++)
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
			float max = 0;
			float min = 0;
			for (int i = 0; i < data.length; i++)
			{
				if (max < data[i])
					max = data[i];
				pixels[i] = data[i];
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

		float x = (params[3] - bounds.x) * scale;
		float y = (params[4] - bounds.y) * scale;

		// Check bounds
		if (x < 0 || x > xlimit || y < 0 || y > ylimit)
			return;

		checkAndUpdateToFrame(peak);

		int x1 = (int) x;
		int y1 = (int) y;

		float[] value = getValue(peak, params, error, x, y, x1, y1);

		int index = y1 * imageWidth + x1;

		// Now add the values to the configured indices
		synchronized (data)
		{
			size++;
			if ((displayFlags & DISPLAY_REPLACE) != 0)
			{
				// Replace the data
				data[index] = value[0];
				data[index + imageWidth] = value[1];
				data[index + 1] = value[2];
				data[index + imageWidth + 1] = value[3];
			}
			else if ((displayFlags & DISPLAY_MAX) != 0)
			{
				// Use the highest value
				data[index] = FastMath.max(data[index], value[0]);
				data[index + imageWidth] = FastMath.max(data[index + imageWidth], value[1]);
				data[index + 1] = FastMath.max(data[index + 1], value[2]);
				data[index + imageWidth + 1] = FastMath.max(data[index + imageWidth + 1], value[3]);
			}
			else
			{
				// Add the data
				data[index] += value[0];
				data[index + imageWidth] += value[1];
				data[index + 1] += value[2];
				data[index + imageWidth + 1] += value[3];
			}
		}

		updateImage();
	}
	
	private float[] getValue(int peak, float[] params, double error, float x, float y, int x1, int y1)
	{
		// Add a count to each adjacent pixel
		float[] value = new float[] { 1, 1, 1, 1 };

		// Use the signal for the count
		if ((displayFlags & DISPLAY_SIGNAL) != 0)
		{
			float signal = PeakResult.getSignal(params);
			for (int i = 0; i < 4; i++)
				value[i] = signal;
		}
		// Use the peak number for the count
		else if ((displayFlags & DISPLAY_PEAK) != 0)
		{
			for (int i = 0; i < 4; i++)
				value[i] = peak;
		}
		// Use the peak number for the count
		else if ((displayFlags & DISPLAY_ERROR) != 0)
		{
			for (int i = 0; i < 4; i++)
				value[i] = (float) error;
		}

		float wx, wy;

		// Use bilinear weighting
		if ((displayFlags & DISPLAY_WEIGHTED) != 0)
		{
			wx = x - x1;
			wy = y - y1;
		}
		else
		{
			// Put the value on the nearest pixel by rounding the weights.
			wx = Math.round(x - x1);
			wy = Math.round(y - y1);
		}

		applyWeights(value, wx, wy);

		return value;
	}

	private void applyWeights(float[] value, float wx, float wy)
	{
		value[0] *= (1f - wx) * (1f - wy);
		value[1] *= (1f - wx) * wy;
		value[2] *= wx * (1f - wy);
		value[3] *= wx * wy;
	}
	
	/**
	 * Simplified method to allow the Image to be reconstructed using just T,X,Y coordinates and a value
	 * @param peak The peak frame
	 * @param x The X coordinate
	 * @param y The Y coordinate
	 * @param v The value
	 */
	public void add(int peak, float x, float y, float v)
	{
		if (!imageActive)
			return;

		x = (x - bounds.x) * scale;
		y = (y - bounds.y) * scale;

		// Check bounds
		if (x < 0 || x > xlimit || y < 0 || y > ylimit)
			return;

		checkAndUpdateToFrame(peak);

		int x1 = (int) x;
		int y1 = (int) y;

		float[] value = getValue(peak, v, x, y, x1, y1);

		int index = y1 * imageWidth + x1;

		// Now add the values to the configured indices
		synchronized (data)
		{
			size++;
			if ((displayFlags & DISPLAY_REPLACE) != 0)
			{
				// Replace the data
				data[index] = value[0];
				data[index + imageWidth] = value[1];
				data[index + 1] = value[2];
				data[index + imageWidth + 1] = value[3];
			}
			else if ((displayFlags & DISPLAY_MAX) != 0)
			{
				// Use the highest value
				data[index] = FastMath.max(data[index], value[0]);
				data[index + imageWidth] = FastMath.max(data[index + imageWidth], value[1]);
				data[index + 1] = FastMath.max(data[index + 1], value[2]);
				data[index + imageWidth + 1] = FastMath.max(data[index + imageWidth + 1], value[3]);
			}
			else
			{
				// Add the data
				data[index] += value[0];
				data[index + imageWidth] += value[1];
				data[index + 1] += value[2];
				data[index + imageWidth + 1] += value[3];
			}
		}

		updateImage();
	}
	
	private float[] getValue(int peak, float v, float x, float y, int x1, int y1)
	{
		// Add a count to each adjacent pixel
		float[] value = new float[] { v, v, v, v };

		float wx, wy;

		// Use bilinear weighting
		if ((displayFlags & DISPLAY_WEIGHTED) != 0)
		{
			wx = x - x1;
			wy = y - y1;
		}
		else
		{
			// Put the value on the nearest pixel by rounding the weights.
			wx = Math.round(x - x1);
			wy = Math.round(y - y1);
		}

		applyWeights(value, wx, wy);

		return value;
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
					nextRepaintSize = 0; // Force repaint 
					updateImage(); // Draw all current data for first frame
				}

				ImageProcessor ip = createNewProcessor();
				stack.addSlice(null, ip);
				currentFrame += rollingWindowSize;
			}

			// Check if any frames were added
			if (i > 0)
			{
				imp.setStack(stack);
				imp.setSlice(stack.getSize());

				// Reset the data
				Arrays.fill(data, 0);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#addAll(java.util.Collection)
	 */
	public void addAll(Collection<PeakResult> results)
	{
		if (!imageActive)
			return;

		int nPoints = 0;

		// Buffer output in batches
		int[] indices = new int[25 * 4];
		float[] values = new float[indices.length];

		for (PeakResult result : results)
		{
			float x = (result.params[3] - bounds.x) * scale;
			float y = (result.params[4] - bounds.y) * scale;

			// Check bounds
			if (x < 0 || x > xlimit || y < 0 || y > ylimit)
				continue;

			if (shouldUpdate(result.peak))
			{
				addData(nPoints, indices, values);
				nPoints = 0;
				updateToFrame(result.peak);
			}

			int x1 = (int) x;
			int y1 = (int) y;

			float[] value = getValue(result.peak, result.params, result.error, x, y, x1, y1);

			int index = y1 * imageWidth + x1;

			indices[nPoints] = index;
			indices[nPoints + 1] = index + imageWidth;
			indices[nPoints + 2] = index + 1;
			indices[nPoints + 3] = index + imageWidth + 1;
			values[nPoints] = value[0];
			values[nPoints + 1] = value[1];
			values[nPoints + 2] = value[2];
			values[nPoints + 3] = value[3];

			nPoints += 4;

			if (nPoints >= indices.length)
			{
				addData(nPoints, indices, values);
				nPoints = 0;
				updateImage();
			}
		}

		// Now add the values to the configured indices
		addData(nPoints, indices, values);

		updateImage();
	}

	private void addData(int nPoints, int[] indices, float[] values)
	{
		// Add the values to the configured indices
		synchronized (data)
		{
			size += nPoints / 4;

			if ((displayFlags & DISPLAY_REPLACE) != 0)
			{
				// Replace the data
				for (int i = 0; i < nPoints; i++)
					data[indices[i]] = values[i];
			}
			else if ((displayFlags & DISPLAY_MAX) != 0)
			{
				// Use the highest value
				for (int i = 0; i < nPoints; i++)
					data[indices[i]] = FastMath.max(data[indices[i]], values[i]);
			}
			else
			{
				// Add the data
				for (int i = 0; i < nPoints; i++)
					data[indices[i]] += values[i];
			}
		}
	}

	protected void updateImage()
	{
		if (size < nextRepaintSize || !imageActive || !displayImage)
			return;

		if (!imp.isVisible())
		{
			//System.out.println("Image has been closed");
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
				nextRepaintSize = (int) (size + size * repaintInterval);

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
	 * Image will be repainted when a fraction of new results have been added.
	 * 
	 * @param repaintInterval
	 *            the repaintInterval to set (range 0.001-1)
	 */
	public void setRepaintInterval(double repaintInterval)
	{
		if (repaintInterval < 0.001)
			repaintInterval = 0.001;
		if (repaintInterval > 1)
			repaintInterval = 1;

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
}
