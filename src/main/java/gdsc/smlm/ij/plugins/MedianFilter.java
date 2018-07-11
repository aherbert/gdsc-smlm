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
package gdsc.smlm.ij.plugins;

import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import gdsc.core.utils.MedianWindowDLLFloat;
import gdsc.core.utils.MedianWindowFloat;
import gdsc.smlm.ij.utils.IJImageConverter;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * Filters each pixel using a sliding median through the time stack. Medians are computed at set intervals and the
 * values interpolated.
 */
public class MedianFilter implements PlugInFilter
{
	private static final String TITLE = "Median Filter";
	private final int FLAGS = DOES_8G | DOES_16 | DOES_32;

	private static int radius = 50;
	private static int interval = 12;
	private static int blockSize = 32;
	private static boolean subtract = false;
	private static float bias = 500;

	ImagePlus imp;
	int counter, size;

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}
		this.imp = imp;
		return showDialog();
	}

	@Override
	public void run(ImageProcessor ip)
	{
		long start = System.currentTimeMillis();

		ImageStack stack = imp.getImageStack();

		final int width = stack.getWidth();
		final int height = stack.getHeight();
		size = width * height;
		float[][] imageStack = new float[stack.getSize()][];
		float[] mean = new float[imageStack.length];

		// Get the mean for each frame and normalise the data using the mean
		ExecutorService threadPool = Executors.newFixedThreadPool(Prefs.getThreads());
		List<Future<?>> futures = new LinkedList<>();

		counter = 0;
		IJ.showStatus("Calculating means...");
		for (int n = 1; n <= stack.getSize(); n++)
		{
			futures.add(threadPool.submit(new ImageNormaliser(stack, imageStack, mean, n)));
		}

		// Finish processing data
		Utils.waitForCompletion(futures);

		futures = new LinkedList<>();

		counter = 0;
		IJ.showStatus("Calculating medians...");
		for (int i = 0; i < size; i += blockSize)
		{
			futures.add(threadPool.submit(new ImageGenerator(imageStack, mean, i, FastMath.min(i + blockSize, size))));
		}

		// Finish processing data
		Utils.waitForCompletion(futures);

		if (Utils.isInterrupted())
			return;

		if (subtract)
		{
			counter = 0;
			IJ.showStatus("Subtracting medians...");
			for (int n = 1; n <= stack.getSize(); n++)
			{
				futures.add(threadPool.submit(new ImageFilter(stack, imageStack, n)));
			}

			// Finish processing data
			Utils.waitForCompletion(futures);
		}

		// Update the image
		ImageStack outputStack = new ImageStack(stack.getWidth(), stack.getHeight(), stack.getSize());
		for (int n = 1; n <= stack.getSize(); n++)
		{
			outputStack.setPixels(imageStack[n - 1], n);
		}

		imp.setStack(outputStack);
		imp.updateAndDraw();

		IJ.showTime(imp, start, "Completed");
		long milliseconds = System.currentTimeMillis() - start;
		Utils.log(TITLE + " : Radius %d, Interval %d, Block size %d = %s, %s / frame", radius, interval, blockSize,
				Utils.timeToString(milliseconds), Utils.timeToString((double) milliseconds / imp.getStackSize()));
	}

	private int showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage(
				"Compute the median using a rolling window at set intervals.\nBlocks of pixels are processed on separate threads.");

		gd.addSlider("Radius", 10, 100, radius);
		gd.addSlider("Interval", 10, 30, interval);
		gd.addSlider("Block_size", 1, 32, blockSize);
		gd.addCheckbox("Subtract", subtract);
		gd.addSlider("Bias", 0, 1000, bias);

		gd.showDialog();

		if (gd.wasCanceled())
			return DONE;

		radius = (int) Math.abs(gd.getNextNumber());
		interval = (int) Math.abs(gd.getNextNumber());
		blockSize = (int) Math.abs(gd.getNextNumber());
		if (blockSize < 1)
			blockSize = 1;
		subtract = gd.getNextBoolean();
		bias = (float) Math.abs(gd.getNextNumber());

		if (gd.invalidNumber() || interval < 1 || radius < 1)
			return DONE;

		// Check the window size is smaller than the stack size
		if (imp.getStackSize() < 2 * radius + 1)
		{
			IJ.error(TITLE,
					"The window size is larger than the stack size.\nThis is equal to a z-stack median projection.");
			return DONE;
		}

		return FLAGS;
	}

	private synchronized void showProgress()
	{
		IJ.showProgress(counter, size);
		counter += blockSize;
	}

	private synchronized void showProgressSingle()
	{
		IJ.showProgress(counter++, size);
	}

	/**
	 * Extract the data for a specified slice, calculate the mean and then normalise by the mean.
	 * <p>
	 * Use a runnable for the image generation to allow multi-threaded operation. Input parameters that are manipulated
	 * should have synchronized methods.
	 */
	private class ImageNormaliser implements Runnable
	{
		final ImageStack inputStack;
		final float[][] imageStack;
		final float[] mean;
		final int n;

		public ImageNormaliser(ImageStack inputStack, float[][] imageStack, float[] mean, int n)
		{
			this.inputStack = inputStack;
			this.imageStack = imageStack;
			this.mean = mean;
			this.n = n;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			showProgressSingle();

			float[] data = imageStack[n - 1] = IJImageConverter.getData(inputStack.getProcessor(n));
			double sum = 0;
			for (float f : data)
				sum += f;
			float av = mean[n - 1] = (float) (sum / data.length);
			for (int i = 0; i < data.length; i++)
				data[i] /= av;
		}
	}

	/**
	 * Compute the rolling median window on a set of pixels in the image stack, interpolating at intervals if necessary.
	 * Convert back into the final image pixel value by multiplying by the mean for the slice.
	 * <p>
	 * Use a runnable for the image generation to allow multi-threaded operation. Input parameters that are manipulated
	 * should have synchronized methods.
	 */
	private class ImageGenerator implements Runnable
	{
		final float[][] imageStack;
		final float[] mean;
		final int start, end;

		public ImageGenerator(float[][] imageStack, float[] mean, int start, int end)
		{
			this.imageStack = imageStack;
			this.mean = mean;
			this.start = start;
			this.end = end;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			if (IJ.escapePressed())
				return;
			showProgress();

			// For each pixel extract the time line of pixel data
			final int nSlices = imageStack.length;
			final int nPixels = end - start;

			if (nPixels == 1)
			{
				if (interval == 1)
				{
					// The rolling window operates effectively in linear time so use this with an interval of 1.
					// There is no need for interpolation and the data can be written directly to the output.
					final int window = 2 * radius + 1;
					float[] data = new float[window];
					for (int slice = 0; slice < window; slice++)
					{
						data[slice] = imageStack[slice][start];
					}

					// Initialise the window with the first n frames.
					MedianWindowDLLFloat mw = new MedianWindowDLLFloat(data);

					// Get the early medians.
					int slice = 0;
					for (; slice < radius; slice++)
					{
						imageStack[slice][start] = mw.getMedianOldest(slice + 1 + radius) * mean[slice];
					}

					// Then increment through the data getting the median when required.
					for (int j = mw.getSize(); j < nSlices; j++, slice++)
					{
						imageStack[slice][start] = mw.getMedian() * mean[slice];
						mw.add(imageStack[j][start]);
					}

					// Then get the later medians as required.
					for (int i = 2 * radius + 1; slice < nSlices; i--, slice++)
					{
						imageStack[slice][start] = mw.getMedianYoungest(i) * mean[slice];
					}
				}
				else
				{
					float[] data = new float[nSlices];
					for (int slice = 0; slice < nSlices; slice++)
					{
						data[slice] = imageStack[slice][start];
					}

					// Create median window filter
					MedianWindowFloat mw = new MedianWindowFloat(data.clone(), radius);

					// Produce the medians
					for (int slice = 0; slice < nSlices; slice += interval)
					{
						data[slice] = mw.getMedian();
						mw.increment(interval);
					}
					// Final position if necessary
					if (mw.getPosition() != nSlices + interval - 1)
					{
						mw.setPosition(nSlices - 1);
						data[nSlices - 1] = mw.getMedian();
					}

					// Interpolate
					for (int slice = 0; slice < nSlices; slice += interval)
					{
						int end = FastMath.min(slice + interval, nSlices - 1);
						final float increment = (data[end] - data[slice]) / (end - slice);
						for (int s = slice + 1, i = 1; s < end; s++, i++)
						{
							data[s] = data[slice] + increment * i;
						}
					}

					// Put back in the image re-scaling using the image mean
					for (int slice = 0; slice < nSlices; slice++)
					{
						imageStack[slice][start] = data[slice] * mean[slice];
					}
				}
			}
			else
			{
				if (interval == 1)
				{
					// The rolling window operates effectively in linear time so use this with an interval of 1.
					// There is no need for interpolation and the data can be written directly to the output.
					final int window = 2 * radius + 1;
					float[][] data = new float[nPixels][window];
					for (int slice = 0; slice < window; slice++)
					{
						float[] sliceData = imageStack[slice];
						for (int pixel = 0, i = start; pixel < nPixels; pixel++, i++)
						{
							data[pixel][slice] = sliceData[i];
						}
					}

					// Initialise the window with the first n frames.
					MedianWindowDLLFloat[] mw = new MedianWindowDLLFloat[nPixels];
					for (int pixel = 0; pixel < nPixels; pixel++)
					{
						mw[pixel] = new MedianWindowDLLFloat(data[pixel]);
					}

					// Get the early medians.
					int slice = 0;
					for (; slice < radius; slice++)
					{
						for (int pixel = 0, i = start; pixel < nPixels; pixel++, i++)
						{
							imageStack[slice][i] = mw[pixel].getMedianOldest(slice + 1 + radius) * mean[slice];
						}
					}

					// Then increment through the data getting the median when required.
					for (int j = mw[0].getSize(); j < nSlices; j++, slice++)
					{
						for (int pixel = 0, i = start; pixel < nPixels; pixel++, i++)
						{
							imageStack[slice][i] = mw[pixel].getMedian() * mean[slice];
							mw[pixel].add(imageStack[j][i]);
						}
					}

					// Then get the later medians as required.
					for (int i = 2 * radius + 1; slice < nSlices; i--, slice++)
					{
						for (int pixel = 0, ii = start; pixel < nPixels; pixel++, ii++)
							imageStack[slice][ii] = mw[pixel].getMedianYoungest(i) * mean[slice];
					}
				}
				else
				{
					float[][] data = new float[nPixels][nSlices];
					for (int slice = 0; slice < nSlices; slice++)
					{
						float[] sliceData = imageStack[slice];
						for (int pixel = 0, i = start; pixel < nPixels; pixel++, i++)
						{
							data[pixel][slice] = sliceData[i];
						}
					}

					// Create median window filter
					MedianWindowFloat[] mw = new MedianWindowFloat[nPixels];
					for (int pixel = 0; pixel < nPixels; pixel++)
					{
						mw[pixel] = new MedianWindowFloat(data[pixel].clone(), radius);
					}

					// Produce the medians
					for (int slice = 0; slice < nSlices; slice += interval)
					{
						for (int pixel = 0; pixel < nPixels; pixel++)
						{
							data[pixel][slice] = mw[pixel].getMedian();
							mw[pixel].increment(interval);
						}
					}
					// Final position if necessary
					if (mw[0].getPosition() != nSlices + interval - 1)
					{
						for (int pixel = 0; pixel < nPixels; pixel++)
						{
							mw[pixel].setPosition(nSlices - 1);
							data[pixel][nSlices - 1] = mw[pixel].getMedian();
						}
					}

					// Interpolate
					float[] increment = new float[nPixels];
					for (int slice = 0; slice < nSlices; slice += interval)
					{
						int end = FastMath.min(slice + interval, nSlices - 1);
						for (int pixel = 0; pixel < nPixels; pixel++)
							increment[pixel] = (data[pixel][end] - data[pixel][slice]) / (end - slice);
						for (int s = slice + 1, i = 1; s < end; s++, i++)
						{
							for (int pixel = 0; pixel < nPixels; pixel++)
								data[pixel][s] = data[pixel][slice] + increment[pixel] * i;
						}
					}

					// Put back in the image re-scaling using the image mean
					for (int slice = 0; slice < nSlices; slice++)
					{
						float[] sliceData = imageStack[slice];
						for (int pixel = 0, i = start; pixel < nPixels; pixel++, i++)
						{
							sliceData[i] = data[pixel][slice] * mean[slice];
						}
					}
				}
			}
		}
	}

	/**
	 * Extract the data for a specified slice, subtract the background median filter and add the bias.
	 * <p>
	 * Use a runnable for the image generation to allow multi-threaded operation. Input parameters that are manipulated
	 * should have synchronized methods.
	 */
	private class ImageFilter implements Runnable
	{
		final ImageStack inputStack;
		final float[][] imageStack;
		final int n;

		public ImageFilter(ImageStack inputStack, float[][] imageStack, int n)
		{
			this.inputStack = inputStack;
			this.imageStack = imageStack;
			this.n = n;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			showProgressSingle();

			final float[] data = IJImageConverter.getData(inputStack.getProcessor(n));
			final float[] filter = imageStack[n - 1];
			final float b = bias;
			for (int i = 0; i < data.length; i++)
				filter[i] = data[i] - filter[i] + b;
		}
	}
}
