package gdsc.smlm.results;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collection;
import java.util.concurrent.atomic.AtomicInteger;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.tsf.TaggedSpotFile.FitMode;
import gdsc.smlm.tsf.TaggedSpotFile.IntensityUnits;
import gdsc.smlm.tsf.TaggedSpotFile.LocationUnits;
import gdsc.smlm.tsf.TaggedSpotFile.Spot;
import gdsc.smlm.tsf.TaggedSpotFile.SpotList;

/**
 * Saves the fit results to file using the Tagged Spot File (TSF) format.
 * 
 * @author Alex Herbert
 */
public class TSFPeakResultsWriter extends AbstractPeakResults
{
	public static final float SD_TO_FWHM_FACTOR = (float) (2.0 * Math.sqrt(2.0 * Math.log(2.0)));

	/**
	 * Application ID assigned to GDSC SMLM ImageJ plugins
	 */
	public static final int APPLICATION_ID = 4;

	private FileOutputStream out = null;

	private String filename = null;

	private int size = 0;
	private AtomicInteger id;

	private FitMode fitMode = FitMode.ONEAXIS;

	private int boxSize = 0;

	public TSFPeakResultsWriter(String filename)
	{
		this.filename = filename;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.AbstractPeakResults#begin()
	 */
	public void begin()
	{
		out = null;
		size = 0;
		id = new AtomicInteger();
		try
		{
			out = new FileOutputStream(filename);
		}
		catch (Exception e)
		{
			System.err.println("Failed to write open TSF file: " + filename);
			e.printStackTrace();
			closeOutput();
			return;
		}

		// Write the offsets used in the TSF format
		try
		{
			DataOutputStream dos = new DataOutputStream(out);
			dos.writeInt(0);
			dos.writeLong(0);
		}
		catch (IOException e)
		{
			System.err.println("Failed to write TSF offset fields");
			e.printStackTrace();
			closeOutput();
		}
	}

	private void closeOutput()
	{
		if (out == null)
			return;

		try
		{
			out.close();
		}
		catch (Exception e)
		{
			// Ignore exception
		}
		finally
		{
			out = null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		return out != null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		if (out == null)
			return;

		Spot.Builder builder = Spot.newBuilder();
		builder.setMolecule(id.incrementAndGet());
		builder.setChannel(1);
		builder.setFrame(peak);
		builder.setXPosition(origX);
		builder.setYPosition(origY);
		builder.setBackground(params[Gaussian2DFunction.BACKGROUND]);
		builder.setIntensity(params[Gaussian2DFunction.SIGNAL]);
		builder.setX(params[Gaussian2DFunction.X_POSITION]);
		builder.setY(params[Gaussian2DFunction.Y_POSITION]);

		setWidth(params, builder);

		if (this.calibration != null)
		{
			double s = (params[Gaussian2DFunction.X_SD] + params[Gaussian2DFunction.Y_SD]) * 0.5 *
					calibration.nmPerPixel;
			float precision = (float) PeakResult.getPrecision(calibration.nmPerPixel, s,
					params[Gaussian2DFunction.SIGNAL] / calibration.gain, noise / calibration.gain, calibration.emCCD);
			builder.setXPrecision(precision);
			builder.setYPrecision(precision);
		}

		// TODO: Extend the Spot with fields for....
		// Note: paramsStdDev for X/Y could be set into the X/YPrecision field.
		//		error;
		//		noise;
		//		origValue;
		//		paramsStdDev;

		Spot spot = builder.build();

		writeResult(1, spot);
	}

	private void setWidth(float[] params, Spot.Builder builder)
	{
		if (params[Gaussian2DFunction.X_SD] == params[Gaussian2DFunction.Y_SD])
		{
			builder.setWidth(SD_TO_FWHM_FACTOR * params[Gaussian2DFunction.X_SD]);
		}
		else
		{
			FitMode newFitMode = FitMode.TWOAXIS;
			builder.setWidth(SD_TO_FWHM_FACTOR *
					(float) Math.sqrt(Math.abs(params[Gaussian2DFunction.X_SD] * params[Gaussian2DFunction.Y_SD])));
			if (params[Gaussian2DFunction.ANGLE] != 0)
			{
				newFitMode = FitMode.TWOAXISANDTHETA;
				builder.setTheta(params[Gaussian2DFunction.ANGLE]);
			}
			builder.setA(params[Gaussian2DFunction.X_SD] / params[Gaussian2DFunction.Y_SD]);

			if (fitMode.getNumber() < newFitMode.getNumber())
				fitMode = newFitMode;
		}
	}

	public void addAll(Collection<PeakResult> results)
	{
		if (out == null)
			return;

		Spot[] spots = new Spot[20];
		int count = 0;
		Spot.Builder builder = Spot.newBuilder();
		for (PeakResult result : results)
		{
			final float[] params = result.params;

			builder.setMolecule(id.incrementAndGet());
			builder.setChannel(1);
			builder.setFrame(result.peak);
			builder.setXPosition(result.origX);
			builder.setYPosition(result.origY);
			builder.setBackground(params[Gaussian2DFunction.BACKGROUND]);
			builder.setIntensity(params[Gaussian2DFunction.SIGNAL]);
			builder.setX(params[Gaussian2DFunction.X_POSITION]);
			builder.setY(params[Gaussian2DFunction.Y_POSITION]);

			setWidth(params, builder);

			if (this.calibration != null)
			{
				double s = (params[Gaussian2DFunction.X_SD] + params[Gaussian2DFunction.Y_SD]) * 0.5 *
						calibration.nmPerPixel;
				float precision = (float) PeakResult.getPrecision(calibration.nmPerPixel, s,
						params[Gaussian2DFunction.SIGNAL] / calibration.gain, result.noise / calibration.gain,
						calibration.emCCD);
				builder.setXPrecision(precision);
				builder.setYPrecision(precision);
			}

			// TODO: Extend the Spot with additional fields for...
			// Note: paramsStdDev for X/Y could be set into the X/YPrecision field.
			//			result.error;
			//			result.noise;
			//			result.getId();
			//			result.getEndFrame();
			//			result.origValue;
			//			result.paramsStdDev;

			spots[count++] = builder.build();

			// Flush the output to allow for very large input lists
			if (count >= spots.length)
			{
				writeResult(count, spots);
				if (!isActive())
					return;
				count = 0;
			}
		}
		writeResult(count, spots);
	}

	private synchronized void writeResult(int count, Spot... spots)
	{
		// In case another thread caused the output to close
		if (out == null)
			return;
		size += count;
		try
		{
			for (int i = 0; i < count; i++)
			{
				spots[i].writeDelimitedTo(out);
			}
		}
		catch (IOException e)
		{
			System.err.println("Failed to write Spot message");
			closeOutput();
		}
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
		// Get the offset to the SpotList message
		long offset = 0;
		try
		{
			offset = out.getChannel().position();
		}
		catch (IOException e)
		{
			// This is bad.
			System.err.println("Failed to determine offset for SpotList message");
			e.printStackTrace();
			closeOutput();
			return;
		}

		// Record the SpotList message
		SpotList.Builder builder = SpotList.newBuilder();

		builder.setApplicationId(APPLICATION_ID);

		builder.setNrSpots(size);

		// Add the standard details the TSF supports
		if (source != null)
		{
			//String.format("#Source %s\n", singleLine(source.toXML()));
			builder.setName(source.getName());
			builder.setNrPixelsX(source.width);
			builder.setNrPixelsY(source.height);
			builder.setNrFrames(source.frames);
		}
		if (name != null)
		{
			// This may have been set already by source
			if (!builder.hasName())
				builder.setName(name);
		}
		if (bounds != null)
		{
			//String.format("#Bounds x%d y%d w%d h%d\n", bounds.x, bounds.y, bounds.width, bounds.height);

			// This may have been set already by source
			if (!builder.hasNrPixelsX())
			{
				// Include the origin so that the number of pixels will cover the range of the data
				builder.setNrPixelsX(bounds.x + bounds.width);
				builder.setNrPixelsY(bounds.y + bounds.height);
			}
		}
		if (calibration != null)
		{
			//String.format("#Calibration %s\n", singleLine(XmlUtils.toXML(calibration)));

			builder.setPixelSize((float) calibration.nmPerPixel);
		}
		if (configuration != null && configuration.length() > 0)
		{
			//String.format("#Configuration %s\n", singleLine(configuration));

		}

		// Have a property so the boxSize can be set
		if (boxSize > 0)
			builder.setBoxSize(boxSize);

		builder.setLocationUnits(LocationUnits.PIXELS);
		builder.setIntensityUnits(IntensityUnits.COUNTS);
		builder.setFitMode(fitMode);

		// TODO: Extend the SpotList with additional fields for ... 
		// Some of these will need to be 'message' types.
		//      name;
		//		source;
		//		bounds;
		//		calibration;
		//		configuration;

		SpotList spotList = builder.build();
		try
		{
			spotList.writeDelimitedTo(out);
		}
		catch (IOException e)
		{
			System.err.println("Failed to write SpotList message");
			e.printStackTrace();
			return;
		}
		finally
		{
			closeOutput();
		}

		// Note: it would be good to be able to use the ability to write to any output stream. However
		// the TSF format requires a seek at the end of writing to record the offset. seek() is not
		// supported by OutputStream. It is supported by: RandomAccessFile, RandomAccessStream (for input). 

		// Write the offset to the SpotList message into the offset position
		RandomAccessFile f = null;
		try
		{
			f = new RandomAccessFile(new File(filename), "rw");
			f.seek(4);
			f.writeLong(offset);
		}
		catch (Exception e)
		{
			System.err.println("Failed to record offset for SpotList message");
			e.printStackTrace();
		}
		finally
		{
			if (f != null)
			{
				try
				{
					f.close();
				}
				catch (IOException e)
				{
				}
			}
		}
	}

	/**
	 * Gets the box size.
	 *
	 * @return the box size (in pixels) of rectangular box used in Gaussian fitting
	 */
	public int getBoxSize()
	{
		return boxSize;
	}

	/**
	 * Sets the box size.
	 *
	 * @param boxSize
	 *            the box size (in pixels) of rectangular box used in Gaussian fitting
	 */
	public void setBoxSize(int boxSize)
	{
		this.boxSize = boxSize;
	}
}
