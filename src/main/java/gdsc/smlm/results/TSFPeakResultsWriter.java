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
import gdsc.smlm.tsf.TaggedSpotFile.FluorophoreType;
import gdsc.smlm.tsf.TaggedSpotFile.IntensityUnits;
import gdsc.smlm.tsf.TaggedSpotFile.LocationUnits;
import gdsc.smlm.tsf.TaggedSpotFile.ROI;
import gdsc.smlm.tsf.TaggedSpotFile.Spot;
import gdsc.smlm.tsf.TaggedSpotFile.Spot.Builder;
import gdsc.smlm.tsf.TaggedSpotFile.SpotList;
import gdsc.smlm.tsf.TaggedSpotFile.ThetaUnits;

/**
 * Saves the fit results to file using the Tagged Spot File (TSF) format.
 * <p>
 * Write out a TSF file assuming the results are in the standard GSDC SMLM format (intensity in counts, angles in
 * degrees).
 * <p>
 * To satisfy the format the calibration must be set including the amplification (electrons/count) and camera bias. The
 * bias is removed from the background. If amplification is not strictly positive then the calibration gain will be
 * written to the TSF 'electron conversion factor' field.
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

	private float bias = 0;

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
		bias = (calibration != null) ? (float) calibration.getBias() : 0;
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
		builder.setFluorophoreType(1); 
		builder.setFrame(peak);
		builder.setXPosition(origX);
		builder.setYPosition(origY);
		setBackground(builder, params[Gaussian2DFunction.BACKGROUND]);
		builder.setIntensity(params[Gaussian2DFunction.SIGNAL]);
		builder.setX(params[Gaussian2DFunction.X_POSITION]);
		builder.setY(params[Gaussian2DFunction.Y_POSITION]);

		setWidth(params, builder);

		if (this.calibration != null)
		{
			double s = (params[Gaussian2DFunction.X_SD] + params[Gaussian2DFunction.Y_SD]) * 0.5 *
					calibration.getNmPerPixel();
			float precision = (float) PeakResult.getPrecision(calibration.getNmPerPixel(), s,
					params[Gaussian2DFunction.SIGNAL] / calibration.getGain(), noise / calibration.getGain(),
					calibration.isEmCCD());
			builder.setXPrecision(precision);
			builder.setYPrecision(precision);
		}

		builder.setError(error);
		builder.setNoise(noise);
		builder.setOriginalValue(origValue);
		if (paramsStdDev != null)
			addNewParamsStdDev(builder, paramsStdDev);

		Spot spot = builder.build();

		writeResult(1, spot);
	}

	/**
	 * Sets the background.
	 *
	 * @param builder
	 *            the builder
	 * @param background
	 *            the background
	 */
	private void setBackground(Builder builder, float background)
	{
		// Q. Should we ensure this is always positive?
		// Since it "should be linearly proportional to the number of photons in the background" 
		// then we assume that it cannot be negative.
		if (background > bias)
			builder.setBackground(background - bias);
		else
			builder.setBackground(0f);
	}

	/**
	 * Sets the width. Convert the X/Y widths used in GDSC SMLM to the single width and shape parameters used in TSF.
	 *
	 * @param params
	 *            the params
	 * @param builder
	 *            the builder
	 */
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
			builder.setA(params[Gaussian2DFunction.X_SD] / params[Gaussian2DFunction.Y_SD]);

			if (params[Gaussian2DFunction.ANGLE] != 0)
			{
				newFitMode = FitMode.TWOAXISANDTHETA;
				builder.setTheta(params[Gaussian2DFunction.ANGLE]);
			}

			if (fitMode.getNumber() < newFitMode.getNumber())
				fitMode = newFitMode;
		}
	}

	/**
	 * Adds the params std dev assuming a new builder.
	 *
	 * @param builder
	 *            the builder
	 * @param paramsStdDev
	 *            the params std dev
	 */
	private void addNewParamsStdDev(Builder builder, float[] paramsStdDev)
	{
		// Note: paramsStdDev for X/Y could be set into the X/Y Precision field.
		for (int i = 0; i < paramsStdDev.length; i++)
			builder.addParamsStdDev(paramsStdDev[i]);
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
			builder.setFluorophoreType(1); 
			builder.setFrame(result.peak);
			builder.setXPosition(result.origX);
			builder.setYPosition(result.origY);
			setBackground(builder, params[Gaussian2DFunction.BACKGROUND]);
			builder.setIntensity(params[Gaussian2DFunction.SIGNAL]);
			builder.setX(params[Gaussian2DFunction.X_POSITION]);
			builder.setY(params[Gaussian2DFunction.Y_POSITION]);

			setWidth(params, builder);

			if (this.calibration != null)
			{
				double s = (params[Gaussian2DFunction.X_SD] + params[Gaussian2DFunction.Y_SD]) * 0.5 *
						calibration.getNmPerPixel();
				float precision = (float) PeakResult.getPrecision(calibration.getNmPerPixel(), s,
						params[Gaussian2DFunction.SIGNAL] / calibration.getGain(), result.noise / calibration.getGain(),
						calibration.isEmCCD());
				builder.setXPrecision(precision);
				builder.setYPrecision(precision);
			}

			builder.setCluster(result.getId());

			builder.setError(result.error);
			builder.setNoise(result.noise);
			builder.setEndFrame(result.getEndFrame());
			builder.setOriginalValue(result.origValue);
			addParamsStdDev(builder, result.paramsStdDev);

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

	/**
	 * Adds the params std dev assuming an existing builder (allowing re-use of the space).
	 *
	 * @param builder
	 *            the builder
	 * @param paramsStdDev
	 *            the params std dev
	 */
	private void addParamsStdDev(Builder builder, float[] paramsStdDev)
	{
		// Note: paramsStdDev for X/Y could be set into the X/Y Precision field.
		if (paramsStdDev == null)
		{
			if (builder.getParamsStdDevCount() != 0)
				builder.clearParamsStdDev();
			return;
		}

		// Reuse the space
		if (builder.getParamsStdDevCount() == paramsStdDev.length)
		{
			for (int i = 0; i < paramsStdDev.length; i++)
				builder.setParamsStdDev(i, paramsStdDev[i]);
		}
		else
		{
			builder.clearParamsStdDev();
			addNewParamsStdDev(builder, paramsStdDev);
		}
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
			// The offset is the amount to skip forward after reading the int 
			// magic number (4 bytes) and long offset (8 bytes)
			//out.flush();
			offset = out.getChannel().position() - 12;
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

		// Add the standard details the TSF supports. We use extensions to add GDSC SMLM data.
		if (name != null)
		{
			builder.setName(name);
		}
		if (source != null)
		{
			builder.setNrPixelsX(source.width);
			builder.setNrPixelsY(source.height);
			builder.setNrFrames(source.frames);

			builder.setSource(singleLine(source.toXML()));
		}
		if (bounds != null)
		{
			ROI.Builder roiBuilder = builder.getRoiBuilder();
			roiBuilder.setX(bounds.x);
			roiBuilder.setY(bounds.y);
			roiBuilder.setXWidth(bounds.width);
			roiBuilder.setYWidth(bounds.height);
			builder.setRoi(roiBuilder.build());
		}
		if (calibration != null)
		{
			builder.setPixelSize((float) calibration.getNmPerPixel());

			builder.setGain(calibration.getGain());
			builder.setExposureTime(calibration.getExposureTime());
			builder.setReadNoise(calibration.getReadNoise());
			builder.setBias(calibration.getBias());
			builder.setEmCCD(calibration.isEmCCD());
			builder.setAmplification(calibration.getAmplification());

			if (calibration.hasGain())
			{
				// Use amplification if present (as this is the correct electrons/count value), otherwise use gain
				if (calibration.hasAmplification())
				{
					double ecf = calibration.getAmplification();
					double qe = calibration.getGain() / ecf;
					builder.addEcf(ecf);
					builder.addQe(qe);
				}
				else
				{
					builder.addEcf(calibration.getGain());
					builder.addQe(1);
				}
			}
		}
		if (configuration != null && configuration.length() > 0)
		{
			builder.setConfiguration(singleLine(configuration));
		}

		// Have a property so the boxSize can be set
		if (boxSize > 0)
			builder.setBoxSize(boxSize);

		builder.setLocationUnits(LocationUnits.PIXELS);
		builder.setIntensityUnits(IntensityUnits.COUNTS);
		builder.setThetaUnits(ThetaUnits.DEGREES);
		builder.setFitMode(fitMode);
		
		FluorophoreType.Builder typeBuilder = FluorophoreType.newBuilder();
		typeBuilder.setId(1);
		typeBuilder.setDescription("Default fluorophore");
		typeBuilder.setIsFiducial(false);
		builder.addFluorophoreTypes(typeBuilder.build());

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

	private String singleLine(String text)
	{
		return text.replaceAll("\n *", "");
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
