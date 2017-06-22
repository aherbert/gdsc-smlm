package gdsc.smlm.results;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;

import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

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

import gdsc.smlm.tsf.TaggedSpotFile.CameraType;
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
	/**
	 * Application ID assigned to GDSC SMLM ImageJ plugins
	 */
	public static final int APPLICATION_ID = 4;

	private FileOutputStream out = null;

	private String filename = null;

	private int size = 0;
	private AtomicInteger id;

	private int isx, isy, ia;
	private FitMode fitMode;

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

		// Only support Gaussian 2D data
		if (psf == null || !PSFHelper.isGaussian2D(psf))
		{
			System.err.println("TSF format requires a Gaussian 2D PSF");
			closeOutput();
			return;
		}
		int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
		isx = indices[0];
		isy = indices[1];
		try
		{
			ia = PSFHelper.getGaussian2DAngleIndex(psf);
			fitMode = FitMode.TWOAXISANDTHETA;
		}
		catch (ConfigurationException e)
		{
			// This is not an angled PSF. Revert to 1/2 axis:
			fitMode = (isx == isy) ? FitMode.ONEAXIS : FitMode.TWOAXIS;
		}

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

		setParams(params, builder);

		builder.setError(error);
		builder.setNoise(noise);
		builder.setOriginalValue(origValue);
		if (paramsStdDev != null)
			addNewParamsStdDev(builder, paramsStdDev);

		Spot spot = builder.build();

		writeResult(1, spot);
	}

	public void add(PeakResult result)
	{
		final float[] params = result.params;

		Spot.Builder builder = Spot.newBuilder();
		builder.setMolecule(id.incrementAndGet());
		builder.setChannel(1);
		builder.setFluorophoreType(1);
		builder.setFrame(result.getFrame());
		builder.setXPosition(result.origX);
		builder.setYPosition(result.origY);

		setParams(params, builder);

		if (result.hasPrecision())
		{
			// Use the actual precision
			float precision = (float) result.getPrecision();
			builder.setXPrecision(precision);
			builder.setYPrecision(precision);
		}

		if (result.hasId())
			builder.setCluster(result.getId());

		builder.setError(result.error);
		builder.setNoise(result.noise);
		if (result.hasEndFrame())
			builder.setEndFrame(result.getEndFrame());
		builder.setOriginalValue(result.origValue);
		if (result.paramsStdDev != null)
			addNewParamsStdDev(builder, result.paramsStdDev);

		Spot spot = builder.build();

		writeResult(1, spot);
	}

	/**
	 * Sets the width. Convert the X/Y widths used in GDSC SMLM to the single width and shape parameters used in TSF.
	 *
	 * @param params
	 *            the params
	 * @param builder
	 *            the builder
	 */
	private void setParams(float[] params, Spot.Builder builder)
	{
		builder.setBackground(params[PeakResult.BACKGROUND]);
		builder.setIntensity(params[PeakResult.INTENSITY]);
		builder.setX(params[PeakResult.X]);
		builder.setY(params[PeakResult.Y]);
		builder.setZ(params[PeakResult.Z]);
		
		switch (fitMode)
		{
			case ONEAXIS:
				builder.setWidth((float) (Gaussian2DFunction.SD_TO_FWHM_FACTOR * params[isx]));
				break;

			case TWOAXIS:
				builder.setWidth((float) (Gaussian2DFunction.SD_TO_FWHM_FACTOR *
						Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy])));
				builder.setA(params[isx] / params[isy]);
				break;

			case TWOAXISANDTHETA:
				builder.setWidth((float) (Gaussian2DFunction.SD_TO_FWHM_FACTOR *
						Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy])));
				builder.setA(params[isx] / params[isy]);
				builder.setTheta(params[ia]);
				break;
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

	public void addAll(PeakResult[] results)
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
			builder.setFrame(result.getFrame());
			builder.setXPosition(result.origX);
			builder.setYPosition(result.origY);

			setParams(params, builder);

			if (result.hasPrecision())
			{
				// Use the actual precision
				float precision = (float) result.getPrecision();
				builder.setXPrecision(precision);
				builder.setYPrecision(precision);
			}

			if (result.hasId())
				builder.setCluster(result.getId());
			else
				builder.clearCluster();

			builder.setError(result.error);
			builder.setNoise(result.noise);
			if (result.hasEndFrame())
				builder.setEndFrame(result.getEndFrame());
			else
				builder.clearEndFrame();
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

	private static CameraType[] cameraTypeMap;
	private static ThetaUnits[] thetaUnitsMap;
	private static LocationUnits[] locationUnitsMap;
	private static IntensityUnits[] intensityUnitsMap;
	static
	{
		// These should have 1:1 mapping. We can extends the TSF proto if necessary.		
		cameraTypeMap = new CameraType[3];
		cameraTypeMap[gdsc.smlm.data.config.SMLMSettings.CameraType.CCD.ordinal()] = CameraType.CCD;
		cameraTypeMap[gdsc.smlm.data.config.SMLMSettings.CameraType.EMCCD.ordinal()] = CameraType.EMCCD;
		cameraTypeMap[gdsc.smlm.data.config.SMLMSettings.CameraType.SCMOS.ordinal()] = CameraType.SCMOS;
		thetaUnitsMap = new ThetaUnits[2];
		thetaUnitsMap[gdsc.smlm.data.config.SMLMSettings.AngleUnit.RADIAN.ordinal()] = ThetaUnits.RADIANS;
		thetaUnitsMap[gdsc.smlm.data.config.SMLMSettings.AngleUnit.DEGREE.ordinal()] = ThetaUnits.DEGREES;
		locationUnitsMap = new LocationUnits[3];
		locationUnitsMap[gdsc.smlm.data.config.SMLMSettings.DistanceUnit.NM.ordinal()] = LocationUnits.NM;
		locationUnitsMap[gdsc.smlm.data.config.SMLMSettings.DistanceUnit.UM.ordinal()] = LocationUnits.UM;
		locationUnitsMap[gdsc.smlm.data.config.SMLMSettings.DistanceUnit.PIXEL.ordinal()] = LocationUnits.PIXELS;
		intensityUnitsMap = new IntensityUnits[2];
		intensityUnitsMap[gdsc.smlm.data.config.SMLMSettings.IntensityUnit.COUNT.ordinal()] = IntensityUnits.COUNTS;
		intensityUnitsMap[gdsc.smlm.data.config.SMLMSettings.IntensityUnit.PHOTON.ordinal()] = IntensityUnits.PHOTONS;
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
			if (calibration.hasNmPerPixel())
				builder.setPixelSize((float) calibration.getNmPerPixel());

			if (calibration.hasExposureTime())
				builder.setExposureTime(calibration.getExposureTime());
			if (calibration.hasReadNoise())
				builder.setReadNoise(calibration.getReadNoise());
			if (calibration.hasBias())
				builder.setBias(calibration.getBias());
			if (calibration.hasCameraType())
				builder.setCameraType(cameraTypeMap[calibration.getCameraType().ordinal()]);
			if (calibration.hasAmplification())
				builder.setAmplification(calibration.getAmplification());

			if (calibration.hasDistanceUnit())
			{
				builder.setLocationUnits(locationUnitsMap[calibration.getDistanceUnit().ordinal()]);
			}
			if (calibration.hasIntensityUnit())
			{
				builder.setIntensityUnits(intensityUnitsMap[calibration.getIntensityUnit().ordinal()]);
			}
			if (calibration.hasAngleUnit())
			{
				builder.setThetaUnits(thetaUnitsMap[calibration.getAngleUnit().ordinal()]);
			}

			// We can use some logic here to get the QE
			if (calibration.hasGain())
			{
				builder.setGain(calibration.getGain());

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
		Printer printer = null;
		if (psf != null)
		{
			try
			{
				if (printer == null)
					printer = JsonFormat.printer().omittingInsignificantWhitespace();
				builder.setPSF(printer.print(psf));
			}
			catch (InvalidProtocolBufferException e)
			{
				// This shouldn't happen so throw it
				throw new NotImplementedException("Unable to serialise the PSF settings", e);
			}
		}

		// Have a property so the boxSize can be set
		if (boxSize > 0)
			builder.setBoxSize(boxSize);

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
