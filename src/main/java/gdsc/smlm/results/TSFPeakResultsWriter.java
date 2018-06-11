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
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import gdsc.smlm.tsf.TSFProtos.CameraType;
import gdsc.smlm.tsf.TSFProtos.FitMode;
import gdsc.smlm.tsf.TSFProtos.FluorophoreType;
import gdsc.smlm.tsf.TSFProtos.IntensityUnits;
import gdsc.smlm.tsf.TSFProtos.LocationUnits;
import gdsc.smlm.tsf.TSFProtos.ROI;
import gdsc.smlm.tsf.TSFProtos.Spot;
import gdsc.smlm.tsf.TSFProtos.Spot.Builder;
import gdsc.smlm.tsf.TSFProtos.SpotList;
import gdsc.smlm.tsf.TSFProtos.ThetaUnits;

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
	@Override
	public void begin()
	{
		out = null;
		size = 0;

		// Only support Gaussian 2D data
		if (getPSF() == null || !PSFHelper.isGaussian2D(getPSF()))
		{
			System.err.println("TSF format requires a Gaussian 2D PSF");
			closeOutput();
			return;
		}
		int[] indices = PSFHelper.getGaussian2DWxWyIndices(getPSF());
		isx = indices[0];
		isy = indices[1];
		try
		{
			ia = PSFHelper.getGaussian2DAngleIndex(getPSF());
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
	@Override
	public boolean isActive()
	{
		return out != null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResults#add(int, int, int, float, double, float, float, float[], float[])
	 */
	@Override
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float meanIntensity,
			float[] params, float[] paramsStdDev)
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

		setParam(params, builder);

		builder.setError(error);
		builder.setNoise(noise);
		builder.setMeanIntensity(meanIntensity);
		builder.setOriginalValue(origValue);
		if (paramsStdDev != null)
			addNewParamStdDevs(builder, paramsStdDev);

		Spot spot = builder.build();

		writeResult(1, spot);
	}

	@Override
	public void add(PeakResult result)
	{
		final float[] params = result.getParameters();

		Spot.Builder builder = Spot.newBuilder();
		builder.setMolecule(id.incrementAndGet());
		builder.setChannel(1);
		builder.setFluorophoreType(1);
		builder.setFrame(result.getFrame());
		builder.setXPosition(result.getOrigX());
		builder.setYPosition(result.getOrigY());

		setParam(params, builder);

		if (result.hasPrecision())
		{
			// Use the actual precision
			float precision = (float) result.getPrecision();
			builder.setXPrecision(precision);
			builder.setYPrecision(precision);
		}

		if (result.hasId())
			builder.setCluster(result.getId());

		builder.setError(result.getError());
		builder.setNoise(result.getNoise());
		builder.setMeanIntensity(result.getMeanIntensity());
		if (result.hasEndFrame())
			builder.setEndFrame(result.getEndFrame());
		builder.setOriginalValue(result.getOrigValue());
		if (result.hasParameterDeviations())
			addNewParamStdDevs(builder, result.getParameterDeviations());

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
	private void setParam(float[] params, Spot.Builder builder)
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
	 * @param paramStdDev
	 *            the params std dev
	 */
	private void addNewParamStdDevs(Builder builder, float[] paramStdDev)
	{
		// Note: paramsStdDev for X/Y could be set into the X/Y Precision field.
		for (int i = 0; i < paramStdDev.length; i++)
			builder.addParamStdDevs(paramStdDev[i]);
	}

	@Override
	public void addAll(PeakResult[] results)
	{
		if (out == null)
			return;

		Spot[] spots = new Spot[20];
		int count = 0;
		Spot.Builder builder = Spot.newBuilder();
		for (PeakResult result : results)
		{
			final float[] params = result.getParameters();

			builder.setMolecule(id.incrementAndGet());
			builder.setChannel(1);
			builder.setFluorophoreType(1);
			builder.setFrame(result.getFrame());
			builder.setXPosition(result.getOrigX());
			builder.setYPosition(result.getOrigY());

			setParam(params, builder);

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

			builder.setError(result.getError());
			builder.setNoise(result.getNoise());
			builder.setMeanIntensity(result.getMeanIntensity());
			if (result.hasEndFrame())
				builder.setEndFrame(result.getEndFrame());
			else
				builder.clearEndFrame();
			builder.setOriginalValue(result.getOrigValue());
			addParamStdDevs(builder, result.getParameterDeviations());

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
	 * @param paramStdDev
	 *            the params std dev
	 */
	private void addParamStdDevs(Builder builder, float[] paramStdDev)
	{
		// Note: paramsStdDev for X/Y could be set into the X/Y Precision field.
		if (paramStdDev == null)
		{
			if (builder.getParamStdDevsCount() != 0)
				builder.clearParamStdDevs();
			return;
		}

		// Reuse the space
		if (builder.getParamStdDevsCount() == paramStdDev.length)
		{
			for (int i = 0; i < paramStdDev.length; i++)
				builder.setParamStdDevs(i, paramStdDev[i]);
		}
		else
		{
			builder.clearParamStdDevs();
			addNewParamStdDevs(builder, paramStdDev);
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
	@Override
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
		cameraTypeMap = new CameraType[gdsc.smlm.data.config.CalibrationProtos.CameraType.values().length];
		cameraTypeMap[gdsc.smlm.data.config.CalibrationProtos.CameraType.CCD.ordinal()] = CameraType.CCD;
		cameraTypeMap[gdsc.smlm.data.config.CalibrationProtos.CameraType.EMCCD.ordinal()] = CameraType.EMCCD;
		cameraTypeMap[gdsc.smlm.data.config.CalibrationProtos.CameraType.SCMOS.ordinal()] = CameraType.SCMOS;
		thetaUnitsMap = new ThetaUnits[gdsc.smlm.data.config.UnitProtos.AngleUnit.values().length];
		thetaUnitsMap[gdsc.smlm.data.config.UnitProtos.AngleUnit.RADIAN.ordinal()] = ThetaUnits.RADIANS;
		thetaUnitsMap[gdsc.smlm.data.config.UnitProtos.AngleUnit.DEGREE.ordinal()] = ThetaUnits.DEGREES;
		locationUnitsMap = new LocationUnits[gdsc.smlm.data.config.UnitProtos.DistanceUnit.values().length];
		locationUnitsMap[gdsc.smlm.data.config.UnitProtos.DistanceUnit.NM.ordinal()] = LocationUnits.NM;
		locationUnitsMap[gdsc.smlm.data.config.UnitProtos.DistanceUnit.UM.ordinal()] = LocationUnits.UM;
		locationUnitsMap[gdsc.smlm.data.config.UnitProtos.DistanceUnit.PIXEL.ordinal()] = LocationUnits.PIXELS;
		intensityUnitsMap = new IntensityUnits[gdsc.smlm.data.config.UnitProtos.IntensityUnit.values().length];
		intensityUnitsMap[gdsc.smlm.data.config.UnitProtos.IntensityUnit.COUNT.ordinal()] = IntensityUnits.COUNTS;
		intensityUnitsMap[gdsc.smlm.data.config.UnitProtos.IntensityUnit.PHOTON.ordinal()] = IntensityUnits.PHOTONS;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#end()
	 */
	@Override
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
		if (!TextUtils.isNullOrEmpty(getName()))
		{
			builder.setName(getName());
		}
		if (getSource() != null)
		{
			builder.setNrPixelsX(getSource().width);
			builder.setNrPixelsY(getSource().height);
			builder.setNrFrames(getSource().frames);

			builder.setSource(singleLine(getSource().toXML()));
		}
		if (getBounds() != null)
		{
			ROI.Builder roiBuilder = builder.getRoiBuilder();
			roiBuilder.setX(getBounds().x);
			roiBuilder.setY(getBounds().y);
			roiBuilder.setXWidth(getBounds().width);
			roiBuilder.setYWidth(getBounds().height);
			builder.setRoi(roiBuilder.build());
		}
		if (hasCalibration())
		{
			CalibrationReader cr = getCalibrationReader();
			if (cr.hasNmPerPixel())
				builder.setPixelSize((float) cr.getNmPerPixel());

			if (cr.hasExposureTime())
				builder.setExposureTime(cr.getExposureTime());
			if (cr.hasReadNoise())
				builder.setReadNoise(cr.getReadNoise());
			if (cr.hasBias())
				builder.setBias(cr.getBias());
			if (cr.hasCameraType())
				builder.setCameraType(cameraTypeMap[cr.getCameraType().ordinal()]);

			if (cr.hasDistanceUnit())
			{
				builder.setLocationUnits(locationUnitsMap[cr.getDistanceUnit().ordinal()]);
			}
			if (cr.hasIntensityUnit())
			{
				builder.setIntensityUnits(intensityUnitsMap[cr.getIntensityUnit().ordinal()]);
			}
			if (cr.hasAngleUnit())
			{
				builder.setThetaUnits(thetaUnitsMap[cr.getAngleUnit().ordinal()]);
			}

			// We can use some logic here to get the QE
			if (cr.hasCountPerPhoton())
			{
				builder.setGain(cr.getCountPerPhoton());

				double qe = (cr.hasQuantumEfficiency()) ? cr.getQuantumEfficiency() : 1;
				// e-/photon / count/photon => e-/count
				double ecf = qe / cr.getCountPerPhoton();

				builder.addEcf(ecf);
				builder.addQe(qe);
			}
		}
		if (!TextUtils.isNullOrEmpty(getConfiguration()))
		{
			builder.setConfiguration(singleLine(getConfiguration()));
		}
		Printer printer = null;
		if (getPSF() != null)
		{
			try
			{
				if (printer == null)
					printer = JsonFormat.printer().omittingInsignificantWhitespace();
				builder.setPSF(printer.print(getPSF()));
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
