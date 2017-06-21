package gdsc.smlm.results;

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

import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.util.JsonFormat;

import gdsc.core.ij.Utils;
import gdsc.core.logging.TrackProgress;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.UnicodeReader;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.PSF;
import gdsc.smlm.data.config.SMLMSettings.PSFType;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
//import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.data.config.UnitHelper;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.utils.XmlUtils;

/**
 * Reads the fit results from file
 */
public class PeakResultsReader
{
	// Set up to read two-axis (and theta) Gaussian 2D data into the current format
	private final static int isx, isy, ia, nTwoAxis, nTwoAxisAndTheta;
	static
	{
		PSF psf = PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
		int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
		isx = indices[0];
		isy = indices[1];
		ia = PSFHelper.getGaussian2DAngleIndex(psf);
		nTwoAxis = PeakResult.STANDARD_PARAMETERS + 2;
		nTwoAxisAndTheta = PeakResult.STANDARD_PARAMETERS + 3;
	}

	/** The columns to recognise in the ImageJ table results header */
	private static String IMAGEJ_TABLE_RESULTS_HEADER = "origX\torigY\torigValue\tError\tNoise";

	/** The space patterm */
	private static Pattern spacePattern = Pattern.compile(" ");
	/** The tab patterm */
	private static Pattern tabPattern = Pattern.compile("\t");
	/** Simple whitespace pattern for tabs of spaces */
	private static Pattern whitespacePattern = Pattern.compile("[\t ]");

	private boolean useScanner = false;
	private boolean rawResults = false;

	private String filename;
	private String header = null;
	private FileFormat format;
	private String version;
	private String name = null;
	private ImageSource source = null;
	private Rectangle bounds = null;
	private Calibration calibration = null;
	private PSF psf = null;
	private String configuration = null;
	private TrackProgress tracker = null;
	private ResultOption[] options = null;

	private boolean deviations, readEndFrame, readId, readSource;
	private int smlmVersion = 3; // Assume the current

	public PeakResultsReader(String filename)
	{
		this.filename = filename;
	}

	/**
	 * @return The header from the results file
	 */
	public String getHeader()
	{
		if (header == null)
		{
			BufferedReader input = null;
			try
			{
				FileInputStream fis = new FileInputStream(filename);
				input = new BufferedReader(new UnicodeReader(fis, null));

				StringBuffer sb = new StringBuffer();
				String line;
				int count = 0;
				while ((line = input.readLine()) != null)
				{
					count++;
					if (count == 1)
					{
						// The NSTORM file format does not have comment characters but does have a single header line
						if (line.startsWith("Channel Name"))
						{
							sb.append(line).append("\n");
							break;
						}
						// User may try and load the text saved directly from the ImageJ Table Results
						if (line.contains(IMAGEJ_TABLE_RESULTS_HEADER))
						{
							sb.append(line).append("\n");
							break;
						}
					}
					if (line.length() == 0)
						continue;
					if (line.charAt(0) == '#')
						sb.append(line).append("\n");
					else
						break;
				}
				header = sb.toString();

				version = getField("FileVersion");

				format = FileFormat.UNKNOWN;
				guessFormatFromVersion();
				if (format == FileFormat.UNKNOWN)
					guessFormat(line);
			}
			catch (IOException e)
			{
				// ignore
			}
			finally
			{
				try
				{
					if (input != null)
						input.close();
				}
				catch (IOException e)
				{
					// Ignore
				}
			}
		}
		return header;
	}

	private void guessFormatFromVersion()
	{
		// Extract information about the file format
		if (version.length() > 0)
		{
			if (version.startsWith("Binary"))
				format = FileFormat.SMLM_BINARY;
			else if (version.startsWith("Text"))
				format = FileFormat.SMLM_TEXT;
			else
				return;

			if (version.contains(".V3"))
			{
				smlmVersion = 3;
			}
			else if (version.contains(".V2"))
			{
				smlmVersion = 2;
			}

			if (smlmVersion == 1)
			{
				// The original files did not have the version tag
				deviations = header.contains((format == FileFormat.SMLM_BINARY) ? "iiiifdfffffffffffffff" : "+/-");
				readEndFrame = readId = false;
			}
			else
			{
				deviations = version.contains(".D1");
				// Extended marker has a bit flag:
				// 1 = showEndFrame
				// 2 = showId
				if (version.contains(".E3"))
				{
					readEndFrame = readId = true;
				}
				else
				{
					readEndFrame = version.contains(".E1");
					readId = version.contains(".E2");
				}
			}
		}
	}

	private void guessFormat(String firstLine)
	{
		if (header.length() == 0)
		{
			// No header. Check non-text formats.
			if (TSFPeakResultsReader.isTSF(filename))
			{
				format = FileFormat.TSF_BINARY;
				return;
			}

			// We cannot continue guess if there is no non-header data
			if (Utils.isNullOrEmpty(firstLine))
				return;
		}

		// Check for Nikon NSTORM header
		if (header.contains("Channel Name"))
		{
			format = FileFormat.NSTORM;
		}
		else if (header.contains(IMAGEJ_TABLE_RESULTS_HEADER))
		{
			format = FileFormat.SMLM_TABLE;
			readId = header.startsWith("#");
			readEndFrame = header.contains("\tEnd ");
			deviations = header.contains("\t+/-\t");
			readSource = header.contains("Source\t");
		}
		// Check for RapidSTORM stuff in the header
		else if (header.contains("<localizations "))
		{
			// RapidSTORM can use the MALK format
			if (isMALKFormat(firstLine))
			{
				format = FileFormat.MALK;
			}
			else
			{
				// There may be many other formats for RapidSTORM. Just support the one we know about.
				format = FileFormat.RAPID_STORM;
			}
		}
		// Check for MALK format: X,Y,T,Signal
		else if (isMALKFormat(firstLine))
		{
			format = FileFormat.MALK;
		}
		else
		{
			// Assume SMLM format
			guessFormatFromVersion();
		}
	}

	/**
	 * @return True if the results file is binary
	 */
	public boolean isBinary()
	{
		getFormat();
		return format == FileFormat.SMLM_BINARY;
	}

	/**
	 * @return The file format
	 */
	public FileFormat getFormat()
	{
		if (format == null)
		{
			getHeader();
		}
		return format;
	}

	/**
	 * @return The bounds specified in the results header
	 */
	public Rectangle getBounds()
	{
		if (bounds == null)
		{
			getHeader();
			bounds = new Rectangle(0, 0, 0, 0);
			if (header != null)
			{
				Pattern pattern = Pattern.compile("x(\\d+) y(\\d+) w(\\d+) h(\\d+)");
				Matcher match = pattern.matcher(header);
				if (match.find())
				{
					bounds.x = Integer.parseInt(match.group(1));
					bounds.y = Integer.parseInt(match.group(2));
					bounds.width = Integer.parseInt(match.group(3));
					bounds.height = Integer.parseInt(match.group(4));
				}
			}
		}
		return bounds;
	}

	/**
	 * @return The name specified in the results header
	 */
	public String getName()
	{
		if (name == null)
		{
			getHeader();
			name = getField("Name");
		}
		return name;
	}

	/**
	 * @return The source specified in the results header
	 */
	public ImageSource getSource()
	{
		if (source == null)
		{
			getHeader();
			if (header != null)
			{
				String xml = getField("Source");
				if (xml != null && xml.length() > 0 && xml.startsWith("<"))
				{
					// Convert the XML back
					source = (ImageSource) ImageSource.fromXML(xml);
				}
			}
		}
		return source;
	}

	/**
	 * @return The calibration specified in the results header
	 */
	public Calibration getCalibration()
	{
		if (calibration == null)
		{
			getHeader();
			if (header != null && header.length() > 0)
			{
				if (format == FileFormat.RAPID_STORM)
				{
					calibration.setDistanceUnit(DistanceUnit.NM);

					// RapidSTORM has a resolution attribute in the header in units of px m^-1
					Pattern pattern = Pattern.compile("resolution=\"([^ ]+) px m");
					Matcher match = pattern.matcher(header);
					if (match.find())
					{
						try
						{
							final float resolution = Float.parseFloat(match.group(1));
							if (Maths.isFinite(resolution) && resolution > 0)
							{
								final double nmPerPixel = (float) (1e9 / resolution);
								calibration = new Calibration();
								calibration.setNmPerPixel(nmPerPixel);
							}
						}
						catch (NumberFormatException e)
						{

						}
					}
				}
				else
				{
					String xml = getField("Calibration");
					if (xml != null && xml.length() > 0 && xml.startsWith("<"))
					{
						// Convert the XML back
						try
						{
							calibration = (Calibration) XmlUtils.fromXML(xml);
							calibration.validate();
							if (smlmVersion < 3)
							{
								// Support reading old objects that had an emCCD boolean
								calibration.setCameraTypeFromEmCCDField();
								// Previous version were always in pixels and counts
								calibration.setDistanceUnit(DistanceUnit.PIXEL);
								calibration.setIntensityUnit(IntensityUnit.COUNT);
							}
						}
						catch (ClassCastException ex)
						{
							ex.printStackTrace();
						}
						catch (Exception ex)
						{
							ex.printStackTrace();
						}
					}
					
					if (format == FileFormat.MALK)
					{
						if (calibration == null)
							calibration = new Calibration();
						calibration.setDistanceUnit(DistanceUnit.NM);
						calibration.setIntensityUnit(IntensityUnit.PHOTON);
					}
				}
			}

			// Calibration is a smart object so we can create an empty one
			if (calibration == null)
				calibration = new Calibration();
			// For debugging we can ensure that the calibration is not used incorrectly 
			calibration.setFieldMissingException(true);
		}
		return calibration;
	}

	public PSF getPSF()
	{
		if (psf == null)
		{
			getHeader();
			if (header != null && header.length() > 0)
			{
				String json = getField("PSF");
				if (json != null && json.length() > 0
				//&& json.startsWith("[")
				)
				{
					// Convert the JSON back
					try
					{
						PSF.Builder psfBuilder = PSF.newBuilder();
						JsonFormat.parser().merge(json, psfBuilder);
						psf = psfBuilder.build();
					}
					catch (InvalidProtocolBufferException e)
					{
						// This should be OK
						System.err.println("Unable to deserialise the PSF settings");
					}
				}

				if (psf == null)
				{
					// Guess base on type
					switch (format)
					{
						case NSTORM:
						case RAPID_STORM:
							// We currently only support two axis data
							psf = PSFHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D);
							break;

						// Note: Older GDSC results were all TwoAxisAndTheta
						case SMLM_BINARY:
						case SMLM_TABLE:
						case SMLM_TEXT:
							if (smlmVersion < 3)
								psf = PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
							break;

						case TSF_BINARY:
							// This is read separately so ignore
							break;

						case MALK:
							// Use the custom type with no additional PSF parameters
							psf = PSFHelper.create(PSFType.CUSTOM);
							break;

						case UNKNOWN:
						default:
							break;
					}
				}
			}
		}
		return psf;
	}

	/**
	 * @return The configuration specified in the results header
	 */
	public String getConfiguration()
	{
		if (configuration == null)
		{
			getHeader();
			configuration = getField("Configuration");
			if (configuration != null && configuration.length() > 0)
			{
				// Format the XML back
				configuration = XmlUtils.formatXml(configuration);
			}
		}
		return configuration;
	}

	private String getField(String name)
	{
		if (header != null)
		{
			Pattern pattern = Pattern.compile(name + " ([^\\n]+)");
			Matcher match = pattern.matcher(header);
			if (match.find())
			{
				return match.group(1);
			}
		}
		return "";
	}

	/**
	 * Read the results from the file. The file is read for each invocation of this method.
	 * 
	 * @return The peak results
	 */
	public MemoryPeakResults getResults()
	{
		getHeader();
		if (header == null || format == null || format == FileFormat.UNKNOWN)
			return null;

		// Always do this as we can guess the PSF based on the file type
		getPSF();

		// Use a switch statement with no break statements to fall through
		switch (format)
		{
			case SMLM_BINARY:
			case SMLM_TEXT:
			case MALK:
				// Read SMLM data. We do this for MALK files because we may have written them.
				getName();
				getSource();
				getBounds();
				getConfiguration();
				// RapidSTORM has calibration too 
			case RAPID_STORM:
				getCalibration();
			default:
				break;
		}

		MemoryPeakResults results = null;
		switch (format)
		{
			case SMLM_BINARY:
				results = readBinary();
				break;
			case SMLM_TEXT:
				results = readText();
				break;
			case SMLM_TABLE:
				results = readTable();
				break;
			case RAPID_STORM:
				results = readRapidSTORM();
				break;
			case NSTORM:
				results = readNSTORM();
				break;
			case MALK:
				results = readMALK();
				break;
			case TSF_BINARY:
				results = readTSF();
				break;
			default:
				break;
		}
		if (results != null)
		{
			results.trimToSize();

			if (!rawResults)
			{
				if (psf != null)
					simplifyPSF(results);

				// Convert to the preferred units if possible
				results.convertToPreferredUnits();
			}
		}
		return results;
	}

	private void simplifyPSF(MemoryPeakResults results)
	{
		switch (psf.getPsfType())
		{
			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				simplifyTwoAxisAndTheta(results);
				break;

			case TWO_AXIS_GAUSSIAN_2D:
				simplifyTwoAxis(results, false);
				break;

			case ASTIGMATIC_GAUSSIAN_2D:
			case ONE_AXIS_GAUSSIAN_2D:
			case UNRECOGNIZED:
			default:
				break;
		}
	}

	/**
	 * Simplify two axis and theta Gaussian 2D data.
	 *
	 * @param results
	 *            the results
	 */
	private void simplifyTwoAxisAndTheta(MemoryPeakResults results)
	{
		int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
		final int isx = indices[0];
		final int isy = indices[1];
		final int ia = PSFHelper.getGaussian2DAngleIndex(psf);

		// Determine if the angle is non-zero with asymmetric widths 
		if (results.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult peakResult)
			{
				return (peakResult.params[ia] != 0 && peakResult.params[isx] != peakResult.params[isy]);
			}
		}))
		{
			// Nothing to simplify
			return;
		}

		simplifyTwoAxis(results, true);
	}

	/**
	 * Simplify two axis Gaussian 2D data.
	 *
	 * @param results
	 *            the results
	 * @param removeTheta
	 *            the remove theta flag
	 */
	private void simplifyTwoAxis(MemoryPeakResults results, boolean removeTheta)
	{
		int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
		final int isx = indices[0];
		final int isy = indices[1];

		// Columns to remove
		int remove = (removeTheta) ? 1 : 0;
		// New PSF type
		PSFType psfType;

		// Determine if sy is redundant
		if (results.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult peakResult)
			{
				return (peakResult.params[isx] != peakResult.params[isy]);
			}
		}))
		{
			if (!removeTheta)
				// Already a TwoAxis Gaussian
				return;

			// Otherwise this was a TwoAxisAndTheta with 1 column to remove 
			// so it should be simplified
			psfType = PSFType.TWO_AXIS_GAUSSIAN_2D;
		}
		else
		{
			// sy is redundant so remove another column
			psfType = PSFType.ONE_AXIS_GAUSSIAN_2D;
			remove++;
		}

		// Update the PSF
		PSF.Builder builder = psf.toBuilder();
		builder.setPsfType(psfType);
		results.setPSF(psf = builder.build());

		// Update the results.
		// We can directly manipulate the params array
		final int newLength = results.getf(0).params.length - remove;
		if (!deviations)
		{
			results.forEach(new PeakResultProcedure()
			{
				public void execute(PeakResult peakResult)
				{
					peakResult.params = Arrays.copyOf(peakResult.params, newLength);
				}
			});
		}
		else
		{
			results.forEach(new PeakResultProcedure()
			{
				public void execute(PeakResult peakResult)
				{
					peakResult.params = Arrays.copyOf(peakResult.params, newLength);
					peakResult.paramsStdDev = Arrays.copyOf(peakResult.paramsStdDev, newLength);
				}
			});
		}
	}

	private int position;

	private MemoryPeakResults readBinary()
	{
		MemoryPeakResults results = createResults();

		// Units were added in version 3
		int nFields;
		boolean gaussian2Dformat;
		if (smlmVersion < 3)
		{
			nFields = 7;
			gaussian2Dformat = true;
			calibration.setIntensityUnit(IntensityUnit.COUNT);
			calibration.setDistanceUnit(DistanceUnit.PIXEL);
			calibration.setAngleUnit(AngleUnit.DEGREE);
		}
		else
		{
			gaussian2Dformat = false;
			// The number of fields should be within the PSF object
			nFields = new PeakResultsHelper(calibration, psf).getNames().length;
		}

		DataInputStream input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new DataInputStream(fis);

			// Seek to the start of the binary data by just reading the header again
			BinaryFilePeakResults.readHeader(input);

			// Format: [i]i[i]iifdf + n*f [+ n*f]
			// where [] are optional and n is the number of fields

			int length = BinaryFilePeakResults.getDataSize(deviations, readEndFrame, readId, nFields);
			byte[] buffer = new byte[length];

			int c = 0;
			final boolean convert = smlmVersion == 1;
			while (true) // Halted by the EOFException
			{
				// Note: Reading single strips seems fast enough at the moment.
				// This could be modified to read larger blocks of data if necessary.

				int bytesRead = input.read(buffer);
				if (bytesRead != length)
					break;

				position = 0;
				final int id = (readId) ? readInt(buffer) : 0;
				int peak = readInt(buffer);
				final int endPeak = (readEndFrame) ? readInt(buffer) : peak;
				int origX = readInt(buffer);
				int origY = readInt(buffer);
				float origValue = readFloat(buffer);
				double error = readDouble(buffer);
				float noise = readFloat(buffer);
				float[] params = readData(buffer, new float[nFields]);
				float[] paramsStdDev = (deviations) ? readData(buffer, new float[nFields]) : null;

				// Convert format which had the full Gaussian2D parameters array
				if (gaussian2Dformat)
				{
					// Convert old binary format with the amplitude to signal
					if (convert)
						params[Gaussian2DFunction.SIGNAL] *= 2 * Math.PI * params[Gaussian2DFunction.X_SD] *
								params[Gaussian2DFunction.Y_SD];
					params = mapGaussian2DFormatParams(params);
					paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
				}

				if (readId || readEndFrame)
					results.add(new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params,
							paramsStdDev, endPeak, id));
				else
					results.add(new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev));

				if (++c % 512 == 0)
					showProgress(channel);
			}
		}
		catch (EOFException e)
		{
			// Ignore
		}
		catch (IOException e)
		{
			// ignore
			e.printStackTrace();
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
		}
		return results;
	}

	private float[] mapGaussian2DFormatParams(float[] params)
	{
		final float[] p = mapGaussian2DFormat(params);
		// New format does not have the bias
		p[PeakResult.BACKGROUND] -= calibration.getBias();
		return p;
	}

	private float[] mapGaussian2DFormatDeviations(float[] params)
	{
		return (params != null) ? mapGaussian2DFormat(params) : null;
	}

	private static float[] mapGaussian2DFormat(float[] params)
	{
		final float[] p = new float[nTwoAxisAndTheta];
		p[PeakResult.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND];
		p[PeakResult.INTENSITY] = params[Gaussian2DFunction.SIGNAL];
		p[PeakResult.X] = params[Gaussian2DFunction.X_POSITION];
		p[PeakResult.Y] = params[Gaussian2DFunction.Y_POSITION];
		p[isx] = params[Gaussian2DFunction.X_SD];
		p[isy] = params[Gaussian2DFunction.Y_SD];
		p[ia] = params[Gaussian2DFunction.SHAPE];
		return p;
	}

	private MemoryPeakResults createResults()
	{
		MemoryPeakResults results = new MemoryPeakResults();
		results.setName(name);
		results.setSource(source);
		results.setBounds(bounds);
		results.setCalibration(calibration);
		results.setConfiguration(configuration);
		results.setPSF(psf);
		return results;
	}

	private final int readInt(byte[] buffer)
	{
		int ch1 = buffer[position] & 0xff;
		int ch2 = buffer[position + 1] & 0xff;
		int ch3 = buffer[position + 2] & 0xff;
		int ch4 = buffer[position + 3] & 0xff;
		position += 4;
		return ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
	}

	private final long readLong(byte[] buffer)
	{
		final long l = (((long) buffer[position + 0] << 56) + ((long) (buffer[position + 1] & 255) << 48) +
				((long) (buffer[position + 2] & 255) << 40) + ((long) (buffer[position + 3] & 255) << 32) +
				((long) (buffer[position + 4] & 255) << 24) + ((buffer[position + 5] & 255) << 16) +
				((buffer[position + 6] & 255) << 8) + ((buffer[position + 7] & 255) << 0));
		position += 8;
		return l;
	}

	private final float readFloat(byte[] buffer)
	{
		return Float.intBitsToFloat(readInt(buffer));
	}

	private final double readDouble(byte[] buffer)
	{
		return Double.longBitsToDouble(readLong(buffer));
	}

	private void showProgress(FileChannel channel) throws IOException
	{
		if (tracker != null)
		{
			tracker.progress(channel.position(), channel.size());
			if (tracker.isEnded())
				// Throw an IOException and it will be caught and ignored by all the file reading methods
				throw new IOException("File read was cancelled");
		}
	}

	private float[] readData(byte[] buffer, float[] params) throws IOException
	{
		for (int i = 0; i < params.length; i++)
			params[i] = readFloat(buffer);
		return params;
	}

	private MemoryPeakResults readText()
	{
		MemoryPeakResults results = createResults();

		// Units were added in version 3
		int nFields;
		if (smlmVersion < 3)
		{
			nFields = 7;
			calibration.setIntensityUnit(IntensityUnit.COUNT);
			calibration.setDistanceUnit(DistanceUnit.PIXEL);
			calibration.setAngleUnit(AngleUnit.DEGREE);

			// Note that in older versions the background included the bias. 
			// The bias is not included in the table and so the user will have to 
			// add this manually. They will also have to add the gain.
		}
		else
		{
			// The number of fields should be within the PSF object
			nFields = new PeakResultsHelper(calibration, psf).getNames().length;
		}

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int errors = 0;

			// Skip the header
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;

				if (line.charAt(0) != '#')
				{
					// This is the first record
					if (!addPeakResult(results, line, smlmVersion, nFields))
						errors = 1;
					break;
				}
			}

			int c = 0;
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				if (line.charAt(0) == '#')
					continue;

				if (!addPeakResult(results, line, smlmVersion, nFields))
				{
					if (++errors >= 10)
					{
						break;
					}
				}

				if (++c % 512 == 0)
					showProgress(channel);
			}
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
		}
		return results;
	}

	private boolean addPeakResult(MemoryPeakResults results, String line, int version, int nFields)
	{
		PeakResult result;
		switch (version)
		{
			case 3:
				// Version 3 has variable PSF parameter fields
				result = (deviations) ? createPeakResultDeviationsV3(line, nFields) : createPeakResultV3(line, nFields);
				break;

			case 2:
				result = (deviations) ? createPeakResultDeviationsV2(line) : createPeakResultV2(line);
				break;

			case 1:
			default:
				result = (deviations) ? createPeakResultDeviationsV1(line) : createPeakResultV1(line);
		}

		if (result != null)
		{
			results.add(result);
			return true;
		}
		return false;
	}

	private PeakResult createPeakResultV1(String line)
	{
		try
		{
			float[] params = new float[7];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(tabPattern);
				scanner.useLocale(Locale.US);
				int id = 0, endPeak = 0;
				if (readId)
					id = scanner.nextInt();
				int peak = scanner.nextInt();
				if (readEndFrame)
					endPeak = scanner.nextInt();
				int origX = scanner.nextInt();
				int origY = scanner.nextInt();
				float origValue = scanner.nextFloat();
				double error = scanner.nextDouble();
				float noise = scanner.nextFloat();
				float signal = scanner.nextFloat(); // Ignored but must be read
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
				}
				scanner.close();
				params[Gaussian2DFunction.SIGNAL] = signal;
				params = mapGaussian2DFormatParams(params);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, null);
			}
			else
			{
				// Code using split and parse
				String[] fields = tabPattern.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double error = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				float signal = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
				}
				params[Gaussian2DFunction.SIGNAL] = signal;
				params = mapGaussian2DFormatParams(params);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, null);
			}
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		catch (IndexOutOfBoundsException e)
		{
		}
		catch (NumberFormatException e)
		{
		}
		return null;
	}

	private PeakResult createPeakResultDeviationsV1(String line)
	{
		try
		{
			float[] params = new float[7];
			float[] paramsStdDev = new float[7];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(tabPattern);
				scanner.useLocale(Locale.US);
				int id = 0, endPeak = 0;
				if (readId)
					id = scanner.nextInt();
				int peak = scanner.nextInt();
				if (readEndFrame)
					endPeak = scanner.nextInt();
				int origX = scanner.nextInt();
				int origY = scanner.nextInt();
				float origValue = scanner.nextFloat();
				double error = scanner.nextDouble();
				float noise = scanner.nextFloat();
				float signal = scanner.nextFloat(); // Ignored but must be read
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
					paramsStdDev[i] = scanner.nextFloat();
				}
				params[Gaussian2DFunction.SIGNAL] = signal;
				scanner.close();
				params = mapGaussian2DFormatParams(params);
				paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
			}
			else
			{
				// JUnit test shows this is faster than the scanner

				// Code using split and parse
				String[] fields = tabPattern.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double error = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				float signal = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
					paramsStdDev[i] = Float.parseFloat(fields[j++]);
				}
				params[Gaussian2DFunction.SIGNAL] = signal;
				params = mapGaussian2DFormatParams(params);
				paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
			}
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		catch (IndexOutOfBoundsException e)
		{
		}
		catch (NumberFormatException e)
		{
		}
		return null;
	}

	private PeakResult createPeakResultV2(String line)
	{
		try
		{
			float[] params = new float[7];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(tabPattern);
				scanner.useLocale(Locale.US);
				int id = 0, endPeak = 0;
				if (readId)
					id = scanner.nextInt();
				int peak = scanner.nextInt();
				if (readEndFrame)
					endPeak = scanner.nextInt();
				int origX = scanner.nextInt();
				int origY = scanner.nextInt();
				float origValue = scanner.nextFloat();
				double error = scanner.nextDouble();
				float noise = scanner.nextFloat();
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
				}
				scanner.close();
				params = mapGaussian2DFormatParams(params);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, null);
			}
			else
			{
				// Code using split and parse
				String[] fields = tabPattern.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double error = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
				}
				params = mapGaussian2DFormatParams(params);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, null);
			}
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		catch (IndexOutOfBoundsException e)
		{
		}
		catch (NumberFormatException e)
		{
		}
		return null;
	}

	private PeakResult createPeakResultDeviationsV2(String line)
	{
		try
		{
			float[] params = new float[7];
			float[] paramsStdDev = new float[7];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(tabPattern);
				scanner.useLocale(Locale.US);
				int id = 0, endPeak = 0;
				if (readId)
					id = scanner.nextInt();
				int peak = scanner.nextInt();
				if (readEndFrame)
					endPeak = scanner.nextInt();
				int origX = scanner.nextInt();
				int origY = scanner.nextInt();
				float origValue = scanner.nextFloat();
				double error = scanner.nextDouble();
				float noise = scanner.nextFloat();
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
					paramsStdDev[i] = scanner.nextFloat();
				}
				scanner.close();
				params = mapGaussian2DFormatParams(params);
				paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
			}
			else
			{
				// JUnit test shows this is faster than the scanner

				// Code using split and parse
				String[] fields = tabPattern.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double error = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
					paramsStdDev[i] = Float.parseFloat(fields[j++]);
				}
				params = mapGaussian2DFormatParams(params);
				paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
			}
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		catch (IndexOutOfBoundsException e)
		{
		}
		catch (NumberFormatException e)
		{
		}
		return null;
	}

	private PeakResult createPeakResultV3(String line, int nFields)
	{
		try
		{
			float[] params = new float[nFields];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(tabPattern);
				scanner.useLocale(Locale.US);
				int id = 0, endPeak = 0;
				if (readId)
					id = scanner.nextInt();
				int peak = scanner.nextInt();
				if (readEndFrame)
					endPeak = scanner.nextInt();
				int origX = scanner.nextInt();
				int origY = scanner.nextInt();
				float origValue = scanner.nextFloat();
				double error = scanner.nextDouble();
				float noise = scanner.nextFloat();
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
				}
				scanner.close();
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, null);
			}
			else
			{
				// Code using split and parse
				String[] fields = tabPattern.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double error = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
				}
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, null);
			}
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		catch (IndexOutOfBoundsException e)
		{
		}
		catch (NumberFormatException e)
		{
		}
		return null;
	}

	private PeakResult createPeakResultDeviationsV3(String line, int nFields)
	{
		try
		{
			float[] params = new float[nFields];
			float[] paramsStdDev = new float[nFields];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(tabPattern);
				scanner.useLocale(Locale.US);
				int id = 0, endPeak = 0;
				if (readId)
					id = scanner.nextInt();
				int peak = scanner.nextInt();
				if (readEndFrame)
					endPeak = scanner.nextInt();
				int origX = scanner.nextInt();
				int origY = scanner.nextInt();
				float origValue = scanner.nextFloat();
				double error = scanner.nextDouble();
				float noise = scanner.nextFloat();
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
					paramsStdDev[i] = scanner.nextFloat();
				}
				scanner.close();
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
			}
			else
			{
				// JUnit test shows this is faster than the scanner

				// Code using split and parse
				String[] fields = tabPattern.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double error = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
					paramsStdDev[i] = Float.parseFloat(fields[j++]);
				}
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
			}
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		catch (IndexOutOfBoundsException e)
		{
		}
		catch (NumberFormatException e)
		{
		}
		return null;
	}

	private MemoryPeakResults readTable()
	{
		MemoryPeakResults results = createResults();
		results.setName(new File(filename).getName());

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int errors = 0;

			// Skip over the single line header
			String header = input.readLine();

			// V1: had the Signal and Amplitude. Parameters 
			// V2: have only the Signal.
			// V3: Has variable columns with units for the PSF parameters. Signal was renamed to Intensity.
			int version;
			int nFields = 0;
			if (header.contains("Signal"))
				version = 1;
			else if (header.contains("Amplitude"))
				version = 2;
			else
			{
				version = 3;
				// Get the number of data fields by counting the standard fields
				String[] columns = header.split("\t");
				int field = 0;
				if (readId)
					field++; // ID #
				if (readSource)
					field++; // Source
				field++; // Frame
				if (readEndFrame)
					field++; // End frame
				field++; // origX
				field++; // origY
				field++; // origValue
				field++; // error
				field++; // noise
				field++; // SNR

				// The remaining fields are PSF parameters with the exception of the final precision field

				nFields = columns.length - field;
				if (columns[columns.length - 1].contains("Precision"))
					nFields--;
				if (deviations)
				{
					nFields /= 2;
				}

				// We can guess part of the calibration.
				if (calibration == null)
					calibration = new Calibration();
				int jump = (deviations) ? 2 : 1;
				// field is currently on Background
				calibration.setIntensityUnit(UnitHelper.guessIntensityUnitFromShortName(extractUnit(columns[field])));
				field += jump; // Move to Intensity
				field += jump; // Move to X
				calibration.setDistanceUnit(UnitHelper.guessDistanceUnitFromShortName(extractUnit(columns[field])));
				field += jump; // Move to Y
				field += jump; // Move to Z
				// The angle may be used in fields above the standard ones
				while (field < columns.length)
				{
					field += jump;
					AngleUnit u = UnitHelper.guessAngleUnitFromShortName(extractUnit(columns[field]));
					if (u != null)
					{
						calibration.setAngleUnit(u);
						break;
					}
				}
			}

			int c = 0;
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;

				if (!addTableResult(results, line, version, nFields))
				{
					if (++errors >= 10)
					{
						break;
					}
				}

				if (++c % 512 == 0)
					showProgress(channel);
			}
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
		}

		return results;
	}

	/**
	 * Extract the unit from a string. The unit is the first string within brackets, e.g. (unit) would return unit.
	 *
	 * @param string
	 *            the string
	 * @return the unit string (or null)
	 */
	public static String extractUnit(String string)
	{
		if (string != null)
		{
			int beginIndex = string.indexOf('(');
			if (beginIndex >= 0)
			{
				int endIndex = string.indexOf(')');
				if (endIndex > beginIndex)
				{
					return string.substring(beginIndex, endIndex);
				}
			}
		}
		return null;
	}

	private boolean addTableResult(MemoryPeakResults results, String line, int version, int nFields)
	{
		final PeakResult result;
		switch (version)
		{
			case 3:
				result = createTableResultV3(line, nFields);
				break;

			case 2:
				result = createTableResultV2(line);
				break;

			case 1:
			default:
				result = createTableResultV1(line);
		}
		if (result != null)
		{
			// Extract the source & bounds from the Source column
			if (results.size() == 0 && readSource)
			{
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(tabPattern);
				scanner.useLocale(Locale.US);
				if (readId)
					scanner.nextInt();
				String source = scanner.next();
				scanner.close();

				if (source.contains(": "))
				{
					String[] fields = source.split(": ");
					results.setName(fields[0]);
					// Bounds is formatted as 'xN yN wN hN'
					Pattern pattern = Pattern.compile("x(\\d+) y(\\d+) w(\\d+) h(\\d+)");
					Matcher match = pattern.matcher(fields[1]);
					if (match.find())
					{
						int x = Integer.parseInt(match.group(1));
						int y = Integer.parseInt(match.group(2));
						int w = Integer.parseInt(match.group(3));
						int h = Integer.parseInt(match.group(4));
						results.setBounds(new Rectangle(x, y, w, h));
					}
				}
				else
				{
					results.setName(source);
				}
			}

			results.add(result);
			return true;
		}
		return false;
	}

	private PeakResult createTableResultV1(String line)
	{
		// Text file with fields:
		// [#]
		// [Source]
		// Peak
		// [End Peak]
		// origX
		// origY
		// origValue
		// Error
		// Noise
		// Signal
		// SNR
		// Background
		// [+/-]
		// Amplitude
		// [+/-]
		// Angle
		// [+/-]
		// X
		// [+/-]
		// Y
		// [+/-]
		// X SD
		// [+/-]
		// Y SD
		// [+/-]
		// [Precision] 
		try
		{
			Scanner scanner = new Scanner(line);
			scanner.useDelimiter(tabPattern);
			scanner.useLocale(Locale.US);
			int id = 0, endPeak = 0;
			if (readId)
				id = scanner.nextInt();
			if (readSource)
				scanner.next();
			int peak = scanner.nextInt();
			if (readEndFrame)
				endPeak = scanner.nextInt();
			int origX = scanner.nextInt();
			int origY = scanner.nextInt();
			float origValue = scanner.nextFloat();
			double error = scanner.nextDouble();
			float noise = scanner.nextFloat();
			@SuppressWarnings("unused")
			float signal = scanner.nextFloat(); // Ignored but must be read
			@SuppressWarnings("unused")
			float snr = scanner.nextFloat(); // Ignored but must be read
			float[] params = new float[7];
			float[] paramsStdDev = (deviations) ? new float[7] : null;
			for (int i = 0; i < params.length; i++)
			{
				params[i] = scanner.nextFloat();
				if (deviations)
					paramsStdDev[i] = scanner.nextFloat();
			}
			scanner.close();

			// For the new format we store the signal (not the amplitude).
			// Convert the amplitude into a signal
			params[Gaussian2DFunction.SIGNAL] *= 2 * Math.PI * params[Gaussian2DFunction.X_SD] *
					params[Gaussian2DFunction.Y_SD];

			params = mapGaussian2DFormatParams(params);
			paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);

			if (readId || readEndFrame)
				return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev,
						endPeak, id);
			else
				return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		return null;
	}

	private PeakResult createTableResultV2(String line)
	{
		// Text file with fields:
		// [#]
		// [Source]
		// Peak
		// [End Peak]
		// origX
		// origY
		// origValue
		// Error
		// Noise
		// SNR
		// Background
		// [+/-]
		// Signal
		// [+/-]
		// Angle
		// [+/-]
		// X
		// [+/-]
		// Y
		// [+/-]
		// X SD
		// [+/-]
		// Y SD
		// [+/-]
		// [Precision] 
		try
		{
			Scanner scanner = new Scanner(line);
			scanner.useDelimiter(tabPattern);
			scanner.useLocale(Locale.US);
			int id = 0, endPeak = 0;
			if (readId)
				id = scanner.nextInt();
			if (readSource)
				scanner.next();
			int peak = scanner.nextInt();
			if (readEndFrame)
				endPeak = scanner.nextInt();
			int origX = scanner.nextInt();
			int origY = scanner.nextInt();
			float origValue = scanner.nextFloat();
			double error = scanner.nextDouble();
			float noise = scanner.nextFloat();
			@SuppressWarnings("unused")
			float snr = scanner.nextFloat(); // Ignored but must be read
			float[] params = new float[7];
			float[] paramsStdDev = (deviations) ? new float[7] : null;
			for (int i = 0; i < params.length; i++)
			{
				params[i] = scanner.nextFloat();
				if (deviations)
					paramsStdDev[i] = scanner.nextFloat();
			}
			scanner.close();

			params = mapGaussian2DFormatParams(params);
			paramsStdDev = mapGaussian2DFormatDeviations(paramsStdDev);

			if (readId || readEndFrame)
				return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev,
						endPeak, id);
			else
				return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		return null;
	}

	private PeakResult createTableResultV3(String line, int nFields)
	{
		// Text file with fields:
		// [#]
		// [Source]
		// Peak
		// [End Peak]
		// origX
		// origY
		// origValue
		// Error
		// Noise[ (units)]
		// SNR
		// Background[ (units)]
		// [+/-]
		// Intensity[ (units)]
		// [+/-]
		// X[ (units)]
		// [+/-]
		// Y[ (units)]
		// [+/-]
		// Z[ (units)]
		// [+/-]
		// Repeated:
		//   Field[ (units)]
		//   [+/-]
		// [Precision] 
		try
		{
			Scanner scanner = new Scanner(line);
			scanner.useDelimiter(tabPattern);
			scanner.useLocale(Locale.US);
			int id = 0, endPeak = 0;
			if (readId)
				id = scanner.nextInt();
			if (readSource)
				scanner.next();
			int peak = scanner.nextInt();
			if (readEndFrame)
				endPeak = scanner.nextInt();
			int origX = scanner.nextInt();
			int origY = scanner.nextInt();
			float origValue = scanner.nextFloat();
			double error = scanner.nextDouble();
			float noise = scanner.nextFloat();
			@SuppressWarnings("unused")
			float snr = scanner.nextFloat(); // Ignored but must be read
			float[] params = new float[nFields];
			float[] paramsStdDev;
			if (deviations)
			{
				paramsStdDev = new float[params.length];
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
					paramsStdDev[i] = scanner.nextFloat();
				}
			}
			else
			{
				paramsStdDev = null;
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
				}
			}
			scanner.close();
			if (readId || readEndFrame)
				return new ExtendedPeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev,
						endPeak, id);
			else
				return new PeakResult(peak, origX, origY, origValue, error, noise, params, paramsStdDev);
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		return null;
	}

	private MemoryPeakResults readRapidSTORM()
	{
		MemoryPeakResults results = createResults();
		results.setName(new File(filename).getName());

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int errors = 0;

			// Skip the header
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				if (line.charAt(0) != '#')
				{
					// This is the first record
					if (!addRapidSTORMResult(results, line))
						errors = 1;
					break;
				}
			}

			int c = 0;
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				if (line.charAt(0) == '#')
					continue;

				if (!addRapidSTORMResult(results, line))
				{
					if (++errors >= 10)
					{
						break;
					}
				}

				if (++c % 512 == 0)
					showProgress(channel);
			}
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
		}
		return results;
	}

	private boolean addRapidSTORMResult(MemoryPeakResults results, String line)
	{
		PeakResult result = createRapidSTORMResult(line);
		if (result != null)
		{
			results.add(result);
			return true;
		}
		return false;
	}

	private PeakResult createRapidSTORMResult(String line)
	{
		// Text file with fields:
		//   X (nm)
		//   Y (nm)
		//   Frame
		//   Amplitude*
		//   sx^2 (pm^2)
		//   sy^2 (pm^2)
		//   2 kernel improvement
		//   Fit residues chi square
		// *Note that the RapidSTORM Amplitude is the signal.
		try
		{
			Scanner scanner = new Scanner(line);
			scanner.useDelimiter(spacePattern);
			scanner.useLocale(Locale.US);
			float x = scanner.nextFloat();
			float y = scanner.nextFloat();
			final int peak = scanner.nextInt();
			final float signal = scanner.nextFloat();
			final float sx2 = scanner.nextFloat();
			final float sy2 = scanner.nextFloat();
			@SuppressWarnings("unused")
			final float kernelImprovement = scanner.nextFloat();
			final double error = scanner.nextDouble();
			scanner.close();

			// Convert from pm^2 to nm
			float sx = (float) (Math.sqrt(sx2) * 1000);
			float sy = (float) (Math.sqrt(sy2) * 1000);

			float[] params = new float[nTwoAxis];
			params[PeakResult.INTENSITY] = signal;
			params[PeakResult.X] = x;
			params[PeakResult.Y] = y;
			params[isx] = sx;
			params[isy] = sy;

			// Store the signal as the original value
			return new PeakResult(peak, (int) x, (int) y, signal, error, 0.0f, params, null);
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		return null;
	}

	private MemoryPeakResults readNSTORM()
	{
		MemoryPeakResults results = createResults();
		results.setName(new File(filename).getName());

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int errors = 0;

			// Skip the single line header
			input.readLine();

			int c = 0;
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;

				if (!addNSTORMResult(results, line))
				{
					if (++errors >= 10)
					{
						break;
					}
				}

				if (++c % 512 == 0)
					showProgress(channel);
			}
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
		}

		// The following relationship holds when length == 1:
		// intensity = height * 2 * pi * sd0 * sd1 / pixel_pitch^2
		// => Pixel_pitch = sqrt(height * 2 * pi * sd0 * sd1 / intensity)
		// Try and create a calibration
		final Statistics pixelPitch = new Statistics();
		results.forEach(new PeakResultProcedureX()
		{
			final double twoPi = 2 * Math.PI;

			public boolean execute(PeakResult p)
			{
				if (p.getFrame() == p.getEndFrame())
				{
					float height = p.origValue;
					float intensity = p.params[PeakResult.INTENSITY];
					float sd0 = p.params[isx];
					float sd1 = p.params[isy];
					pixelPitch.add(Math.sqrt(height * twoPi * sd0 * sd1 / intensity));
					// Stop when we have enough for a good guess
					return (pixelPitch.getN() > 100);
				}
				return false;
			}
		});

		// TODO - Support all the NSTORM formats: one-axis, two-axis, rotated, 3D.
		// Is this information in the header?
		// We could support setting the PSF as a Gaussian2D with one/two axis SD.
		// This would mean updating all the result params if it is a one axis PSF.
		// For now just record it as a 2 axis PSF.

		// Create a calibration
		calibration = new Calibration();

		// Q. Is NSTORM in photons?
		calibration.setIntensityUnit(IntensityUnit.COUNT);
		calibration.setDistanceUnit(DistanceUnit.NM);

		if (pixelPitch.getN() > 0)
		{
			final double nmPerPixel = pixelPitch.getMean();
			calibration.setNmPerPixel(nmPerPixel);
		}

		results.setCalibration(calibration);

		return results;
	}

	private boolean addNSTORMResult(MemoryPeakResults results, String line)
	{
		PeakResult result = createNSTORMResult(line);
		if (result != null)
		{
			results.add(result);
			return true;
		}
		return false;
	}

	// So that the fields can be named  
	@SuppressWarnings("unused")
	private PeakResult createNSTORMResult(String line)
	{
		// Note that the NSTORM file contains traced molecules hence the Frame
		// and Length fields.

		// Text file with tab delimited fields:
		//Channel Name: This is the name of the Channel where the molecule was
		//              detected.
		//X: The X position of the centroid of the molecule in
		//   nanometers. Similar to the conventional image, molecules
		//   positions in the image are relative to the upper left corner of
		//   the image.
		//Y: The Y position of the centroid of the molecule in
		//   nanometers. Similar to the conventional image, molecules
		//   positions in the image are relative to the upper left corner of
		//   the image.
		//Xc: The X position of the centroid of the molecule (in
		//    nanometers) with drift correction applied. If no drift
		//    correction was applied to this data then Xc= X.
		//Yc: The Y position of the centroid of the molecule (in
		//    nanometers) with drift correction applied. If no drift
		//    correction was applied to this data then Yc= Y.
		//Height: Intensity of the peak height in the detection frame (after the
		//        detection process)
		//Area: Volume under the peak. Units are intensity * pixel^2
		//Width: Square Root(Wx*Wy)
		//Phi: The Angle of the molecule. This is the axial angle for 2D
		//     and distance from Z calibration curve in nm in Wx,Wy
		//     space for 3D)
		//Ax: Axial ratio of Wy/Wx
		//BG: The local background for the molecule
		//I: Accumulated intensity
		//Frame: The sequential frame number where the molecule was first
		//       detected.
		//Length: The number of consecutive frames the molecule was
		//        detected.
		//Link: For Internal Use only
		//Valid: For Internal Use only
		//Z: The Z coordinate of the molecule in Nanometers (origin is
		//   the cover glass).
		//Zc: The Z position of the molecule (in nanometers) with drift
		//    correction applied. If no drift correction was applied to this
		//    data then Zc= Z.		
		try
		{
			Scanner scanner = new Scanner(line);
			scanner.useDelimiter(tabPattern);
			scanner.useLocale(Locale.US);
			String channelName = scanner.next();
			float x = scanner.nextFloat();
			float y = scanner.nextFloat();
			float xc = scanner.nextFloat();
			float yc = scanner.nextFloat();
			float height = scanner.nextFloat();
			float area = scanner.nextFloat();
			float width = scanner.nextFloat();
			float phi = scanner.nextFloat();
			float ax = scanner.nextFloat();
			float bg = scanner.nextFloat();
			float i = scanner.nextFloat();
			int frame = scanner.nextInt();
			int length = scanner.nextInt();
			// These are not needed
			//float link = scanner.nextFloat();
			//float valid = scanner.nextFloat();
			//float z = scanner.nextFloat();
			//float zc = scanner.nextFloat();
			scanner.close();

			// The coordinates are in nm
			// The values are in ADUs. The area value is the signal.

			// The following relationship holds when length == 1:
			// Area = Height * 2 * pi * (Width / (pixel_pitch*2) )^2
			// => Pixel_pitch = 0.5 * Width / sqrt(Area / (Height * 2 * pi))

			float[] params = new float[nTwoAxis];
			params[PeakResult.BACKGROUND] = bg;
			//params[Gaussian2DFunction.ANGLE] = ax;
			params[PeakResult.INTENSITY] = area;
			params[PeakResult.X] = xc;
			params[PeakResult.Y] = yc;

			// Convert width (2*SD) to SD
			width /= 2f;

			// Convert to separate XY widths using the axial ratio
			if (ax == 1)
			{
				params[isx] = width;
				params[isy] = width;
			}
			else
			{
				// Ensure the axial ratio is long/short
				if (ax < 1)
					ax = 1.0f / ax;
				double a = Math.sqrt(ax);

				params[isx] = (float) (width * a);
				params[isy] = (float) (width / a);
			}

			// Store the signal as the original value
			return new ExtendedPeakResult(frame, (int) xc, (int) yc, height, 0.0, 0.0f, params, null,
					frame + length - 1, 0);
		}
		catch (InputMismatchException e)
		{
			// Ignore
		}
		catch (NoSuchElementException e)
		{
			// Ignore
		}
		return null;
	}

	private MemoryPeakResults readMALK()
	{
		MemoryPeakResults results = createResults();
		if (Utils.isNullOrEmpty(name))
			results.setName(new File(filename).getName());

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int errors = 0;

			// Skip the header
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;

				if (line.charAt(0) != '#')
				{
					// This is the first record
					if (!addMALKResult(results, line))
						errors = 1;
					break;
				}
			}

			int c = 0;
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;
				if (line.charAt(0) == '#')
					continue;

				if (!addMALKResult(results, line))
				{
					if (++errors >= 10)
					{
						break;
					}
				}

				if (++c % 512 == 0)
					showProgress(channel);
			}
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			try
			{
				if (input != null)
					input.close();
			}
			catch (IOException e)
			{
				// Ignore
			}
		}

		// Set default calibration for MALK format.
		// The calibration may not be null if this was a GDSC MALK file since that has a header.
		if (calibration == null)
		{
			calibration = new Calibration();
			// Default assumption is nm
			calibration.setDistanceUnit(DistanceUnit.NM);
			// MALK uses photons
			calibration.setIntensityUnit(IntensityUnit.PHOTON);
		}

		return results;
	}

	private boolean addMALKResult(MemoryPeakResults results, String line)
	{
		PeakResult result = createMALKResult(line);
		if (result != null)
		{
			results.add(result);
			return true;
		}
		return false;
	}

	private boolean isMALKFormat(String firstLine)
	{
		// The MALK file format is very simple: X,Y,T,Signal
		String[] fields = whitespacePattern.split(firstLine);
		if (fields.length != 4)
			return false;
		if (createMALKResult(firstLine) == null)
			return false;
		return true;
	}

	private PeakResult createMALKResult(String line)
	{
		try
		{
			float[] params = new float[PeakResult.STANDARD_PARAMETERS];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter(whitespacePattern);
				scanner.useLocale(Locale.US);
				params[PeakResult.X] = scanner.nextFloat();
				params[PeakResult.Y] = scanner.nextFloat();
				int peak = scanner.nextInt();
				params[PeakResult.INTENSITY] = scanner.nextFloat();
				scanner.close();

				return new PeakResult(peak, 0, 0, 0, 0, 0, params, null);
			}
			else
			{
				// Code using split and parse
				String[] fields = whitespacePattern.split(line);

				params[PeakResult.X] = Float.parseFloat(fields[0]);
				params[PeakResult.Y] = Float.parseFloat(fields[1]);
				int peak = Integer.parseInt(fields[2]);
				params[PeakResult.INTENSITY] = Float.parseFloat(fields[3]);

				return new PeakResult(peak, 0, 0, 0, 0, 0, params, null);
			}
		}
		catch (InputMismatchException e)
		{
		}
		catch (NoSuchElementException e)
		{
		}
		catch (IndexOutOfBoundsException e)
		{
		}
		catch (NumberFormatException e)
		{
		}
		return null;
	}

	private MemoryPeakResults readTSF()
	{
		TSFPeakResultsReader reader = new TSFPeakResultsReader(filename);
		reader.setOptions(options);
		return reader.read();
	}

	/**
	 * @return the tracker
	 */
	public TrackProgress getTracker()
	{
		return tracker;
	}

	/**
	 * @param tracker
	 *            the tracker to set
	 */
	public void setTracker(TrackProgress tracker)
	{
		this.tracker = tracker;
	}

	public boolean isUseScanner()
	{
		return useScanner;
	}

	public void setUseScanner(boolean useScanner)
	{
		this.useScanner = useScanner;
	}

	public boolean isRawResults()
	{
		return rawResults;
	}

	public void setRawResults(boolean rawResults)
	{
		this.rawResults = rawResults;
	}

	/**
	 * Gets the options for reading the results. Allows specific file formats to provide options for how to read the
	 * data.
	 *
	 * @return the options
	 */
	public ResultOption[] getOptions()
	{
		getHeader();
		if (header == null || format == null)
			return null;
		if (format == FileFormat.TSF_BINARY)
		{
			TSFPeakResultsReader reader = new TSFPeakResultsReader(filename);
			return reader.getOptions();
		}
		return null;
	}

	/**
	 * Sets the options for reading the results.
	 *
	 * @param options
	 *            the new options for reading the results
	 */
	public void setOptions(ResultOption[] options)
	{
		this.options = options;
	}
}
