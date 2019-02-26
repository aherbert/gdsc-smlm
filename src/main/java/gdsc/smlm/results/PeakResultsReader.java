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
import java.util.InputMismatchException;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.UnicodeReader;
import uk.ac.sussex.gdsc.core.utils.XmlUtils;

import gdsc.smlm.utils.XStreamXmlUtils;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Reads the fit results from file
 */
public class PeakResultsReader
{
	private boolean useScanner = false;

	private String filename;
	private String header = null;
	private FileFormat format;
	private String version;
	private String name = null;
	private ImageSource source = null;
	private Rectangle bounds = null;
	private Calibration calibration = null;
	private String configuration = null;
	private TrackProgress tracker = null;

	private boolean deviations, readEndFrame, readId, readSource;
	// Original file data contains signal and amplitude
	private boolean readAmplitude = true;

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
						if (line.contains("origX\torigY\torigValue\tError\tNoise\tSignal\tSNR\tBackground"))
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

				guessFormat();
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

	private void guessFormat()
	{
		format = FileFormat.UNKNOWN;

		// This cannot be done for empty header
		if (header.length() == 0)
			return;

		// Check for Nikon NSTORM header
		if (header.contains("Channel Name"))
		{
			format = FileFormat.NSTORM;
		}
		else if (header.contains("origX\torigY\torigValue\tError\tNoise\tSignal\tSNR\tBackground"))
		{
			format = FileFormat.SMLM_TABLE;
			// Header contains:
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
			readId = header.startsWith("#");
			readEndFrame = header.contains("\tEnd ");
			deviations = header.contains("\t+/-\t");
			readSource = header.contains("Source\t");
		}
		// Check for RapidSTORM stuff in the header
		else if (header.contains("<localizations "))
		{
			format = FileFormat.RAPID_STORM;
		}
		else
		{
			// Assume SMLM format

			// Extract information about the file format
			if (version.length() > 0)
			{
				format = (version.contains("Binary")) ? FileFormat.SMLM_BINARY : FileFormat.SMLM_TEXT;
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
				if (version.contains(".V2"))
				{
					readAmplitude = false;
				}
			}
			else
			{
				// The original files did not have the version tag
				format = (version.contains("Binary")) ? FileFormat.SMLM_BINARY : FileFormat.SMLM_TEXT;
				deviations = header.contains((format == FileFormat.SMLM_BINARY) ? "iiiifdfffffffffffffff" : "+/-");
				readEndFrame = readId = false;
			}
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
					// RapidSTORM has a resolution attribute in the header in units of px m^-1
					Pattern pattern = Pattern.compile("resolution=\"([^ ]+) px m");
					Matcher match = pattern.matcher(header);
					if (match.find())
					{
						try
						{
							final float resolution = Float.parseFloat(match.group(1));
							final double nmPerPixel = (float) (1e9 / resolution);
							calibration = new Calibration(nmPerPixel, 1, 0);
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
						calibration = (Calibration) XStreamXmlUtils.fromXML(xml);
					}
				}
			}
		}
		return calibration;
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

		// Use a switch statement with no break statements to fall through
		switch (format)
		{
			case SMLM_BINARY:
			case SMLM_TEXT:
				// Read SMLM data
				getName();
				getSource();
				getBounds();
				getConfiguration();
			case RAPID_STORM:
				// RapidSTORM has calibration too 
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
			default:
				break;
		}
		if (results != null)
			results.trimToSize();
		return results;
	}

	private int position;

	private MemoryPeakResults readBinary()
	{
		MemoryPeakResults results = createResults();

		DataInputStream input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new DataInputStream(fis);

			// Seek to the start of the binary data by just reading the header again
			BinaryFilePeakResults.readHeader(input);

			// Format: [i]i[i]iifdffffffff[fffffff]
			// where [] are optional

			int length = BinaryFilePeakResults.getDataSize(deviations, readEndFrame, readId);
			byte[] buffer = new byte[length];

			int c = 0;
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
				double chiSquared = readDouble(buffer);
				float noise = readFloat(buffer);
				float[] params = readData(buffer, new float[7]);
				float[] paramsStdDev = (deviations) ? readData(buffer, new float[7]) : null;

				// Convert old binary format with the amplitude to signal
				if (readAmplitude)
					params[Gaussian2DFunction.SIGNAL] *= 2 * Math.PI * params[Gaussian2DFunction.X_SD] *
							params[Gaussian2DFunction.Y_SD];

				if (readId || readEndFrame)
					results.add(new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params,
							paramsStdDev, endPeak, id));
				else
					results.add(new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev));

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

	private MemoryPeakResults createResults()
	{
		MemoryPeakResults results = new MemoryPeakResults();
		results.setName(name);
		results.setSource(source);
		results.setBounds(bounds);
		results.setCalibration(calibration);
		results.setConfiguration(configuration);
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

		BufferedReader input = null;
		try
		{
			FileInputStream fis = new FileInputStream(filename);
			FileChannel channel = fis.getChannel();
			input = new BufferedReader(new UnicodeReader(fis, null));

			String line;
			int errors = 0;

			// Read different versions
			int version = (readAmplitude) ? 1 : 2;

			// Skip the header
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;

				if (line.charAt(0) != '#')
				{
					// This is the first record
					if (!addPeakResult(results, line, version))
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

				if (!addPeakResult(results, line, version))
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

	private boolean addPeakResult(MemoryPeakResults results, String line, int version)
	{
		PeakResult result;
		switch (version)
		{
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

	private static Pattern p = Pattern.compile("\t");

	private PeakResult createPeakResultV1(String line)
	{
		try
		{
			float[] params = new float[7];

			if (isUseScanner())
			{
				// Code using a Scanner
				Scanner scanner = new Scanner(line);
				scanner.useDelimiter("\t");
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
				double chiSquared = scanner.nextDouble();
				float noise = scanner.nextFloat();
				float signal = scanner.nextFloat(); // Ignored but must be read
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
				}
				scanner.close();
				params[Gaussian2DFunction.SIGNAL] = signal;
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null);
			}
			else
			{
				// Code using split and parse
				String[] fields = p.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double chiSquared = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				float signal = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
				}
				params[Gaussian2DFunction.SIGNAL] = signal;
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null);
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
				scanner.useDelimiter("\t");
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
				double chiSquared = scanner.nextDouble();
				float noise = scanner.nextFloat();
				float signal = scanner.nextFloat(); // Ignored but must be read
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
					paramsStdDev[i] = scanner.nextFloat();
				}
				params[Gaussian2DFunction.SIGNAL] = signal;
				scanner.close();
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev);
			}
			else
			{
				// Tried using:
				// -String.split()
				// -Pattern.split()
				// -StringTokenizer
				// -Custom tokenizer routine the iterates the line
				// All are slower than the Scanner

				// Code using split and parse
				String[] fields = p.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double chiSquared = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				float signal = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
					paramsStdDev[i] = Float.parseFloat(fields[j++]);
				}
				params[Gaussian2DFunction.SIGNAL] = signal;
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev);
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
				scanner.useDelimiter("\t");
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
				double chiSquared = scanner.nextDouble();
				float noise = scanner.nextFloat();
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
				}
				scanner.close();
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null);
			}
			else
			{
				// Code using split and parse
				String[] fields = p.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double chiSquared = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
				}
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null,
							endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, null);
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
				scanner.useDelimiter("\t");
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
				double chiSquared = scanner.nextDouble();
				float noise = scanner.nextFloat();
				for (int i = 0; i < params.length; i++)
				{
					params[i] = scanner.nextFloat();
					paramsStdDev[i] = scanner.nextFloat();
				}
				scanner.close();
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev);
			}
			else
			{
				// Tried using:
				// -String.split()
				// -Pattern.split()
				// -StringTokenizer
				// -Custom tokenizer routine the iterates the line
				// All are slower than the Scanner

				// Code using split and parse
				String[] fields = p.split(line);
				int j = 0;
				int id = (readId) ? Integer.parseInt(fields[j++]) : 0;
				int peak = Integer.parseInt(fields[j++]);
				int endPeak = (readEndFrame) ? Integer.parseInt(fields[j++]) : 0;
				int origX = Integer.parseInt(fields[j++]);
				int origY = Integer.parseInt(fields[j++]);
				float origValue = Float.parseFloat(fields[j++]);
				double chiSquared = Double.parseDouble(fields[j++]);
				float noise = Float.parseFloat(fields[j++]);
				for (int i = 0; i < params.length; i++)
				{
					params[i] = Float.parseFloat(fields[j++]);
					paramsStdDev[i] = Float.parseFloat(fields[j++]);
				}
				if (readId || readEndFrame)
					return new ExtendedPeakResult(peak, origX, origY, origValue, chiSquared, noise, params,
							paramsStdDev, endPeak, id);
				else
					return new PeakResult(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev);
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

			// Old table results had the Signal and Amplitude.
			// New table results have only the Signal.
			int version = 2;
			if (header.contains("Amplitude"))
				version = 1;

			int c = 0;
			while ((line = input.readLine()) != null)
			{
				if (line.length() == 0)
					continue;

				if (!addTableResult(results, line, version))
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

	private boolean addTableResult(MemoryPeakResults results, String line, int version)
	{
		final PeakResult result;
		switch (version)
		{
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
				scanner.useDelimiter("\t");
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
			scanner.useDelimiter("\t");
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
			scanner.useDelimiter("\t");
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
		// *Note that the RapidSTORM Amplitude is the signal. To get the Amplitude we must divide by the 2*pi*sx*sy
		try
		{
			Scanner scanner = new Scanner(line);
			scanner.useDelimiter(" ");
			scanner.useLocale(Locale.US);
			float x = scanner.nextFloat();
			float y = scanner.nextFloat();
			final int peak = scanner.nextInt();
			final float signal = scanner.nextFloat();
			final float sx2 = scanner.nextFloat();
			final float sy2 = scanner.nextFloat();
			@SuppressWarnings("unused")
			final float kernelImprovement = scanner.nextFloat();
			final double chiSquared = scanner.nextDouble();
			scanner.close();

			// Convert from pm^2 to nm
			float sx = (float) (Math.sqrt(sx2) * 1000);
			float sy = (float) (Math.sqrt(sy2) * 1000);

			// If calibration was found convert to pixels
			if (calibration != null)
			{
				x /= calibration.nmPerPixel;
				y /= calibration.nmPerPixel;
				sx /= calibration.nmPerPixel;
				sy /= calibration.nmPerPixel;
			}

			float[] params = new float[7];
			params[Gaussian2DFunction.SIGNAL] = signal;
			params[Gaussian2DFunction.X_POSITION] = x;
			params[Gaussian2DFunction.Y_POSITION] = y;
			params[Gaussian2DFunction.X_SD] = sx;
			params[Gaussian2DFunction.Y_SD] = sy;

			// Store the signal as the original value
			return new PeakResult(peak, (int) x, (int) y, signal, chiSquared, 0.0f, params, null);
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
		// Area = Height * 2 * pi * (Width / (pixel_pitch*2) )^2
		// => Pixel_pitch = 0.5 * Width / sqrt(Area / (Height * 2 * pi))
		// Try and create a calibration
		Statistics pixelPitch = new Statistics();
		for (PeakResult p : results.getResults())
		{
			if (p.peak == p.getEndFrame())
			{
				float width = p.params[Gaussian2DFunction.X_SD];
				float height = p.params[Gaussian2DFunction.SIGNAL];
				float area = p.origValue;
				pixelPitch.add(0.5 * width / Math.sqrt(area / (height * 2 * Math.PI)));
				if (pixelPitch.getN() > 100)
					break;
			}
		}

		if (pixelPitch.getN() > 0)
		{
			final double nmPerPixel = pixelPitch.getMean();
			final double widthConversion = 1.0 / (2 * nmPerPixel);

			// Create a calibration
			calibration = new Calibration(nmPerPixel, 1, 0);
			results.setCalibration(calibration);

			// Convert data
			for (PeakResult p : results.getResults())
			{
				p.params[Gaussian2DFunction.X_POSITION] /= nmPerPixel;
				p.params[Gaussian2DFunction.Y_POSITION] /= nmPerPixel;
				// Since the width is 2*pixel pitch
				p.params[Gaussian2DFunction.X_SD] *= widthConversion;
				p.params[Gaussian2DFunction.Y_SD] *= widthConversion;
			}
		}

		// We initially stored the height of the peak in the signal field. 
		// Swap to the intensity stored in the origValue field.
		for (PeakResult p : results.getResults())
		{
			final float origValue = p.params[Gaussian2DFunction.SIGNAL];
			p.params[Gaussian2DFunction.SIGNAL] = p.origValue;
			p.origValue = origValue;
		}

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

	@SuppressWarnings("unused")
	// So that the fields can be named  
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
			scanner.useDelimiter("\t");
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

			float[] params = new float[7];
			params[Gaussian2DFunction.BACKGROUND] = bg;
			//params[Gaussian2DFunction.ANGLE] = ax;
			params[Gaussian2DFunction.SIGNAL] = height;
			params[Gaussian2DFunction.X_POSITION] = xc;
			params[Gaussian2DFunction.Y_POSITION] = yc;
			params[Gaussian2DFunction.X_SD] = width;
			params[Gaussian2DFunction.Y_SD] = width;

			// Store the signal as the original value
			return new ExtendedPeakResult(frame, (int) xc, (int) yc, area, 0.0, 0.0f, params, null, frame + length - 1,
					0);
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
}
