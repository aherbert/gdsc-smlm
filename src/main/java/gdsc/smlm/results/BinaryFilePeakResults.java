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

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;

/**
 * Saves the fit results to a binary file format
 */
public class BinaryFilePeakResults extends SMLMFilePeakResults
{
	public static final String END_HEADER = "END_HEADER";

	public BinaryFilePeakResults(String filename)
	{
		super(filename);
	}

	public BinaryFilePeakResults(String filename, boolean showDeviations)
	{
		super(filename, showDeviations);
	}

	public BinaryFilePeakResults(String resultsFilename, boolean showDeviations, boolean showEndFrame)
	{
		super(resultsFilename, showDeviations, showEndFrame);
	}

	public BinaryFilePeakResults(String resultsFilename, boolean showDeviations, boolean showEndFrame, boolean showId)
	{
		super(resultsFilename, showDeviations, showEndFrame, showId);
	}

	@Override
	protected void openOutput()
	{
		// Nothing to do as we write direct to the file output stream
	}
	
	@Override
	protected void write(String data)
	{
		try
		{
			fos.write(data.getBytes());
		}
		catch (IOException e)
		{
			closeOutput();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#getHeaderEnd()
	 */
	protected String getHeaderEnd()
	{
		return END_HEADER;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#getHeaderComments()
	 */
	protected String[] getHeaderComments()
	{
		String[] comments = new String[2];
		comments[0] = "Records start after the final comment line";
		String format = "iiifdffffffff";
		if (isShowDeviations())
			format += "fffffff";
		if (isShowEndFrame())
			format = "i" + format;
		if (isShowId())
			format = "i" + format;
		comments[1] = "Binary Format (raw Java bytes) = " + format;
		return comments;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#getFieldNames()
	 */
	protected String[] getFieldNames()
	{
		ArrayList<String> names = new ArrayList<String>(20);
		if (isShowId())
			names.add("Id");
		names.add(peakIdColumnName);
		if (isShowEndFrame())
			names.add("End " + peakIdColumnName);
		names.add("origX");
		names.add("origY");
		names.add("origValue");
		names.add("Error");
		names.add("Noise");
		names.add("Signal");
		String[] fieldNames = new String[] { "Background", "Amplitude", "Angle", "X", "Y", "X SD", "Y SD" };
		for (String field : fieldNames)
		{
			names.add(field);
		}
		if (isShowDeviations())
		{
			for (String field : fieldNames)
			{
				names.add(field + " StdDev");
			}
		}
		return names.toArray(new String[names.size()]);
	}

	protected void closeOutput()
	{
		super.closeOutput();

		if (fos == null)
			return;

		try
		{
			fos.close();
		}
		catch (Exception e)
		{
			// Ignore exception
		}
		finally
		{
			fos = null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		if (fos == null)
			return;

		// Buffer the output for the synchronized write method
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		DataOutputStream buffer = new DataOutputStream(bytes);

		try
		{
			addResult(buffer, 0, peak, peak, origX, origY, origValue, error, noise, params, paramsStdDev);
			buffer.flush();
		}
		catch (IOException e)
		{
			// Do nothing - This result will not be added to the file
			return;
		}

		writeResult(1, bytes.toByteArray());
	}

	private void addResult(DataOutputStream buffer, final int id, final int peak, final int endPeak, final int origX,
			final int origY, final float origValue, final double error, final float noise, final float[] params,
			float[] paramsStdDev) throws IOException
	{
		if (isShowId())
			buffer.writeInt(id);
		buffer.writeInt(peak);
		if (isShowEndFrame())
			buffer.writeInt(endPeak);
		buffer.writeInt(origX);
		buffer.writeInt(origY);
		buffer.writeFloat(origValue);
		buffer.writeDouble(error);
		buffer.writeFloat(noise);

		for (int i = 0; i < 7; i++)
			buffer.writeFloat(params[i]);
		if (isShowDeviations())
		{
			if (paramsStdDev == null)
				paramsStdDev = new float[7];
			for (int i = 0; i < 7; i++)
				buffer.writeFloat(paramsStdDev[i]);
		}
	}

	public void addAll(Collection<PeakResult> results)
	{
		if (fos == null)
			return;

		int count = 0;

		// Buffer the output for the synchronized write method
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		DataOutputStream buffer = new DataOutputStream(bytes);

		for (PeakResult result : results)
		{
			try
			{
				addResult(buffer, result.getId(), result.getFrame(), result.getEndFrame(), result.origX, result.origY,
						result.origValue, result.error, result.noise, result.params, result.paramsStdDev);
			}
			catch (IOException e)
			{
				// Do nothing - This result will not be added to the file
				return;
			}

			// Flush the output to allow for very large input lists
			if (++count >= 20)
			{
				try
				{
					buffer.flush();
				}
				catch (IOException e)
				{
					// Do nothing - This result will not be added to the file
					return;
				}
				writeResult(count, bytes.toByteArray());
				if (!isActive())
					return;
				bytes.reset();
				count = 0;
			}
		}

		try
		{
			buffer.flush();
		}
		catch (IOException e)
		{
			// Do nothing - This result will not be added to the file
			return;
		}
		writeResult(count, bytes.toByteArray());
	}

	private synchronized void writeResult(int count, byte[] bytes)
	{
		// In case another thread caused the output to close
		if (fos == null)
			return;
		size += count;
		try
		{
			fos.write(bytes);
		}
		catch (IOException ioe)
		{
			closeOutput();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#sort()
	 */
	protected void sort() throws IOException
	{
		try
		{
			ArrayList<Result> results = new ArrayList<Result>(size);

			DataInputStream input = new DataInputStream(new FileInputStream(filename));

			String header;
			try
			{
				header = readHeader(input);

				byte[] line = new byte[getDataSize(isShowDeviations(), isShowEndFrame(), isShowId())];
				while (input.read(line) == line.length)
				{
					results.add(new Result(line));
				}
			}
			finally
			{
				input.close();
			}

			Collections.sort(results);

			// Must write using the same method as the main code so use a FileOutputStream again
			FileOutputStream output = new FileOutputStream(filename);
			try
			{
				output.write(header.getBytes());
				for (Result result : results)
				{
					output.write(result.line);
				}
			}
			finally
			{
				output.close();
			}
		}
		catch (IOException e)
		{
			throw e;
		}
		finally
		{
			fos = null;
		}
	}

	/**
	 * Read all lines from the input stream that begin with '#' and collates them into a header.
	 * Stops reading if a line contains {@value #END_HEADER}.
	 * <p>
	 * Lines are defined by the '\n' character. The input stream will read the first non-header line unless the header
	 * is terminated by the {@value #END_HEADER} tag.
	 * 
	 * @param input
	 * @return The header
	 * @throws IOException
	 */
	public static String readHeader(DataInputStream input) throws IOException
	{
		StringBuffer sb = new StringBuffer();
		String line;
		do
		{
			line = readLine(input);
			if (line.charAt(0) == '#')
				sb.append(line);
			else
				break;
		} while (!line.startsWith("#" + END_HEADER));
		return sb.toString();
	}

	private static String readLine(DataInputStream input) throws IOException
	{
		StringBuffer sb = new StringBuffer();
		byte b;
		do
		{
			b = input.readByte();
			sb.append((char) b);
		} while (b != '\n');
		return sb.toString();
	}

	public static int getDataSize(boolean deviations, boolean endFrame, boolean id)
	{
		final int BYTES_INT = 4;
		final int BYTES_FLOAT = 4;
		final int BYTES_DOUBLE = 8;

		// iiifdffffffff
		//	or
		// iiifdfffffffffffffff
		// + Extra i added for the end frame after the first integer
		// + Extra i added for the id in the first field
		int size = 3 * BYTES_INT + BYTES_FLOAT + BYTES_DOUBLE + 8 * BYTES_FLOAT;
		if (deviations)
			size += 7 * BYTES_FLOAT;
		if (endFrame)
			size += BYTES_INT;
		if (id)
			size += BYTES_INT;
		return size;
	}

	private class Result implements Comparable<Result>
	{
		byte[] line;
		int slice = 0;

		public Result(byte[] line)
		{
			this.line = Arrays.copyOf(line, line.length);
			extractSlice();
		}

		private void extractSlice()
		{
			int offset = (isShowId()) ? 4 : 0;
			int ch1 = line[offset + 0] & 0xff;
			int ch2 = line[offset + 1] & 0xff;
			int ch3 = line[offset + 2] & 0xff;
			int ch4 = line[offset + 3] & 0xff;
			slice = ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
		}

		public int compareTo(Result o)
		{
			// Sort by slice number
			// (Note: peak height is already done in the run(...) method)
			return slice - o.slice;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#isBinary()
	 */
	@Override
	public boolean isBinary()
	{
		return true;
	}
}
