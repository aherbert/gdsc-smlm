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

import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import gdsc.core.utils.BitFlags;

/**
 * Saves the fit results to a binary file format
 */
public class BinaryFilePeakResults extends SMLMFilePeakResults
{
	public static final String END_HEADER = "END_HEADER";

	private String[] fieldNames;
	private int nFields;

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

	public BinaryFilePeakResults(String resultsFilename, boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean showPrecision)
	{
		super(resultsFilename, showDeviations, showEndFrame, showId, showPrecision);
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
	@Override
	protected String getHeaderEnd()
	{
		return END_HEADER;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FilePeakResults#getHeaderComments()
	 */
	@Override
	protected String[] getHeaderComments()
	{
		fieldNames = new PeakResultConversionHelper(null, getPSF()).getNames();
		nFields = fieldNames.length;

		String[] comments = new String[2];
		comments[0] = "Records start after the final comment line";
		StringBuilder sb = new StringBuilder();
		if (isShowEndFrame())
			sb.append('i');
		if (isShowId())
			sb.append('i');
		sb.append("iiifdff");
		if (isShowDeviations())
		{
			for (int i = 0; i < nFields; i++)
				sb.append("ff");
		}
		else
		{
			for (int i = 0; i < nFields; i++)
				sb.append('f');
		}
		if (isShowPrecision())
			sb.append("f");
		comments[1] = "Binary Format (raw Java bytes) = " + sb.toString();
		return comments;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FilePeakResults#getFieldNames()
	 */
	@Override
	protected String[] getFieldNames()
	{
		ArrayList<String> names = new ArrayList<>(20);
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
		for (int i = 0; i < nFields; i++)
		{
			names.add(fieldNames[i]);
		}
		if (isShowDeviations())
		{
			for (int i = 0; i < nFields; i++)
			{
				names.add(fieldNames[i] + " StdDev");
			}
		}
		if (isShowPrecision())
			names.add("Precision (nm)");
		return names.toArray(new String[names.size()]);
	}

	@Override
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
	 * @see gdsc.smlm.results.PeakResults#add(int, int, int, float, double, float, float, float[], float[])
	 */
	@Override
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float meanIntensity,
			float[] params, float[] paramsStdDev)
	{
		if (fos == null)
			return;

		// Buffer the output for the synchronized write method
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		DataOutputStream buffer = new DataOutputStream(bytes);

		try
		{
			addResult(buffer, 0, peak, peak, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev,
					0.0);
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
			final int origY, final float origValue, final double error, final float noise, float meanIntensity,
			final float[] params, float[] paramsStdDev, double precision) throws IOException
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
		buffer.writeFloat(meanIntensity);

		checkSize(nFields, params);
		for (int i = 0; i < nFields; i++)
			buffer.writeFloat(params[i]);
		if (isShowDeviations())
		{
			if (paramsStdDev == null)
			{
				for (int i = 0; i < nFields; i++)
					buffer.writeInt(0); // An empty int is the same size as an empty float
			}
			else
			{
				checkSize(nFields, paramsStdDev);
				for (int i = 0; i < nFields; i++)
					buffer.writeFloat(paramsStdDev[i]);
			}
		}
		if (isShowPrecision())
			buffer.writeFloat((float) precision);
	}

	@Override
	public void add(PeakResult result)
	{
		if (fos == null)
			return;

		// Buffer the output for the synchronized write method
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		DataOutputStream buffer = new DataOutputStream(bytes);

		try
		{
			addResult(buffer, result.getId(), result.getFrame(), result.getEndFrame(), result.getOrigX(),
					result.getOrigY(), result.getOrigValue(), result.getError(), result.getNoise(),
					result.getMeanIntensity(), result.getParameters(), result.getParameterDeviations(),
					result.getPrecision());
			buffer.flush();
		}
		catch (IOException e)
		{
			// Do nothing - This result will not be added to the file
			return;
		}

		writeResult(1, bytes.toByteArray());
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.PeakResults#addAll(gdsc.smlm.results.PeakResult[])
	 */
	@Override
	public void addAll(PeakResult[] results)
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
				addResult(buffer, result.getId(), result.getFrame(), result.getEndFrame(), result.getOrigX(),
						result.getOrigY(), result.getOrigValue(), result.getError(), result.getNoise(),
						result.getMeanIntensity(), result.getParameters(), result.getParameterDeviations(),
						result.getPrecision());
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
	@Override
	protected void sort() throws IOException
	{
		try
		{
			ArrayList<Result> results = new ArrayList<>(size);

			DataInputStream input = new DataInputStream(new FileInputStream(filename));

			String header;
			try
			{
				header = readHeader(input);

				int flags = 0;
				if (isShowEndFrame())
					flags += FLAG_END_FRAME;
				if (isShowId())
					flags += FLAG_ID;
				if (isShowPrecision())
					flags += FLAG_PRECISION;
				byte[] line = new byte[getDataSize(isShowDeviations(), flags, nFields)];
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
	 *            the input
	 * @return The header
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public static String readHeader(DataInputStream input) throws IOException
	{
		StringBuilder sb = new StringBuilder();
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
		StringBuilder sb = new StringBuilder();
		byte b;
		do
		{
			b = input.readByte();
			sb.append((char) b);
		} while (b != '\n');
		return sb.toString();
	}

	static int getDataSize(boolean deviations, int flags, int nFields)
	{
		final int BYTES_INT = 4;
		final int BYTES_FLOAT = 4;
		final int BYTES_DOUBLE = 8;

		// iiifdff + n*f
		//	or
		// iiifdff + 2*n*f
		// + Extra i added for the end frame after the first integer
		// + Extra i added for the id in the first field
		int size = 3 * BYTES_INT + BYTES_FLOAT + BYTES_DOUBLE + BYTES_FLOAT + BYTES_FLOAT + nFields * BYTES_FLOAT;
		if (deviations)
			size += nFields * BYTES_FLOAT;
		if (BitFlags.areSet(flags, FLAG_END_FRAME))
			size += BYTES_INT;
		if (BitFlags.areSet(flags, FLAG_ID))
			size += BYTES_INT;
		if (BitFlags.areSet(flags, FLAG_PRECISION))
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

		@Override
		public int compareTo(Result o)
		{
			// Sort by slice number
			// (Note: peak height is already done in the run(...) method)
			if (slice < o.slice)
				return -1;
			if (slice > o.slice)
				return 1;
			return 0;
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
