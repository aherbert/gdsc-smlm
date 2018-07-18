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
package uk.ac.sussex.gdsc.smlm.results;

import java.io.FileOutputStream;
import java.io.IOException;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.MessageOrBuilder;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Printer;

import uk.ac.sussex.gdsc.core.utils.NotImplementedException;
import uk.ac.sussex.gdsc.core.utils.TextUtils;

/**
 * Saves the fit results to file
 */
public abstract class FilePeakResults extends AbstractPeakResults implements ThreadSafePeakResults
{
	/** Only write to a single results file */
	protected FileOutputStream fos = null;

	/** The filename. */
	protected String filename;
	private boolean sortAfterEnd = false;

	/** The size of the results. */
	protected int size = 0;

	/**
	 * Instantiates a new file peak results.
	 *
	 * @param filename
	 *            the filename
	 */
	public FilePeakResults(String filename)
	{
		this.filename = filename;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.utils.fitting.PeakResults#begin()
	 */
	@Override
	public void begin()
	{
		fos = null;
		size = 0;
		try
		{
			fos = new FileOutputStream(filename);
			openOutput();
			write(createResultsHeader());
		}
		catch (final Exception e)
		{
			// TODO - Add better handling of errors
			e.printStackTrace();
			closeOutput();
		}
	}

	/**
	 * Open the required output from the open file output stream
	 */
	protected abstract void openOutput();

	/**
	 * Write the data to the output.
	 *
	 * @param data
	 *            the data
	 */
	protected abstract void write(String data);

	/**
	 * Write the result and increment the size by the count.
	 * <p>
	 * This method is synchronised to ensure that the change to the size or the output file are thread safe.
	 *
	 * @param count
	 *            the count
	 * @param result
	 *            the result
	 */
	protected synchronized void writeResult(int count, String result)
	{
		// In case another thread caused the output to close
		if (fos == null)
			return;
		size += count;
		write(result);
	}

	/**
	 * Creates the results header.
	 *
	 * @return the header
	 */
	protected String createResultsHeader()
	{
		final StringBuilder sb = new StringBuilder();

		addComment(sb, getHeaderTitle());
		sb.append(String.format("#FileVersion %s\n", getVersion()));

		Printer printer = null;

		// Add the standard details
		if (!TextUtils.isNullOrEmpty(getName()))
			sb.append(String.format("#Name %s\n", singleLine(getName())));
		if (getSource() != null)
			sb.append(String.format("#Source %s\n", singleLine(getSource().toXML())));
		if (getBounds() != null)
			sb.append(String.format("#Bounds x%d y%d w%d h%d\n", getBounds().x, getBounds().y, getBounds().width,
					getBounds().height));
		if (getCalibration() != null)
			printer = addMessage(sb, printer, "Calibration", getCalibration());
		if (!TextUtils.isNullOrEmpty(getConfiguration()))
			sb.append(String.format("#Configuration %s\n", singleLine(getConfiguration())));
		if (getPSF() != null)
			printer = addMessage(sb, printer, "PSF", getPSF());

		// Add any extra comments
		final String[] comments = getHeaderComments();
		if (comments != null)
			for (final String comment : comments)
				addComment(sb, comment);

		// Output the field names
		final String[] fields = getFieldNames();
		if (fields != null)
		{
			sb.append('#');
			for (int i = 0; i < fields.length; i++)
			{
				if (i != 0)
					sb.append('\t');
				sb.append(fields[i]);
			}
			sb.append('\n');
		}

		addComment(sb, getHeaderEnd());

		return sb.toString();
	}

	private static Printer addMessage(StringBuilder sb, Printer printer, String name, MessageOrBuilder msg)
	{
		try
		{
			if (printer == null)
				printer = JsonFormat.printer().omittingInsignificantWhitespace()
				//.includingDefaultValueFields()
				;
			sb.append(String.format("#%s %s\n", name, printer.print(msg)));
		}
		catch (final InvalidProtocolBufferException e)
		{
			// This shouldn't happen so throw it
			throw new NotImplementedException("Unable to serialise the " + name + " settings", e);
		}
		return printer;
	}

	private static void addComment(StringBuilder sb, String comment)
	{
		if (comment != null)
			sb.append("#").append(comment).append('\n');
	}

	/**
	 * @return The first line added to the header
	 */
	protected String getHeaderTitle()
	{
		return "Localisation Results File";
	}

	/**
	 * @return The last line added to the header (e.g. a header end tag)
	 */
	protected String getHeaderEnd()
	{
		return null;
	}

	/**
	 * @return A line containing the file format version
	 */
	protected abstract String getVersion();

	/**
	 * @return Any comment lines to add to the header after the standard output of source, name, bounds, etc.
	 */
	protected String[] getHeaderComments()
	{
		return null;
	}

	/**
	 * @return The names of the fields in each record. Will be the last comment of the header
	 */
	protected abstract String[] getFieldNames();

	/**
	 * Convert the text to a single line.
	 *
	 * @param text
	 *            the text
	 * @return the new text
	 */
	protected static String singleLine(String text)
	{
		return text.replaceAll("\n *", "");
	}

	/**
	 * Close the output.
	 */
	protected void closeOutput()
	{
		if (fos == null)
			return;

		try
		{
			fos.close();
		}
		catch (final Exception e)
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
	 * @see gdsc.utils.fitting.PeakResults#size()
	 */
	@Override
	public int size()
	{
		return size;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.utils.fitting.PeakResults#end()
	 */
	@Override
	public void end()
	{
		if (fos == null)
			return;

		// Close the file.
		try
		{
			closeOutput();

			if (!isSortAfterEnd())
				return;

			sort();
		}
		catch (final IOException e)
		{
			// ignore
		}
		finally
		{
			fos = null;
		}
	}

	/**
	 * Sort the data file records. This is called once the file has been closed for input.
	 *
	 * @throws IOException
	 *             if an IO error occurs
	 */
	protected abstract void sort() throws IOException;

	/**
	 * @param sortAfterEnd
	 *            True if the results should be sorted after the {@link #end()} method
	 */
	public void setSortAfterEnd(boolean sortAfterEnd)
	{
		this.sortAfterEnd = sortAfterEnd;
	}

	/**
	 * @return True if the results should be sorted after the {@link #end()} method
	 */
	public boolean isSortAfterEnd()
	{
		return sortAfterEnd;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	@Override
	public boolean isActive()
	{
		return fos != null;
	}

	/**
	 * @return true if the records are stored as binary data
	 */
	public boolean isBinary()
	{
		return false;
	}
}
