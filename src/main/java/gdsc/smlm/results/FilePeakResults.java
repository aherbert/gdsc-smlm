package gdsc.smlm.results;

import gdsc.smlm.utils.XmlUtils;

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

import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;

/**
 * Saves the fit results to file
 */
public class FilePeakResults extends AbstractPeakResults
{
	// Only write to a single results file
	protected OutputStreamWriter out = null;

	protected String filename;
	private boolean showDeviations = true;
	private boolean showEndFrame = false;
	private boolean showId = false;
	protected boolean sortAfterEnd = false;
	protected String peakIdColumnName = "Peak";

	protected int size = 0;

	public FilePeakResults(String filename)
	{
		this.filename = filename;
	}

	public FilePeakResults(String filename, boolean showDeviations)
	{
		this.filename = filename;
		this.showDeviations = showDeviations;
	}

	public FilePeakResults(String filename, boolean showDeviations, boolean showEndFrame)
	{
		this.filename = filename;
		this.showDeviations = showDeviations;
		this.showEndFrame = showEndFrame;
	}

	public FilePeakResults(String filename, boolean showDeviations, boolean showEndFrame, boolean showId)
	{
		this.filename = filename;
		this.showDeviations = showDeviations;
		this.showEndFrame = showEndFrame;
		this.showId = showId;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#begin()
	 */
	public void begin()
	{
		out = null;
		size = 0;
		try
		{
			FileOutputStream fos = new FileOutputStream(filename);
			out = new OutputStreamWriter(fos, "UTF-8");
			out.write(createResultsHeader());
		}
		catch (Exception e)
		{
			// TODO - Add better handling of errors
			e.printStackTrace();
			closeOutput();
		}
	}

	protected String createResultsHeader()
	{
		StringBuilder sb = new StringBuilder();

		addComment(sb, getHeaderTitle());
		sb.append(String.format("#FileVersion %s\n", getVersion()));

		// Add the standard details
		if (name != null)
			sb.append(String.format("#Name %s\n", singleLine(name)));
		if (source != null)
			sb.append(String.format("#Source %s\n", singleLine(source.toXML())));
		if (bounds != null)
			sb.append(String.format("#Bounds x%d y%d w%d h%d\n", bounds.x, bounds.y, bounds.width, bounds.height));
		if (calibration != null)
			sb.append(String.format("#Calibration %s\n", singleLine(XmlUtils.toXML(calibration))));
		if (configuration != null && configuration.length() > 0)
			sb.append(String.format("#Configuration %s\n", singleLine(configuration)));

		// Add any extra comments
		String[] comments = getHeaderComments();
		if (comments != null)
		{
			for (String comment : comments)
				addComment(sb, comment);
		}

		// Output the field names
		String[] fields = getFieldNames();
		if (fields != null)
		{
			sb.append("#");
			for (int i = 0; i < fields.length; i++)
			{
				if (i != 0)
					sb.append("\t");
				sb.append(fields[i]);
			}
			sb.append('\n');
		}

		addComment(sb, getHeaderEnd());

		return sb.toString();
	}

	private void addComment(StringBuilder sb, String comment)
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
	protected String getVersion()
	{
		StringBuilder sb = new StringBuilder();
		sb.append(isBinary() ? "Binary" : "Text");
		sb.append(".");
		sb.append(isShowDeviations() ? "D1" : "D0");
		sb.append(".E");
		int extended = isShowEndFrame() ? 1 : 0;
		extended += isShowId() ? 2 : 0;
		sb.append(extended);
		// Version 1 had signal and amplitude in the results
		// Version 2 has only signal in the results
		sb.append(".V2");
		return sb.toString();
	}

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
	protected String[] getFieldNames()
	{
		ArrayList<String> names = new ArrayList<String>(20);
		if (showId)
			names.add("Id");
		names.add(peakIdColumnName);
		if (showEndFrame)
			names.add("End " + peakIdColumnName);
		names.add("origX");
		names.add("origY");
		names.add("origValue");
		names.add("Error");
		names.add("Noise");
		for (String field : new String[] { "Background", "Signal", "Angle", "X", "Y", "X SD", "Y SD" })
		{
			names.add(field);
			if (showDeviations)
				names.add("+/-");
		}
		if (this.calibration != null)
			names.add("Precision");
		return names.toArray(new String[names.size()]);
	}

	protected String singleLine(String text)
	{
		return text.replaceAll("\n *", "");
	}

	protected void closeOutput()
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
	 * @see gdsc.utils.fitting.results.PeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise, float[] params,
			float[] paramsStdDev)
	{
		if (out == null)
			return;

		StringBuilder sb = new StringBuilder();

		addStandardData(sb, 0, peak, peak, origX, origY, origValue, chiSquared, noise);

		// Add the parameters		
		if (showDeviations)
		{
			if (paramsStdDev != null)
				paramsStdDev = new float[7];
			addResult(sb, params[0], paramsStdDev[0], params[1], paramsStdDev[1], params[2], paramsStdDev[2],
					params[3], paramsStdDev[3], params[4], paramsStdDev[4], params[5], paramsStdDev[5], params[6],
					paramsStdDev[6]);
		}
		else
		{
			addResult(sb, params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
		}

		if (this.calibration != null)
		{
			double s = (params[Gaussian2DFunction.X_SD] + params[Gaussian2DFunction.Y_SD]) * 0.5 *
					calibration.nmPerPixel;
			float precision = (float) PeakResult.getPrecision(calibration.nmPerPixel, s,
					params[Gaussian2DFunction.SIGNAL] / calibration.gain, noise / calibration.gain, calibration.emCCD);
			addResult(sb, precision);
		}

		sb.append('\n');
		writeResult(1, sb.toString());
	}

	private void addStandardData(StringBuilder sb, final int id, final int peak, final int endPeak, final int origX,
			final int origY, final float origValue, final double chiSquared, final float noise)
	{
		if (showId)
		{
			sb.append(id);
			sb.append('\t');
		}
		sb.append(peak);
		sb.append('\t');
		if (showEndFrame)
		{
			sb.append(endPeak);
			sb.append('\t');
		}
		sb.append(origX);
		sb.append('\t');
		sb.append(origY);
		sb.append('\t');
		sb.append(origValue);
		sb.append('\t');
		sb.append(chiSquared);
		sb.append('\t');
		sb.append(noise);
	}

	private void addResult(StringBuilder sb, float... args)
	{
		for (float f : args)
			//sb.append(String.format("\t%g", f));
			sb.append("\t").append(f);
	}

	public void addAll(Collection<PeakResult> results)
	{
		if (out == null)
			return;

		int count = 0;

		StringBuilder sb = new StringBuilder();
		for (PeakResult result : results)
		{
			// Add the standard data
			addStandardData(sb, result.getId(), result.peak, result.getEndFrame(), result.origX, result.origY,
					result.origValue, result.error, result.noise);

			// Add the parameters		
			if (showDeviations)
			{
				final float[] paramsStdDev = (result.paramsStdDev != null) ? result.paramsStdDev : new float[7];
				addResult(sb, result.params[0], paramsStdDev[0], result.params[1], paramsStdDev[1], result.params[2],
						paramsStdDev[2], result.params[3], paramsStdDev[3], result.params[4], paramsStdDev[4],
						result.params[5], paramsStdDev[5], result.params[6], paramsStdDev[6]);
			}
			else
			{
				addResult(sb, result.params[0], result.params[1], result.params[2], result.params[3], result.params[4],
						result.params[5], result.params[6]);
			}

			if (this.calibration != null)
			{
				double s = result.getSD() * calibration.nmPerPixel;
				float precision = (float) PeakResult.getPrecision(calibration.nmPerPixel, s,
						result.params[Gaussian2DFunction.SIGNAL] / calibration.gain, result.noise / calibration.gain,
						calibration.emCCD);
				addResult(sb, precision);
			}
			sb.append('\n');

			// Flush the output to allow for very large input lists
			if (++count >= 20)
			{
				writeResult(count, sb.toString());
				if (!isActive())
					return;
				sb.setLength(0);
				count = 0;
			}
		}
		writeResult(count, sb.toString());
	}

	/**
	 * Output a cluster to the results file.
	 * <p>
	 * Note: This is not synchronised
	 * 
	 * @param cluster
	 */
	public void addCluster(Cluster cluster)
	{
		if (out == null)
			return;
		if (cluster.size() > 0)
		{
			float[] centroid = cluster.getCentroid();
			writeResult(
					0,
					String.format("#Cluster %f %f (+/-%f) n=%d\n", centroid[0], centroid[1],
							cluster.getStandardDeviation(), cluster.size()));
			addAll(cluster);
		}
	}

	private void addAll(Cluster cluster)
	{
		if (!showId || cluster.getId() == 0)
		{
			addAll(cluster.getPoints());
		}
		else
		{
			// Store the ID from the trace
			final int id = cluster.getId();
			ArrayList<PeakResult> results = cluster.getPoints();
			ArrayList<PeakResult> results2 = new ArrayList<PeakResult>(results.size());
			for (PeakResult result : results)
			{
				if (result.getId() == id)
					results2.add(result);
				else
				{
					results2.add(new ExtendedPeakResult(result.peak, result.origX, result.origY, result.origValue,
							result.error, result.noise, result.params, result.paramsStdDev, result.getEndFrame(), id));
				}
			}
			addAll(results2);
		}
	}

	/**
	 * Output a trace to the results file.
	 * <p>
	 * Note: This is not synchronised
	 * 
	 * @param trace
	 */
	public void addTrace(Trace trace)
	{
		if (out == null)
			return;
		if (trace.size() > 0)
		{
			float[] centroid = trace.getCentroid();
			writeResult(0, String.format("#Trace %f %f (+/-%f) n=%d, b=%d, on=%f, off=%f, signal= %f\n", centroid[0],
					centroid[1], trace.getStandardDeviation(), trace.size(), trace.getNBlinks(), trace.getOnTime(),
					trace.getOffTime(), trace.getSignal()));
			addAll(trace);
		}
	}

	/**
	 * Output a comment to the results file.
	 * <p>
	 * Note: This is not synchronised
	 * 
	 * @param text
	 */
	public void addComment(String text)
	{
		if (out == null)
			return;
		// Ensure comments are preceded by the comment character
		if (!text.startsWith("#"))
			text = "#" + text;
		if (text.contains("\n"))
			text.replace("\n", "\n#");
		if (!text.endsWith("\n"))
			text += "\n";
		writeResult(0, text);
	}

	protected synchronized void writeResult(int count, String result)
	{
		// In case another thread caused the output to close
		if (out == null)
			return;
		size += count;
		try
		{
			out.write(result);
		}
		catch (IOException ioe)
		{
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
		if (out == null)
			return;

		// Close the file.
		try
		{
			closeOutput();

			if (!isSortAfterEnd())
				return;

			ArrayList<Result> results = new ArrayList<Result>(size);

			StringBuffer header = new StringBuffer();
			BufferedReader input = new BufferedReader(new FileReader(filename));
			try
			{
				String line;
				// Skip the header
				while ((line = input.readLine()) != null)
				{
					if (line.charAt(0) != '#')
					{
						// This is the first record
						results.add(new Result(line));
						break;
					}
					else
						header.append(line).append('\n');
				}

				while ((line = input.readLine()) != null)
				{
					results.add(new Result(line));
				}
			}
			finally
			{
				input.close();
			}

			Collections.sort(results);

			BufferedWriter output = new BufferedWriter(new FileWriter(filename));
			try
			{
				output.write(header.toString());
				for (Result result : results)
				{
					output.write(result.line);
					output.write("\n");
				}
			}
			finally
			{
				output.close();
			}
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			out = null;
		}
	}

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

	private class Result implements Comparable<Result>
	{
		String line;
		int slice = 0;

		public Result(String line)
		{
			this.line = line;
			extractSlice();
		}

		private void extractSlice()
		{
			Scanner scanner = new Scanner(line);
			scanner.useDelimiter("\t");

			try
			{
				slice = scanner.nextInt();
				if (showId)
					// The peak is the second column
					slice = scanner.nextInt();
				scanner.close();
			}
			catch (InputMismatchException e)
			{
			}
			catch (NoSuchElementException e)
			{
			}
		}

		public int compareTo(Result o)
		{
			// Sort by slice number
			// (Note: peak height is already done in the run(...) method)
			return slice - o.slice;
		}
	}

	/**
	 * @return the name of the peak column
	 */
	public String getPeakIdColumnName()
	{
		return peakIdColumnName;
	}

	/**
	 * @param peakIdColumnName
	 *            the name of the peak column
	 */
	public void setPeakIdColumnName(String peakIdColumnName)
	{
		this.peakIdColumnName = peakIdColumnName;
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

	/**
	 * @return True if the records contain the parameter deviations
	 */
	public boolean isShowDeviations()
	{
		return showDeviations;
	}

	/**
	 * @return True if the records contain the result end frame
	 */
	public boolean isShowEndFrame()
	{
		return showEndFrame;
	}

	/**
	 * @return true if the records are stored as binary data
	 */
	public boolean isBinary()
	{
		return false;
	}

	/**
	 * @return True if the records contain a result Id
	 */
	public boolean isShowId()
	{
		return showId;
	}
}
