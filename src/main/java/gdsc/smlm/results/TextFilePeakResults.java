package gdsc.smlm.results;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.results.Calibration.CameraType;

/**
 * Saves the fit results to file
 */
public class TextFilePeakResults extends SMLMFilePeakResults
{
	private OutputStreamWriter out;
	
	private boolean canComputePrecision, emCCD;

	public TextFilePeakResults(String filename)
	{
		super(filename);
	}

	public TextFilePeakResults(String filename, boolean showDeviations)
	{
		super(filename, showDeviations);
	}

	public TextFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame)
	{
		super(filename, showDeviations, showEndFrame);
	}

	public TextFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame, boolean showId)
	{
		super(filename, showDeviations, showEndFrame, showId);
	}


	@Override
	protected void openOutput()
	{
		try
		{
			out = new OutputStreamWriter(fos, "UTF-8");
		}
		catch (UnsupportedEncodingException e)
		{
			throw new RuntimeException(e);
		}
	}

	@Override
	protected void write(String data)
	{
		try
		{
			out.write(data);
		}
		catch (IOException e)
		{
			closeOutput();
		}
	}
	
	@Override
	protected void closeOutput()
	{
		if (fos == null)
			return;

		try
		{
			// Make sure we close the writer since it may be buffered
			out.close();
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
	 * @see gdsc.smlm.results.FilePeakResults#createStandardConverters()
	 */
	@Override
	protected void createStandardConverters()
	{
		super.createStandardConverters();
		
		// TODO - make this optional
		canComputePrecision = false;

		if (calibration != null)
		{
			if (isCCD() && toNMConverter != null && toPhotonConverter != null)
			{
				emCCD = calibration.getCameraType() == CameraType.EM_CCD;
				canComputePrecision = true;
			}

			// TODO - Add ability to write output in selected units

			// Clone the calibration as we may change it
			//this.calibration = calibration.clone();

			// Check the units. If not correct then set-up a conversion and change the units
			// in the calibration
			// This will require a conversion multiplier for distance
			// It will require a conversion for signal 
			// It will require a conversion plus bias for background:
			// (b-bias1)*conversion + bias2

		}
	}

	private boolean isCCD()
	{
		if (calibration.hasCameraType())
		{
			switch (calibration.getCameraType())
			{
				case CCD:
				case EM_CCD:
					return true;
				default:
					break;
			}
		}
		return false;
	}

	/**
	 * @return The names of the fields in each record. Will be the last comment of the header
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
		String[] fields = new String[] { "Background", "Signal", "Angle", "X", "Y", "X SD", "Y SD" };
		// Add units
		if (calibration != null)
		{
			if (calibration.hasIntensityUnit())
			{
				String unit = String.format(" (%s)", calibration.getIntensityUnit());
				fields[0] += unit;
				fields[1] += unit;
			}
			if (calibration.hasDistanceUnit())
			{
				String unit = String.format(" (%s)", calibration.getDistanceUnit());
				fields[3] += unit;
				fields[4] += unit;
				fields[5] += unit;
				fields[6] += unit;
			}
		}
		for (String field : fields)
		{
			names.add(field);
			if (isShowDeviations())
				names.add("+/-");
		}
		if (canComputePrecision)
			names.add("Precision (nm)");
		return names.toArray(new String[names.size()]);
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

		StringBuilder sb = new StringBuilder();

		addStandardData(sb, 0, peak, peak, origX, origY, origValue, error, noise);

		// Add the parameters		
		if (isShowDeviations())
		{
			if (paramsStdDev != null)
				paramsStdDev = new float[7];
			addResult(sb, params[0], paramsStdDev[0], params[1], paramsStdDev[1], params[2], paramsStdDev[2], params[3],
					paramsStdDev[3], params[4], paramsStdDev[4], params[5], paramsStdDev[5], params[6],
					paramsStdDev[6]);
		}
		else
		{
			addResult(sb, params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
		}

		if (canComputePrecision)
		{
			double s = toNMConverter
					.convert(PeakResult.getSD(params[Gaussian2DFunction.X_SD], params[Gaussian2DFunction.Y_SD]));
			float precision = (float) PeakResult.getPrecision(nmPerPixel, s,
					toPhotonConverter.convert(params[Gaussian2DFunction.SIGNAL]), toPhotonConverter.convert(noise),
					emCCD);
			addResult(sb, precision);
		}

		sb.append('\n');
		writeResult(1, sb.toString());
	}

	private void addStandardData(StringBuilder sb, final int id, final int peak, final int endPeak, final int origX,
			final int origY, final float origValue, final double chiSquared, final float noise)
	{
		if (isShowId())
		{
			sb.append(id);
			sb.append('\t');
		}
		sb.append(peak);
		sb.append('\t');
		if (isShowEndFrame())
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
			addStandardData(sb, result.getId(), result.getFrame(), result.getEndFrame(), result.origX, result.origY,
					result.origValue, result.error, result.noise);

			// Add the parameters		
			if (isShowDeviations())
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

			if (canComputePrecision)
			{
				double s = toNMConverter.convert(result.getSD());
				float precision = (float) PeakResult.getPrecision(nmPerPixel, s,
						toPhotonConverter.convert(result.params[Gaussian2DFunction.SIGNAL]),
						toPhotonConverter.convert(result.noise), emCCD);
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
			writeResult(0, String.format("#Cluster %f %f (+/-%f) n=%d\n", centroid[0], centroid[1],
					cluster.getStandardDeviation(), cluster.size()));
			addAll(cluster);
		}
	}

	protected void addAll(Cluster cluster)
	{
		if (!isShowId() || cluster.getId() == 0)
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
					results2.add(new ExtendedPeakResult(result.getFrame(), result.origX, result.origY, result.origValue,
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
			writeResult(0,
					String.format("#Trace %f %f (+/-%f) n=%d, b=%d, on=%f, off=%f, signal= %f\n", centroid[0],
							centroid[1], trace.getStandardDeviation(), trace.size(), trace.getNBlinks(),
							trace.getOnTime(), trace.getOffTime(), trace.getSignal()));
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
	 * @see gdsc.smlm.results.FilePeakResults#sort()
	 */
	protected void sort() throws IOException
	{
		try
		{
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
			throw e;
		}
		finally
		{
			out = null;
		}
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
				if (isShowId())
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
}
