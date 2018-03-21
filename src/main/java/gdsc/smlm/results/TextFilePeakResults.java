package gdsc.smlm.results;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;

import gdsc.core.data.utils.ConversionException;

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

import gdsc.core.data.utils.Converter;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Saves the fit results to file
 */
public class TextFilePeakResults extends SMLMFilePeakResults
{
	private Gaussian2DPeakResultCalculator calculator;

	private PeakResultConversionHelper helper;
	private Converter[] converters;

	private DistanceUnit distanceUnit = null;
	private IntensityUnit intensityUnit = null;
	private AngleUnit angleUnit = null;
	private boolean computePrecision = false;
	private boolean canComputePrecision = false;

	private OutputStreamWriter out;

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

	public TextFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean showPrecision)
	{
		super(filename, showDeviations, showEndFrame, showId, showPrecision);
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
	 * @see gdsc.smlm.results.FilePeakResults#begin()
	 */
	@Override
	public void begin()
	{
		calculator = null;
		canComputePrecision = false;

		if (isShowPrecision() && hasCalibration())
		{
			// Determine if we can compute the precision using the current settings
			if (computePrecision)
			{
				try
				{
					calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibrationReader(),
							Gaussian2DPeakResultHelper.LSE_PRECISION);
					canComputePrecision = true;
				}
				catch (ConfigurationException e)
				{
				}
				catch (ConversionException e)
				{
				}
			}
		}

		// We must correctly convert all the PSF parameter types
		helper = new PeakResultConversionHelper(getCalibration(), getPSF());
		helper.setIntensityUnit(intensityUnit);
		helper.setDistanceUnit(distanceUnit);
		helper.setAngleUnit(angleUnit);
		converters = helper.getConverters();
		// Update the calibration if converters were created
		if (helper.isCalibrationChanged())
			setCalibration(helper.getCalibration());

		super.begin();
	}

	/**
	 * @return The names of the fields in each record. Will be the last comment of the header
	 */
	protected String[] getFieldNames()
	{
		String[] unitNames = helper.getUnitNames();

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
		String noiseField = "Noise";
		if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY]))
			noiseField += " (" + (unitNames[PeakResult.INTENSITY] + ")");
		names.add(noiseField);

		String[] fields = helper.getNames();

		for (int i = 0; i < fields.length; i++)
		{
			String f = fields[i];
			// Add units
			if (!TextUtils.isNullOrEmpty(unitNames[i]))
				f += " (" + unitNames[i] + ")";
			names.add(f);
			if (isShowDeviations())
				names.add("+/-");
		}
		if (isShowPrecision())
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
		if (fos == null)
			return;

		StringBuilder sb = new StringBuilder();
		addStandardData(sb, 0, peak, peak, origX, origY, origValue, error, noise);

		// Add the parameters
		if (isShowDeviations())
		{
			if (paramsStdDev != null)
			{
				checkSize(converters.length, params, paramsStdDev);
				for (int i = 0; i < converters.length; i++)
				{
					add(sb, converters[i].convert(params[i]));
					add(sb, converters[i].convert(paramsStdDev[i]));
				}
			}
			else
			{
				checkSize(converters.length, params);
				for (int i = 0; i < converters.length; i++)
				{
					add(sb, converters[i].convert(params[i]));
					sb.append("\t0");
				}
			}
		}
		else
		{
			checkSize(converters.length, params);
			for (int i = 0; i < converters.length; i++)
			{
				add(sb, converters[i].convert(params[i]));
			}
		}
		if (isShowPrecision())
		{
			if (canComputePrecision)
				addPrecision(sb, calculator.getLSEPrecision(params, noise), true);
			else
				sb.append("\t0");
		}
		sb.append('\n');
		writeResult(1, sb.toString());
	}

	private void addStandardData(StringBuilder sb, final int id, final int peak, final int endPeak, final int origX,
			final int origY, final float origValue, final double error, final float noise)
	{
		if (isShowId())
		{
			sb.append(id);
			sb.append('\t');
		}
		sb.append(peak);
		if (isShowEndFrame())
			sb.append('\t').append(endPeak);
		sb.append('\t').append(origX).append('\t').append(origY);
		add(sb, origValue);
		add(sb, error);
		add(sb, converters[PeakResult.INTENSITY].convert(noise));
	}

	private void add(StringBuilder sb, float value)
	{
		sb.append('\t').append(value);
	}

	private void add(StringBuilder sb, double value)
	{
		sb.append('\t').append(value);
	}

	private void addPrecision(StringBuilder sb, double value, boolean computed)
	{
		// Cast to a float as the precision is probably limited in significant figures 
		sb.append('\t').append((float) value);
		if (computed)
			sb.append('*');
	}

	public void add(PeakResult result)
	{
		if (fos == null)
			return;

		StringBuilder sb = new StringBuilder();
		add(sb, result);
		writeResult(1, sb.toString());
	}

	private void add(StringBuilder sb, PeakResult result)
	{
		addStandardData(sb, result.getId(), result.getFrame(), result.getEndFrame(), result.getOrigX(),
				result.getOrigY(), result.getOrigValue(), result.getError(), result.getNoise());

		// Add the parameters		
		final float[] params = result.getParameters();
		if (isShowDeviations())
		{
			final float[] paramsStdDev = result.getParameterDeviations();
			if (paramsStdDev != null)
			{
				checkSize(converters.length, params, paramsStdDev);
				for (int i = 0; i < converters.length; i++)
				{
					add(sb, converters[i].convert(params[i]));
					add(sb, converters[i].convert(paramsStdDev[i]));
				}
			}
			else
			{
				checkSize(converters.length, params);
				for (int i = 0; i < converters.length; i++)
				{
					add(sb, converters[i].convert(params[i]));
					sb.append("\t0");
				}
			}
		}
		else
		{
			checkSize(converters.length, params);
			for (int i = 0; i < converters.length; i++)
			{
				add(sb, converters[i].convert(params[i]));
			}
		}
		if (isShowPrecision())
		{
			if (result.hasPrecision())
				addPrecision(sb, result.getPrecision(), false);
			else if (canComputePrecision)
				addPrecision(sb, calculator.getLSEPrecision(params, result.getNoise()), true);
			else
				sb.append("\t0");
		}
		sb.append('\n');
	}

	public void addAll(PeakResult[] results)
	{
		if (fos == null)
			return;

		int count = 0;

		StringBuilder sb = new StringBuilder();
		for (PeakResult result : results)
		{
			add(sb, result);

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
		if (fos == null)
			return;
		if (cluster.size() > 0)
		{
			float[] centroid = cluster.getCentroid();
			writeResult(0,
					String.format("#Cluster %f %f (+/-%f) n=%d\n", converters[PeakResult.X].convert(centroid[0]),
							converters[PeakResult.X].convert(centroid[1]),
							converters[PeakResult.X].convert(cluster.getStandardDeviation()), cluster.size()));
			addAll(cluster);
		}
	}

	protected void addAll(Cluster cluster)
	{
		if (!isShowId() || cluster.getId() == 0)
		{
			addAll(cluster.getPoints().toArray());
		}
		else
		{
			// Store the ID from the trace
			final int id = cluster.getId();
			final ArrayPeakResultStore results2 = new ArrayPeakResultStore(cluster.size());
			cluster.getPoints().forEach(new PeakResultProcedure()
			{
				public void execute(PeakResult result)
				{
					if (result.getId() == id)
						results2.add(result);
					else
					{
						results2.add(new ExtendedPeakResult(result.getFrame(), result.getOrigX(), result.getOrigY(),
								result.getOrigValue(), result.getError(), result.getNoise(), result.getParameters(),
								result.getParameterDeviations(), result.getEndFrame(), id));
					}
				}
			});
			addAll(results2.toArray());
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
		if (fos == null)
			return;
		if (trace.size() > 0)
		{
			float[] centroid = trace.getCentroid();
			writeResult(0, String.format("#Trace %f %f (+/-%f) n=%d, b=%d, on=%f, off=%f, signal= %f\n",
					converters[PeakResult.X].convert(centroid[0]), converters[PeakResult.X].convert(centroid[1]),
					converters[PeakResult.X].convert(trace.getStandardDeviation()), trace.size(), trace.getNBlinks(),
					trace.getOnTime(), trace.getOffTime(),
					converters[PeakResult.INTENSITY].convert(trace.getSignal())));
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
		if (fos == null)
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
		if (fos == null)
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
			TurboList<Result> results = new TurboList<Result>(size);

			StringBuilder header = new StringBuilder();
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
				for (int i = 0; i < results.size(); i++)
				{
					output.write(results.getf(i).line);
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

	/**
	 * Gets the distance unit.
	 *
	 * @return the distance unit
	 */
	public DistanceUnit getDistanceUnit()
	{
		return distanceUnit;
	}

	/**
	 * Sets the distance unit.
	 *
	 * @param distanceUnit
	 *            the new distance unit
	 */
	public void setDistanceUnit(DistanceUnit distanceUnit)
	{
		this.distanceUnit = distanceUnit;
	}

	/**
	 * Gets the intensity unit.
	 *
	 * @return the intensity unit
	 */
	public IntensityUnit getIntensityUnit()
	{
		return intensityUnit;
	}

	/**
	 * Sets the intensity unit.
	 *
	 * @param intensityUnit
	 *            the new intensity unit
	 */
	public void setIntensityUnit(IntensityUnit intensityUnit)
	{
		this.intensityUnit = intensityUnit;
	}

	/**
	 * Gets the angle unit.
	 *
	 * @return the angle unit
	 */
	public AngleUnit getAngleUnit()
	{
		return angleUnit;
	}

	/**
	 * Sets the angle unit.
	 *
	 * @param angleUnit
	 *            the new angle unit
	 */
	public void setAngleUnit(AngleUnit angleUnit)
	{
		this.angleUnit = angleUnit;
	}

	/**
	 * Checks if the precision will be computed if needed. This is only relevant if show precision is true (see
	 * {@link #isShowPrecision()}).
	 *
	 * @return true, if the precision will be computed
	 */
	public boolean isComputePrecision()
	{
		return computePrecision;
	}

	/**
	 * Sets the compute precision flag. This is only relevant if show precision is true (see
	 * {@link #isShowPrecision()}).
	 *
	 * @param computePrecision
	 *            set to true to compute the precision
	 */
	public void setComputePrecision(boolean computePrecision)
	{
		this.computePrecision = computePrecision;
	}
}
