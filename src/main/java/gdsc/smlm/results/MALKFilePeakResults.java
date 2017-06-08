package gdsc.smlm.results;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.InputMismatchException;
import java.util.NoSuchElementException;
import java.util.Scanner;

import gdsc.core.ij.Utils;

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
import gdsc.smlm.units.DistanceUnit;
import gdsc.smlm.units.IdentityUnitConverter;
import gdsc.smlm.units.IntensityUnit;

/**
 * Saves the fit results to file using the simple MALK file format (Molecular Accuracy Localisation Keep). This consists
 * of [X,Y,T,Signal] data in a white-space separated format. Comments are allowed with the # character.
 */
public class MALKFilePeakResults extends FilePeakResults
{
	private OutputStreamWriter out;

	public MALKFilePeakResults(String filename)
	{
		super(filename);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#getHeaderEnd()
	 */
	protected String getHeaderEnd()
	{
		return null;
	}

	@Override
	protected String getVersion()
	{
		return "MALK";
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

	@Override
	public void begin()
	{
		super.begin();

		// Ensure we write out in nm and photons if possible.
		// If converters were not created then use dummy converters.

		// Copy it so it can be modified
		setCalibration(calibration.clone());

		if (toNMConverter == null)
			toNMConverter = new IdentityUnitConverter<DistanceUnit>(DistanceUnit.NM);
		else
			calibration.setDistanceUnit(DistanceUnit.NM);

		if (toPhotonConverter == null)
			toPhotonConverter = new IdentityUnitConverter<IntensityUnit>(IntensityUnit.PHOTON);
		else
			calibration.setIntensityUnit(IntensityUnit.PHOTON);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#getHeaderComments()
	 */
	protected String[] getHeaderComments()
	{
		String[] comments = new String[3];
		int count = 0;
		if (calibration != null)
		{
			if (calibration.hasNmPerPixel())
			{
				comments[count++] = String.format("Pixel pitch %s (nm)", Utils.rounded(calibration.getNmPerPixel()));
			}
			if (calibration.hasGain())
			{
				comments[count++] = String.format("Gain %s (Count/photon)", Utils.rounded(calibration.getGain()));
			}
			if (calibration.hasExposureTime())
			{
				comments[count++] = String.format("Exposure time %s (seconds)",
						Utils.rounded(calibration.getExposureTime() * 1e-3));
			}
		}
		return Arrays.copyOf(comments, count);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FilePeakResults#getFieldNames()
	 */
	protected String[] getFieldNames()
	{
		String[] names = new String[] { "X", "Y", "Frame", "Signal" };
		if (toNMConverter != null)
		{
			names[0] += " (nm)";
			names[1] += " (nm)";
		}
		if (toPhotonConverter != null)
		{
			names[3] += " (photon)";
		}
		return names;
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

		StringBuilder sb = new StringBuilder(100);

		addStandardData(sb, params[Gaussian2DFunction.X_POSITION], params[Gaussian2DFunction.Y_POSITION], peak,
				params[Gaussian2DFunction.SIGNAL]);

		writeResult(1, sb.toString());
	}

	private void addStandardData(StringBuilder sb, final float x, final float y, final int frame, final float signal)
	{
		sb.append(toNMConverter.convert(x));
		sb.append('\t');
		sb.append(toNMConverter.convert(y));
		sb.append('\t');
		sb.append(frame);
		sb.append('\t');
		sb.append(toPhotonConverter.convert(signal));
		sb.append('\n');
	}

	public void addAll(Collection<PeakResult> results)
	{
		if (out == null)
			return;

		int count = 0;

		StringBuilder sb = new StringBuilder(2000);
		for (PeakResult result : results)
		{
			// Add the standard data
			addStandardData(sb, result.getXPosition(), result.getYPosition(), result.getFrame(), result.getSignal());

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

	protected void addAll(Cluster cluster)
	{
		addAll(cluster.getPoints());
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
						header.append(line).append("\n");
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
				scanner.nextFloat(); // X
				scanner.nextFloat(); // Y
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
