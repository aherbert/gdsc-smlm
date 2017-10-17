package gdsc.smlm.function.cspline;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.DataOutput;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import gdsc.core.logging.Ticker;
import gdsc.core.logging.TrackProgress;
import gdsc.core.math.interpolation.CustomTricubicFunction;
import gdsc.core.math.interpolation.FloatCustomTricubicFunction;

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

/**
 * Stores a cubic spline data.
 */
public class CubicSplineData
{
	final int maxx;
	final int maxy;
	final CustomTricubicFunction[][] splines;

	/**
	 * Instantiates a new cubic spline data.
	 *
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param splines
	 *            the splines
	 */
	public CubicSplineData(int maxx, int maxy, CustomTricubicFunction[][] splines)
	{
		if (maxx < 1 || maxy < 1 || splines.length < 1)
			throw new IllegalArgumentException("No splines");
		int size = maxx * maxy;
		for (int z = 0; z < splines.length; z++)
			if (splines[z].length != size)
				throw new IllegalArgumentException("Incorrect XY splines size");
		this.maxx = maxx;
		this.maxy = maxy;
		this.splines = splines;
	}

	/**
	 * Instantiates a new cubic spline data.
	 *
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param splines
	 *            the splines
	 * @param dummy
	 *            the dummy flag
	 */
	private CubicSplineData(int maxx, int maxy, CustomTricubicFunction[][] splines, boolean dummy)
	{
		this.maxx = maxx;
		this.maxy = maxy;
		this.splines = splines;
	}

	/**
	 * Checks if is single precision.
	 *
	 * @return true, if is single precision
	 */
	public boolean isSinglePrecision()
	{
		return splines[0][0] instanceof FloatCustomTricubicFunction;
	}

	private static interface SplineWriter
	{
		void write(DataOutput out, CustomTricubicFunction f) throws IOException;
	}

	private static class FloatSplineWriter implements SplineWriter
	{
		public void write(DataOutput out, CustomTricubicFunction f) throws IOException
		{
			for (int i = 0; i < 64; i++)
				out.writeFloat(f.getf(i));
		}
	}

	private static class DoubleSplineWriter implements SplineWriter
	{
		public void write(DataOutput out, CustomTricubicFunction f) throws IOException
		{
			for (int i = 0; i < 64; i++)
				out.writeDouble(f.get(i));
		}
	}

	private static interface SplineReader
	{
		CustomTricubicFunction read(DataInput in) throws IOException;
	}

	private static class FloatSplineReader implements SplineReader
	{
		float[] splines = new float[64];

		public CustomTricubicFunction read(DataInput in) throws IOException
		{
			for (int i = 0; i < 64; i++)
				splines[i] = in.readFloat();
			return CustomTricubicFunction.create(splines.clone());
		}
	}

	private static class DoubleSplineReader implements SplineReader
	{
		double[] splines = new double[64];

		public CustomTricubicFunction read(DataInput in) throws IOException
		{
			for (int i = 0; i < 64; i++)
				splines[i] = in.readDouble();
			return CustomTricubicFunction.create(splines.clone());
		}
	}

	/**
	 * Write a tricubic splines to the output stream. The output will be buffered for performance.
	 *
	 * @param outputStream
	 *            the output stream
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void write(OutputStream outputStream) throws IOException
	{
		write(outputStream, null);
	}

	/**
	 * Write a tricubic splines to the output stream. The output will be buffered for performance.
	 *
	 * @param outputStream
	 *            the output stream
	 * @param progress
	 *            the progress
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void write(OutputStream outputStream, TrackProgress progress) throws IOException
	{
		// Write dimensions
		int maxz = splines.length;
		Ticker ticker = Ticker.create(progress, (long) maxx * maxy * maxz, false);
		ticker.start();
		BufferedOutputStream buffer = new BufferedOutputStream(outputStream);
		DataOutput out = new DataOutputStream(buffer);
		out.writeInt(maxx);
		out.writeInt(maxy);
		out.writeInt(maxz);
		// Write precision
		boolean singlePrecision = isSinglePrecision();
		out.writeBoolean(singlePrecision);
		SplineWriter writer = (singlePrecision) ? new FloatSplineWriter() : new DoubleSplineWriter();
		int size = maxx * maxy;
		for (int z = 0; z < maxz; z++)
		{
			for (int i = 0; i < size; i++)
			{
				ticker.tick();
				writer.write(out, splines[z][i]);
			}
		}
		ticker.stop();
		buffer.flush();
	}

	/**
	 * Read a tricubic splines from the input stream.
	 * <p>
	 * Note: For best performance a buffered input stream should be used.
	 *
	 * @param inputStream
	 *            the input stream
	 * @return the custom tricubic interpolating function
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public static CubicSplineData read(InputStream inputStream) throws IOException
	{
		return read(inputStream, null);
	}

	/**
	 * Read a tricubic splines from the input stream.
	 * <p>
	 * Note: For best performance a buffered input stream should be used.
	 *
	 * @param inputStream
	 *            the input stream
	 * @param progress
	 *            the progress
	 * @return the custom tricubic interpolating function
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public static CubicSplineData read(InputStream inputStream, TrackProgress progress) throws IOException
	{
		// Read dimensions
		BufferedInputStream buffer = new BufferedInputStream(inputStream);
		DataInput in = new DataInputStream(buffer);
		int maxx = in.readInt();
		int maxy = in.readInt();
		int maxz = in.readInt();
		Ticker ticker = Ticker.create(progress, (long) maxx * maxy * maxz, false);
		ticker.start();
		// Read precision
		boolean singlePrecision = in.readBoolean();
		SplineReader reader = (singlePrecision) ? new FloatSplineReader() : new DoubleSplineReader();
		int size = maxx * maxy;
		CustomTricubicFunction[][] splines = new CustomTricubicFunction[maxz][maxx * maxy];
		for (int z = 0; z < maxz; z++)
		{
			for (int i = 0; i < size; i++)
			{
				ticker.tick();
				splines[z][i] = reader.read(in);
			}
		}
		ticker.stop();
		// Skip validation
		return new CubicSplineData(maxx, maxy, splines, false);
	}
}