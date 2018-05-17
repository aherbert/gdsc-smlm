package gdsc.smlm.function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Implement the FastLog methods using Math.log
 */
public class NonFastLog extends FastLog
{
	public static final NonFastLog INSTANCE = new NonFastLog();

	@Override
	public double getBase()
	{
		return Math.E;
	}

	@Override
	public double getScale()
	{
		return LN2;
	}

	@Override
	public int getN()
	{
		return 52;
	}

	@Override
	public float log2(float x)
	{
		return (float) (Math.log(x) / LN2);
	}

	@Override
	public float fastLog2(float x)
	{
		return (float) (Math.log(x) / LN2);
	}

	@Override
	public float log(float x)
	{
		return (float) Math.log(x);
	}

	@Override
	public float fastLog(float x)
	{
		return (float) Math.log(x);
	}

	@Override
	public float log2(double x)
	{
		return (float) (Math.log(x) / LN2);
	}

	@Override
	public float fastLog2(double x)
	{
		return (float) (Math.log(x) / LN2);
	}

	@Override
	public float log(double x)
	{
		return (float) Math.log(x);
	}

	@Override
	public float fastLog(double x)
	{
		return (float) Math.log(x);
	}

	@Override
	public double log2D(double x)
	{
		return Math.log(x) / LN2;
	}

	@Override
	public double fastLog2D(double x)
	{
		return Math.log(x) / LN2;
	}

	@Override
	public double logD(double x)
	{
		return Math.log(x);
	}

	@Override
	public double fastLogD(double x)
	{
		return Math.log(x);
	}
}
