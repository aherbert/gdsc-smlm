package gdsc.smlm.function;

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

import java.util.Arrays;

import gnu.trove.list.array.TDoubleArrayList;

/**
 * Allow optimisation using Apache Commons Math 3 Optimiser
 */
public abstract class OptimiserFunction
{
	protected TDoubleArrayList x = null;
	protected TDoubleArrayList y = null;

	public void addPoint(double x, double y)
	{
		if (this.x == null)
		{
			this.x = new TDoubleArrayList();
			this.y = new TDoubleArrayList();
		}
		this.x.add(x);
		this.y.add(y);
	}

	public void addData(float[] x, float[] y)
	{
		this.x = new TDoubleArrayList();
		this.y = new TDoubleArrayList();
		for (int i = 0; i < x.length; i++)
		{
			this.x.add((double) x[i]);
			this.y.add((double) y[i]);
		}
	}

	public void addData(double[] x, double[] y)
	{
		this.x = new TDoubleArrayList();
		this.y = new TDoubleArrayList();
		for (int i = 0; i < x.length; i++)
		{
			this.x.add(x[i]);
			this.y.add(y[i]);
		}
	}

	public double[] getX()
	{
		return x.toArray();
	}

	public double[] getY()
	{
		return y.toArray();
	}

	public double[] getWeights()
	{
		double[] w = new double[y.size()];
		Arrays.fill(w, 1);
		return w;
	}

	public int size()
	{
		return x.size();
	}
}