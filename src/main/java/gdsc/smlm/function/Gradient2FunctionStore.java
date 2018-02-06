package gdsc.smlm.function;

// TODO: Auto-generated Javadoc
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
 * Wrap a function and store the values from the procedure.
 */
public class Gradient2FunctionStore extends Gradient1FunctionStore implements Gradient2Function, Gradient2Procedure
{
	private Gradient2Function f;
	private Gradient2Procedure procedure;

	/** The gradients from the last call to {@link #forEach(Gradient2Procedure)}. */
	public double[][] d2yda2;

	/**
	 * Instantiates a new gradient 2 function store.
	 *
	 * @param f
	 *            the f
	 */
	public Gradient2FunctionStore(Gradient2Function f)
	{
		this(f, null, null, null);
	}

	/**
	 * Instantiates a new gradient 2 function store with storage.
	 *
	 * @param f
	 *            the f
	 * @param values
	 *            the values
	 * @param dyda
	 *            the dyda
	 * @param d2yda2
	 *            the d2yda2
	 */
	public Gradient2FunctionStore(Gradient2Function f, double[] values, double[][] dyda, double[][] d2yda2)
	{
		super(f, values, dyda);
		this.f = f;
		this.d2yda2 = d2yda2;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient2Function#initialise2(double[])
	 */
	public void initialise2(double[] a)
	{
		f.initialise2(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient2Function#forEach(gdsc.smlm.function.Gradient2Procedure)
	 */
	public void forEach(Gradient2Procedure procedure)
	{
		i = 0;
		if (values == null || values.length != f.size())
			values = new double[f.size()];
		if (dyda == null || dyda.length != f.size())
			dyda = new double[values.length][length];
		if (d2yda2 == null || d2yda2.length != f.size())
			d2yda2 = new double[values.length][length];
		this.procedure = procedure;
		f.forEach((Gradient2Procedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient2Procedure#execute(double, double[], double[])
	 */
	public void execute(double value, double[] dy_da, double[] d2y_da2)
	{
		values[i] = value;
		System.arraycopy(dy_da[i], 0, dyda[i], 0, length);
		System.arraycopy(d2y_da2[i], 0, d2yda2[i], 0, length);
		i++;
		procedure.execute(value, dy_da, d2y_da2);
	}
}