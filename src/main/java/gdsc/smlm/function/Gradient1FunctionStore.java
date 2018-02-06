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
public class Gradient1FunctionStore extends ValueFunctionStore implements Gradient1Function, Gradient1Procedure
{
	private Gradient1Function f;
	private Gradient1Procedure procedure;
	protected final int length;

	/** The gradients from the last call to {@link #forEach(Gradient1Procedure)}. */
	public double[][] dyda;

	/**
	 * Instantiates a new gradient 1 function store.
	 *
	 * @param f
	 *            the f
	 */
	public Gradient1FunctionStore(Gradient1Function f)
	{
		this(f, null, null);
	}

	/**
	 * Instantiates a new gradient 1 function store with storage.
	 *
	 * @param f
	 *            the f
	 * @param values
	 *            the values
	 * @param dyda
	 *            the dyda
	 */
	public Gradient1FunctionStore(Gradient1Function f, double[] values, double[][] dyda)
	{
		super(f, values);
		this.f = f;
		this.dyda = dyda;
		length = f.getNumberOfGradients();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#initialise(double[])
	 */
	public void initialise(double[] a)
	{
		f.initialise(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Function#initialise1(double[])
	 */
	public void initialise1(double[] a)
	{
		f.initialise(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return f.gradientIndices();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#getNumberOfGradients()
	 */
	public int getNumberOfGradients()
	{
		return f.getNumberOfGradients();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Function#forEach(gdsc.smlm.function.Gradient1Procedure)
	 */
	public void forEach(Gradient1Procedure procedure)
	{
		i = 0;
		if (values ==null || values.length != f.size())
			values = new double[f.size()];
		if (dyda==null || dyda.length != f.size())
			dyda = new double[values.length][length];
		this.procedure = procedure;
		f.forEach((Gradient1Procedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		values[i] = value;
		System.arraycopy(dy_da[i], 0, dyda[i], 0, length);
		i++;
		procedure.execute(value, dy_da);
	}
}