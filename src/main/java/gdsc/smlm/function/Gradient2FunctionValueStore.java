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
 * Wrap a function and store the only the values from the procedure.
 */
public class Gradient2FunctionValueStore extends ValueFunctionStore
		implements Gradient1Function, Gradient1Procedure, Gradient2Function, Gradient2Procedure
{
	private Gradient1Function f1;
	private Gradient1Procedure p1;
	private Gradient2Function f2;
	private Gradient2Procedure p2;

	/**
	 * Instantiates a new gradient 2 function store.
	 *
	 * @param f
	 *            the f
	 */
	public Gradient2FunctionValueStore(ValueFunction f)
	{
		super(f);
	}

	/**
	 * Instantiates a new gradient 2 function store.
	 *
	 * @param f
	 *            the f
	 */
	public Gradient2FunctionValueStore(Gradient1Function f)
	{
		super(f);
		this.f1 = f;
	}

	/**
	 * Instantiates a new gradient 2 function store.
	 *
	 * @param f
	 *            the f
	 */
	public Gradient2FunctionValueStore(Gradient2Function f)
	{
		super(f);
		this.f1 = f;
		this.f2 = f;
	}

	/**
	 * Instantiates a new gradient 2 function store.
	 *
	 * @param f
	 *            the f
	 * @param values
	 *            the values
	 */
	public Gradient2FunctionValueStore(ValueFunction f, double[] values)
	{
		super(f, values);
	}

	/**
	 * Instantiates a new gradient 2 function store.
	 *
	 * @param f
	 *            the f
	 * @param values
	 *            the values
	 */
	public Gradient2FunctionValueStore(Gradient1Function f, double[] values)
	{
		super(f, values);
		this.f1 = f;
	}

	/**
	 * Instantiates a new gradient 2 function store.
	 *
	 * @param f
	 *            the f
	 * @param values
	 *            the values
	 */
	public Gradient2FunctionValueStore(Gradient2Function f, double[] values)
	{
		super(f, values);
		this.f1 = f;
		this.f2 = f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#initialise(double[])
	 */
	public void initialise(double[] a)
	{
		f1.initialise(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Function#initialise1(double[])
	 */
	public void initialise1(double[] a)
	{
		f1.initialise(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#gradientIndices()
	 */
	public int[] gradientIndices()
	{
		return f1.gradientIndices();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#getNumberOfGradients()
	 */
	public int getNumberOfGradients()
	{
		return f1.getNumberOfGradients();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Function#forEach(gdsc.smlm.function.Gradient1Procedure)
	 */
	public void forEach(Gradient1Procedure procedure)
	{
		i = 0;
		createValues();
		this.p1 = procedure;
		f1.forEach((Gradient1Procedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		values[i++] = value;
		p1.execute(value, dy_da);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient2Function#initialise2(double[])
	 */
	public void initialise2(double[] a)
	{
		f2.initialise2(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient2Function#forEach(gdsc.smlm.function.Gradient2Procedure)
	 */
	public void forEach(Gradient2Procedure procedure)
	{
		i = 0;
		createValues();
		this.p2 = procedure;
		f2.forEach((Gradient2Procedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient2Procedure#execute(double, double[], double[])
	 */
	public void execute(double value, double[] dy_da, double[] d2y_da2)
	{
		values[i++] = value;
		p2.execute(value, dy_da, d2y_da2);
	}
}