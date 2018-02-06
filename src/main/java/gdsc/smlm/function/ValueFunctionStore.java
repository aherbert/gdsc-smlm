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
 * Wrap a function and store the values from the procedure
 */
public class ValueFunctionStore implements ValueFunction, ValueProcedure
{
	private ValueFunction f;
	private ValueProcedure procedure;

	protected int i;
	/** The values from the last call to {@link #forEach(ValueProcedure)} */
	public double[] values;

	/**
	 * Instantiates a new value function store.
	 *
	 * @param f
	 *            the f
	 */
	public ValueFunctionStore(ValueFunction f)
	{
		this(f, null);
	}

	/**
	 * Instantiates a new value function store with storage.
	 *
	 * @param f
	 *            the f
	 * @param values
	 *            the values
	 */
	public ValueFunctionStore(ValueFunction f, double[] values)
	{
		this.f = f;
		this.values = values;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueFunction#size()
	 */
	public int size()
	{
		return f.size();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueFunction#initialise0(double[])
	 */
	public void initialise0(double[] a)
	{
		f.initialise0(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueFunction#forEach(gdsc.smlm.function.ValueProcedure)
	 */
	public void forEach(ValueProcedure procedure)
	{
		i = 0;
		if (values == null || values.length != f.size())
			values = new double[f.size()];
		this.procedure = procedure;
		f.forEach(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(double)
	 */
	public void execute(double value)
	{
		values[i++] = value;
		procedure.execute(value);
	}
}