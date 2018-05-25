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

/**
 * Wraps a value function to add a pre-computed offset to the value during the forEach procedure
 */
public class OffsetGradient1Function extends OffsetValueFunction
		implements Gradient1Function, Gradient1Procedure, NonLinearFunction
{
	protected final Gradient1Function f1;
	protected Gradient1Procedure procedure;

	/**
	 * Class for evaluating a function and storing the values and gradients
	 */
	protected class FunctionStore implements Gradient1Procedure
	{
		private int i;

		public final double[] values;
		public final double[][] dyda;
		public final int length;

		public FunctionStore(double[] values, double[][] dyda)
		{
			length = f1.getNumberOfGradients();
			if (values == null)
			{
				values = new double[f1.size()];
				dyda = new double[values.length][length];
			}
			this.values = values;
			this.dyda = dyda;
		}

		/**
		 * Gets the values.
		 */
		public void getValues()
		{
			i = 0;
			f1.forEach(this);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
		 */
		public void execute(double value, double[] dy_da)
		{
			values[i] = value;
			System.arraycopy(dy_da, 0, dyda[i], 0, length);
			i++;
		}
	}

	// Used to store all the values and gradients for the NonLinearFunction interface
	protected FunctionStore store = null;
	protected double[] all_values;
	protected double[][] all_dyda;

	/**
	 * Instantiates a new offset gradient1 function.
	 *
	 * @param f
	 *            the function
	 * @param values
	 *            the precomputed values
	 * @throws IllegalArgumentException
	 *             if the values length does not match the function size
	 */
	protected OffsetGradient1Function(Gradient1Function f, double[] values)
	{
		super(f, values);
		f1 = f;
	}

	protected OffsetGradient1Function(OffsetGradient1Function pre, double[] values2)
	{
		super(pre, values2);
		f1 = (Gradient1Function) f;
	}

	public Gradient1Function getGradient1Function()
	{
		return f1;
	}

	public void initialise(double[] a)
	{
		store = null;
		f1.initialise(a);
		i = 0;
	}

	public int[] gradientIndices()
	{
		return f1.gradientIndices();
	}

	public int getNumberOfGradients()
	{
		return f1.getNumberOfGradients();
	}

	public void initialise1(double[] a)
	{
		f1.initialise1(a);
	}

	public void forEach(Gradient1Procedure procedure)
	{
		this.procedure = procedure;
		i = 0;
		f1.forEach((Gradient1Procedure) this);
	}

	public void execute(double value, double[] dy_da)
	{
		procedure.execute(value + values[i++], dy_da);
	}

	/**
	 * Wrap a function with pre-computed values.
	 *
	 * @param func
	 *            the function
	 * @param b
	 *            Baseline pre-computed y-values
	 * @return the wrapped function (or the original if pre-computed values are null or wrong length)
	 */
	public static Gradient1Function wrapGradient1Function(final Gradient1Function func, final double[] b)
	{
		if (b != null && b.length == func.size())
		{
			// Avoid multiple wrapping
			if (func instanceof OffsetGradient1Function)
			{
				return new OffsetGradient1Function((OffsetGradient1Function) func, b);
			}
			return new OffsetGradient1Function(func, b);
		}
		return func;
	}

	public double eval(int x, double[] dyda)
	{
		createStore();
		System.arraycopy(all_dyda[i], 0, dyda, 0, store.length);
		return all_values[x];
	}

	private void createStore()
	{
		if (store == null)
		{
			store = new FunctionStore(all_values, all_dyda);
			store.getValues();
			// Re-use space
			all_values = store.values;
			all_dyda = store.dyda;
		}
	}

	public double eval(int x)
	{
		createStore();
		return store.values[x];
	}

	public double eval(int x, double[] dyda, double[] w)
	{
		w[0] = 1;
		return eval(x, dyda);
	}

	public double evalw(int x, double[] w)
	{
		w[0] = 1;
		return eval(x);
	}

	public boolean canComputeWeights()
	{
		return false;
	}
}
