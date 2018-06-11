/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.function;


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
		createValues();
		createDYDA();
		this.procedure = procedure;
		f.forEach((Gradient1Procedure) this);
	}

	protected void createDYDA()
	{
		if (dyda==null || dyda.length != f.size())
			dyda = new double[values.length][length];
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
