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
		createValues();
		createDYDA();
		createD2YDA2();
		this.procedure = procedure;
		f.forEach((Gradient2Procedure) this);
	}

	protected void createD2YDA2()
	{
		if (d2yda2 == null || d2yda2.length != f.size())
			d2yda2 = new double[values.length][length];
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
