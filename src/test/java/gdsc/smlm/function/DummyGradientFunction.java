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

public class DummyGradientFunction implements Gradient1Function, Gradient2Function
{
	int n;

	public DummyGradientFunction(int n)
	{
		this.n = n;
	}

	@Override
	public int size()
	{
		return 0;
	}

	@Override
	public void initialise0(double[] a)
	{
	}

	@Override
	public void initialise1(double[] a)
	{
	}

	@Override
	public void initialise2(double[] a)
	{
	}

	@Override
	public int[] gradientIndices()
	{
		return null;
	}

	@Override
	public int getNumberOfGradients()
	{
		return n;
	}

	@Override
	public void forEach(ValueProcedure procedure)
	{
	}

	@Override
	public void forEach(Gradient1Procedure procedure)
	{
	}

	@Override
	public void forEach(Gradient2Procedure procedure)
	{
	}

	@Override
	public void initialise(double[] a)
	{
	}
}
