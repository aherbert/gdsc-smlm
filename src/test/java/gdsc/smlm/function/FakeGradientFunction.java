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

import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.PseudoRandomSequence;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class FakeGradientFunction
		implements ExtendedGradient2Function, Gradient2Function, Gradient1Function, NonLinearFunction
{
	private final int maxx, n, nparams;
	private final PseudoRandomSequence r;
	private final double[] dy_da;

	public FakeGradientFunction(int maxx, int nparams)
	{
		this(maxx, nparams, 1000, 30051977, 10.0);
	}

	public FakeGradientFunction(int maxx, int nparams, double scale)
	{
		this(maxx, nparams, 1000, 30051977, scale);
	}

	public FakeGradientFunction(int maxx, int nparams, int randomSize, int randomSeed, double scale)
	{
		this.maxx = maxx;
		this.n = maxx * maxx;
		this.nparams = nparams;
		this.r = new PseudoRandomSequence(randomSize, TestSettings.getRandomGenerator(randomSeed), scale);
		this.dy_da = new double[nparams];
	}

	@Override
	public int size()
	{
		return n;
	}

	@Override
	public void initialise(double[] a)
	{
		int seed = 0;
		for (int i = a.length; i-- > 0;)
			seed += hashCode(a[i]);
		//System.out.printf("seed = %d\n", seed);
		r.setSeed(seed);
	}

	@Override
	public void initialise0(double[] a)
	{
		initialise(a);
	}

	@Override
	public void initialise1(double[] a)
	{
		initialise(a);
	}

	@Override
	public void initialise2(double[] a)
	{
		initialise(a);
	}

	@Override
	public void initialiseExtended2(double[] a)
	{
		initialise(a);
	}

	private int hashCode(double d)
	{
		long bits = Double.doubleToLongBits(d);
		return (int) (bits ^ (bits >>> 32));
	}

	@Override
	public int[] gradientIndices()
	{
		return SimpleArrayUtils.newArray(nparams, 0, 1);
	}

	@Override
	public int getNumberOfGradients()
	{
		return nparams;
	}

	@Override
	public void forEach(ValueProcedure procedure)
	{
		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
			for (int x = 0; x < maxx; x++)
				procedure.execute(r.nextDouble());
	}

	@Override
	public void forEach(Gradient1Procedure procedure)
	{
		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				for (int j = nparams; j-- > 0;)
					dy_da[j] = r.nextDouble() * x + y;
				//System.out.println(Arrays.toString(dy_da));
				procedure.execute(r.nextDouble(), dy_da);
			}
		}
	}

	@Override
	public void forEach(Gradient2Procedure procedure)
	{
		final double[] d2y_da2 = new double[nparams];

		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				for (int j = nparams; j-- > 0;)
				{
					dy_da[j] = r.nextDouble() * x + y;
					d2y_da2[j] = r.nextDouble() * x + y;
				}
				//System.out.println(Arrays.toString(dy_da));
				procedure.execute(r.nextDouble(), dy_da, d2y_da2);
			}
		}
	}

	@Override
	public void forEach(ExtendedGradient2Procedure procedure)
	{
		final double[] d2y_dadb = new double[nparams * nparams];

		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				for (int j = nparams; j-- > 0;)
					dy_da[j] = r.nextDouble() * x + y;
				for (int j = d2y_dadb.length; j-- > 0;)
					d2y_dadb[j] = r.nextDouble() * x + y;

				//System.out.println(Arrays.toString(dy_da));
				procedure.executeExtended(r.nextDouble(), dy_da, d2y_dadb);
			}
		}
	}

	@Override
	public double eval(int i, double[] dy_da)
	{
		// Unpack the predictor to the 2D coordinates
		int y = i / maxx;
		int x = i % maxx;
		for (int j = nparams; j-- > 0;)
			dy_da[j] = r.nextDouble() * x + y;
		//System.out.println(Arrays.toString(dy_da));
		return r.nextDouble();
	}

	@Override
	public double eval(int x)
	{
		return r.nextDouble();
	}

	@Override
	public double eval(int x, double[] dyda, double[] w)
	{
		throw new NotImplementedException();
	}

	@Override
	public double evalw(int x, double[] w)
	{
		throw new NotImplementedException();
	}

	@Override
	public boolean canComputeWeights()
	{
		return false;
	}
}
