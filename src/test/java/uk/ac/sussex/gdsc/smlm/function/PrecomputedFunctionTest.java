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
package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;

import uk.ac.sussex.gdsc.core.utils.PseudoRandomGenerator;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class PrecomputedFunctionTest
{
	@SeededTest
	public void precomputedValueFunctionWrapsPrecomputedValues(RandomSeed seed)
	{
		final UniformRandomProvider r = TestSettings.getRandomGenerator(seed.getSeed());
		final int size = 100;
		final double[] v = new PseudoRandomGenerator(size, r).getSequence();
		final ValueFunction f = new PrecomputedValueFunction(v);
		final double[] vo = evaluateValueFunction(f);
		Assertions.assertArrayEquals(v, vo, "values");
	}

	private static double[] evaluateValueFunction(ValueFunction f)
	{
		final double[] v = new double[f.size()];
		f.initialise0(null);
		f.forEach(new ValueProcedure()
		{
			int i = 0;

			@Override
			public void execute(double value)
			{
				v[i++] = value;
			}
		});
		return v;
	}

	@SeededTest
	public void precomputedGradient1FunctionWrapsPrecomputedValues(RandomSeed seed)
	{
		final int n = 3;
		final UniformRandomProvider r = TestSettings.getRandomGenerator(seed.getSeed());
		final int size = 100;
		final double[] v = new PseudoRandomGenerator(size, r).getSequence();
		final double[][] g1 = new double[size][];
		for (int i = 0; i < g1.length; i++)
			g1[i] = new PseudoRandomGenerator(n, r).getSequence();
		final Gradient1Function f = new PrecomputedGradient1Function(v, g1);
		final double[][] g1o = new double[size][];
		final double[] vo = evaluateGradient1Function(f, g1o);
		Assertions.assertArrayEquals(v, vo, "values");
		Assertions.assertArrayEquals(g1, g1o, "g1");
	}

	private static double[] evaluateGradient1Function(Gradient1Function f, final double[][] g1)
	{
		final double[] v = new double[f.size()];
		f.initialise1(null);
		f.forEach(new Gradient1Procedure()
		{
			int i = 0;

			@Override
			public void execute(double value, double[] dy_da)
			{
				v[i] = value;
				g1[i] = dy_da;
				i++;
			}
		});
		return v;
	}

	@SeededTest
	public void precomputedGradient2FunctionWrapsPrecomputedValues(RandomSeed seed)
	{
		final int n = 3;
		final UniformRandomProvider r = TestSettings.getRandomGenerator(seed.getSeed());
		final int size = 100;
		final double[] v = new PseudoRandomGenerator(size, r).getSequence();
		final double[][] g1 = new double[size][];
		final double[][] g2 = new double[size][];
		for (int i = 0; i < g1.length; i++)
		{
			g1[i] = new PseudoRandomGenerator(n, r).getSequence();
			g2[i] = new PseudoRandomGenerator(n, r).getSequence();
		}
		final Gradient2Function f = new PrecomputedGradient2Function(v, g1, g2);
		final double[][] g1o = new double[size][];
		final double[][] g2o = new double[size][];
		final double[] vo = evaluateGradient2Function(f, g1o, g2o);
		Assertions.assertArrayEquals(v, vo, "values");
		Assertions.assertArrayEquals(g1, g1o, "g1");
		Assertions.assertArrayEquals(g2, g2o, "g2");
	}

	private static double[] evaluateGradient2Function(Gradient2Function f, final double[][] g1, final double[][] g2)
	{
		final double[] v = new double[f.size()];
		f.initialise2(null);
		f.forEach(new Gradient2Procedure()
		{
			int i = 0;

			@Override
			public void execute(double value, double[] dy_da, double[] d2y_da2)
			{
				v[i] = value;
				g1[i] = dy_da;
				g2[i] = d2y_da2;
				i++;
			}
		});
		return v;
	}
}
