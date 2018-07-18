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
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.TestSettings;

/**
 * Test that a bounded fitter can return the same results with and without bounds.
 */
@SuppressWarnings({ "javadoc" })
public class BoundedFunctionSolverTest extends BaseFunctionSolverTest
{
	// The following tests ensure that the LVM can fit data without
	// requiring a bias (i.e. an offset to the background).
	// In a previous version the LVM fitter was stable only if a bias existed.
	// The exact source of this instability is unknown as it could be due to
	// how the data was processed before or after fitting, or within the LVM
	// fitter itself. However the process should be the same without a bias
	// and these tests ensure that is true.

	@Test
	public void fitSingleGaussianLVMWithoutBias()
	{
		fitSingleGaussianLVMWithoutBias(false, 0);
	}

	@Test
	public void fitSingleGaussianCLVMWithoutBias()
	{
		fitSingleGaussianLVMWithoutBias(false, 1);
	}

	@Test
	public void fitSingleGaussianDCLVMWithoutBias()
	{
		fitSingleGaussianLVMWithoutBias(false, 2);
	}

	@Test
	public void fitSingleGaussianBLVMWithoutBias()
	{
		fitSingleGaussianLVMWithoutBias(true, 0);
	}

	@Test
	public void fitSingleGaussianBCLVMWithoutBias()
	{
		fitSingleGaussianLVMWithoutBias(true, 1);
	}

	@Test
	public void fitSingleGaussianBDCLVMWithoutBias()
	{
		fitSingleGaussianLVMWithoutBias(true, 2);
	}

	private void fitSingleGaussianLVMWithoutBias(boolean applyBounds, int clamping)
	{
		final double bias = 100;

		final NonLinearFit solver = getLVM((applyBounds) ? 2 : 1, clamping, false);
		final NonLinearFit solver2 = getLVM((applyBounds) ? 2 : 1, clamping, false);

		final String name = getLVMName(applyBounds, clamping, false);

		final int LOOPS = 5;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final StoredDataStatistics[] stats = new StoredDataStatistics[6];

		for (final double s : signal)
		{
			final double[] expected = createParams(1, s, 0, 0, 1);
			double[] lower = null, upper = null;
			if (applyBounds)
			{
				lower = createParams(0, s * 0.5, -0.2, -0.2, 0.8);
				upper = createParams(3, s * 2, 0.2, 0.2, 1.2);
				solver.setBounds(lower, upper);
			}

			final double[] expected2 = addBiasToParams(expected, bias);
			if (applyBounds)
			{
				final double[] lower2 = addBiasToParams(lower, bias);
				final double[] upper2 = addBiasToParams(upper, bias);
				solver2.setBounds(lower2, upper2);
			}

			for (int loop = LOOPS; loop-- > 0;)
			{
				final double[] data = drawGaussian(expected, rg);
				final double[] data2 = data.clone();
				for (int i = 0; i < data.length; i++)
					data2[i] += bias;

				for (int i = 0; i < stats.length; i++)
					stats[i] = new StoredDataStatistics();

				for (final double db : base)
					for (final double dx : shift)
						for (final double dy : shift)
							for (final double dsx : factor)
							{
								final double[] p = createParams(db, s, dx, dy, dsx);
								final double[] p2 = addBiasToParams(p, bias);

								final double[] fp = fitGaussian(solver, data, p, expected);
								final double[] fp2 = fitGaussian(solver2, data2, p2, expected2);

								// The result should be the same without a bias
								Assert.assertEquals(name + " Iterations", solver.getEvaluations(),
										solver2.getEvaluations());
								fp2[0] -= bias;
								Assert.assertArrayEquals(name + " Solution", fp, fp2, 1e-6);
							}
			}
		}
	}

	// Standard LVM
	@Test
	public void canFitSingleGaussianLVM()
	{
		fitSingleGaussianLVM(0, 0, false);
	}

	// Bounded/Clamped LVM

	@Test
	public void canFitSingleGaussianBLVMNoBounds()
	{
		fitSingleGaussianLVM(1, 0, false);
	}

	@Test
	public void canFitSingleGaussianBLVM()
	{
		fitSingleGaussianLVM(2, 0, false);
	}

	@Test
	public void canFitSingleGaussianCLVM()
	{
		fitSingleGaussianLVM(0, 1, false);
	}

	@Test
	public void canFitSingleGaussianDCLVM()
	{
		fitSingleGaussianLVM(0, 2, false);
	}

	@Test
	public void canFitSingleGaussianBCLVM()
	{
		fitSingleGaussianLVM(2, 1, false);
	}

	@Test
	public void canFitSingleGaussianBDCLVM()
	{
		fitSingleGaussianLVM(2, 2, false);
	}

	// MLE LVM

	@Test
	public void canFitSingleGaussianLVMMLE()
	{
		fitSingleGaussianLVM(0, 0, true);
	}

	@Test
	public void canFitSingleGaussianBLVMMLENoBounds()
	{
		fitSingleGaussianLVM(1, 0, true);
	}

	@Test
	public void canFitSingleGaussianBLVMMLE()
	{
		fitSingleGaussianLVM(2, 0, true);
	}

	private void fitSingleGaussianLVM(int bounded, int clamping, boolean mle)
	{
		canFitSingleGaussian(getLVM(bounded, clamping, mle), bounded == 2);
	}

	// Is Bounded/Clamped LVM better?

	@Test
	public void fitSingleGaussianBLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 0, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 1, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 1, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianDCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 2, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBDCLVMBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 2, false, false, 0, false);
	}

	@Test
	public void fitSingleGaussianLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 0, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 0, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 1, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 1, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianDCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(false, 2, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBDCLVMMLEBetterThanLVM()
	{
		fitSingleGaussianBetterLVM(true, 2, true, false, 0, false);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(true, 0, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianCLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(false, 1, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianDCLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(false, 2, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianBDCLVMMLEBetterThanLVMMLE()
	{
		fitSingleGaussianBetterLVM(true, 2, true, false, 0, true);
	}

	@Test
	public void fitSingleGaussianBLVMMLEBetterThanBLVM()
	{
		fitSingleGaussianBetterLVM(true, 0, true, true, 0, false);
	}

	@Test
	public void fitSingleGaussianBCLVMMLEBetterThanBCLVM()
	{
		fitSingleGaussianBetterLVM(true, 1, true, true, 1, false);
	}

	@Test
	public void fitSingleGaussianBDCLVMMLEBetterThanBDCLVM()
	{
		fitSingleGaussianBetterLVM(true, 2, true, true, 2, false);
	}

	private void fitSingleGaussianBetterLVM(boolean bounded2, int clamping2, boolean mle2, boolean bounded,
			int clamping, boolean mle)
	{
		TestSettings.assumeMediumComplexity();
		final NonLinearFit solver = getLVM((bounded) ? 2 : 1, clamping, mle);
		final NonLinearFit solver2 = getLVM((bounded2) ? 2 : 1, clamping2, mle2);
		canFitSingleGaussianBetter(solver, bounded, solver2, bounded2, getLVMName(bounded, clamping, mle),
				getLVMName(bounded2, clamping2, mle2));
	}

	private static NonLinearFit getLVM(int bounded, int clamping, boolean mle)
	{
		final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, flags, null);
		final StoppingCriteria sc = new ErrorStoppingCriteria(5);
		sc.setMaximumIterations(100);
		final NonLinearFit solver = (bounded != 0 || clamping != 0) ? new BoundedNonLinearFit(f, sc, null)
				: new NonLinearFit(f, sc);
		if (clamping != 0)
		{
			final BoundedNonLinearFit bsolver = (BoundedNonLinearFit) solver;
			final ParameterBounds bounds = new ParameterBounds(f);
			bounds.setClampValues(defaultClampValues);
			bounds.setDynamicClamp(clamping == 2);
			bsolver.setBounds(bounds);
		}
		solver.setMLE(mle);
		solver.setInitialLambda(1);
		return solver;
	}

	private static String getLVMName(boolean bounded, int clamping, boolean mle)
	{
		return ((bounded) ? "B" : "") + ((clamping == 0) ? "" : ((clamping == 1) ? "C" : "DC")) + "LVM" +
				((mle) ? " MLE" : "");
	}

	private class CheatingStoppingCriteria extends StoppingCriteria
	{
		@Override
		public void evaluate(double oldError, double newError, double[] a)
		{
			// Do nothing
		}

		@Override
		public boolean areAchieved()
		{
			return true;
		}

		@Override
		public boolean areNotSatisfied()
		{
			return false;
		}
	}

	@Test
	public void canFitAndComputeDeviationsLVM()
	{
		canFitAndComputeDeviationsLVM(false);
	}

	@Test
	public void canFitAndComputeDeviationsLVMMLE()
	{
		canFitAndComputeDeviationsLVM(true);
	}

	private void canFitAndComputeDeviationsLVM(boolean mle)
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final NonLinearFit solver1 = getLVM(0, 0, mle);
		final NonLinearFit solver2 = getLVM(0, 0, mle);
		solver1.setStoppingCriteria(new CheatingStoppingCriteria());
		fitAndComputeDeviationsMatch(rg, solver1, solver2, NoiseModel.EMCCD, false);
	}
}
