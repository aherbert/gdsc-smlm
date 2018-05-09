package gdsc.smlm.fitting;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureFactory;
import gdsc.smlm.function.FisherInformation;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.PoissonFisherInformation;
import gdsc.smlm.function.PoissonGaussianApproximationFisherInformation;
import gdsc.smlm.function.PrecomputedFunctionFactory;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;

public class UnivariateLikelihoodFisherInformationCalculatorTest
{
	@Test
	public void canComputePoissonFisherInformation()
	{
		RandomGenerator r = new Well19937c(30051977);
		for (int n = 1; n < 10; n++)
		{
			canComputePoissonFisherInformation(r, false);
		}
	}

	@Test
	public void canComputePoissonGaussianApproximationFisherInformation()
	{
		RandomGenerator r = new Well19937c(30051977);
		for (int n = 1; n < 10; n++)
		{
			canComputePoissonFisherInformation(r, true);
		}
	}

	private void canComputePoissonFisherInformation(RandomGenerator r, boolean modelGaussian)
	{
		RandomDataGenerator rdg = new RandomDataGenerator(r);

		// Create function
		Gaussian2DFunction func = GaussianFunctionFactory.create2D(1, 10, 10, GaussianFunctionFactory.FIT_ERF_CIRCLE,
				null);
		double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		params[Gaussian2DFunction.BACKGROUND] = rdg.nextUniform(0.1, 0.3);
		params[Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
		params[Gaussian2DFunction.X_POSITION] = rdg.nextUniform(4, 6);
		params[Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(4, 6);
		params[Gaussian2DFunction.X_SD] = rdg.nextUniform(1, 1.3);

		Gradient1Function f1 = func;
		FisherInformation fi;

		if (modelGaussian)
		{
			// Get a variance
			double var = 0.9 + 0.2 * r.nextDouble();
			fi = new PoissonGaussianApproximationFisherInformation(Math.sqrt(var));
			f1 = (Gradient1Function) PrecomputedFunctionFactory.wrapFunction(func,
					SimpleArrayUtils.newDoubleArray(func.size(), var));
		}
		else
		{
			fi = new PoissonFisherInformation();
		}

		// This introduces a dependency on a different package, and relies on that 
		// computing the correct answer. However that code predates this and so the
		// test ensures that the FisherInformationCalculator functions correctly.
		PoissonGradientProcedure p1 = PoissonGradientProcedureFactory.create(f1);
		p1.computeFisherInformation(params);
		double[] e = p1.getLinear();

		FisherInformationCalculator calc = new UnivariateLikelihoodFisherInformationCalculator(func, fi);
		FisherInformationMatrix I = calc.compute(params);
		double[] o = I.getMatrix().data;

		Assert.assertArrayEquals(e, o, 1e-6);

		if (!modelGaussian)
		{
			// Get the Mortensen approximation for fitting Poisson data with a Gaussian. 
			// Set a to 100 for the square pixel adjustment.
			double a = 100;
			double s = params[Gaussian2DFunction.X_SD] * a;
			double N = params[Gaussian2DFunction.SIGNAL];
			double b2 = params[Gaussian2DFunction.BACKGROUND];
			double var = Gaussian2DPeakResultHelper.getMLVarianceX(a, s, N, b2, false);

			// Convert expected variance to pixels
			var /= (a * a);

			// Get the limits by inverting the Fisher information
			double[] crlb = I.crlb();

			Assert.assertEquals(var, crlb[2], var * 1e-2);
			Assert.assertEquals(var, crlb[3], var * 1e-2);
		}
	}

	@Test
	public void canComputePerPixelPoissonGaussianApproximationFisherInformation()
	{
		RandomGenerator r = new Well19937c(30051977);
		for (int n = 1; n < 10; n++)
		{
			canComputePerPixelPoissonGaussianApproximationFisherInformation(r);
		}
	}

	private void canComputePerPixelPoissonGaussianApproximationFisherInformation(RandomGenerator r)
	{
		RandomDataGenerator rdg = new RandomDataGenerator(r);

		// Create function
		Gaussian2DFunction func = GaussianFunctionFactory.create2D(1, 10, 10, GaussianFunctionFactory.FIT_ERF_CIRCLE,
				null);
		double[] params = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		params[Gaussian2DFunction.BACKGROUND] = rdg.nextUniform(0.1, 0.3);
		params[Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
		params[Gaussian2DFunction.X_POSITION] = rdg.nextUniform(4, 6);
		params[Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(4, 6);
		params[Gaussian2DFunction.X_SD] = rdg.nextUniform(1, 1.3);

		Gradient1Function f1 = func;
		FisherInformation[] fi;

		// Get a per-pixel variance
		double[] var = new double[func.size()];
		
		fi = new FisherInformation[var.length]; 
		for (int i = var.length; i-- > 0;)
		{
			var[i] = 0.9 + 0.2 * r.nextDouble();
			fi[i] = new PoissonGaussianApproximationFisherInformation(Math.sqrt(var[i]));
		}

		f1 = (Gradient1Function) PrecomputedFunctionFactory.wrapFunction(func, var);

		// This introduces a dependency on a different package, and relies on that 
		// computing the correct answer. However that code predates this and so the
		// test ensures that the FisherInformationCalculator functions correctly.
		PoissonGradientProcedure p1 = PoissonGradientProcedureFactory.create(f1);
		p1.computeFisherInformation(params);
		double[] e = p1.getLinear();

		FisherInformationCalculator calc = new UnivariateLikelihoodFisherInformationCalculator(func, fi);
		FisherInformationMatrix I = calc.compute(params);
		double[] o = I.getMatrix().data;

		Assert.assertArrayEquals(e, o, 1e-6);
	}
}
