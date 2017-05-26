package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.function.FakeGradientFunction;
import gdsc.smlm.function.gaussian.HoltzerAstimatismZModel;
import gdsc.smlm.function.gaussian.erf.ErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;

/**
 * Contains speed tests for the methods for calculating the Hessian and gradient vector
 * for use in the LVM algorithm.
 */
public class FastMLEJacobianGradient2ProcedureTest extends FastMLEGradient2ProcedureTest
{
	// Skip super-class tests ...
	@Test
	public void gradientProcedureFactoryCreatesOptimisedProcedures()
	{
		Assume.assumeTrue(false);
	}

	@Test
	public void gradientProcedureComputesSameLogLikelihoodAsMLEGradientCalculator()
	{
		Assume.assumeTrue(false);
	}

	@Test
	public void gradientProcedureIsNotSlowerThanGradientCalculator()
	{
		Assume.assumeTrue(false);
	}

	@Test
	public void gradientProcedureComputesSameWithPrecomputed()
	{
		Assume.assumeTrue(false);
	}

	@Test
	public void gradientProcedureUnrolledComputesSameAsGradientProcedure()
	{
		Assume.assumeTrue(false);
	}

	@Test
	public void gradientProcedureIsFasterUnrolledThanGradientProcedure()
	{
		Assume.assumeTrue(false);
	}

	@Test
	public void gradientProcedureComputesSameAsBaseGradientProcedure()
	{
		// Test the base functionality of computing the partial derivatives is the same
		gradientProcedureComputesSameAsBaseGradientProcedure(4);
		gradientProcedureComputesSameAsBaseGradientProcedure(5);
		gradientProcedureComputesSameAsBaseGradientProcedure(6);
		gradientProcedureComputesSameAsBaseGradientProcedure(11);
		gradientProcedureComputesSameAsBaseGradientProcedure(21);
	}

	private void gradientProcedureComputesSameAsBaseGradientProcedure(int nparams)
	{
		int iter = 10;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createFakeData(nparams, iter, paramsList, yList);
		FakeGradientFunction func = new FakeGradientFunction(blockWidth, nparams);

		for (int i = 0; i < paramsList.size(); i++)
		{
			FastMLEGradient2Procedure p = FastMLEGradient2ProcedureFactory.createUnrolled(yList.get(i), func);
			FastMLEJacobianGradient2Procedure p2 = new FastMLEJacobianGradient2Procedure(yList.get(i), func);
			p.computeSecondDerivative(paramsList.get(i));
			p2.computeSecondDerivative(paramsList.get(i));
			// Virtually the same ...
			Assert.assertArrayEquals(p.d1, p2.d1, 1e-5);
			Assert.assertArrayEquals(p.d2, p2.d2, 1e-5);
		}
	}

	@Test
	public void gradientCalculatorComputesGradient()
	{
		gradientCalculatorComputesGradient(new SingleFreeCircularErfGaussian2DFunction(blockWidth, blockWidth));

		// Use a reasonable z-depth function from the Smith, et al (2010) paper (page 377)
		double gamma = 0.389;
		double d = 0.531;
		double Ax = -0.0708;
		double Bx = -0.073;
		double Ay = 0.164;
		double By = 0.0417;
		HoltzerAstimatismZModel zModel = HoltzerAstimatismZModel.create(gamma, d, Ax, Bx, Ay, By);
		gradientCalculatorComputesGradient(new SingleAstigmatismErfGaussian2DFunction(blockWidth, blockWidth, zModel));
	}

	private void gradientCalculatorComputesGradient(ErfGaussian2DFunction func)
	{
		// Check the first and second derivatives
		int nparams = func.getNumberOfGradients();
		int[] indices = func.gradientIndices();

		int iter = 100;
		rdg = new RandomDataGenerator(new Well19937c(30051977));

		ArrayList<double[]> paramsList = new ArrayList<double[]>(iter);
		ArrayList<double[]> yList = new ArrayList<double[]>(iter);

		createData(1, iter, paramsList, yList, true);

		double delta = 1e-5;
		DoubleEquality eq = new DoubleEquality(1e-4, 1e-3);

		for (int i = 0; i < paramsList.size(); i++)
		{
			double[] y = yList.get(i);
			double[] a = paramsList.get(i);
			double[] a2 = a.clone();
			FastMLEJacobianGradient2Procedure p = new FastMLEJacobianGradient2Procedure(y, func);
			//double ll = p.computeLogLikelihood(a);
			p.computeSecondDerivative(a);
			double[] d1 = p.d1.clone();
			double[] d2 = p.d2.clone();
			DenseMatrix64F J = DenseMatrix64F.wrap(nparams, nparams, p.J.clone());
			for (int j = 0; j < nparams; j++)
			{
				int k = indices[j];
				double d = Precision.representableDelta(a[k], (a[k] == 0) ? delta : a[k] * delta);
				a2[k] = a[k] + d;
				double llh = p.computeLogLikelihood(a2);
				p.computeFirstDerivative(a2);
				double[] d1h = p.d1.clone();
				a2[k] = a[k] - d;
				double lll = p.computeLogLikelihood(a2);
				p.computeFirstDerivative(a2);
				double[] d1l = p.d1.clone();
				a2[k] = a[k];

				double gradient1 = (llh - lll) / (2 * d);
				double gradient2 = (d1h[j] - d1l[j]) / (2 * d);
				//System.out.printf("[%d,%d] ll - %f  (%s %f+/-%f) d1 %f ?= %f : d2 %f ?= %f\n", i, k, ll, func.getName(k), a[k], d, 
				//		gradient1, d1[j], gradient2, d2[j]);
				Assert.assertTrue("Not same gradient1 @ " + j, eq.almostEqualRelativeOrAbsolute(gradient1, d1[j]));
				Assert.assertTrue("Not same gradient2 @ " + j, eq.almostEqualRelativeOrAbsolute(gradient2, d2[j]));

				// Test the Jacobian ...

				for (int jj = 0; jj < nparams; jj++)
				{
					if (j == jj)
					{
						// This is done above
						// Check it anyway to ensure the Jacobian is correct
						//continue;
					}

					int kk = indices[jj];
					double dd = Precision.representableDelta(a[kk], (a[kk] == 0) ? delta : a[kk] * delta);
					//a2[k] = a[k] + d;
					a2[kk] = a[kk] + dd;
					//llh = p.computeLogLikelihood(a2);
					p.computeFirstDerivative(a2);
					d1h = p.d1.clone();
					//a2[k] = a[k] - d;
					a2[kk] = a[kk] - dd;
					//lll = p.computeLogLikelihood(a2);
					p.computeFirstDerivative(a2);
					d1l = p.d1.clone();
					//a2[k] = a[k];
					a2[kk] = a[kk];

					//gradient1 = (llh - lll) / (2 * d) / (2 * dd);
					gradient2 = (d1h[jj] - d1l[jj]) / (2 * dd);
					//System.out.printf("[%d,%d,%d] (%s %f  %s %f+/-%f) J %f ?= %f\n", 
					//		i, k, kk, func.getName(k), a[k], func.getName(kk), a[kk], dd,  
					//		gradient2, J.get(j, jj));
					Assert.assertTrue(String.format("Not same gradientJ @ [%d,%d]", j, jj), 
							eq.almostEqualRelativeOrAbsolute(gradient2, J.get(j, jj)));
				}
			}
		}
	}
}
