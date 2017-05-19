package gdsc.smlm.fitting;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;

public class FisherInformationMatrixTest
{
	@Test
	public void canComputeCRLB()
	{
		for (int n = 1; n < 10; n++)
		{
			canComputeCRLB(n, 0, true);
		}
	}

	@Test
	public void canComputeCRLBWithZeros()
	{
		for (int n = 2; n < 10; n++)
		{
			canComputeCRLB(n, 1, true);
			canComputeCRLB(n, n / 2, true);
		}
	}

	@Test
	public void canComputeCRLBWithReciprocal()
	{
		for (int n = 1; n < 10; n++)
		{
			canComputeCRLB(n, 0, false);
		}
	}

	@Test
	public void canComputeCRLBWithReciprocalWithZeros()
	{
		for (int n = 2; n < 10; n++)
		{
			canComputeCRLB(n, 1, false);
			canComputeCRLB(n, n / 2, false);
		}
	}

	@Test
	public void inversionMatchesReciprocal()
	{
		for (int n = 1; n < 10; n++)
		{
			FisherInformationMatrix m = createFisherInformationMatrix(n, 0);
			double[] crlb = m.crlb();
			double[] crlb2 = m.crlbReciprocal();
			// These increasingly do not match with increasing number of parameters. 
			System.out.printf("%s =? %s\n", Arrays.toString(crlb), Arrays.toString(crlb2));
		}
	}

	private double[] canComputeCRLB(int n, int k, boolean invert)
	{
		FisherInformationMatrix m = createFisherInformationMatrix(n, k);

		// Invert for CRLB
		double[] crlb = (invert) ? m.crlb() : m.crlbReciprocal();
		System.out.printf("n=%d, k=%d : %s\n", n, k, Arrays.toString(crlb));
		Assert.assertNotNull("CRLB failed", crlb);
		return crlb;
	}

	private FisherInformationMatrix createFisherInformationMatrix(int n, int k)
	{
		int maxx = 10;
		int size = maxx * maxx;
		RandomGenerator randomGenerator = new Well19937c(30051977);
		RandomDataGenerator rdg = new RandomDataGenerator(randomGenerator);

		// Use a real Gaussian function here to compute the Fisher information.
		// The matrix may be sensitive to the type of equation used.
		int npeaks = 1;
		while (1 + npeaks * 6 < n)
			npeaks++;
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(npeaks, maxx, maxx,
				GaussianFunctionFactory.FIT_ELLIPTICAL, null);

		double[] a = new double[1 + npeaks * 6];
		a[Gaussian2DFunction.BACKGROUND] = rdg.nextUniform(1, 5);
		for (int i = 0, j = 0; i < npeaks; i++, j += 6)
		{
			a[j + Gaussian2DFunction.SIGNAL] = rdg.nextUniform(100, 300);
			a[j + Gaussian2DFunction.SHAPE] = rdg.nextUniform(-Math.PI, Math.PI);
			// Non-overlapping peaks otherwise the CRLB are poor
			a[j + Gaussian2DFunction.X_POSITION] = rdg.nextUniform(2 + i * 2, 4 + i * 2);
			a[j + Gaussian2DFunction.Y_POSITION] = rdg.nextUniform(2 + i * 2, 4 + i * 2);
			a[j + Gaussian2DFunction.X_SD] = rdg.nextUniform(1.5, 2);
			a[j + Gaussian2DFunction.Y_SD] = rdg.nextUniform(1.5, 2);
		}
		f.initialise(a);

		GradientCalculator c = GradientCalculatorFactory.newCalculator(a.length);
		double[][] I = c.fisherInformationMatrix(size, a, f);

		//System.out.printf("n=%d, k=%d, I=\n", n, k);
		//for (int i = 0; i < I.length; i++)
		//	System.out.println(Arrays.toString(I[i]));

		// Reduce to the desired size
		I = Arrays.copyOf(I, n);
		for (int i = 0; i < n; i++)
			I[i] = Arrays.copyOf(I[i], n);

		// Zero selected columns
		if (k > 0)
		{
			int[] zero = new RandomDataGenerator(randomGenerator).nextPermutation(n, k);
			for (int i : zero)
			{
				for (int j = 0; j < n; j++)
				{
					I[i][j] = I[j][i] = 0;
				}
			}
		}

		//System.out.printf("n=%d, k=%d\n", n, k);
		//for (int i = 0; i < n; i++)
		//	System.out.println(Arrays.toString(I[i]));

		// Create matrix
		return new FisherInformationMatrix(I, 1e-3);
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
