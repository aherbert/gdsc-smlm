package gdsc.smlm.fitting;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;

public class FisherInformationMatrixTest
{
	@Test
	public void canComputeCRLB()
	{
		for (int n = 1; n < 10; n++)
		{
			canComputeCRLB(n, 0);
		}
	}

	@Test
	public void canComputeCRLBWithZeros()
	{
		canComputeCRLB(4, 1);

		for (int n = 2; n < 10; n++)
		{
			canComputeCRLB(n, 1);
			canComputeCRLB(n, n / 2);
		}
	}

	private void canComputeCRLB(int n, int k)
	{
		RandomGenerator randomGenerator = new Well19937c(30051977);

		// TODO - Use a real Gaussian function here to compute the Fisher information.
		// The matrix may be sensitive to the type of equation used.
		
		double[] dx = new double[n];
		for (int i = 0; i < n; i++)
			dx[i] = randomGenerator.nextDouble() * 100;
		double[][] I = new double[n][n];
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				I[i][j] = dx[i] * dx[j];
			}
		}

		// Symmetric
		for (int i = 0; i < n - 1; i++)
			for (int j = i + 1; j < n; j++)
				I[i][j] = I[j][i];

		// Zero
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

		FisherInformationMatrix m = new FisherInformationMatrix(I);
		DoubleEquality eq = new DoubleEquality(3, 1e-6);
		m.setEqual(eq);

		double[] crlb = m.crlb(true);
		System.out.printf("n=%d, k=%d : %s\n", n, k, Arrays.toString(crlb));
		Assert.assertNotNull("CRLB failed", crlb);
	}

	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
