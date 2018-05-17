package gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;

public class InterpolatedPoissonFisherInformationTest
{
	@Test
	public void canInterpolateFisherInformation()
	{
		for (int i = 0; i < 4; i++)
		{
			double s = (1 << i) * 0.25;
			canInterpolateFisherInformation(s);
		}
	}

	private void canInterpolateFisherInformation(double s)
	{
		canComputeFisherInformation(new PoissonGaussianApproximationFisherInformation(s));
	}

	private void canComputeFisherInformation(BasePoissonFisherInformation f)
	{
		// Build a range for the Fisher information
		int min = -100;
		int max = 20;
		double[] logU = new double[max - min + 1];
		double[] alpha = new double[logU.length];

		for (int exp = min, i = 0; exp <= max; exp++, i++)
		{
			logU[i] = exp;
			alpha[i] = f.getAlpha(FastMath.exp(exp));
		}

		InterpolatedPoissonFisherInformation fi = new InterpolatedPoissonFisherInformation(logU, alpha, f);

		// Lower bound
		check(f, fi, min - 1, 1e-10);
		check(f, fi, min, 1e-2);

		// Within
		for (int exp = min; exp < max; exp++)
		{
			check(f, fi, exp + 0.5, 1e-2);
		}

		// Upper bound
		check(f, fi, max, 1e-2);
		check(f, fi, max + 1, 1e-2);
	}

	private void check(BasePoissonFisherInformation f, InterpolatedPoissonFisherInformation fi, double logU, double tol)
	{
		double u = FastMath.exp(logU);
		double e = f.getAlpha(u);
		double o = fi.getAlpha(u);
		//System.out.printf("logU=%g  u=%g  e=%g  o=%g  error=%g\n", logU, u, e, o, DoubleEquality.relativeError(o, e));

		// Small numbers may have a large relative error but the absolute error is small
		Assert.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e, o, 5e-3, 1e-20));
	}
}