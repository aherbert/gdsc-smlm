package gdsc.smlm.function.gaussian;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.DoubleEquality;

public class HoltzerAstigmatismZModelTest
{
	protected DoubleEquality eq = new DoubleEquality(1e-5, 1e-7);

	// Compute as per Numerical Recipes 5.7.
	// Approximate error accuracy in single precision: Ef
	// Step size for derivatives:
	// h ~ (Ef)^(1/3) * xc
	// xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x is close to zero)
	protected double h_ = 0.0001; //(double) (Math.pow(1e-3f, 1.0 / 3));

	@Test
	public void canStaticComputeGradient()
	{
		double d = 0.531;
		double Ax = -0.0708;
		double Bx = -0.073;
		double Ay = 0.164;
		double By = 0.0417;
		
		canStaticComputeGradient(d, Ax, Bx);
		canStaticComputeGradient(d, Ay, By);
	}
	
	private void canStaticComputeGradient(double d, double Ax, double Bx)
	{
		double one_d2 = 1.0 / d / d;

		double[] ds_dz = new double[1];
		double[] ds_dz2 = new double[2];
		double[] ds_duz = new double[1];
		double[] ds_dlz = new double[1];
		for (double z = -0.5; z < 0.5; z += 0.01)
		{
			double s0 = HoltzerAstimatismZModel.getS(z, one_d2, Ax, Bx);
			double s1 = HoltzerAstimatismZModel.getS1(z, one_d2, Ax, Bx, ds_dz);
			double s2 = HoltzerAstimatismZModel.getS2(z, one_d2, Ax, Bx, ds_dz2);

			Assert.assertEquals(s0, s1, 0);
			Assert.assertEquals(s0, s2, 0);
			Assert.assertEquals(ds_dz[0], ds_dz2[0], 0);

			double uz = z + h_;
			double lz = z - h_;
			double upper = HoltzerAstimatismZModel.getS1(uz, one_d2, Ax, Bx, ds_duz);
			double lower = HoltzerAstimatismZModel.getS1(lz, one_d2, Ax, Bx, ds_dlz);

			double e1 = (upper - lower) / (uz - lz);
			double o1 = ds_dz[0];

			// Second gradient
			double e2 = (ds_duz[0] - ds_dlz[0]) / (uz - lz);
			double o2 = ds_dz2[1];

			System.out.printf("z=%f s=%f : ds_dz=%g  %g  (%g): d2s_dz2=%g   %g  (%g)\n", z, s0, e1, o1,
					DoubleEquality.relativeError(o1, e1), e2, o2, DoubleEquality.relativeError(o2, e2));

			//double error = DoubleEquality.relativeError(o, e);
			if (Math.abs(z) > 0.02)
				Assert.assertTrue(e1 + " sign != " + o1, (e1 * o1) >= 0);
			Assert.assertTrue(e1 + " != " + o1, eq.almostEqualRelativeOrAbsolute(e1, o1));

			if (Math.abs(z) > 0.02)
				Assert.assertTrue(e2 + " sign != " + o2, (e2 * o2) >= 0);
			Assert.assertTrue(e2 + " != " + o2, eq.almostEqualRelativeOrAbsolute(e2, o2));
		}
	}
}
