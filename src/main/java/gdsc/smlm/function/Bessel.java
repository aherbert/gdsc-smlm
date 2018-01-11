package gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;

/**
 * Class for computing various Bessel functions
 * <p>
 * The implementation is based upon that presented in: Numerical Recipes in C++, The Art of Scientific Computing, Second
 * Edition, W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery (Cambridge University Press, Cambridge, 2002).
 */
public class Bessel
{
	/**
	 * Compute the zero th order Bessel function of the first kind.
	 * 
	 * @param x
	 * 
	 * @return the Bessel function J0
	 */
	public static double J0(final double x)
	{
		double ax;

		if ((ax = Math.abs(x)) < 8.0)
		{
			final double y = x * x;
			final double ans1 = 57568490574.0 + y *
					(-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
			final double ans2 = 57568490411.0 + y *
					(1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y * 1.0))));

			return ans1 / ans2;
		}

		final double z = 8.0 / ax;
		final double y = z * z;
		final double xx = ax - 0.785398164;
		final double ans1 = 1.0 + y *
				(-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
		final double ans2 = -0.1562499995e-1 + y *
				(0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)));

		return Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
	}

	/**
	 * Compute the first order Bessel function of the first kind.
	 * 
	 * @param x
	 * 
	 * @return the Bessel function J1
	 */
	public static double J1(final double x)
	{
		double ax;

		if ((ax = Math.abs(x)) < 8.0)
		{
			final double y = x * x;
			final double ans1 = x *
					(72362614232.0 + y *
							(-7895059235.0 + y *
									(242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
			final double ans2 = 144725228442.0 + y *
					(2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
			return ans1 / ans2;
		}

		final double z = 8.0 / ax;
		final double xx = ax - 2.356194491;
		final double y = z * z;

		final double ans1 = 1.0 + y *
				(0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
		final double ans2 = 0.04687499995 + y *
				(-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
		final double ans = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
		return (x < 0.0) ? -ans : ans;
	}

	/**
	 * Compute the second order Bessel function of the first kind.
	 * <p>
	 * This is stable when abs(x) > n.
	 * 
	 * @param x
	 * 
	 * @return the Bessel function J2
	 */
	public static double J2(double x)
	{
		if (x == 0.0)
		{
			return 0.0;
		}
		double value0 = J0(x);
		double value1 = J1(x);
		return 2.0 * value1 / x + value0;
	}

	/**
	 * Compute the modified Bessel function of the first kind.
	 * 
	 * @param x
	 * 
	 * @return the modified Bessel function I1
	 */
	public static double I1(final double x)
	{
		double ax, ans, y;

		if ((ax = Math.abs(x)) < 3.75)
		{
			y = x / 3.75;
			y *= y;
			ans = ax *
					(0.5 + y *
							(0.87890594 + y *
									(0.51498869 + y *
											(0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
		}
		else
		{
			y = 3.75 / ax;
			ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
			ans = 0.39894228 + y *
					(-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
			ans *= (FastMath.exp(ax) / Math.sqrt(ax));
		}
		return x < 0.0 ? -ans : ans;
	}
}
