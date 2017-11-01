package gdsc.smlm.function.cspline;

public class SingleCubicSplineFunctionTest extends CubicSplineFunctionTest
{
	protected void init()
	{
		f1 = new SingleCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, scale);
		f1f = new SingleCubicSplineFunction(splineDataFloat, maxx, maxy, cx, cy, cz, scale);
	}
}
