package gdsc.smlm.function.cspline;

public class MultiCubicSplineFunctionTest extends CubicSplineFunctionTest
{
	protected void init()
	{
		f1 = new MultiCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, scale);
		f1f = new MultiCubicSplineFunction(splineDataFloat, maxx, maxy, cx, cy, cz, scale);
		MultiCubicSplineFunction f = new MultiCubicSplineFunction(splineData, maxx, maxy, cx, cy, cz, scale);
		f.setN(2);
		f2 = f;
		f = new MultiCubicSplineFunction(splineDataFloat, maxx, maxy, cx, cy, cz, scale);
		f.setN(2);
		f2f = f;
	}
}
