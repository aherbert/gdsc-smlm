package gdsc.smlm.ij.utils;

import gdsc.core.utils.SimpleArrayUtils;

public class FloatImage2DTest extends Image2DTest
{
	protected FloatImage2D createData(int w, int h)
	{
		float[] data = SimpleArrayUtils.newArray(w * h, 1f, 1f);
		return new FloatImage2D(w, h, data);
	}
	
	protected FloatImage2D createEmptyData(int w, int h)
	{
		return new FloatImage2D(w, h);
	}
}
