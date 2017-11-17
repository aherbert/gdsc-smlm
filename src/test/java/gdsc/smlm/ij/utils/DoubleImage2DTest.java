package gdsc.smlm.ij.utils;

import gdsc.core.utils.SimpleArrayUtils;

public class DoubleImage2DTest extends Image2DTest
{
	protected DoubleImage2D createData(int w, int h)
	{
		double[] data = SimpleArrayUtils.newArray(w * h, 1.0, 1.0);
		return new DoubleImage2D(w, h, data);
	}
	
	protected DoubleImage2D createEmptyData(int w, int h)
	{
		return new DoubleImage2D(w, h);
	}
}
