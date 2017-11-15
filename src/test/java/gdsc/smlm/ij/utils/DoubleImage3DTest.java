package gdsc.smlm.ij.utils;

import gdsc.core.utils.SimpleArrayUtils;

public class DoubleImage3DTest extends Image3DTest 
{
	protected DoubleImage3D createData(int w, int h, int d)
	{
		double[] data = SimpleArrayUtils.newArray(w * h * d, 1.0, 1.0);
		return new DoubleImage3D(w, h, d, data);
	}
}
