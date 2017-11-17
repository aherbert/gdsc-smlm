package gdsc.smlm.ij.utils;

import gdsc.core.utils.SimpleArrayUtils;

public class FloatImage3DTest extends Image3DTest
{
	protected FloatImage3D createData(int w, int h, int d)
	{
		float[] data = SimpleArrayUtils.newArray(w * h * d, 1f, 1f);
		return new FloatImage3D(w, h, d, data);
	}

	protected FloatImage3D createEmptyData(int w, int h, int d)
	{
		return new FloatImage3D(w, h, d);
	}
}
