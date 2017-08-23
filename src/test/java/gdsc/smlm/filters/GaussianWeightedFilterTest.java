package gdsc.smlm.filters;

public class GaussianWeightedFilterTest extends WeightedFilterTest
{
	@Override
	DataFilter createDataFilter()
	{
		return new DataFilter("gaussian", true)
		{
			GaussianFilter f = new GaussianFilter();

			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.convolve(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.convolveInternal(data, width, height, boxSize);
			}

			@Override
			public void setWeights(float[] w, int width, int height)
			{
				f.setWeights(w, width, height);
			}
		};
	}
}
