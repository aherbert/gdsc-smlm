package gdsc.smlm.filters;

public class RollingBlockMeanWeightedFilterTest extends WeightedFilterTest
{
	@Override
	DataFilter createDataFilter()
	{
		return new DataFilter("rollingBlockMean", false)
		{
			BlockMeanFilter f = new BlockMeanFilter();

			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilterInternal(data, width, height, (int) boxSize);
			}

			@Override
			public void setWeights(float[] w, int width, int height)
			{
				f.setWeights(w, width, height);
			}
		};
	}
}
