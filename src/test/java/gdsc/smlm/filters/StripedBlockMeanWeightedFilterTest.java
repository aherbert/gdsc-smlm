package gdsc.smlm.filters;

public class StripedBlockMeanWeightedFilterTest extends WeightedFilterTest
{
	@Override
	DataFilter createDataFilter()
	{
		return new DataFilter("stripedBlockMean", true)
		{
			BlockMeanFilter f = new BlockMeanFilter();

			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterInternal(data, width, height, boxSize);
			}

			@Override
			public void setWeights(float[] w, int width, int height)
			{
				f.setWeights(w, width, height);
			}
		};
	}
}
