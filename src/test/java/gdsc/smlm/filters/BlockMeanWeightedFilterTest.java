package gdsc.smlm.filters;

public class BlockMeanWeightedFilterTest extends WeightedFilterTest
{
	@Override
	DataFilter createDataFilter()
	{
		return new DataFilter("blockMean", true)
		{
			BlockMeanFilter f = new BlockMeanFilter();

			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, boxSize);
			}

			@Override
			public void setWeights(float[] w, int width, int height)
			{
				f.setWeights(w, width, height);
			}
		};
	}
}
