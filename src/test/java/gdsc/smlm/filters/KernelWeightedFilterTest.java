package gdsc.smlm.filters;

public class KernelWeightedFilterTest extends WeightedFilterTest
{
	@Override
	DataFilter createDataFilter()
	{
		// Do not support non-integer box sizes
		return new DataFilter("kernel", false, 3)
		{
			float[] w;
			int width, height;

			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				int k = (int) boxSize;
				// Only do odd box sizes
				if ((k & 1) != 1)
					return;
				KernelFilter f = createKernelFilter(k);
				f.convolve(data, width, height);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				int k = (int) boxSize;
				// Only do odd box sizes
				if ((k & 1) != 1)
					return;
				KernelFilter f = createKernelFilter(k);
				f.convolve(data, width, height, k / 2);
			}

			private KernelFilter createKernelFilter(int k)
			{
				KernelFilter f = new KernelFilter(KernelFilterTest.createKernel(k, k), k, k);
				if (w != null)
					f.setWeights(w, width, height);
				return f;
			}

			@Override
			public void setWeights(float[] w, int width, int height)
			{
				this.w = w;
				this.width = width;
				this.height = height;
			}
		};
	}
}
