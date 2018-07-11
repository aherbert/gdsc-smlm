/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.filters;

@SuppressWarnings({ "javadoc" })
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
