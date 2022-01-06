/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.filters;

@SuppressWarnings({"javadoc"})
class KernelWeightedFilterTest extends WeightedKernelFilterTest {
  @Override
  DataFilter createDataFilter() {
    // Do not support non-integer box sizes
    return new DataFilter("kernel", false, 3) {
      float[] weights;
      int width;
      int height;
      KernelFilter filter;
      int kernelWidth = 0;

      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        final int k = (int) boxSize;
        // Only do odd box sizes
        if ((k & 1) != 1) {
          return;
        }
        updateKernelFilter(k);
        filter.convolve(data, width, height);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        final int kernelWidth = (int) boxSize;
        // Only do odd box sizes
        if ((kernelWidth & 1) != 1) {
          return;
        }
        updateKernelFilter(kernelWidth);
        filter.convolve(data, width, height, kernelWidth / 2);
      }

      private void updateKernelFilter(int kernelWidth) {
        if (this.kernelWidth != kernelWidth || filter == null) {
          filter = createKernelFilter(kernelWidth);
        }
        if (weights != null) {
          filter.setWeights(weights, this.width, this.height);
        }
      }

      private KernelFilter createKernelFilter(int kernelWidth) {
        return new KernelFilter(KernelFilterTest.createKernel(kernelWidth, kernelWidth),
            kernelWidth, kernelWidth);
      }

      @Override
      public void setWeights(float[] weights, int width, int height) {
        this.weights = weights;
        this.width = width;
        this.height = height;
      }
    };
  }
}
