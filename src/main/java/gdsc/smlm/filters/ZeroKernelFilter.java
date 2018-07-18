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

/**
 * Computes a convolution in the spatial domain for each point within the array. Pixels outside the array are assumed to
 * be zero.
 * <p>
 * Adapted from {@link ij.plugin.filter.Convolver}
 */
public class ZeroKernelFilter extends KernelFilter
{
	/**
	 * Instantiates a new kernel filter.
	 *
	 * @param kernel
	 *            the kernel
	 * @param kw
	 *            the kernel width (must be odd)
	 * @param kh
	 *            the kernel height (must be odd)
	 */
	public ZeroKernelFilter(float[] kernel, int kw, int kh)
	{
		super(kernel, kw, kh);
	}

	@Override
	protected void convolveData(float[] in, float[] out, final int width, final int height, int border)
	{
		final int x1 = border;
		final int y1 = border;
		final int x2 = width - border;
		final int y2 = height - border;
		final int uc = kw / 2;
		final int vc = kh / 2;
		final int xedge = width - uc;
		final int yedge = height - vc;
		for (int y = y1; y < y2; y++)
		{
			final boolean edgeY = y < vc || y >= yedge;
			for (int x = x1, c = x1 + y * width; x < x2; x++)
			{
				double sum = 0.0;
				int i = 0;
				// Determine if at the edge
				if (edgeY || x < uc || x >= xedge)
					for (int v = -vc; v <= vc; v++)
					{
						// Create a safe y-index
						int yIndex = y + v;
						if (yIndex < 0 || yIndex >= height)
						{
							// Nothing to convolve so skip forward
							i += kw;
							continue;
						}
						yIndex *= width;

						for (int u = -uc; u <= uc; u++)
							//if (i >= kernel.length) // work around for JIT compiler bug on Linux
							//	IJ.log("kernel index error: " + i);
							sum += getPixel(x + u, yIndex, in, width) * kernel[i++];
					}
				else
					// Internal
					for (int v = -vc; v <= vc; v++)
						for (int u = -uc, offset = x - uc + (y + v) * width; u++ <= uc;)
							sum += in[offset++] * kernel[i++];
				out[c++] = (float) (sum * scale);
			}
		}
	}

	/**
	 * Gets the pixel respecting the image boundaries.
	 *
	 * @param x
	 *            the x
	 * @param yIndex
	 *            the y index in the 2D array
	 * @param pixels
	 *            the pixels
	 * @param width
	 *            the width
	 * @return the pixel
	 */
	private static float getPixel(int x, int yIndex, float[] pixels, int width)
	{
		if (x < 0 || x >= width)
			return 0f;
		return pixels[x + yIndex];
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#clone()
	 */
	@Override
	public ZeroKernelFilter clone()
	{
		final ZeroKernelFilter o = (ZeroKernelFilter) super.clone();
		return o;
	}
}
