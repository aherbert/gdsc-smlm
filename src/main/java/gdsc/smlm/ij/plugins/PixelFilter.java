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
package gdsc.smlm.ij.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.util.Tools;

import java.awt.AWTEvent;
import java.awt.Label;

import org.apache.commons.math3.util.FastMath;

/**
 * Filters pixels using the surrounding region.
 */
public class PixelFilter implements ExtendedPlugInFilter, DialogListener
{
	private static final String TITLE = "Pixel Filter";
	private final int FLAGS = DOES_8G | DOES_16 | DOES_32 | PARALLELIZE_STACKS;

	private static int radius = 1;
	private static double error = 3;

	private PlugInFilterRunner pfr = null;
	private double[] cachedS = null;
	private double[] cachedSS = null;
	private boolean preview = false;
	private Label label = null;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}
		return FLAGS;
	}

	@Override
	public void run(ImageProcessor ip)
	{
		// Compute rolling sums
		FloatProcessor fp = ip.toFloat(0, null);
		float[] data = (float[]) ip.toFloat(0, null).getPixels();
		double[] s = null, ss = null;
		if (preview && cachedS != null)
		{
			s = cachedS;
			ss = cachedSS;
		}
		if (s == null)
		{
			s = new double[ip.getPixelCount()];
			ss = new double[s.length];
			calculateRollingSums(fp, s, ss);
		}

		int count = 0;
		final int maxx = ip.getWidth();
		final int maxy = ip.getHeight();
		for (int y = 0, i = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, i++)
			{
				double sum = 0;
				double sumSquares = 0;

				int minU = x - radius - 1;
				int maxU = FastMath.min(x + radius, maxx - 1);
				int minV = y - radius - 1;
				int maxV = FastMath.min(y + radius, maxy - 1);

				// Compute sum from rolling sum using:
				// sum(u,v) = 
				// + s(u+N,v+N) 
				// - s(u-N-1,v+N)
				// - s(u+N,v-N-1)
				// + s(u-N-1,v-N-1)
				// Note: 
				// s(u,v) = 0 when either u,v < 0
				// s(u,v) = s(umax,v) when u>umax
				// s(u,v) = s(u,vmax) when v>vmax
				// s(u,v) = s(umax,vmax) when u>umax,v>vmax
				// Likewise for ss

				// + s(u+N-1,v+N-1) 
				int index = maxV * maxx + maxU;
				sum += s[index];
				sumSquares += ss[index];

				if (minU >= 0)
				{
					// - s(u-1,v+N-1)
					index = maxV * maxx + minU;
					sum -= s[index];
					sumSquares -= ss[index];
				}
				if (minV >= 0)
				{
					// - s(u+N-1,v-1)
					index = minV * maxx + maxU;
					sum -= s[index];
					sumSquares -= ss[index];

					if (minU >= 0)
					{
						// + s(u-1,v-1)
						index = minV * maxx + minU;
						sum += s[index];
						sumSquares += ss[index];
					}
				}

				// Reset to bounds to calculate the number of pixels
				if (minU < 0)
					minU = -1;
				if (minV < 0)
					minV = -1;

				int n = (maxU - minU) * (maxV - minV);

				if (n < 2)
					continue;

				// Get the sum of squared differences
				double residuals = sumSquares - (sum * sum) / n;

				//// -----------------------------
				//// Check using the original data
				//// -----------------------------
				//double sx = 0;
				//double ssx = 0;
				//int nn = 0;
				//for (int yy = y - radius; yy <= y + radius; yy++)
				//	for (int xx = x - radius; xx <= x + radius; xx++)
				//	{
				//		if (xx >= 0 && xx < maxx && yy >= 0 && yy < maxy)
				//		{
				//			float value = fp.getf(xx, yy);
				//			sx += value;
				//			ssx += value * value;
				//			nn++;
				//		}
				//	}
				//DoubleEquality eq = new DoubleEquality(8, 1e-16);
				//if (n != nn)
				//{
				//	System.out.printf("Wrong @ %d,%d %d <> %d\n", x, y, n, nn);
				//	residuals = ssx - sx * sx / nn;
				//}
				//else if (!eq.almostEqualComplement(sx, sum) || !eq.almostEqualComplement(ssx, sumSquares))
				//{
				//	System.out.printf("Wrong @ %d,%d %g <> %g : %g <> %g\n", x, y, sx, sum, ssx, sumSquares);
				//	residuals = ssx - sx * sx / nn;
				//}
				//// -----------------------------

				if (residuals > 0.0)
				{
					double stdDev = Math.sqrt(residuals / (n - 1.0));
					double mean = sum / n;

					if (Math.abs(data[i] - mean) / stdDev > error)
					{
						ip.setf(i, (float) mean);
						count++;
					}
				}
			}
		}
		if (preview)
		{
			cachedS = s;
			cachedSS = ss;
			label.setText("Replaced " + count);
		}
		else if (pfr != null && count > 0)
			IJ.log(String.format("Slice %d : Replaced %d pixels", pfr.getSliceNumber(), count));
	}

	private void calculateRollingSums(FloatProcessor ip, double[] s_, double[] ss)
	{
		// Compute the rolling sum and sum of squares
		// s(u,v) = f(u,v) + s(u-1,v) + s(u,v-1) - s(u-1,v-1) 
		// ss(u,v) = f(u,v) * f(u,v) + ss(u-1,v) + ss(u,v-1) - ss(u-1,v-1)
		// where s(u,v) = ss(u,v) = 0 when either u,v < 0

		int maxx = ip.getWidth();
		int maxy = ip.getHeight();
		float[] originalData = (float[]) ip.getPixels();
		double[] data = Tools.toDouble(originalData);

		// First row
		double cs_ = 0; // Column sum
		double css = 0; // Column sum-squares
		for (int i = 0; i < maxx; i++)
		{
			cs_ += data[i];
			css += data[i] * data[i];
			s_[i] = cs_;
			ss[i] = css;
		}

		// Remaining rows:
		// sum = rolling sum of row + sum of row above
		for (int y = 1; y < maxy; y++)
		{
			int i = y * maxx;
			cs_ = 0;
			css = 0;

			// Remaining columns
			for (int x = 0; x < maxx; x++, i++)
			{
				cs_ += data[i];
				css += data[i] * data[i];

				s_[i] = s_[i - maxx] + cs_;
				ss[i] = ss[i - maxx] + css;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#showDialog(ij.ImagePlus, java.lang.String,
	 * ij.plugin.filter.PlugInFilterRunner)
	 */
	@Override
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		this.pfr = pfr;
		preview = true;

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Replace pixels with mean if they are N StdDevs from the mean");

		gd.addSlider("Radius", 1, 5, radius);
		gd.addSlider("Error (SD units)", 2.5, 7, error);

		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);

		gd.addMessage("");
		label = (Label) gd.getMessage();

		gd.showDialog();

		if (gd.wasCanceled() || !dialogItemChanged(gd, null))
			return DONE;

		preview = false;
		cachedS = cachedSS = null;
		label = null;
		return IJ.setupDialog(imp, FLAGS);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		label.setText("");
		radius = (int) gd.getNextNumber();
		error = gd.getNextNumber();
		if (gd.invalidNumber() || radius < 1 || error < 0)
			return false;
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	@Override
	public void setNPasses(int nPasses)
	{
		// Ignore		
	}
}
