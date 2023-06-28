/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.engine;

import java.awt.Rectangle;
import java.util.Objects;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration.PeakResultValidationData;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.utils.ImageConverter;

/**
 * Returns a simple estimate of local background and noise for fitting a Gaussian 2D function. The
 * estimate is valid as long as the centre is within the data range.
 */
public class SimplePeakResultValidationData implements PeakResultValidationData {
  private final ImageConverter ic = new ImageConverter();
  private final GaussianFunctionFactory factory;
  private final Object data;
  private final int ox;
  private final int oy;
  private final int maxx;
  private final int maxy;

  private int peakNumber;
  private double[] params;
  private double localBackground;
  private double noise = -1;

  /**
   * Instantiates a new simple peak result validation data using function (for 1 peak) that is being
   * used to fit the data (may be multiple peaks).
   *
   * <p>Input data is any array type supported by the ImageConverter class.
   *
   * @param factory the function factory
   * @param ox the x-origin of the fit region within the data
   * @param oy the y-origin of the fit region within the data
   * @param data the data
   * @param maxx the maxx of the data
   * @param maxy the maxy of the data
   */
  public SimplePeakResultValidationData(GaussianFunctionFactory factory, int ox, int oy,
      Object data, int maxx, int maxy) {
    Objects.requireNonNull(factory, "Factory is null");
    SimpleArrayUtils.check2DSize(maxx, maxy);
    if (!ic.isSupported(data)) {
      throw new IllegalArgumentException("Data is not supported");
    }
    this.factory = factory;
    this.ox = ox;
    this.oy = oy;
    this.data = data;
    this.maxx = maxx;
    this.maxy = maxy;
  }

  @Override
  public void setResult(int peakNumber, double[] initialParams, double[] params,
      double[] paramDevs) {
    noise = -1;
    this.peakNumber = peakNumber;
    this.params = params;
  }

  @Override
  public double getLocalBackground() {
    compute();
    return localBackground;
  }

  @Override
  public double getNoise() {
    compute();
    return noise;
  }

  private void compute() {
    if (noise != -1) {
      return;
    }

    final double[] spotParams = extractSpotParams(params, peakNumber);

    // Adjust to the data frame
    spotParams[Gaussian2DFunction.X_POSITION] += ox;
    spotParams[Gaussian2DFunction.Y_POSITION] += oy;

    // Add the 0.5 pixel offset to get the centre pixel
    final int x = (int) (spotParams[Gaussian2DFunction.X_POSITION] + 0.5);
    final int y = (int) (spotParams[Gaussian2DFunction.Y_POSITION] + 0.5);
    // Do not evaluate over a large region for speed.
    // Use only 50% of the Gaussian volume or 3 pixels.
    int nx = getRange(spotParams[Gaussian2DFunction.X_SD] * Gaussian2DPeakResultHelper.R_2D_50, 3);
    int ny = getRange(spotParams[Gaussian2DFunction.Y_SD] * Gaussian2DPeakResultHelper.R_2D_50, 3);

    final Rectangle r1 = new Rectangle(x - nx, y - ny, 2 * nx + 1, 2 * ny + 1);
    final Rectangle r2 = r1.intersection(new Rectangle(0, 0, maxx, maxy));

    if (r2.width * r2.height <= 1) {
      localBackground = params[Gaussian2DFunction.BACKGROUND];
      noise = Math.sqrt(localBackground); // Assume photon shot noise.
      return;
    }

    // Get the region of the data
    final double[] region = ic.getDoubleData(data, maxx, maxy, r2, null);

    // Compute the function in the same region.
    // Adjust the coordinates for clipping (r2 will be >= r1).
    // This effectively makes nx/ny represent the number of pixels before the centre pixel
    nx -= r2.x - r1.x;
    ny -= r2.y - r1.y;
    // Put the spot in the centre of the region
    spotParams[Gaussian2DFunction.X_POSITION] += nx - x;
    spotParams[Gaussian2DFunction.Y_POSITION] += ny - y;
    final Gaussian2DFunction f = factory.create2D(1, r2.width, r2.height);
    final double[] v = f.computeValues(spotParams);

    final Statistics stats = new Statistics();
    for (int i = 0; i < v.length; i++) {
      stats.add(region[i] - v[i]);
    }

    localBackground = stats.getMean();
    noise = stats.getStandardDeviation();
  }

  /**
   * Extract parameters for the specified peak. The background is ignored.
   *
   * @param params the params
   * @param n the peak
   * @return the extracted params
   */
  private static double[] extractSpotParams(double[] params, int n) {
    final double[] newParams = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    System.arraycopy(params, n * Gaussian2DFunction.PARAMETERS_PER_PEAK + 1, newParams, 1,
        Gaussian2DFunction.PARAMETERS_PER_PEAK);
    return newParams;
  }

  /**
   * Gets the range over which to evaluate a Gaussian using a factor of the standard deviation.
   *
   * <p>The range is clipped to 1 to max.
   *
   * @param range the range factor
   * @param max the max value to return
   * @return the range
   */
  private static int getRange(double range, int max) {
    final double l = Math.ceil(range);
    if (l < 1) {
      return 1;
    }
    if (l >= max) {
      return max;
    }
    return (int) l;
  }
}
