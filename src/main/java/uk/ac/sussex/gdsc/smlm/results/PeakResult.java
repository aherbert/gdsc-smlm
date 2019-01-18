/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.results;

import uk.ac.sussex.gdsc.core.utils.MathUtils;

import org.apache.commons.lang3.ArrayUtils;

import java.util.Arrays;

/**
 * Specifies a peak fitting result.
 */
public class PeakResult {
  /** Index of the background in the parameters array. */
  public static final int BACKGROUND = 0;
  /** Index of the intensity in the parameters array. */
  public static final int INTENSITY = 1;
  /** Index of the x-position in the parameters array. */
  public static final int X = 2;
  /** Index of the y-position in the parameters array. */
  public static final int Y = 3;
  /** Index of the z-position in the parameters array. */
  public static final int Z = 4;
  /** Number of standard parameters. */
  public static final int STANDARD_PARAMETERS = 5;

  private static final String[] NAMES = {"Background", "Intensity", "X", "Y", "Z"};

  private int frame;
  private int origX;
  private int origY;
  private float origValue;
  private double error;
  private float noise;
  private float meanIntensity;

  /**
   * The parameters (for the standard parameters plus any PSF specific parameters). This is never
   * null.
   */
  private float[] params;
  /**
   * The parameter standard deviations (for the standard parameters plus any PSF specific
   * parameters). This may be null or the same length as {@link #params}.
   */
  private float[] paramStdDevs;

  /**
   * Gets the parameter name.
   *
   * <p>Valid for the standard parameters ({@code index <} {@link #STANDARD_PARAMETERS}).
   *
   * @param index the index
   * @return the parameter name
   */
  public static String getParameterName(int index) {
    return NAMES[index];
  }

  /**
   * Instantiates a new peak result.
   *
   * @param frame the frame
   * @param origX the original X position
   * @param origY the original Y position
   * @param origValue the original value
   * @param error the error
   * @param noise the noise
   * @param meanIntensity the mean intensity
   * @param params the params (must not be null and must have at least {@value #STANDARD_PARAMETERS}
   *        parameters)
   * @param paramsStdDev the params standard deviations (if not null must match the length of the
   *        params array)
   * @throws IllegalArgumentException the illegal argument exception if the parameters are invalid
   */
  public PeakResult(int frame, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    if (params == null) {
      throw new IllegalArgumentException("Parameters must not be null");
    }
    if (params.length < STANDARD_PARAMETERS) {
      throw new IllegalArgumentException("Parameters must contain all standard parameters");
    }
    if (paramsStdDev != null && paramsStdDev.length != params.length) {
      throw new IllegalArgumentException("Parameter deviations length must match parameters");
    }
    this.frame = frame;
    this.origX = origX;
    this.origY = origY;
    this.origValue = origValue;
    this.error = error;
    this.noise = noise;
    this.meanIntensity = meanIntensity;
    this.params = params;
    this.paramStdDevs = paramsStdDev;
  }

  /**
   * Simple constructor to create a result with frame, location, and intensity.
   *
   * @param frame the frame
   * @param x the x position
   * @param y the y position
   * @param intensity the intensity
   */
  public PeakResult(int frame, float x, float y, float intensity) {
    setFrame(frame);
    origX = (int) x;
    origY = (int) y;
    params = new float[5];
    params[X] = x;
    params[Y] = y;
    params[INTENSITY] = intensity;
  }

  /**
   * Simple constructor to create a result with location and intensity.
   *
   * @param x the x position
   * @param y the y position
   * @param intensity the intensity
   */
  public PeakResult(float x, float y, float intensity) {
    this(0, x, y, intensity);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected PeakResult(PeakResult source) {
    frame = source.frame;
    origX = source.origX;
    origY = source.origY;
    origValue = source.origValue;
    error = source.error;
    noise = source.noise;
    meanIntensity = source.meanIntensity;
    params = source.params.clone();
    paramStdDevs = ArrayUtils.clone(source.paramStdDevs);
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public PeakResult copy() {
    return new PeakResult(this);
  }

  /**
   * Creates the params array for a peak result.
   *
   * @param background the background
   * @param intensity the intensity
   * @param x the x
   * @param y the y
   * @param z the z
   * @return the params array
   */
  public static float[] createParams(float background, float intensity, float x, float y, float z) {
    return new float[] {background, intensity, x, y, z};
  }

  /**
   * Utility function to check for equality.
   *
   * @param r1 the first result
   * @param r2 the second result
   * @return true, if equal
   */
  public static boolean equals(PeakResult r1, PeakResult r2) {
    if (r1 == r2) {
      // The same or both null
      return true;
    }
    if (r1 == null || r2 == null) {
      // At least one is not null so this is not equal
      return false;
    }

    // Check parameters
    if (r1.getNumberOfParameters() != r2.getNumberOfParameters()) {
      return false;
    }
    for (int i = 0; i < r1.params.length; i++) {
      if (r1.params[i] != r2.params[i]) {
        return false;
      }
    }

    // Check primitive fields
    if (r1.frame != r2.frame) {
      return false;
    }
    if (r1.origX != r2.origX) {
      return false;
    }
    if (r1.origY != r2.origY) {
      return false;
    }
    if (r1.origValue != r2.origValue) {
      return false;
    }
    if (r1.noise != r2.noise) {
      return false;
    }
    if (r1.meanIntensity != r2.meanIntensity) {
      return false;
    }

    // Check optional properties
    if (r1.hasId()) {
      if (!r2.hasId() || r1.getId() != r2.getId()) {
        return false;
      }
    } else if (r2.hasId()) {
      return false;
    }
    if (r1.hasEndFrame()) {
      if (!r2.hasEndFrame() || r1.getEndFrame() != r2.getEndFrame()) {
        return false;
      }
    } else if (r2.hasEndFrame()) {
      return false;
    }
    if (r1.hasPrecision()) {
      if (!r2.hasPrecision() || r1.getPrecision() != r2.getPrecision()) {
        return false;
      }
    } else if (r2.hasPrecision()) {
      return false;
    }

    // Check parameters deviations. Do this last as they are not often used.
    if (r1.paramStdDevs != null) {
      if (r2.paramStdDevs == null) {
        return false;
      }
      for (int i = 0; i < r1.paramStdDevs.length; i++) {
        if (r1.paramStdDevs[i] != r2.paramStdDevs[i]) {
          return false;
        }
      }
    } else if (r2.paramStdDevs != null) {
      return false;
    }

    return true;
  }

  /**
   * Gets the background for the first peak.
   *
   * @return The background for the first peak.
   */
  public float getBackground() {
    return params[BACKGROUND];
  }

  /**
   * Sets the background for the first peak.
   *
   * @param background the new background
   */
  public void setBackground(float background) {
    params[BACKGROUND] = background;
  }

  /**
   * Get the intensity for the first peak.
   *
   * @return The intensity of the first peak
   */
  public float getIntensity() {
    return params[INTENSITY];
  }

  /**
   * Sets the intensity for the first peak.
   *
   * @param intensity the new intensity
   */
  public void setIntensity(float intensity) {
    params[INTENSITY] = intensity;
  }

  /**
   * Gets the x position for the first peak.
   *
   * @return The x position for the first peak
   */
  public float getXPosition() {
    return params[X];
  }

  /**
   * Sets the x position for the first peak.
   *
   * @param x the new x position
   */
  public void setXPosition(float x) {
    params[X] = x;
  }

  /**
   * Gets the y position for the first peak.
   *
   * @return The y position for the first peak
   */
  public float getYPosition() {
    return params[Y];
  }

  /**
   * Sets the y position for the first peak.
   *
   * @param y the new y position
   */
  public void setYPosition(float y) {
    params[Y] = y;
  }

  /**
   * Gets the z position for the first peak.
   *
   * @return The z position for the first peak
   */
  public float getZPosition() {
    return params[Z];
  }

  /**
   * Sets the z position for the first peak.
   *
   * @param z the new z position
   */
  public void setZPosition(float z) {
    params[Z] = z;
  }

  /**
   * Gets the time frame.
   *
   * @return The time frame
   */
  public int getFrame() {
    return frame;
  }

  /**
   * Sets the time frame.
   *
   * @param frame The time frame
   */
  public void setFrame(int frame) {
    this.frame = frame;
  }

  /**
   * Gets the original X position in the fitted frame.
   *
   * @return the original X position
   */
  public int getOrigX() {
    return origX;
  }

  /**
   * Sets the original X position in the fitted frame.
   *
   * @param origX the original X position in the fitted frame.
   */
  public void setOrigX(int origX) {
    this.origX = origX;
  }

  /**
   * Gets the original Y position in the fitted frame.
   *
   * @return the original Y position in the fitted frame.
   */
  public int getOrigY() {
    return origY;
  }

  /**
   * Sets the original Y position in the fitted frame.
   *
   * @param origY the original Y position in the fitted frame.
   */
  public void setOrigY(int origY) {
    this.origY = origY;
  }

  /**
   * Gets the original value in the fitted frame.
   *
   * @return the original value in the fitted frame.
   */
  public float getOrigValue() {
    return origValue;
  }

  /**
   * Sets the original value in the fitted frame.
   *
   * @param origValue the original value in the fitted frame.
   */
  public void setOrigValue(float origValue) {
    this.origValue = origValue;
  }

  /**
   * Gets the error.
   *
   * @return the error
   */
  public double getError() {
    return error;
  }

  /**
   * Sets the error.
   *
   * @param error the new error
   */
  public void setError(double error) {
    this.error = error;
  }

  /**
   * Gets the noise.
   *
   * @return the noise
   */
  public float getNoise() {
    return noise;
  }

  /**
   * Sets the noise.
   *
   * @param noise the new noise
   */
  public void setNoise(float noise) {
    this.noise = noise;
  }

  /**
   * Checks for noise.
   *
   * @return true, if successful
   */
  public boolean hasNoise() {
    return noise > 0;
  }

  /**
   * Gets the mean intensity.
   *
   * <p>This requires a knowledge of the PSF used to create the result. It could be the peak signal
   * in the PSF or the average signal over a range of the PSF, e.g. the area covered from the maxima
   * to half-maxima for spots.
   *
   * @return the mean intensity
   */
  public float getMeanIntensity() {
    return meanIntensity;
  }

  /**
   * Sets the mean intensity.
   *
   * <p>This requires a knowledge of the PSF used to create the result. It could be the peak signal
   * in the PSF or the average signal over a range of the PSF, e.g. the area covered from the maxima
   * to half-maxima for spots.
   *
   * @param meanIntensity the new mean intensity
   */
  public void setMeanIntensity(float meanIntensity) {
    this.meanIntensity = meanIntensity;
  }

  /**
   * Checks for mean intensity.
   *
   * @return true, if successful
   */
  public boolean hasMeanIntensity() {
    return meanIntensity > 0;
  }

  /**
   * Gets the Signal-to-Noise Ratio (SNR). This is the ratio of the average signal value to the
   * background standard deviation.
   *
   * @return the snr
   * @see #getMeanIntensity()
   * @see #getNoise()
   */
  public float getSnr() {
    return MathUtils.div0(meanIntensity, noise);
  }

  /**
   * Checks for end frame. Derived classes can override this.
   *
   * @return true, if successful
   */
  public boolean hasEndFrame() {
    return false;
  }

  /**
   * Gets the end frame. Default = {@link #getFrame()}. Derived classes can override this.
   *
   * @return The last time frame that this result corresponds to
   */
  public int getEndFrame() {
    return frame;
  }

  /**
   * Checks for id. Derived classes can override this.
   *
   * @return true, if successful
   */
  public boolean hasId() {
    return false;
  }

  /**
   * Gets the id. Default = 0. Derived classes can override this.
   *
   * @return The results identifier
   */
  public int getId() {
    return 0;
  }

  /**
   * Checks for precision. Derived classes can override this.
   *
   * @return true, if successful
   */
  public boolean hasPrecision() {
    return false;
  }

  /**
   * Gets the localisation precision. Default = Double.NaN. Derived classes can override this.
   *
   * <p>This is provided so that results can be loaded from external sources.
   *
   * @return the precision (in nm)
   */
  public double getPrecision() {
    return Double.NaN;
  }

  /**
   * Return the true positive score for use in classification analysis.
   *
   * @return The true positive score
   */
  public double getTruePositiveScore() {
    return (origValue != 0) ? 1 : 0;
  }

  /**
   * Return the false positive score for use in classification analysis.
   *
   * @return The false positive score
   */
  public double getFalsePositiveScore() {
    return 1 - getTruePositiveScore();
  }

  /**
   * Return the true negative score for use in classification analysis.
   *
   * @return The true negative score
   */
  public double getTrueNegativeScore() {
    return (origValue != 0) ? 0 : 1;
  }

  /**
   * Return the false negative score for use in classification analysis.
   *
   * @return The false negative score
   */
  public double getFalseNegativeScore() {
    return 1 - getTrueNegativeScore();
  }

  /**
   * Return the squared distance to the other peak result.
   *
   * @param result The result
   * @return The squared distance
   */
  public double distance2(PeakResult result) {
    final double dx = getXPosition() - result.getXPosition();
    final double dy = getYPosition() - result.getYPosition();
    return dx * dx + dy * dy;
  }

  /**
   * Return the squared distance to the other coordinate.
   *
   * @param x the x
   * @param y the y
   * @return The squared distance
   */
  public double distance2(final double x, final double y) {
    final double dx = getXPosition() - x;
    final double dy = getYPosition() - y;
    return dx * dx + dy * dy;
  }

  /**
   * Return the distance to the other peak result.
   *
   * @param result The result
   * @return The distance
   */
  public double distance(PeakResult result) {
    return Math.sqrt(distance2(result));
  }

  /**
   * Return the distance to the other coordinate.
   *
   * @param x the x
   * @param y the y
   * @return The distance
   */
  public double distance(final double x, final double y) {
    return Math.sqrt(distance2(x, y));
  }

  /**
   * This methods return the x-position. To allow filters to use the actual shift requires either
   * off-setting the position with the initial fit position, or extending this class so the shift
   * can be stored.
   *
   * @return the x shift
   */
  public float getXShift() {
    return getXPosition();
  }

  /**
   * This methods return the y-position. To allow filters to use the actual shift requires either
   * off-setting the position with the initial fit position, or extending this class so the shift
   * can be stored.
   *
   * @return the y shift
   */
  public float getYShift() {
    return getYPosition();
  }

  /**
   * Gets the parameters. This is a direct reference to the instance parameter array so use with
   * caution.
   *
   * @return the parameters
   */
  public float[] getParameters() {
    return params;
  }

  /**
   * Checks for parameter deviations.
   *
   * @return true, if successful
   */
  public boolean hasParameterDeviations() {
    return paramStdDevs != null;
  }

  /**
   * Gets the parameter deviations. This is a direct reference to the instance parameter array so
   * use with caution.
   *
   * @return the parameter deviations
   */
  public float[] getParameterDeviations() {
    return paramStdDevs;
  }

  /**
   * Gets the number of parameters.
   *
   * @return the number of parameters
   */
  public int getNumberOfParameters() {
    return params.length;
  }

  /**
   * Gets the parameter for the given index.
   *
   * @param index the index
   * @return the parameter
   */
  public float getParameter(int index) {
    return params[index];
  }

  /**
   * Sets the parameter for the given index.
   *
   * @param index the index
   * @param value the value
   */
  public void setParameter(int index, float value) {
    params[index] = value;
  }

  /**
   * Gets the parameter deviation for the given index.
   *
   * @param index the index
   * @return the parameter deviation
   */
  public float getParameterDeviation(int index) {
    return paramStdDevs[index];
  }

  /**
   * Sets the parameter deviation for the given index.
   *
   * @param index the index
   * @param value the value
   */
  public void setParameterDeviation(int index, float value) {
    paramStdDevs[index] = value;
  }

  /**
   * Resize the parameters.
   *
   * @param length the new length
   */
  void resizeParameters(int length) {
    params = Arrays.copyOf(params, length);
  }

  /**
   * Resize parameter deviations.
   *
   * @param length the length
   */
  void resizeParameterDeviations(int length) {
    paramStdDevs = Arrays.copyOf(paramStdDevs, length);
  }
}
