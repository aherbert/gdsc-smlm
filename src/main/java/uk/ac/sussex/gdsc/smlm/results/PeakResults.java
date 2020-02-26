/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import java.awt.Rectangle;
import java.util.Collection;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;

/**
 * Specifies the interface for saving peak fitting results.
 */
public interface PeakResults {
  /**
   * Should be called at the start of fitting to prepare the output.
   */
  void begin();

  /**
   * Add a fitted peak result.
   *
   * @param peak The peak number
   * @param origX The original X value
   * @param origY The original Y value
   * @param origValue The original value
   * @param error The error of the fit
   * @param noise Estimate of the noise in the signal
   * @param meanSignal Estimate of the mean signal
   * @param params The peak parameters
   * @param paramsStdDev The peak parameters standard deviations
   */
  void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanSignal, float[] params, float[] paramsStdDev);

  /**
   * Add a fitted peak result.
   *
   * @param result the result
   */
  void add(PeakResult result);

  /**
   * Add a series of fitted peak results.
   *
   * @param results the results
   */
  void addAll(Collection<PeakResult> results);

  /**
   * Add a series of fitted peak results.
   *
   * @param results the results
   */
  void addAll(PeakResult[] results);

  /**
   * Add a series of fitted peak results.
   *
   * @param results the results
   */
  void addAll(PeakResultStore results);

  /**
   * Get the number of results added since {@link #begin()}.
   *
   * @return The number of results added
   */
  int size();

  /**
   * Called at the end of fitting to finalise the output.
   */
  void end();

  /**
   * Checks if still accepting results using the add methods.
   *
   * @return True if still accepting results using the add methods.
   */
  boolean isActive();

  /**
   * Sets the source used to create the results.
   *
   * @param source The source used to create the results
   */
  void setSource(ImageSource source);

  /**
   * Gets the source used to create the results.
   *
   * @return The source used to create the results.
   */
  ImageSource getSource();

  /**
   * Set the bounds of the results. All fitting results are expected to be within the bounds, i.e.
   * results were created by fitting a rectangle taken from the image source.
   *
   * <p>Note that the bounds are relative to the width and height of the image source. They do not
   * include the (x,y) origin of the image source.
   *
   * @param bounds The bounds of the image source used to create the results
   */
  void setBounds(Rectangle bounds);

  /**
   * Get the bounds of the rectangle taken from the image source that encapsulates all the fitting
   * results.
   *
   * <p>Note that the bounds are relative to the width and height of the image source. They do not
   * include the (x,y) origin of the image source.
   *
   * @return The bounds used to create the results
   */
  Rectangle getBounds();

  /**
   * Sets the calibration used to obtain the results.
   *
   * @param calibration The calibration used to obtain the results
   */
  void setCalibration(Calibration calibration);

  /**
   * Gets the calibration used to obtain the results.
   *
   * @return The calibration used to obtain the results.
   */
  Calibration getCalibration();

  /**
   * Gets the Point Spread Function (PSF) used when fitting the results.
   *
   * @return the PSF
   */
  PSF getPsf();

  /**
   * Sets the Point Spread Function (PSF) used when fitting the results.
   *
   * @param psf the new PSF
   */
  void setPsf(PSF psf);

  /**
   * Sets the configuration used to obtain the results.
   *
   * @param configuration The configuration used to create the results
   */
  void setConfiguration(String configuration);

  /**
   * Gets the configuration used to obtain the results.
   *
   * @return The configuration used to create the results.
   */
  String getConfiguration();

  /**
   * Gets the name of the results set.
   *
   * @return The name of the results set.
   */
  String getName();

  /**
   * Sets the name of the results set.
   *
   * @param name The name of the results set
   */
  void setName(String name);

  /**
   * Copy the settings (source, bounds, configuration) from the given results.
   *
   * @param peakResults the peak results
   */
  void copySettings(PeakResults peakResults);
}
