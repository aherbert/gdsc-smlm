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
package gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.Collection;

import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;


/**
 * Specifies the interface for saving peak fitting results
 */
public interface PeakResults
{
	/**
	 * Should be called at the start of fitting to prepare the output.
	 */
	public void begin();

	/**
	 * Add a fitted peak result.
	 *
	 * @param peak
	 *            The peak number
	 * @param origX
	 *            The original X value
	 * @param origY
	 *            The original Y value
	 * @param origValue
	 *            The original value
	 * @param error
	 *            The error of the fit
	 * @param noise
	 *            Estimate of the noise in the signal
	 * @param meanSignal
	 *            Estimate of the mean signal
	 * @param params
	 *            The peak parameters
	 * @param paramsStdDev
	 *            The peak parameters standard deviations
	 */
	public void add(int peak, int origX, int origY, float origValue, double error, float noise, float meanSignal,
			float[] params, float[] paramsStdDev);

	/**
	 * Add a fitted peak result
	 *
	 * @param result
	 *            the result
	 */
	public void add(PeakResult result);

	/**
	 * Add a series of fitted peak results.
	 *
	 * @param results
	 *            the results
	 */
	public void addAll(Collection<PeakResult> results);

	/**
	 * Add a series of fitted peak results.
	 *
	 * @param results
	 *            the results
	 */
	public void addAll(PeakResult[] results);

	/**
	 * Add a series of fitted peak results.
	 *
	 * @param results
	 *            the results
	 */
	public void addAll(PeakResultStore results);

	/**
	 * @return The number of results added since begin()
	 */
	public int size();

	/**
	 * Called at the end of fitting to finalise the output.
	 */
	public void end();

	/**
	 * @return True if still accepting results using the add methods
	 */
	public boolean isActive();

	/**
	 * @param source
	 *            The source used to create the results
	 */
	public void setSource(ImageSource source);

	/**
	 * @return The source used to create the results
	 */
	public ImageSource getSource();

	/**
	 * Set the bounds of the results. All fitting results are expected to be within the bounds, i.e. results were
	 * created by fitting a rectangle taken from the image source.
	 * <p>
	 * Note that the bounds are relative to the width and height of the image source. They do not include the (x,y)
	 * origin of the image source.
	 * 
	 * @param bounds
	 *            The bounds of the image source used to create the results
	 */
	public void setBounds(Rectangle bounds);

	/**
	 * Get the bounds of the rectangle taken from the image source that encapsulates all the fitting results.
	 * <p>
	 * Note that the bounds are relative to the width and height of the image source. They do not include the (x,y)
	 * origin of the image source.
	 * 
	 * @return The bounds used to create the results
	 */
	public Rectangle getBounds();

	/**
	 * @param calibration
	 *            The calibration used to obtain the results
	 */
	public void setCalibration(Calibration calibration);

	/**
	 * @return The calibration used to obtain the results
	 */
	public Calibration getCalibration();

	/**
	 * Gets the Point Spread Function (PSF) used when fitting the results.
	 *
	 * @return the psf
	 */
	public PSF getPSF();

	/**
	 * Sets the Point Spread Function (PSF) used when fitting the results.
	 *
	 * @param psf
	 *            the new psf
	 */
	public void setPSF(PSF psf);

	/**
	 * @param configuration
	 *            The configuration used to create the results
	 */
	public void setConfiguration(String configuration);

	/**
	 * @return The configuration used to create the results
	 */
	public String getConfiguration();

	/**
	 * @return The name of the results set
	 */
	public String getName();

	/**
	 * @param name
	 *            The name of the results set
	 */
	public void setName(String name);

	/**
	 * Copy the settings (source, bounds, configuration) from the given results
	 * 
	 * @param peakResults
	 */
	public void copySettings(PeakResults peakResults);
}
