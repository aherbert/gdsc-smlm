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
package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.core.match.FractionalAssignment;

/**
 * Specifies a peak fitting result for use in filtering. Any result implementing this interface can be directly filtered
 * without requiring the filter to be initialised with calibration data.
 */
public interface PreprocessedPeakResult
{
	/**
	 * Get the frame containing the result
	 *
	 * @return The frame
	 */
	public int getFrame();

	/**
	 * Gets the unique id. This should be unique within the entire set of results.
	 *
	 * @return the unique id
	 */
	public int getUniqueId();

	/**
	 * Gets the id.
	 *
	 * @return the id
	 */
	public int getId();

	/**
	 * Return the candidate Id of this result (i.e. the candidate used to identify this position for fitting)
	 *
	 * @return the candidate Id (-1 for no candidate)
	 */
	public int getCandidateId();

	/**
	 * Get the signal strength (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 *
	 * @return The signal (in photons)
	 */
	public float getSignal();

	/**
	 * Get the mean signal.
	 * <p>
	 * This requires a knowledge of the PSF used to create the result. It could be the peak signal in the PSF or the
	 * average signal over a range of the PSF, e.g. the area covered from the maxima to half-maxima for spots.
	 *
	 * @return The mean signal (in photons)
	 */
	public float getMeanSignal();

	/**
	 * Get the signal-to-noise ratio (SNR). This is ratio of the average signal value to the standard deviation of the
	 * background. Ideally the standard deviation of the background is computed in the region around the centre.
	 * <p>
	 * The default implementation is {@link #getMeanSignal()} / {@link #getNoise()}.
	 *
	 * @return The SNR
	 */
	public default float getSNR()
	{
		return getMeanSignal() / getNoise();
	}

	/**
	 * Get the background noise
	 *
	 * @return The noise (in photons)
	 */
	public float getNoise();

	/**
	 * The localisation variance using the noise estimate for the data to approximate the local noise
	 *
	 * @return The location variance in nm
	 */
	public double getLocationVariance();

	/**
	 * The localisation variance using the local background estimate for the local noise
	 *
	 * @return The location variance in nm
	 */
	public double getLocationVariance2();

	/**
	 * The localisation variance using the fitted parameter variance. This is computed using the Cramér-Rao lower bound
	 * (CRLB) on the variance of the estimators. The variance for the fitted X and Y position is averaged to produce a
	 * localisation precision.
	 *
	 * @return The location variance in nm
	 */
	public double getLocationVarianceCRLB();

	/**
	 * @return The average peak standard deviation in the X and Y dimension
	 */
	public float getSD();

	/**
	 * @return The background
	 */
	public float getBackground();

	/**
	 * Get the amplitude. Amplitude = Signal / (2*pi*sd0*sd1).
	 *
	 * @return The amplitude
	 */
	public float getAmplitude();

	/**
	 * @return The angle (for an elliptical Gaussian peak)
	 */
	public float getAngle();

	/**
	 * @return The x position
	 */
	public float getX();

	/**
	 * @return The y position
	 */
	public float getY();

	/**
	 * @return The z position
	 */
	public float getZ();

	/**
	 * Return the squared shift between the initial and the final x-position, i.e. how much the position moved during
	 * fitting, relative to the initial peak SD
	 *
	 * @return The relative x position shift squared
	 */
	public float getXRelativeShift2();

	/**
	 * Return the squared shift between the initial and the final y-position, i.e. how much the position moved during
	 * fitting, relative to the initial peak SD.
	 *
	 * @return The relative y position shift squared
	 */
	public float getYRelativeShift2();

	/**
	 * @return The x-dimension standard deviation
	 */
	public float getXSD();

	/**
	 * @return The y-dimension standard deviation
	 */
	public float getYSD();

	/**
	 * Return the ratio between the initial and the final x-dimension standard deviation, i.e. how much the width
	 * changed during fitting
	 *
	 * @return The x-dimension width factor
	 */
	public float getXSDFactor();

	/**
	 * Return the ratio between the initial and the final y-dimension standard deviation, i.e. how much the width
	 * changed during fitting
	 *
	 * @return The y-dimension width factor
	 */
	public float getYSDFactor();

	/**
	 * Return true if this is a result that has previously been fitted.
	 *
	 * @return True if this result is an existing fit result
	 */
	public boolean isExistingResult();

	/**
	 * Return true if this is a new result. It is expected that this can then be classified for use in
	 * scoring analysis. The {@link #getAssignments(int)} method can then be called to return the assignments.
	 *
	 * @return True if this result is a new fit result
	 */
	public boolean isNewResult();

	/**
	 * Get the assignments between this result and the true data.
	 * <p>
	 * The assignments should all have the same predicted Id. The actual Id should be a value starting from 0 and
	 * incrementing for each actual result in the frame that is scored.
	 *
	 * @param predictedId
	 *            The predicted Id
	 * @return The assignments
	 */
	public FractionalAssignment[] getAssignments(int predictedId);

	/**
	 * Ignore this result during assignment scoring. It is expected that the result will return null from
	 * {@link #getAssignments(int)}.
	 *
	 * @return true, if this should be ignored (i.e. not counted as a false positive)
	 */
	public boolean ignore();

	/**
	 * Convert this to the parameters for a Gaussian2DFunction
	 *
	 * @return the parameters
	 */
	public double[] toGaussian2DParameters();

	/**
	 * Sets the validation result.
	 *
	 * @param result
	 *            the new validation result
	 */
	public void setValidationResult(int result);

	/**
	 * Gets the validation result.
	 *
	 * @return the validation result
	 */
	public int getValidationResult();

	/**
	 * Returns true if this result is not a duplicate. The default value should be false.
	 * <p>
	 * Implementations can preprocess a results set to check if this is close to any preceeding results. If it is
	 * impossible to be a duplicate then the return value is true.
	 *
	 * @return true, if is not duplicate.
	 */
	public boolean isNotDuplicate();
}
