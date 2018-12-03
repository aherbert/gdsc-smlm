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
package uk.ac.sussex.gdsc.smlm.results;

/**
 * Specifies a peak fitting result that spans multiple frames.
 */
public class ExtendedPeakResult extends IdPeakResult
{
    private int endFrame;

    /**
     * Instantiates a new peak result.
     *
     * @param startFrame
     *            the start frame
     * @param origX
     *            the original X position
     * @param origY
     *            the original Y position
     * @param origValue
     *            the original value
     * @param error
     *            the error
     * @param noise
     *            the noise
     * @param meanIntensity
     *            the mean intensity
     * @param params
     *            the params (must not be null and must have at least
     *            {@value uk.ac.sussex.gdsc.smlm.results.PeakResult#STANDARD_PARAMETERS} parameters)
     * @param paramsStdDev
     *            the params standard deviations (if not null must match the length of the params array)
     * @param endFrame
     *            the end frame
     * @param id
     *            the id
     * @throws IllegalArgumentException
     *             the illegal argument exception if the parameters are invalid
     */
    public ExtendedPeakResult(int startFrame, int origX, int origY, float origValue, double error, float noise,
            float meanIntensity, float[] params, float[] paramsStdDev, int endFrame, int id)
            throws IllegalArgumentException
    {
        super(startFrame, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev, id);
        setEndFrame(endFrame);
    }

    /**
     * Simple constructor to create a result with frame, location, width, strength, and id.
     *
     * @param frame
     *            the frame
     * @param x
     *            the x position
     * @param y
     *            the y position
     * @param intensity
     *            the intensity
     * @param id
     *            the id
     */
    public ExtendedPeakResult(int frame, float x, float y, float intensity, int id)
    {
        super(frame, x, y, intensity, id);
    }

    /**
     * Simple constructor to create a result with location, width, strength, and id.
     *
     * @param x
     *            the x position
     * @param y
     *            the y position
     * @param intensity
     *            the intensity
     * @param id
     *            the id
     */
    public ExtendedPeakResult(float x, float y, float intensity, int id)
    {
        super(x, y, intensity, id);
    }

    /** {@inheritDoc} */
    @Override
    public boolean hasEndFrame()
    {
        return endFrame > getFrame();
    }

    /** {@inheritDoc} */
    @Override
    public int getEndFrame()
    {
        return endFrame;
    }

    /**
     * Sets the end frame.
     *
     * @param endFrame
     *            the new end frame
     */
    public void setEndFrame(int endFrame)
    {
        // Ensure that the end frame is valid
        this.endFrame = Math.max(endFrame, super.getFrame());
    }

    @Override
    public void setFrame(int frame)
    {
        // Set the new start frame
        super.setFrame(frame);
        // Validate the current end frame
        setEndFrame(endFrame);
    }
}
