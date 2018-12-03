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
package uk.ac.sussex.gdsc.smlm.search;

/**
 * Specify the dimensions for a search.
 */
public interface Dimension
{
    /**
     * Gets the current lower bound of the range.
     *
     * @return the current lower bound of the range
     */
    public double getLower();

    /**
     * Gets the current upper bound of the range.
     *
     * @return the current upper bound of the range
     */
    public double getUpper();

    /**
     * Gets the current centre of the range.
     *
     * @return the current centre of the range
     */
    public double getCentre();

    /**
     * Gets the minimum allowed value of the range.
     *
     * @return the minimum
     */
    public double getMin();

    /**
     * Gets the maximum allowed value of the range.
     *
     * @return the maximum
     */
    public double getMax();

    /**
     * Checks if is active.
     *
     * @return true, if is active
     */
    public boolean isActive();

    /**
     * Checks if the value is at (or beyond) the lower/upper bounds of the current dimension range.
     *
     * @param v
     *            the value
     * @return true, if is at bounds
     */
    public boolean isAtBounds(double v);

    /**
     * Creates a new dimension with the given bounds. The current min/max and other settings should be respected.
     *
     * @param lower
     *            the lower
     * @param upper
     *            the upper
     * @return the dimension
     */
    public Dimension create(double lower, double upper);

    /**
     * Round the value to the working resolution of the dimension. The resolution defines the minimum delta between
     * values in the dimension; it can be zero in which case no rounding is performed.
     *
     * @param value
     *            the value
     * @return the rounded value
     */
    public double round(double value);

    /**
     * True if the dimension can round values.
     *
     * @return true, if can round
     */
    public boolean canRound();
}
