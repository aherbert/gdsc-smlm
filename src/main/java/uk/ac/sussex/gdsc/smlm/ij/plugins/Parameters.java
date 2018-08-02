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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.util.EnumSet;

/**
 * Utility methods for checking parameters
 */
public class Parameters
{
    /**
     * Define the requirements for the parameter
     */
    public enum Requirement
    {
    /**
     * Value is greater than zero
     */
    ABOVE_ZERO,
    /**
     * Value is zero or above
     */
    POSITIVE
    }

    /**
     * Check if the named parameter value is greater than zero.
     *
     * @param name
     *            the name
     * @param value
     *            the value
     * @throws IllegalArgumentException
     *             the illegal argument exception
     */
    public static void isAboveZero(String name, double value)
    {
        if (value <= 0)
            throw new IllegalArgumentException(name + " should be above zero");
    }

    /**
     * Check if the named parameter value is greater than the given limit.
     *
     * @param name
     *            the name
     * @param value
     *            the value
     * @param limit
     *            the limit
     * @throws IllegalArgumentException
     *             the illegal argument exception
     */
    public static void isAbove(String name, double value, double limit)
    {
        if (value <= limit)
            throw new IllegalArgumentException(name + " should be > " + limit);
    }

    /**
     * Check if the named parameter value is greater than or equal to the given limit.
     *
     * @param name
     *            the name
     * @param value
     *            the value
     * @param limit
     *            the limit
     * @throws IllegalArgumentException
     *             the illegal argument exception
     */
    public static void isEqualOrAbove(String name, double value, double limit)
    {
        if (value < limit)
            throw new IllegalArgumentException(name + " should be >= " + limit);
    }

    /**
     * Check if the named parameter value is less than the given limit.
     *
     * @param name
     *            the name
     * @param value
     *            the value
     * @param limit
     *            the limit
     * @throws IllegalArgumentException
     *             the illegal argument exception
     */
    public static void isBelow(String name, double value, double limit)
    {
        if (value >= limit)
            throw new IllegalArgumentException(name + " should be < " + limit);
    }

    /**
     * Check if the named parameter value is less then or equal to the given limit.
     *
     * @param name
     *            the name
     * @param value
     *            the value
     * @param limit
     *            the limit
     * @throws IllegalArgumentException
     *             the illegal argument exception
     */
    public static void isEqualOrBelow(String name, double value, double limit)
    {
        if (value > limit)
            throw new IllegalArgumentException(name + " should be <= " + limit);
    }

    /**
     * Check if the named parameter value is zero or greater.
     *
     * @param name
     *            the name
     * @param value
     *            the value
     * @throws IllegalArgumentException
     *             the illegal argument exception
     */
    public static void isPositive(String name, double value)
    {
        if (value < 0)
            throw new IllegalArgumentException(name + " should be positive");
    }

    /**
     * Check if the named parameter meets the requirements.
     *
     * @param name
     *            the name
     * @param value
     *            the value
     * @param requirements
     *            the requirements
     */
    public static void isValid(String name, double value, EnumSet<Requirement> requirements)
    {
        if (requirements.contains(Requirement.ABOVE_ZERO))
            isAboveZero(name, value);
        if (requirements.contains(Requirement.POSITIVE))
            isPositive(name, value);
    }
}
