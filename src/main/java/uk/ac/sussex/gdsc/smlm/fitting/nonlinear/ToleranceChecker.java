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
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

/**
 * Check if converged using a tolerance on the value, parameters, and/or the number of iterations
 */
public class ToleranceChecker
{
    /** The constant to use for a tolerance that is ignored */
    public static final double IGNORE_TOLERANCE = -1.0;
    /** The constant to use for a max iterations that is ignored */
    public static final int IGNORE_MAX_ITERATIONS = 0;

    /** Flag to indicate the max iterations have been reached. This is a failure to converge. */
    public static final int STATUS_MAX_ITERATIONS = 0x00000001;
    /** Flag to indicate convergence on the value. */
    public static final int STATUS_VALUE = 0x00000002;
    /** Flag to indicate convergence on the parameters. */
    public static final int STATUS_PARAMETERS = 0x00000004;
    /** Flag to indicate convergence on the target number of iterations. */
    public static final int STATUS_TARGET_ITERATIONS = 0x00000008;
    /** Flag to indicate convergence was set manually. */
    public static final int STATUS_MANUAL_CONVERGENCE = 0x00000010;
    /** Flag to indicate all valid convergence flags. */
    public static final int STATUS_CONVERGED = STATUS_VALUE | STATUS_PARAMETERS | STATUS_TARGET_ITERATIONS |
            STATUS_MANUAL_CONVERGENCE;

    /** The relative tolerance threshold for the value. Set to negative to disable. */
    public final double relativeValue;
    /** The absolute tolerance threshold for the value. Set to negative to disable. */
    public final double absoluteValue;
    /** The relative tolerance threshold for the parameters. Set to negative to disable. */
    public final double relativeParameters;
    /** The absolute tolerance threshold for the parameters. Set to negative to disable. */
    public final double absoluteParameters;
    /**
     * Flag indicating if the value will be checked for convergence. Either {@link #relativeValue} or
     * {@link #absoluteValue} will be above zero.
     */
    public final boolean checkValue;
    /**
     * Flag indicating if the parameters will be checked for convergence. Either {@link #relativeParameters} or
     * {@link #absoluteParameters} will be above zero.
     */
    public final boolean checkParameters;
    /**
     * Flag indicating that the value must be minimised. Convergence is only signalled on the value with a marginal
     * improvement in the value with the direction being either minimal or maximal, i.e. if the value becomes worse but
     * the change is within the tolerance then convergence will not be signalled.
     */
    private boolean minimiseValue;
    /**
     * The maximum number of allowed iterations.
     * Set above zero to limit the iterations. Set below zero to define a specific number of iterations until
     * convergence is signalled. Set to zero to disable.
     */
    public final int maxIterations;

    private int iterations = 0;
    private boolean manualConvergence = false;

    /**
     * Build an instance with specified thresholds. This only checks convergence using the parameters.
     * <p>
     * In order to perform only relative checks, the absolute tolerance must be set to a negative value. In order to
     * perform only absolute checks, the relative tolerance must be set to a negative value.
     * <p>
     * Note: If a tolerance is set to zero then an exact match will achieve equality. Only negative values disable the
     * tolerance.
     *
     * @param relativeParameters
     *            relative tolerance threshold on the parameters. Set to negative to disable.
     * @param absoluteParameters
     *            absolute tolerance threshold on the parameters. Set to negative to disable.
     * @throws IllegalArgumentException
     *             if none of the convergence criteria are valid (i.e. convergence is not possible)
     */
    public ToleranceChecker(double relativeParameters, double absoluteParameters)
    {
        this(false, IGNORE_TOLERANCE, IGNORE_TOLERANCE, relativeParameters, absoluteParameters, IGNORE_MAX_ITERATIONS);
    }

    /**
     * Build an instance with specified thresholds.
     * <p>
     * In order to perform only relative checks, the absolute tolerance must be set to a negative value. In order to
     * perform only absolute checks, the relative tolerance must be set to a negative value.
     * <p>
     * Note: If a tolerance is set to zero then an exact match will achieve equality. Only negative values disable the
     * tolerance.
     *
     * @param relativeValue
     *            relative tolerance threshold on the value. Set to negative to disable.
     * @param absoluteValue
     *            absolute tolerance threshold on the value. Set to negative to disable.
     * @param relativeParameters
     *            relative tolerance threshold on the parameters. Set to negative to disable.
     * @param absoluteParameters
     *            absolute tolerance threshold on the parameters. Set to negative to disable.
     * @param maxIterations
     *            The maximum number of allowed iterations.
     *            Set above zero to limit the iterations. Set below zero to define a specific number of iterations until
     *            convergence is signalled. Set to zero to disable.
     * @throws IllegalArgumentException
     *             if none of the convergence criteria are valid (i.e. convergence is not possible)
     */
    public ToleranceChecker(double relativeValue, double absoluteValue, double relativeParameters,
            double absoluteParameters, int maxIterations)
    {
        this(true, relativeValue, absoluteValue, relativeParameters, absoluteParameters, maxIterations);
    }

    /**
     * Build an instance with specified thresholds.
     * <p>
     * In order to perform only relative checks, the absolute tolerance must be set to a negative value. In order to
     * perform only absolute checks, the relative tolerance must be set to a negative value.
     * <p>
     * Note: If a tolerance is set to zero then an exact match will achieve equality. Only negative values disable the
     * tolerance.
     *
     * @param minimiseValue
     *            Set to true to ensure the value is minimised at converge (otherwise it is maximised)
     * @param relativeValue
     *            relative tolerance threshold on the value. Set to negative to disable.
     * @param absoluteValue
     *            absolute tolerance threshold on the value. Set to negative to disable.
     * @param relativeParameters
     *            relative tolerance threshold on the parameters. Set to negative to disable.
     * @param absoluteParameters
     *            absolute tolerance threshold on the parameters. Set to negative to disable.
     * @param maxIterations
     *            The maximum number of allowed iterations.
     *            Set above zero to limit the iterations. Set below zero to define a specific number of iterations until
     *            convergence is signalled. Set to zero to disable.
     * @throws IllegalArgumentException
     *             if none of the convergence criteria are valid (i.e. convergence is not possible)
     */
    public ToleranceChecker(boolean minimiseValue, double relativeValue, double absoluteValue,
            double relativeParameters, double absoluteParameters, int maxIterations)
    {
        checkValue = (relativeValue >= 0 || absoluteValue >= 0);
        checkParameters = (relativeParameters >= 0 || absoluteParameters >= 0);

        if (!(checkValue || checkParameters || maxIterations != 0))
            throw new IllegalArgumentException("No valid convergence criteria");

        this.minimiseValue = minimiseValue;
        this.relativeValue = relativeValue;
        this.absoluteValue = absoluteValue;
        this.relativeParameters = relativeParameters;
        this.absoluteParameters = absoluteParameters;
        this.maxIterations = maxIterations;
    }

    /**
     * Check if all the pairs of values are equal
     *
     * @param p
     *            Previous
     * @param c
     *            Current
     * @param relative
     *            relative tolerance threshold (set to negative to ignore)
     * @param absolute
     *            absolute tolerance threshold (set to negative to ignore)
     * @return True if equal
     */
    public static boolean areEqual(final double[] p, final double[] c, double absolute, double relative)
    {
        for (int i = 0; i < p.length; ++i)
            if (!areEqual(p[i], c[i], absolute, relative))
                return false;
        return true;
    }

    /**
     * Check if the pair of values are equal
     *
     * @param p
     *            Previous
     * @param c
     *            Current
     * @param relative
     *            relative tolerance threshold (set to negative to ignore)
     * @param absolute
     *            absolute tolerance threshold (set to negative to ignore)
     * @return True if equal
     */
    public static boolean areEqual(final double p, final double c, double absolute, double relative)
    {
        final double difference = Math.abs(p - c);
        if (difference <= absolute)
            return true;
        final double size = max(Math.abs(p), Math.abs(c));
        return (difference <= size * relative);
    }

    private static double max(final double a, final double b)
    {
        // Ignore NaN
        return (a > b) ? a : b;
    }

    /**
     * Check if converged. All the conditions for convergence will be set in the returned status flag.
     *
     * @param previousValue
     *            the previous value
     * @param previousParameters
     *            the previous parameters
     * @param currentValue
     *            the current value
     * @param currentParameters
     *            the current parameters
     * @return The status flag. Non-zero for convergence.
     */
    public int converged(double previousValue, double[] previousParameters, double currentValue,
            double[] currentParameters)
    {
        iterations++;
        int status = 0;
        if (checkValue && correctDirection(previousValue, currentValue) &&
                areEqual(previousValue, currentValue, absoluteValue, relativeValue))
            status |= STATUS_VALUE;
        if (checkParameters && areEqual(previousParameters, currentParameters, absoluteParameters, relativeParameters))
            status |= STATUS_PARAMETERS;
        if (maxIterations != 0 && iterations >= Math.abs(maxIterations))
            status |= (maxIterations < 0) ? STATUS_TARGET_ITERATIONS : STATUS_MAX_ITERATIONS;
        if (manualConvergence)
            status |= STATUS_MANUAL_CONVERGENCE;
        return status;
    }

    private boolean correctDirection(double previousValue, double currentValue)
    {
        return (minimiseValue) ? currentValue <= previousValue : currentValue >= previousValue;
    }

    /**
     * Gets the iterations. The iterations count is incremented on each call to
     * {@link #converged(double, double[], double, double[])}.
     *
     * @return the iterations
     */
    public int getIterations()
    {
        return iterations;
    }

    /**
     * Checks if the value must be minimised.
     *
     * @return true, if the value must be minimised
     */
    public boolean isMinimiseValue()
    {
        return minimiseValue;
    }

    /**
     * Sets the flag indicating that the value must be minimised. Convergence is only signalled on the value with a
     * marginal improvement in the value with the direction being either minimal or maximal, i.e. if the value becomes
     * worse but the change is within the tolerance then convergence will not be signalled.
     *
     * @param minimiseValue
     *            true, if the value must be minimised
     */
    public void setMinimiseValue(boolean minimiseValue)
    {
        this.minimiseValue = minimiseValue;
    }

    /**
     * Reset the iterations allowing the checker to be reused.
     */
    public void reset()
    {
        iterations = 0;
        manualConvergence = false;
    }

    /**
     * Sets the converged flag to true. This can be used to manually set convergence, for example in the event that
     * further computation is not needed.
     */
    public void setConverged()
    {
        manualConvergence = true;
    }
}
