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
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stop;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.StoppingCriteria;

/**
 * Defines the stopping criteria for the {@link uk.ac.sussex.gdsc.smlm.fitting.nonlinear.NonLinearFit } class.
 * <p>
 * Stop when N successive iterations reduce error by a negligible amount.
 */
public class ErrorStoppingCriteria extends StoppingCriteria
{
	private int iterationCount;
	private int iterationLimit = 1;
	private int significantDigits;
	private double maxRelativeError;
	private boolean avoidPlateau = false;

	private int insignificantImprovmentIteration;
	private int improvementExponent;
	private int insignificantImprovmentCounter;
	private static final int CONSECUTIVE_INSIGNIFICANT_IMPROVEMENT_LIMIT = 3;

	/**
	 * Instantiates a new error stopping criteria.
	 */
	public ErrorStoppingCriteria()
	{
		this(5);
	}

	/**
	 * Instantiates a new error stopping criteria.
	 *
	 * @param significantDigits
	 *            the significant digits for a negligible change in error
	 */
	public ErrorStoppingCriteria(int significantDigits)
	{
		setSignificantDigits(significantDigits);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.fitting.model.StoppingCriteria#initialise(double[])
	 */
	@Override
	public void initialise(double[] a)
	{
		super.initialise(a);
		iterationCount = 0;

		// Used to avoid flats of insignificant improvement later in the routine
		if (avoidPlateau)
		{
			insignificantImprovmentIteration = getMaximumIterations() / 2;
			improvementExponent = 0;
			insignificantImprovmentCounter = 0;
		}
		else
			insignificantImprovmentIteration = getMaximumIterations();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.fitting.model.StoppingCriteria#copyCoefficients(double[])
	 */
	@Override
	protected void copyCoefficients(double[] a)
	{
		// Do nothing
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stoppingCriteria#evaluate(double, double, double[])
	 */
	@Override
	public void evaluate(double oldError, double newError, double[] a)
	{
		// Use a comparison of error to a set number of significant digits

		int result;

		if (newError > oldError)
		{
			// Fit is worse
			// Setting the iteration count to zero forces the negligible improvements to be sequential.
			//
			// Note: The NonLinearFit algorithm may jump around a low minimum finding marginally
			// higher error fits. If the count is reset each time then N-sequential improvements
			// may not be found. In practice this reset is not necessary.
			//iterationCount = 0;

			result = 1;
			increment(a, false);
		}
		else
		{
			//final double diff = oldError - newError;
			// Allow debugging of the change in error
			// Record how many significant digits the improvement was
			//System.out.printf("Iter %d : Abs %g : Rel %g (%d)\n", getIteration(), diff,
			//		DoubleEquality.relativeError(newError, oldError),
			//		getExponent(DoubleEquality.relativeError(newError, oldError)));

			// Check if equal or the fit is near perfect
			boolean negligable = false;
			if (DoubleEquality.almostEqualRelativeOrAbsolute(oldError, newError, maxRelativeError, 0) ||
					newError < 0.001 || false)
				negligable = true;
			else // The desire is to avoid a lot of computation if the error is not really moving anywhere.
			// Once we have reached a set number of iterations, e.g. maximumIterations/2
			// then allow stopping if the absolute change is the same.
			if (getIteration() > insignificantImprovmentIteration)
			{
				// Get the exponent of the improvement
				final int exp = getExponent(DoubleEquality.relativeError(newError, oldError));

				// All we know is that there is an improvement but not enough to signal convergence.
				// If the improvements get smaller then it may converge anyway.
				// If the improvements get larger then the routine should be allowed to run.
				// So check to see if the improvement is the same.
				// Q. Should a check be made to see if the exponent is small?
				if (improvementExponent == exp)
					negligable = (++insignificantImprovmentCounter >= CONSECUTIVE_INSIGNIFICANT_IMPROVEMENT_LIMIT);
					//System.out.printf("Tiny improvement %f = %f @ %d\n", oldError - newError,
					//		DoubleEquality.relativeError(newError, oldError), getIteration());
				else
				{
					insignificantImprovmentCounter = 0;
					improvementExponent = exp;
				}
			}

			if (negligable)
			{
				// Improvement is negligible
				result = 0;
				iterationCount++;
				if (iterationCount >= getIterationLimit() && getIteration() > getMinimumIterations())
				{
					areAchieved = true;
					notSatisfied = false;
				}
			}
			else
			{
				// Improvement is significant
				// Setting the iteration count to zero forces negligible improvements to be successive
				// (i.e. no significant jumps between insignificant jumps)
				result = -1;
				iterationCount = 0;
			}

			increment(a, true);
		}

		if (log != null)
			log.info("iter = %d, error = %f -> %f : %s : Continue = %b", getIteration(), oldError, newError,
					(result == 1) ? "worse"
							: "Delta = " + DoubleEquality.relativeError(oldError, newError) +
									((result == 0) ? " (negligible)" : ""),
					notSatisfied);
	}

	/**
	 * Bit mask to isolate the exponent field of a <code>double</code>.
	 */
	private static final long EXP_BIT_MASK = 0x7FF0000000000000L;

	/**
	 * The number of logical bits in the significand of a <code>double</code> number, including the implicit bit.
	 */
	private static final int SIGNIFICAND_WIDTH = 53;

	/**
	 * Bias used in representing a <code>double</code> exponent.
	 */
	private static final int EXP_BIAS = 1023;

	/**
	 * Returns unbiased exponent of a <code>double</code>.
	 *
	 * @param d
	 *            the double
	 * @return the exponent
	 */
	public static int getExponent(double d)
	{
		// Copied from java.lang.StrictMath.getExponent(double d)
		// This is only available in Java 1.6

		/*
		 * Bitwise convert d to long, mask out exponent bits, shift
		 * to the right and then subtract out double's bias adjust to
		 * get true exponent value.
		 */
		return (int) (((Double.doubleToRawLongBits(d) & EXP_BIT_MASK) >> (SIGNIFICAND_WIDTH - 1)) - EXP_BIAS);
	}

	/**
	 * Set the number of iterations that the fit has to improve by a negligible amount
	 *
	 * @param iterationLimit
	 *            the iterationLimit to set
	 */
	public void setIterationLimit(int iterationLimit)
	{
		this.iterationLimit = iterationLimit;
	}

	/**
	 * @return the iterationLimit
	 */
	public int getIterationLimit()
	{
		return iterationLimit;
	}

	/**
	 * @param significantDigits
	 *            the significant digits for a negligible change in error
	 */
	public void setSignificantDigits(int significantDigits)
	{
		this.significantDigits = significantDigits;
		maxRelativeError = DoubleEquality.getMaxRelativeError(significantDigits);
	}

	/**
	 * @return the significantDigits
	 */
	public int getSignificantDigits()
	{
		return significantDigits;
	}

	/**
	 * @return true if avoiding plateaus
	 */
	public boolean isAvoidPlateau()
	{
		return avoidPlateau;
	}

	/**
	 * During an optimisation the error can asymptote below the improvement required to achieve convergence. This
	 * can be avoided by setting the avoid plateau flag to true. In this case an additional check will be made for
	 * slowly improving error when more than half-way towards the maximum number of iterations. If 3 consecutive
	 * improvements in error are at the same level of floating point accuracy then convergence is achieved.
	 * <p>
	 * Note that no check is made that the amount of improvement in error is small (as is done using the standard
	 * significant digits check). Thus it is possible that a plateau will be detected when the error value is still
	 * significantly improving (but by the same amount each step). Chances of this will be minimised by using a maximum
	 * number of iterations approximately twice that which allows reasonable fits to converge. E.g. If typical fitted
	 * data converges within 10 iterations to 6 significant digits then the maximum iterations should be set to 20.
	 *
	 * @param avoidPlateau
	 *            Set to true to avoid plateaus
	 */
	public void setAvoidPlateau(boolean avoidPlateau)
	{
		this.avoidPlateau = avoidPlateau;
	}
}
