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
package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.FastLog;
import gdsc.smlm.function.Gradient1Function;

/**
 * Create a gradient procedure for use in the Levenberg–Marquardt (LVM) algorithm
 */
public class LVMGradientProcedureFactory
{
	/**
	 * The type of LVM gradient procedure
	 */
	public enum Type
	{
		//@formatter:off
		/** Least-squares */
		LSQ,
		/** Maximum Likelihood Estimation (using LVM) */
		MLE { @Override public boolean isMLE()	{return true;} },
		/** Weighted least-squares */
		WLSQ,
		/** Fast Maximum Likelihood Estimation (using Newton iteration) */
		FastLogMLE { @Override public boolean isMLE()	{return true;} };
		//@formatter:on

		/**
		 * Checks if is MLE.
		 *
		 * @return true, if is MLE
		 */
		public boolean isMLE()
		{
			return false;
		}
	}

	/**
	 * Create a new gradient calculator.
	 *
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 * @param type
	 *            the type
	 * @param fastLog
	 *            the fast log
	 * @return the gradient procedure
	 */
	public static LVMGradientProcedure create(final double[] y, final Gradient1Function func, Type type,
			FastLog fastLog)
	{
		switch (type)
		{
			case WLSQ:
				// Do not support per observation weights
				return WLSQLVMGradientProcedureFactory.create(y, null, func);
			case MLE:
				return MLELVMGradientProcedureFactory.create(y, func);
			case LSQ:
				return LSQLVMGradientProcedureFactory.create(y, func);
			case FastLogMLE:
				return MLELVMGradientProcedureFactory.create(y, func, fastLog);
			default:
				throw new IllegalStateException("Unknown type: " + type);
		}
	}
}
