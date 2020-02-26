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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.function.FastLog;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Create a gradient procedure for use in the Levenberg–Marquardt (LVM) algorithm.
 */
public final class LvmGradientProcedureUtils {
  /** No public constructor. */
  private LvmGradientProcedureUtils() {}

  /**
   * The type of LVM gradient procedure.
   */
  public enum Type {
    /** Least-squares. */
    LSQ,
    /** Maximum Likelihood Estimation (using LVM). */
    MLE {
      @Override
      public boolean isMle() {
        return true;
      }
    },

    /** Weighted least-squares. */
    WLSQ,
    /** Fast Maximum Likelihood Estimation (using Newton iteration). */
    FAST_LOG_MLE {
      @Override
      public boolean isMle() {
        return true;
      }
    };

    /**
     * Checks if is MLE.
     *
     * @return true, if is MLE
     */
    public boolean isMle() {
      return false;
    }
  }

  /**
   * Create a new gradient calculator.
   *
   * @param y Data to fit
   * @param func Gradient function
   * @param type the type
   * @param fastLog the fast log
   * @return the gradient procedure
   */
  public static LvmGradientProcedure create(final double[] y, final Gradient1Function func,
      Type type, FastLog fastLog) {
    switch (type) {
      case WLSQ:
        // Do not support per observation weights
        return WLsqLvmGradientProcedureUtils.create(y, null, func);
      case MLE:
        return MleLvmGradientProcedureUtils.create(y, func);
      case LSQ:
        return LsqLvmGradientProcedureUtils.create(y, func);
      case FAST_LOG_MLE:
        return MleLvmGradientProcedureUtils.create(y, func, fastLog);
      default:
        throw new IllegalArgumentException("Unknown type: " + type);
    }
  }
}
