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
package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;

import org.apache.commons.math3.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public class InterpolatedPoissonFisherInformationTest {
  @Test
  public void canInterpolateFisherInformation() {
    for (int i = 0; i < 4; i++) {
      final double s = (1 << i) * 0.25;
      canInterpolateFisherInformation(s);
    }
  }

  private static void canInterpolateFisherInformation(double s) {
    canComputeFisherInformation(new PoissonGaussianApproximationFisherInformation(s));
  }

  private static void canComputeFisherInformation(BasePoissonFisherInformation f) {
    // Build a range for the Fisher information
    final int min = -100;
    final int max = 20;
    final double[] logU = new double[max - min + 1];
    final double[] alpha = new double[logU.length];

    for (int exp = min, i = 0; exp <= max; exp++, i++) {
      logU[i] = exp;
      alpha[i] = f.getAlpha(FastMath.exp(exp));
    }

    final InterpolatedPoissonFisherInformation fi =
        new InterpolatedPoissonFisherInformation(logU, alpha, false, f);

    // No check for beyond the range since a separate test does this.
    check(f, fi, min, 5e-3);

    // Within
    for (int exp = min; exp < max; exp++) {
      check(f, fi, exp + 0.5, 5e-3);
    }

    // Upper bound
    check(f, fi, max, 5e-3);

    // Beyond upper bound
    check(f, fi, max + 1, 5e-3);
  }

  private static void check(BasePoissonFisherInformation f, InterpolatedPoissonFisherInformation fi,
      double logU, double tol) {
    final double u = FastMath.exp(logU);
    final double e = f.getAlpha(u);
    final double o = fi.getAlpha(u);
    // logger.fine(FunctionUtils.getSupplier("logU=%g u=%g e=%g o=%g error=%g", logU, u, e, o,
    // DoubleEquality.relativeError(o, e));

    // Small numbers may have a large relative error but the absolute error is small
    // Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e, o, 5e-3, 1e-20));
    Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e, o, tol, 1e-20));
  }

  @Test
  public void canInterpolateLowerFisherInformation() {
    for (int i = 0; i < 4; i++) {
      final double s = (1 << i) * 0.25;
      canInterpolateLowerFisherInformation(s);
    }
  }

  private static void canInterpolateLowerFisherInformation(double s) {
    canComputeLowerFisherInformation(new PoissonGaussianApproximationFisherInformation(s));
  }

  private static void canComputeLowerFisherInformation(BasePoissonFisherInformation f) {
    // Build a range for the Fisher information where it should plateau
    final int min = -500;
    final int max = min + 4;
    final double[] logU = new double[max - min + 1];
    final double[] alpha = new double[logU.length];

    for (int exp = min, i = 0; exp <= max; exp++, i++) {
      logU[i] = exp;
      alpha[i] = f.getAlpha(FastMath.exp(exp));
    }

    // Lower fixed I
    InterpolatedPoissonFisherInformation fi =
        new InterpolatedPoissonFisherInformation(logU, alpha, true, f);

    final double I = f.getFisherInformation(FastMath.exp(min));
    final BasePoissonFisherInformation fixedI = new BasePoissonFisherInformation() {
      @Override
      public double getFisherInformation(double t) throws IllegalArgumentException {
        return I;
      }

      @Override
      public double getAlpha(double t) {
        return t * I;
      }

      @Override
      protected void postClone() {
        // Do nothing
      }
    };
    check(fixedI, fi, min - 1, 0);
    check(fixedI, fi, min - 2, 0);
    check(fixedI, fi, min - 20, 0);

    // Lower fixed alpha
    fi = new InterpolatedPoissonFisherInformation(logU, alpha, false, f);

    final double A = alpha[0];
    final BasePoissonFisherInformation fixedA = new BasePoissonFisherInformation() {
      @Override
      public double getFisherInformation(double t) throws IllegalArgumentException {
        return t / A;
      }

      @Override
      public double getAlpha(double t) {
        return A;
      }

      @Override
      protected void postClone() {
        // Do nothing
      }
    };

    check(fixedA, fi, min - 1, 0);
    check(fixedA, fi, min - 2, 0);
    check(fixedA, fi, min - 20, 0);
  }
}
