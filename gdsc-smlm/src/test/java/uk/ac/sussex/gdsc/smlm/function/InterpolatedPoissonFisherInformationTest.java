/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

@SuppressWarnings({"javadoc"})
class InterpolatedPoissonFisherInformationTest {
  @Test
  void canInterpolateFisherInformation() {
    for (int i = 0; i < 4; i++) {
      final double s = (1 << i) * 0.25;
      canInterpolateFisherInformation(s);
    }
  }

  private static void canInterpolateFisherInformation(double sd) {
    canComputeFisherInformation(new PoissonGaussianApproximationFisherInformation(sd));
  }

  private static void canComputeFisherInformation(BasePoissonFisherInformation fi) {
    // Build a range for the Fisher information
    final int min = -100;
    final int max = 20;
    final double[] logU = new double[max - min + 1];
    final double[] alpha = new double[logU.length];

    for (int exp = min, i = 0; exp <= max; exp++, i++) {
      logU[i] = exp;
      alpha[i] = fi.getAlpha(StdMath.exp(exp));
    }

    final InterpolatedPoissonFisherInformation intFi =
        new InterpolatedPoissonFisherInformation(logU, alpha, false, fi);

    // No check for beyond the range since a separate test does this.
    check(intFi, intFi, min, 5e-3);

    // Within
    for (int exp = min; exp < max; exp++) {
      check(intFi, intFi, exp + 0.5, 5e-3);
    }

    // Upper bound
    check(intFi, intFi, max, 5e-3);

    // Beyond upper bound
    check(intFi, intFi, max + 1, 5e-3);
  }

  private static void check(BasePoissonFisherInformation fi,
      InterpolatedPoissonFisherInformation intFi, double logU, double tol) {
    final double u = StdMath.exp(logU);
    final double e = fi.getAlpha(u);
    final double o = intFi.getAlpha(u);
    // logger.fine(FormatSupplier.getSupplier("logU=%g u=%g e=%g o=%g error=%g", logU, u, e, o,
    // DoubleEquality.relativeError(o, e));

    // Small numbers may have a large relative error but the absolute error is small
    // Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e, o, 5e-3, 1e-20));
    Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e, o, tol, 1e-20));
  }

  @Test
  void canInterpolateLowerFisherInformation() {
    for (int i = 0; i < 4; i++) {
      final double s = (1 << i) * 0.25;
      canInterpolateLowerFisherInformation(s);
    }
  }

  private static void canInterpolateLowerFisherInformation(double sd) {
    canComputeLowerFisherInformation(new PoissonGaussianApproximationFisherInformation(sd));
  }

  private static void canComputeLowerFisherInformation(BasePoissonFisherInformation fi) {
    // Build a range for the Fisher information where it should plateau
    final int min = -500;
    final int max = min + 4;
    final double[] logU = new double[max - min + 1];
    final double[] alpha = new double[logU.length];

    for (int exp = min, i = 0; exp <= max; exp++, i++) {
      logU[i] = exp;
      alpha[i] = fi.getAlpha(StdMath.exp(exp));
    }

    // Lower fixed I
    InterpolatedPoissonFisherInformation intFi =
        new InterpolatedPoissonFisherInformation(logU, alpha, true, fi);

    final double I = fi.getFisherInformation(StdMath.exp(min));
    final BasePoissonFisherInformation fixedI = new BasePoissonFisherInformation() {
      @Override
      public double getFisherInformation(double theta) {
        return I;
      }

      @Override
      public double getAlpha(double theta) {
        return theta * I;
      }

      @Override
      public BasePoissonFisherInformation copy() {
        return this;
      }
    };
    check(fixedI, intFi, min - 1, 0);
    check(fixedI, intFi, min - 2, 0);
    check(fixedI, intFi, min - 20, 0);

    // Lower fixed alpha
    intFi = new InterpolatedPoissonFisherInformation(logU, alpha, false, fi);

    final double A = alpha[0];
    final BasePoissonFisherInformation fixedA = new BasePoissonFisherInformation() {
      @Override
      public double getFisherInformation(double theta) {
        return theta / A;
      }

      @Override
      public double getAlpha(double theta) {
        return A;
      }

      @Override
      public BasePoissonFisherInformation copy() {
        return this;
      }
    };

    check(fixedA, intFi, min - 1, 0);
    check(fixedA, intFi, min - 2, 0);
    check(fixedA, intFi, min - 20, 0);
  }
}
