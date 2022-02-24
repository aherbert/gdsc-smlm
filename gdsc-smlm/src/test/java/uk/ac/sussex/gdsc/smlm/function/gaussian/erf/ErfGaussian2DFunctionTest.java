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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.IntegralValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunctionTest;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

@SuppressWarnings({"javadoc"})
public abstract class ErfGaussian2DFunctionTest extends Gaussian2DFunctionTest {

  /**
   * Instantiates a new erf gaussian 2 D function test.
   */
  public ErfGaussian2DFunctionTest() {
    super();
    // The derivative check can be tighter with the ERF since it is a true integration
    stepH = 0.0001;
    eq3 = new DoubleEquality(5e-3, 1e-3); // For the Gaussian integral
  }

  @Test
  void factoryDefaultsToErfGaussian2DFunction() {
    Gaussian2DFunction f1;
    Gaussian2DFunction f2;

    final int flags2 = BitFlagUtils.unset(flags, GaussianFunctionFactory.FIT_ERF);

    if (this.f2 != null) {
      f1 = GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
      f2 = GaussianFunctionFactory.create2D(2, maxx, maxy, flags2, zModel);
      Assertions.assertTrue(f1.getClass() == f2.getClass(), "Incorrect function2");
    } else {
      f1 = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
      f2 = GaussianFunctionFactory.create2D(1, maxx, maxy, flags2, zModel);
      Assertions.assertTrue(f1.getClass() == f2.getClass(), "Incorrect function1");
    }
  }

  @Test
  void functionComputesSecondBackgroundGradient() {
    if (f1.evaluatesBackground()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.BACKGROUND);
    }
  }

  @Test
  void functionComputesSecondSignalGradient() {
    if (f1.evaluatesSignal()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.SIGNAL);
    }
  }

  @Test
  void functionComputesSecondXGradient() {
    functionComputesSecondTargetGradient(Gaussian2DFunction.X_POSITION);
  }

  @Test
  void functionComputesSecondYGradient() {
    functionComputesSecondTargetGradient(Gaussian2DFunction.Y_POSITION);
  }

  @Test
  void functionComputesSecondZGradient() {
    if (f1.evaluatesZ()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.Z_POSITION);
    }
  }

  @Test
  void functionComputesSecondXWidthGradient() {
    if (f1.evaluatesSD0()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.X_SD);
    }
  }

  @Test
  void functionComputesSecondYWidthGradient() {
    if (f1.evaluatesSD1()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.Y_SD);
    }
  }

  @Test
  void functionComputesSecondAngleGradient() {
    if (f1.evaluatesAngle()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.ANGLE);
    }
  }

  private void functionComputesSecondTargetGradient(int targetParameter) {
    final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;
    final int gradientIndex = findGradientIndex(f1, targetParameter);
    final double[] dyda = new double[f1.getNumberOfGradients()];
    final double[] dyda2 = new double[dyda.length];
    double[] params;

    // Test fitting of second derivatives
    final ErfGaussian2DFunction f1a =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
    final ErfGaussian2DFunction f1b =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);
    final Statistics s = new Statistics();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  params =
                      createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
                  f1.initialise2(params);

                  // Numerically solve gradient.
                  // Calculate the step size h to be an exact numerical representation
                  final double xx = params[targetParameter];

                  // Get h to minimise roundoff error
                  final double h = Precision.representableDelta(xx, stepH);

                  // Evaluate at (x+h) and (x-h)
                  params[targetParameter] = xx + h;
                  f1a.initialise1(params.clone());

                  params[targetParameter] = xx - h;
                  f1b.initialise1(params.clone());

                  for (final int x : testx) {
                    for (final int y : testy) {
                      final int i = y * maxx + x;
                      f1a.eval(i, dyda);
                      final double value2 = dyda[gradientIndex];
                      f1b.eval(i, dyda);
                      final double value3 = dyda[gradientIndex];
                      f1.eval2(i, dyda, dyda2);

                      final double gradient = (value2 - value3) / (2 * h);
                      final double error =
                          DoubleEquality.relativeError(gradient, dyda2[gradientIndex]);
                      s.add(error);
                      Assertions.assertTrue((gradient * dyda2[gradientIndex]) >= 0,
                          () -> gradient + " sign != " + dyda2[gradientIndex]);
                      // TestLog.fine(logger,"[%d,%d] %f == [%d] %f? (%g)", x, y, gradient,
                      // gradientIndex, dyda2[gradientIndex], error);
                      Assertions.assertTrue(
                          eq.almostEqualRelativeOrAbsolute(gradient, dyda2[gradientIndex]),
                          () -> gradient + " != " + dyda2[gradientIndex]);

                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    logger.log(TestLevel.TEST_INFO, () -> {
      return String.format("functionComputesSecondTargetGradient %s %s (error %s +/- %s)",
          f1.getClass().getSimpleName(), Gaussian2DFunction.getName(targetParameter),
          MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  @Test
  void functionComputesSecondBackgroundGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesBackground()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.BACKGROUND);
    }
  }

  @Test
  void functionComputesSecondSignalGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesSignal()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.SIGNAL + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  void functionComputesSecondXGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION);
    functionComputesSecondTargetGradientWith2Peaks(
        Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesSecondYGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION);
    functionComputesSecondTargetGradientWith2Peaks(
        Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesSecondZGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesZ()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Z_POSITION);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.Z_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  void functionComputesSecondXWidthGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesSD0()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.X_SD);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.X_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  void functionComputesSecondYWidthGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesSD1()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.Y_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  void functionComputesSecondAngleGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesAngle()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.ANGLE);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.ANGLE + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  private void functionComputesSecondTargetGradientWith2Peaks(int targetParameter) {
    final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;
    final int gradientIndex = findGradientIndex(f2, targetParameter);
    final double[] dyda = new double[f2.getNumberOfGradients()];
    final double[] dyda2 = new double[dyda.length];
    double[] params;

    // Test fitting of second derivatives
    final ErfGaussian2DFunction f2a =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
    final ErfGaussian2DFunction f2b =
        (ErfGaussian2DFunction) GaussianFunctionFactory.create2D(2, maxx, maxy, flags, zModel);
    final Statistics s = new Statistics();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  // Peak 2
                  for (final double signal2 : testsignal2) {
                    for (final double cx2 : testcx2) {
                      for (final double cy2 : testcy2) {
                        for (final double cz2 : testcz2) {
                          for (final double[] w2 : testw2) {
                            for (final double angle2 : testangle2) {
                              params = createParameters(background, signal1, cx1, cy1, cz1, w1[0],
                                  w1[1], angle1, signal2, cx2, cy2, cz2, w2[0], w2[1], angle2);
                              f2.initialise2(params);

                              // Numerically solve gradient.
                              // Calculate the step size h to be an exact numerical representation
                              final double xx = params[targetParameter];

                              // Get h to minimise roundoff error
                              final double h = Precision.representableDelta(xx, stepH);

                              // Evaluate at (x+h) and (x-h)
                              params[targetParameter] = xx + h;
                              f2a.initialise1(params.clone());

                              params[targetParameter] = xx - h;
                              f2b.initialise1(params.clone());

                              for (final int x : testx) {
                                for (final int y : testy) {
                                  final int i = y * maxx + x;
                                  f2a.eval(i, dyda);
                                  final double value2 = dyda[gradientIndex];
                                  f2b.eval(i, dyda);
                                  final double value3 = dyda[gradientIndex];
                                  f2.eval2(i, dyda, dyda2);

                                  final double gradient = (value2 - value3) / (2 * h);
                                  final double error =
                                      DoubleEquality.relativeError(gradient, dyda2[gradientIndex]);
                                  s.add(error);
                                  Assertions.assertTrue(

                                      (gradient * dyda2[gradientIndex]) >= 0,
                                      () -> gradient + " sign != " + dyda2[gradientIndex]);
                                  // TestLog.fine(logger,"[%d,%d] %f == [%d] %f? (%g)", x, y,
                                  // gradient,
                                  // gradientIndex, dyda2[gradientIndex], error);
                                  Assertions.assertTrue(

                                      eq.almostEqualRelativeOrAbsolute(gradient,
                                          dyda2[gradientIndex]),
                                      () -> gradient + " != " + dyda2[gradientIndex]);
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    logger.log(TestLevel.TEST_INFO, () -> {
      return String.format("functionComputesSecondTargetGradient %s [%d] %s (error %s +/- %s)",
          f2.getClass().getSimpleName(), Gaussian2DFunction.getPeak(targetParameter),
          Gaussian2DFunction.getName(targetParameter), MathUtils.rounded(s.getMean()),
          MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  @Test
  void functionComputesGradientForEach() {
    final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

    final int n = f1.size();
    final double[] du_da = new double[f1.getNumberOfGradients()];
    final double[] du_db = new double[f1.getNumberOfGradients()];
    final double[] d2u_da2 = new double[f1.getNumberOfGradients()];

    final double[] values = new double[n];
    final double[][] jacobian = new double[n][];
    final double[][] jacobian2 = new double[n][];

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  final double[] a =
                      createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
                  f1.initialiseExtended2(a);

                  // Compute single
                  for (int i = 0; i < n; i++) {
                    final double o1 = f1.eval(i, du_da);
                    final double o2 = f1.eval2(i, du_db, d2u_da2);
                    Assertions.assertEquals(o1, o2, 1e-10, "Value");
                    Assertions.assertArrayEquals(du_da, du_db, 1e-10, "Jacobian!=Jacobian");
                    values[i] = o1;
                    jacobian[i] = du_da.clone();
                    jacobian2[i] = d2u_da2.clone();
                  }

                  // Use procedures
                  f1.forEach(new ValueProcedure() {
                    int index = 0;

                    @Override
                    public void execute(double value) {
                      Assertions.assertEquals(values[index], value, 1e-10, "Value ValueProcedure");
                      index++;
                    }
                  });

                  f1.forEach(new Gradient1Procedure() {
                    int index = 0;

                    @Override
                    public void execute(double value, double[] dyDa) {
                      Assertions.assertEquals(values[index], value, 1e-10,
                          "Value Gradient1Procedure");
                      Assertions.assertArrayEquals(jacobian[index], dyDa, 1e-10,
                          "du_da Gradient1Procedure");
                      index++;
                    }
                  });

                  f1.forEach(new Gradient2Procedure() {
                    int index = 0;

                    @Override
                    public void execute(double value, double[] dyDa, double[] d2yDa2) {
                      Assertions.assertEquals(values[index], value, 1e-10,
                          "Value Gradient2Procedure");
                      Assertions.assertArrayEquals(jacobian[index], dyDa, 1e-10,
                          "du_da Gradient2Procedure");
                      Assertions.assertArrayEquals(jacobian2[index], d2yDa2, 1e-10,
                          "d2u_da2 Gradient2Procedure");
                      index++;
                    }
                  });

                  f1.forEach(new ExtendedGradient2Procedure() {
                    int index = 0;

                    @Override
                    public void executeExtended(double value, double[] dyDa, double[] d2yDaDb) {
                      Assertions.assertEquals(values[index], value, 1e-10,
                          "Value ExtendedGradient2Procedure");
                      Assertions.assertArrayEquals(jacobian[index], dyDa, 1e-10,
                          "du_da ExtendedGradient2Procedure");
                      for (int j = 0, k = 0; j < d2u_da2.length; j++, k += d2u_da2.length + 1) {
                        d2u_da2[j] = d2yDaDb[k];
                      }
                      Assertions.assertArrayEquals(jacobian2[index], d2u_da2, 1e-10,
                          "d2u_da2 Gradient2Procedure");
                      index++;
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }

  @Test
  void functionComputesExtendedGradientForEach() {
    final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

    final int nparams = f1.getNumberOfGradients();
    final int[] gradientIndices = f1.gradientIndices();

    final ErfGaussian2DFunction[] fHigh = new ErfGaussian2DFunction[nparams];
    final ErfGaussian2DFunction[] fLow = new ErfGaussian2DFunction[nparams];
    final double[] delta = new double[nparams];
    for (int j = 0; j < nparams; j++) {
      fHigh[j] = f1.copy();
      fLow[j] = f1.copy();
    }

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  final double[] a =
                      createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
                  f1.initialiseExtended2(a);

                  // Create a set of functions initialised +/- delta in each parameter
                  for (int j = 0; j < nparams; j++) {
                    final int targetParameter = gradientIndices[j];
                    // Numerically solve gradient.
                    // Calculate the step size h to be an exact numerical representation
                    final double xx = a[targetParameter];

                    // Get h to minimise roundoff error
                    final double h = Precision.representableDelta(xx, stepH);

                    // Evaluate at (x+h) and (x-h)
                    a[targetParameter] = xx + h;
                    fHigh[j].initialise1(a.clone());

                    a[targetParameter] = xx - h;
                    fLow[j].initialise1(a.clone());

                    a[targetParameter] = xx;
                    delta[j] = 2 * h;
                  }

                  f1.forEach(new ExtendedGradient2Procedure() {
                    int index = -1;
                    final double[] duDa = new double[f1.getNumberOfGradients()];
                    final double[] duDb = new double[f1.getNumberOfGradients()];

                    @Override
                    public void executeExtended(double value, double[] dyDa, double[] d2yDaDb) {
                      index++;
                      final DenseMatrix64F m = DenseMatrix64F.wrap(nparams, nparams, d2yDaDb);
                      for (int j = 0; j < nparams; j++) {
                        // Evaluate the function +/- delta for parameter j
                        fHigh[j].eval(index, duDa);
                        fLow[j].eval(index, duDb);
                        // Check the gradient with respect to parameter k
                        for (int k = 0; k < nparams; k++) {
                          final double gradient = (duDa[k] - duDb[k]) / delta[j];
                          final boolean ok =
                              eq.almostEqualRelativeOrAbsolute(gradient, m.get(j, k));
                          if (!ok) {
                            logger.log(TestLogUtils.getRecord(TestLevel.TEST_INFO,
                                "%d [%d,%d] %f ?= %f", index, j, k, gradient, m.get(j, k)));
                            Assertions.fail(String.format("%d [%d,%d] %f != %f", index, j, k,
                                gradient, m.get(j, k)));
                          }
                        }
                      }
                    }
                  });
                }
              }
            }
          }
        }
      }
    }
  }

  @Test
  void functionComputesGradientForEachWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;

    final int n = f2.size();
    final double[] du_da = new double[f2.getNumberOfGradients()];
    final double[] du_db = new double[f2.getNumberOfGradients()];
    final double[] d2u_da2 = new double[f2.getNumberOfGradients()];

    final double[] values = new double[n];
    final double[][] jacobian = new double[n][];
    final double[][] jacobian2 = new double[n][];

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  // Peak 2
                  for (final double signal2 : testsignal2) {
                    for (final double cx2 : testcx2) {
                      for (final double cy2 : testcy2) {
                        for (final double cz2 : testcz2) {
                          for (final double[] w2 : testw2) {
                            for (final double angle2 : testangle2) {
                              final double[] a =
                                  createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1],
                                      angle1, signal2, cx2, cy2, cz2, w2[0], w2[1], angle2);
                              f2.initialiseExtended2(a);

                              // Compute single
                              for (int i = 0; i < n; i++) {
                                final double o1 = f2.eval(i, du_da);
                                final double o2 = f2.eval2(i, du_db, d2u_da2);
                                Assertions.assertEquals(o1, o2, 1e-10, "Value");
                                Assertions.assertArrayEquals(du_da, du_db, 1e-10,
                                    "Jacobian!=Jacobian");
                                values[i] = o1;
                                jacobian[i] = du_da.clone();
                                jacobian2[i] = d2u_da2.clone();
                              }

                              // Use procedures
                              f2.forEach(new ValueProcedure() {
                                int index = 0;

                                @Override
                                public void execute(double value) {
                                  Assertions.assertEquals(values[index], value, 1e-10,
                                      "Value ValueProcedure");
                                  index++;
                                }
                              });

                              f2.forEach(new Gradient1Procedure() {
                                int index = 0;

                                @Override
                                public void execute(double value, double[] dyDa) {
                                  Assertions.assertEquals(values[index], value, 1e-10,
                                      "Value Gradient1Procedure");
                                  Assertions.assertArrayEquals(jacobian[index], dyDa, 1e-10,
                                      "du_da Gradient1Procedure");
                                  index++;
                                }
                              });

                              f2.forEach(new Gradient2Procedure() {
                                int index = 0;

                                @Override
                                public void execute(double value, double[] dyDa, double[] d2yDa2) {
                                  Assertions.assertEquals(values[index], value, 1e-10,
                                      "Value Gradient2Procedure");
                                  Assertions.assertArrayEquals(jacobian[index], dyDa, 1e-10,
                                      "du_da Gradient2Procedure");
                                  Assertions.assertArrayEquals(jacobian2[index], d2yDa2, 1e-10,
                                      "d2u_da2 Gradient2Procedure");
                                  index++;
                                }
                              });

                              f2.forEach(new ExtendedGradient2Procedure() {
                                int index = 0;

                                @Override
                                public void executeExtended(double value, double[] dyDa,
                                    double[] d2yDaDb) {
                                  Assertions.assertEquals(

                                      values[index], value, 1e-10,
                                      "Value ExtendedGradient2Procedure");
                                  Assertions.assertArrayEquals(

                                      jacobian[index], dyDa, 1e-10,
                                      "du_da ExtendedGradient2Procedure");
                                  for (int j = 0, k = 0; j < d2u_da2.length;
                                      j++, k += d2u_da2.length + 1) {
                                    d2u_da2[j] = d2yDaDb[k];
                                  }
                                  Assertions.assertArrayEquals(jacobian2[index], d2u_da2, 1e-10,
                                      "d2u_da2 Gradient2Procedure");
                                  index++;
                                }
                              });
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  @Test
  void functionComputesExtendedGradientForEachWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    final ErfGaussian2DFunction f2 = (ErfGaussian2DFunction) this.f2;

    final int nparams = f2.getNumberOfGradients();
    final int[] gradientIndices = f2.gradientIndices();

    final ErfGaussian2DFunction[] fHigh = new ErfGaussian2DFunction[nparams];
    final ErfGaussian2DFunction[] fLow = new ErfGaussian2DFunction[nparams];
    final double[] delta = new double[nparams];
    for (int j = 0; j < nparams; j++) {
      fHigh[j] = f2.copy();
      fLow[j] = f2.copy();
    }

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  // Peak 2
                  for (final double signal2 : testsignal2) {
                    for (final double cx2 : testcx2) {
                      for (final double cy2 : testcy2) {
                        for (final double cz2 : testcz2) {
                          for (final double[] w2 : testw2) {
                            for (final double angle2 : testangle2) {
                              final double[] a =
                                  createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1],
                                      angle1, signal2, cx2, cy2, cz2, w2[0], w2[1], angle2);

                              f2.initialiseExtended2(a);

                              // Create a set of functions initialised +/- delta in each parameter
                              for (int j = 0; j < nparams; j++) {
                                final int targetParameter = gradientIndices[j];
                                // Numerically solve gradient.
                                // Calculate the step size h to be an exact numerical representation
                                final double xx = a[targetParameter];

                                // Get h to minimise roundoff error
                                final double h = Precision.representableDelta(xx, stepH);

                                // Evaluate at (x+h) and (x-h)
                                a[targetParameter] = xx + h;
                                fHigh[j].initialise1(a.clone());

                                a[targetParameter] = xx - h;
                                fLow[j].initialise1(a.clone());

                                a[targetParameter] = xx;
                                delta[j] = 2 * h;
                              }

                              f2.forEach(new ExtendedGradient2Procedure() {
                                int index = -1;
                                final double[] duDa = new double[f2.getNumberOfGradients()];
                                final double[] duDb = new double[f2.getNumberOfGradients()];

                                @Override
                                public void executeExtended(double value, double[] dyDa,
                                    double[] d2yDaDb) {
                                  index++;
                                  // if (i!=f2.size()/2) return;
                                  final DenseMatrix64F m =
                                      DenseMatrix64F.wrap(nparams, nparams, d2yDaDb);
                                  for (int j = 0; j < nparams; j++) {
                                    // Evaluate the function +/- delta for parameter j
                                    fHigh[j].eval(index, duDa);
                                    fLow[j].eval(index, duDb);
                                    // Check the gradient with respect to parameter k
                                    for (int k = 0; k < nparams; k++) {
                                      final double gradient = (duDa[k] - duDb[k]) / delta[j];
                                      final boolean ok =
                                          eq.almostEqualRelativeOrAbsolute(gradient, m.get(j, k));
                                      // logger.log(TestLog.getRecord(TestLevel.TEST_DEBUG,
                                      // "%d [%d,%d] %f ?= %f", i, j, k,
                                      // gradient, m.get(j, k)));
                                      if (!ok) {
                                        Assertions.fail(String.format("%d [%d,%d] %f != %f", index,
                                            j, k, gradient, m.get(j, k)));
                                      }
                                    }
                                  }
                                }
                              });
                              // return;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  abstract static class SimpleProcedure {
    ErfGaussian2DFunction func;
    double sum = 0;

    SimpleProcedure(ErfGaussian2DFunction func) {
      this.func = func;
    }

    void reset() {
      sum = 0;
    }

    void run(double[] a) {
      func = func.copy();
      initialise(a);
      forEach();
    }

    abstract void initialise(double[] a);

    abstract void forEach();
  }

  static class Procedure0 extends SimpleProcedure implements ValueProcedure {
    Procedure0(ErfGaussian2DFunction func) {
      super(func);
    }

    @Override
    void initialise(double[] a) {
      func.initialise0(a);
    }

    @Override
    void forEach() {
      func.forEach(this);
    }

    @Override
    public void execute(double value) {
      sum += value;
    }
  }

  static class Procedure1 extends SimpleProcedure implements Gradient1Procedure {
    Procedure1(ErfGaussian2DFunction func) {
      super(func);
    }

    @Override
    void initialise(double[] a) {
      func.initialise1(a);
    }

    @Override
    void forEach() {
      func.forEach(this);
    }

    @Override
    public void execute(double value, double[] dyDa) {
      sum += value;
    }
  }

  static class Procedure2 extends SimpleProcedure implements Gradient2Procedure {
    Procedure2(ErfGaussian2DFunction func) {
      super(func);
    }

    @Override
    void initialise(double[] a) {
      func.initialise2(a);
    }

    @Override
    void forEach() {
      func.forEach(this);
    }

    @Override
    public void execute(double value, double[] dyDa, double[] d2yDa2) {
      sum += value;
    }
  }

  @Test
  void functionCanComputeIntegral() {
    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-8, 0);
    double[] params;
    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  params =
                      createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
                  final double e = new IntegralValueProcedure().getIntegral(f1, params);
                  final double o = f1.integral(params);
                  TestAssertions.assertTest(e, o, predicate);
                }
              }
            }
          }
        }
      }
    }
  }

  @Test
  void computeIntegralIsFaster() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final LocalList<double[]> p = new LocalList<>();
    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  p.add(createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1));
                }
              }
            }
          }
        }
      }
    }
    final int n = (int) Math.ceil(10000.0 / p.size());
    double s1 = 0;
    double s2 = 0;
    long t1 = System.nanoTime();
    for (int i = n; i-- > 0;) {
      for (int j = p.size(); j-- > 0;) {
        s1 += new IntegralValueProcedure().getIntegral(f1, p.unsafeGet(j));
      }
    }
    long t2 = System.nanoTime();
    for (int i = n; i-- > 0;) {
      for (int j = p.size(); j-- > 0;) {
        s2 += f1.integral(p.unsafeGet(j));
      }
    }
    final long t3 = System.nanoTime();
    t1 = t2 - t1;
    t2 = t3 - t2;
    logger.log(
        TestLogUtils.getRecord(TestLevel.TEST_INFO, "computeIntegralIsFaster %s %d vs %d (%gx)",
            f1.getClass().getSimpleName(), t1, t2, (double) t1 / t2));
    TestAssertions.assertTest(s1, s2, TestHelper.doublesAreClose(1e-3, 0));
    Assertions.assertTrue(t2 < t1);
  }

  @Test
  void functionCanComputeIntegralWith2Peaks() {
    Assumptions.assumeTrue(null != f2);

    final DoubleDoubleBiPredicate predicate = TestHelper.doublesAreClose(1e-8, 0);
    double[] params;
    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  // Peak 2
                  for (final double signal2 : testsignal2) {
                    for (final double cx2 : testcx2) {
                      for (final double cy2 : testcy2) {
                        for (final double cz2 : testcz2) {
                          for (final double[] w2 : testw2) {
                            for (final double angle2 : testangle2) {
                              params = createParameters(background, signal1, cx1, cy1, cz1, w1[0],
                                  w1[1], angle1, signal2, cx2, cy2, cz2, w2[0], w2[1], angle2);
                              final double e = new IntegralValueProcedure().getIntegral(f2, params);
                              final double o = f2.integral(params);
                              TestAssertions.assertTest(e, o, predicate);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
