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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import java.util.logging.Level;
import org.apache.commons.math3.util.Precision;
import org.ejml.data.DenseMatrix64F;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.TurboList;
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
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingResult;
import uk.ac.sussex.gdsc.test.utils.TimingService;

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
  public void factoryDefaultsToErfGaussian2DFunction() {
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
  public void functionComputesSecondBackgroundGradient() {
    if (f1.evaluatesBackground()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.BACKGROUND);
    }
  }

  @Test
  public void functionComputesSecondSignalGradient() {
    if (f1.evaluatesSignal()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.SIGNAL);
    }
  }

  @Test
  public void functionComputesSecondXGradient() {
    functionComputesSecondTargetGradient(Gaussian2DFunction.X_POSITION);
  }

  @Test
  public void functionComputesSecondYGradient() {
    functionComputesSecondTargetGradient(Gaussian2DFunction.Y_POSITION);
  }

  @Test
  public void functionComputesSecondZGradient() {
    if (f1.evaluatesZ()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.Z_POSITION);
    }
  }

  @Test
  public void functionComputesSecondXWidthGradient() {
    if (f1.evaluatesSD0()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.X_SD);
    }
  }

  @Test
  public void functionComputesSecondYWidthGradient() {
    if (f1.evaluatesSD1()) {
      functionComputesSecondTargetGradient(Gaussian2DFunction.Y_SD);
    }
  }

  @Test
  public void functionComputesSecondAngleGradient() {
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
    logger.info(() -> {
      return String.format("functionComputesSecondTargetGradient %s %s (error %s +/- %s)",
          f1.getClass().getSimpleName(), Gaussian2DFunction.getName(targetParameter),
          MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  @Test
  public void functionComputesSecondBackgroundGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesBackground()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.BACKGROUND);
    }
  }

  @Test
  public void functionComputesSecondSignalGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesSignal()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.SIGNAL);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.SIGNAL + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  public void functionComputesSecondXGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.X_POSITION);
    functionComputesSecondTargetGradientWith2Peaks(
        Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  public void functionComputesSecondYGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Y_POSITION);
    functionComputesSecondTargetGradientWith2Peaks(
        Gaussian2DFunction.Y_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  public void functionComputesSecondZGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesZ()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Z_POSITION);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.Z_POSITION + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  public void functionComputesSecondXWidthGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesSD0()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.X_SD);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.X_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  public void functionComputesSecondYWidthGradientWith2Peaks() {
    Assumptions.assumeTrue(null != f2);
    if (f2.evaluatesSD1()) {
      functionComputesSecondTargetGradientWith2Peaks(Gaussian2DFunction.Y_SD);
      functionComputesSecondTargetGradientWith2Peaks(
          Gaussian2DFunction.Y_SD + Gaussian2DFunction.PARAMETERS_PER_PEAK);
    }
  }

  @Test
  public void functionComputesSecondAngleGradientWith2Peaks() {
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
    logger.info(() -> {
      return String.format("functionComputesSecondTargetGradient %s [%d] %s (error %s +/- %s)",
          f2.getClass().getSimpleName(), Gaussian2DFunction.getPeak(targetParameter),
          Gaussian2DFunction.getName(targetParameter), MathUtils.rounded(s.getMean()),
          MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  private class FunctionTimingTask extends BaseTimingTask {
    Gaussian2DFunction f1;
    ErfGaussian2DFunction f2;
    double[][] x;
    int order;
    final double[] dyda;
    final double[] d2yda2;
    final int size = f1.size();

    public FunctionTimingTask(Gaussian2DFunction f1, double[][] x, int order) {
      super(f1.getClass().getSimpleName() + " " + order + " eval");
      this.f1 = f1;
      if (order == 2) {
        f2 = (ErfGaussian2DFunction) f1;
      }
      this.x = x;
      this.order = order;
      dyda = new double[f1.getNumberOfGradients()];
      d2yda2 = new double[f1.getNumberOfGradients()];
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int index) {
      return null;
    }

    @Override
    public Object run(Object data) {
      double sum = 0;
      f1 = f1.copy();
      if (order == 0) {
        for (int i = 0; i < x.length; i++) {
          f1.initialise0(x[i]);
          for (int j = 0; j < size; j++) {
            sum += f1.eval(j);
          }
        }
      } else if (order == 1) {
        for (int i = 0; i < x.length; i++) {
          f1.initialise1(x[i]);
          for (int j = 0; j < size; j++) {
            sum += f1.eval(j, dyda);
          }
        }
      } else {
        for (int i = 0; i < x.length; i++) {
          f2.initialise2(x[i]);
          for (int j = 0; j < size; j++) {
            sum += f2.eval2(j, dyda, d2yda2);
          }
        }
      }
      return sum;
    }
  }

  // Speed test verses equivalent Gaussian2DFunction
  @SpeedTag
  @Test
  public void functionIsFasterThanEquivalentGaussian2DFunction() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final int flags = this.flags & ~GaussianFunctionFactory.FIT_ERF;
    final Gaussian2DFunction gf = GaussianFunctionFactory.create2D(1, maxx, maxy, flags, zModel);

    final boolean zDepth = (flags & GaussianFunctionFactory.FIT_Z) != 0;

    final TurboList<double[]> params1 = new TurboList<>();
    final TurboList<double[]> params2 = new TurboList<>();
    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              for (final double[] w1 : testw1) {
                for (final double angle1 : testangle1) {
                  double[] params =
                      createParameters(background, signal1, cx1, cy1, cz1, w1[0], w1[1], angle1);
                  params1.add(params);
                  if (zDepth) {
                    // Change to a standard free circular function
                    params = params.clone();
                    params[Gaussian2DFunction.X_SD] *=
                        zModel.getSx(params[Gaussian2DFunction.Z_POSITION]);
                    params[Gaussian2DFunction.Y_SD] *=
                        zModel.getSy(params[Gaussian2DFunction.Z_POSITION]);
                    params[Gaussian2DFunction.Z_POSITION] = 0;
                    params2.add(params);
                  }
                }
              }
            }
          }
        }
      }
    }
    final double[][] x = params1.toArray(new double[params1.size()][]);
    final double[][] x2 = (zDepth) ? params2.toArray(new double[params2.size()][]) : x;

    final int runs = 10000 / x.length;
    final TimingService ts = new TimingService(runs);
    ts.execute(new FunctionTimingTask(gf, x2, 1));
    ts.execute(new FunctionTimingTask(gf, x2, 0));
    ts.execute(new FunctionTimingTask(f1, x, 2));
    ts.execute(new FunctionTimingTask(f1, x, 1));
    ts.execute(new FunctionTimingTask(f1, x, 0));

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport());
    }

    for (int i = 1; i <= 2; i++) {
      final TimingResult slow = ts.get(-i - 3);
      final TimingResult fast = ts.get(-i);
      logger.log(TestLogUtils.getTimingRecord(slow, fast));
    }
  }

  @Test
  public void functionComputesGradientForEach() {
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
  public void functionComputesExtendedGradientForEach() {
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
                            logger.log(TestLogUtils.getRecord(Level.INFO, "%d [%d,%d] %f ?= %f",
                                index, j, k, gradient, m.get(j, k)));
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
  public void functionComputesGradientForEachWith2Peaks() {
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
  public void functionComputesExtendedGradientForEachWith2Peaks() {
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
                                      // logger.log(TestLog.getRecord(Level.FINE,
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

  private static class ForEachTimingTask extends BaseTimingTask {
    double[][] x;
    SimpleProcedure procedure;

    public ForEachTimingTask(ErfGaussian2DFunction func, double[][] x, int order) {
      super(func.getClass().getSimpleName() + " " + order + " forEach");
      this.x = x;
      if (order == 0) {
        procedure = new Procedure0(func);
      } else if (order == 1) {
        procedure = new Procedure1(func);
      } else {
        procedure = new Procedure2(func);
      }
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int index) {
      return null;
    }

    @Override
    public Object run(Object data) {
      procedure.reset();
      for (int i = 0; i < x.length; i++) {
        procedure.run(x[i]);
      }
      return procedure.sum;
    }
  }

  // Speed test forEach verses equivalent eval() function calls
  @SpeedTag
  @Test
  public void functionIsFasterUsingForEach() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final ErfGaussian2DFunction f1 = (ErfGaussian2DFunction) this.f1;

    final TurboList<double[]> params = new TurboList<>();
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
                  params.add(a);
                }
              }
            }
          }
        }
      }
    }
    final double[][] x = params.toArray(new double[params.size()][]);

    final int runs = 10000 / x.length;
    final TimingService ts = new TimingService(runs);
    ts.execute(new FunctionTimingTask(f1, x, 2));
    ts.execute(new FunctionTimingTask(f1, x, 1));
    ts.execute(new FunctionTimingTask(f1, x, 0));
    ts.execute(new ForEachTimingTask(f1, x, 2));
    ts.execute(new ForEachTimingTask(f1, x, 1));
    ts.execute(new ForEachTimingTask(f1, x, 0));

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport());
    }

    for (int i = 1; i <= 3; i++) {
      final TimingResult slow = ts.get(-i - 3);
      final TimingResult fast = ts.get(-i);
      logger.log(TestLogUtils.getTimingRecord(slow, fast));
    }
  }

  @Test
  public void functionCanComputeIntegral() {
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
  public void computeIntegralIsFaster() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final TurboList<double[]> p = new TurboList<>();
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
        s1 += new IntegralValueProcedure().getIntegral(f1, p.getf(j));
      }
    }
    long t2 = System.nanoTime();
    for (int i = n; i-- > 0;) {
      for (int j = p.size(); j-- > 0;) {
        s2 += f1.integral(p.getf(j));
      }
    }
    final long t3 = System.nanoTime();
    t1 = t2 - t1;
    t2 = t3 - t2;
    logger.log(TestLogUtils.getRecord(Level.INFO, "computeIntegralIsFaster %s %d vs %d (%gx)",
        f1.getClass().getSimpleName(), t1, t2, (double) t1 / t2));
    TestAssertions.assertTest(s1, s2, TestHelper.doublesAreClose(1e-3, 0));
    Assertions.assertTrue(t2 < t1);
  }

  @Test
  public void functionCanComputeIntegralWith2Peaks() {
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

  @Test
  public void computeIntegralIsFasterWith2Peaks() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    Assumptions.assumeTrue(null != f2);

    final TurboList<double[]> p = new TurboList<>();
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
                              p.add(createParameters(background, signal1, cx1, cy1, cz1, w1[0],
                                  w1[1], angle1, signal2, cx2, cy2, cz2, w2[0], w2[1], angle2));
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
    final int n = (int) Math.ceil(10000.0 / p.size());
    double s1 = 0;
    double s2 = 0;
    long t1 = System.nanoTime();
    for (int i = n; i-- > 0;) {
      for (int j = p.size(); j-- > 0;) {
        s1 += new IntegralValueProcedure().getIntegral(f2, p.getf(j));
      }
    }
    long t2 = System.nanoTime();
    for (int i = n; i-- > 0;) {
      for (int j = p.size(); j-- > 0;) {
        s2 += f2.integral(p.getf(j));
      }
    }
    final long t3 = System.nanoTime();
    t1 = t2 - t1;
    t2 = t3 - t2;
    logger.log(
        TestLogUtils.getRecord(Level.INFO, "computeIntegralIsFasterWith2Peaks %s %d vs %d (%gx)",
            f1.getClass().getSimpleName(), t1, t2, (double) t1 / t2));
    TestAssertions.assertTest(s1, s2, TestHelper.doublesAreClose(1e-3, 0));
    Assertions.assertTrue(t2 < t1);
  }
}
