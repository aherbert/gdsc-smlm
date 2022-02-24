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

package uk.ac.sussex.gdsc.smlm.function.cspline;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.Precision;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.data.DoubleStackTrivalueProvider;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicInterpolatingFunction;
import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicInterpolator;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.smlm.function.StandardGradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.StandardGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.StandardValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;

@SuppressWarnings({"javadoc"})
public abstract class CubicSplineFunctionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(CubicSplineFunctionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  protected DoubleEquality eq = new DoubleEquality(1e-2, 1e-3);
  protected DoubleEquality eq2 = new DoubleEquality(1e-5, 1e-8);
  protected DoubleEquality eq3 = new DoubleEquality(1e-1, 1e-3); // For the Gaussian integral

  // Compute as per Numerical Recipes 5.7.
  // Approximate error accuracy in single precision: Ef
  // Step size for derivatives:
  // h ~ (Ef)^(1/3) * xc
  // xc is the characteristic scale over which x changes, assumed to be 1 (not x as per NR since x
  // is close to zero)
  protected double stepH = 0.0001; // (double) (Math.pow(1e-3f, 1.0 / 3));

  protected int[] testx = new int[] {4, 5, 6};
  protected int[] testy = new int[] {4, 5, 6};
  protected double[] testbackground = new double[] {0, 400};
  protected double[] testsignal1 = new double[] {15, 105};
  // Pick some to fall on the node boundaries as the second order
  // numerical gradients evaluate poorly on the node boundaries.
  protected double[] testcx1 = new double[] {5.3, 5.0};
  protected double[] testcy1 = new double[] {4.5, 5.2};
  protected double[] testcz1 = new double[] {-1.5, 1.1};
  protected double[] testsignal2 = new double[] {20, 50};
  protected double[] testcx2 = new double[] {4.8, 5.3};
  protected double[] testcy2 = new double[] {5.1, 4.9};
  protected double[] testcz2 = new double[] {-1.9, 0.7};

  // Different widths to test for non-square function evaluation
  protected int maxx = 8;
  protected int maxy = 9;
  protected double background = 50;
  protected CubicSplineFunction f1;
  protected CubicSplineFunction f1f;
  protected CubicSplineFunction f2 = null;
  protected CubicSplineFunction f2f = null;

  // Test Astigmatic Gaussian
  static final double gamma = 2;
  static final int zDepth = 5;
  protected QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

  static final CubicSplineData splineData;
  static final CubicSplineData splineDataFloat;
  static final double cx;
  static final double cy;
  static final double cz;
  static final int scale;

  static {
    // Create a Gaussian PSF twice the size of the test Gaussian for interpolation
    scale = 2;
    final QuadraticAstigmatismZModel zModel =
        new QuadraticAstigmatismZModel(scale * gamma, scale * zDepth);
    final int size = 40;
    final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
        GaussianFunctionFactory.FIT_ASTIGMATISM, zModel);
    final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    a[Gaussian2DFunction.SIGNAL] = 1;
    a[Gaussian2DFunction.X_POSITION] = size / scale;
    a[Gaussian2DFunction.Y_POSITION] = size / scale;
    a[Gaussian2DFunction.X_SD] = scale;
    a[Gaussian2DFunction.Y_SD] = scale;

    // Create the Gaussian data for different z-depths
    final int minz = -scale * zDepth;
    final int maxz = -minz;
    final double[][] val = new double[maxz - minz + 1][];
    final StandardValueProcedure p = new StandardValueProcedure();
    for (int z = minz, i = 0; z <= maxz; z++, i++) {
      a[CubicSplineFunction.Z_POSITION] = z;
      val[i] = p.getValues(f, a);
    }
    final DoubleStackTrivalueProvider fval = new DoubleStackTrivalueProvider(val, size, size);

    // Utils.display("CubicSplineData", val, size, size);

    //@formatter:off
    final CustomTricubicInterpolatingFunction function = new CustomTricubicInterpolator.Builder()
        // The axis value are ignored ...
        .setIntegerAxisValues(true)
        .setFValue(fval)
        .interpolate();
    //@formatter:on
    splineData = new CubicSplineData(function);
    cx = a[CubicSplineFunction.X_POSITION];
    cy = a[CubicSplineFunction.Y_POSITION];
    cz = splineData.getMaxZ() / scale;

    function.toSinglePrecision();
    splineDataFloat = new CubicSplineData(function);
  }

  /**
   * Instantiates a new cubic spline function test.
   */
  public CubicSplineFunctionTest() {
    init();

    // Setup Tests
    if (!f1.evaluatesBackground()) {
      testbackground = new double[] {testbackground[0]};
    }
    if (!f1.evaluatesSignal()) {
      testsignal1 = new double[] {testsignal1[0]};
      testsignal2 = new double[] {testsignal2[0]};
    }
    // XY Position is always evaluated
    if (!f1.evaluatesZ()) {
      testcz1 = new double[] {0};
      testcz2 = new double[] {0};
    }

    postInit();
  }

  /**
   * Create the CubicSplineFunction for 1 and 2 peaks. Creates the flags for the factory
   */
  protected abstract void init();

  protected void postInit() {
    // To be overridden
  }

  @Test
  void functionCreatesCorrectGradientIndices() {
    checkGradientIndices(1, f1);
    checkGradientIndices(2, f2);
  }

  private static void checkGradientIndices(int npeaks, CubicSplineFunction cf) {
    if (cf == null) {
      return;
    }

    final int[] gradientIndices = cf.gradientIndices();
    logger.log(TestLogUtils.getRecord(Level.INFO, "Function%d %s %s", npeaks,
        cf.getClass().getName(), Arrays.toString(gradientIndices)));

    Assertions.assertEquals(cf.getN(), npeaks, "Incorrect number of peaks");

    int count = 0;
    if (cf.evaluatesBackground()) {
      Assertions.assertEquals(0, gradientIndices[count++], "Background");
    }
    for (int peak = 1, i = 1; peak <= npeaks;
        peak++, i += CubicSplineFunction.PARAMETERS_PER_PEAK) {
      final int ii = i;
      if (cf.evaluatesSignal()) {
        Assertions.assertEquals(i, gradientIndices[count++], () -> CubicSplineFunction.getName(ii));
      }
      if (cf.evaluatesPosition()) {
        Assertions.assertEquals(i + 1, gradientIndices[count++],
            () -> CubicSplineFunction.getName(ii + 1));
        Assertions.assertEquals(i + 2, gradientIndices[count++],
            () -> CubicSplineFunction.getName(ii + 2));
      }
      if (cf.evaluatesZ()) {
        Assertions.assertEquals(i + 3, gradientIndices[count++],
            () -> CubicSplineFunction.getName(ii + 3));
      }
    }
  }

  @Test
  void factoryCreatesCorrectFunction() {
    CubicSplineFunction func;

    if (f2 != null) {
      func = CubicSplineFunctionFactory.createCubicSplineFunction(splineData, maxx, maxy, cx, cy,
          cz, 2, 2);
      Assertions.assertTrue(func.getClass() == f2.getClass(), "Incorrect function2");
    } else {
      func = CubicSplineFunctionFactory.createCubicSplineFunction(splineData, maxx, maxy, cx, cy,
          cz, 2, 1);
      Assertions.assertTrue(func.getClass() == f1.getClass(), "Incorrect function1");
    }
  }

  @Test
  void functionComputesTargetWithAndWithoutGradient() {
    final StandardValueProcedure p0 = new StandardValueProcedure();
    final StandardGradient1Procedure p1 = new StandardGradient1Procedure();
    final StandardGradient2Procedure p2 = new StandardGradient2Procedure();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              final double[] a = createParameters(background, signal1, cx1, cy1, cz1);

              final double[] e = p0.getValues(f1, a);
              final double[] o1 = p1.getValues(f1, a);
              final double[] o2 = p2.getValues(f1, a);

              Assertions.assertArrayEquals(e, o1);
              Assertions.assertArrayEquals(e, o2);
              for (int i = e.length; i-- > 0;) {
                Assertions.assertArrayEquals(p1.gradients[i], p2.gradients1[i]);
              }
            }
          }
        }
      }
    }
  }

  @Test
  void functionComputesBackgroundGradient1() {
    Assumptions.assumeTrue(f1.evaluatesBackground());
    functionComputesTargetGradient1(CubicSplineFunction.BACKGROUND);
  }

  @Test
  void functionComputesSignalGradient1() {
    Assumptions.assumeTrue(f1.evaluatesSignal());
    functionComputesTargetGradient1(CubicSplineFunction.SIGNAL);
  }

  @Test
  void functionComputesXGradient1() {
    functionComputesTargetGradient1(CubicSplineFunction.X_POSITION);
  }

  @Test
  void functionComputesYGradient1() {
    functionComputesTargetGradient1(CubicSplineFunction.Y_POSITION);
  }

  @Test
  void functionComputesZGradient1() {
    Assumptions.assumeTrue(f1.evaluatesZ());
    functionComputesTargetGradient1(CubicSplineFunction.Z_POSITION);
  }

  private void functionComputesTargetGradient1(int targetParameter) {
    final int gradientIndex = findGradientIndex(f1, targetParameter);

    final Statistics s = new Statistics();

    final StandardValueProcedure p1a = new StandardValueProcedure();
    final StandardValueProcedure p1b = new StandardValueProcedure();
    final StandardGradient1Procedure p2 = new StandardGradient1Procedure();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              final double[] a = createParameters(background, signal1, cx1, cy1, cz1);

              // System.out.println(java.util.Arrays.toString(a));

              // Evaluate all gradients
              p2.getValues(f1, a);

              // Numerically solve gradient.
              // Calculate the step size h to be an exact numerical representation
              final double xx = a[targetParameter];

              // Get h to minimise roundoff error
              final double h = Precision.representableDelta(xx, stepH);

              // Evaluate at (x+h) and (x-h)
              a[targetParameter] = xx + h;
              p1a.getValues(f1, a);

              a[targetParameter] = xx - h;
              p1b.getValues(f1, a);

              // Only test close to the XY centre
              for (final int x : testx) {
                for (final int y : testy) {
                  final int i = y * maxx + x;
                  final double high = p1a.values[i];
                  final double low = p1b.values[i];

                  final double gradient = (high - low) / (2 * h);
                  final double dyda = p2.gradients[i][gradientIndex];
                  final double error = DoubleEquality.relativeError(gradient, dyda);
                  s.add(error);
                  if ((gradient * dyda) < 0) {
                    Assertions.fail(String.format("%s sign != %s", gradient, dyda));
                  }
                  // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f == [%d] %f? (%g)", x, y,
                  // gradient, gradientIndex, dyda, error);
                  if (!eq.almostEqualRelativeOrAbsolute(gradient, dyda)) {
                    Assertions.fail(String.format("%s != %s", gradient, dyda));
                  }
                }
              }
            }
          }
        }
      }
    }
    logger.info(() -> {
      return String.format("functionComputesTargetGradient1 %s %s (error %s +/- %s)",
          f1.getClass().getSimpleName(), CubicSplineFunction.getName(targetParameter),
          MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  protected int findGradientIndex(CubicSplineFunction func, int targetParameter) {
    final int index = func.findGradientIndex(targetParameter);
    Assertions.assertTrue(index >= 0, "Cannot find gradient index");
    return index;
  }

  @Test
  void functionComputesBackgroundGradient2() {
    Assumptions.assumeTrue(f1.evaluatesBackground());
    functionComputesTargetGradient2(CubicSplineFunction.BACKGROUND);
  }

  @Test
  void functionComputesSignalGradient2() {
    Assumptions.assumeTrue(f1.evaluatesSignal());
    functionComputesTargetGradient2(CubicSplineFunction.SIGNAL);
  }

  @Test
  void functionComputesXGradient2() {
    functionComputesTargetGradient2(CubicSplineFunction.X_POSITION);
  }

  @Test
  void functionComputesYGradient2() {
    functionComputesTargetGradient2(CubicSplineFunction.Y_POSITION);
  }

  @Test
  void functionComputesZGradient2() {
    Assumptions.assumeTrue(f1.evaluatesZ());
    functionComputesTargetGradient2(CubicSplineFunction.Z_POSITION);
  }

  private void functionComputesTargetGradient2(int targetParameter) {
    final int gradientIndex = findGradientIndex(f1, targetParameter);

    final Statistics s = new Statistics();

    final StandardGradient1Procedure p1a = new StandardGradient1Procedure();
    final StandardGradient1Procedure p1b = new StandardGradient1Procedure();
    final StandardGradient2Procedure p2 = new StandardGradient2Procedure();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              final double[] a = createParameters(background, signal1, cx1, cy1, cz1);

              // System.out.println(java.util.Arrays.toString(a));

              f1.initialise2(a);
              final boolean test = !f1.isNodeBoundary(gradientIndex);
              // Comment out when printing errors
              if (!test) {
                continue;
              }

              // Evaluate all gradients
              p2.getValues(f1, a);

              // Numerically solve gradient.
              // Calculate the step size h to be an exact numerical representation
              final double xx = a[targetParameter];

              // Get h to minimise roundoff error
              final double h = Precision.representableDelta(xx, stepH);

              // Evaluate at (x+h) and (x-h)
              a[targetParameter] = xx + h;
              p1a.getValues(f1, a);

              a[targetParameter] = xx - h;
              p1b.getValues(f1, a);

              // Only test close to the XY centre
              for (final int x : testx) {
                for (final int y : testy) {
                  final int i = y * maxx + x;
                  final double high = p1a.gradients[i][gradientIndex];
                  final double low = p1b.gradients[i][gradientIndex];

                  final double gradient = (high - low) / (2 * h);
                  final double d2yda2 = p2.gradients2[i][gradientIndex];
                  final double error = DoubleEquality.relativeError(gradient, d2yda2);
                  // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f == [%d] %f? (%g)", x, y,
                  // gradient, gradientIndex, d2yda2, error);
                  if (test) {
                    s.add(error);
                    Assertions.assertTrue((gradient * d2yda2) >= 0,
                        () -> String.format("%s sign != %s", gradient, d2yda2));
                    // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f == [%d] %f? (%g)", x, y,
                    // gradient, gradientIndex, d2yda2, error);
                    Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(gradient, d2yda2),
                        () -> String.format("%s != %s", gradient, d2yda2));
                  }
                }
              }
            }
          }
        }
      }
    }
    logger.info(() -> {
      return String.format("functionComputesTargetGradient2 %s %s (error %s +/- %s)",
          f1.getClass().getSimpleName(), CubicSplineFunction.getName(targetParameter),
          MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  @Test
  void functionComputesTargetWithAndWithoutGradientWith2Peaks() {
    if (f2 == null) {
      return;
    }

    final StandardValueProcedure p0 = new StandardValueProcedure();
    final StandardGradient1Procedure p1 = new StandardGradient1Procedure();
    final StandardGradient2Procedure p2 = new StandardGradient2Procedure();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              // Peak 2
              for (final double signal2 : testsignal2) {
                for (final double cx2 : testcx2) {
                  for (final double cy2 : testcy2) {
                    for (final double cz2 : testcz2) {
                      final double[] a = createParameters(background, signal1, cx1, cy1, cz1,
                          signal2, cx2, cy2, cz2);

                      final double[] e = p0.getValues(f1, a);
                      final double[] o1 = p1.getValues(f1, a);
                      final double[] o2 = p2.getValues(f1, a);

                      Assertions.assertArrayEquals(e, o1);
                      Assertions.assertArrayEquals(e, o2);
                      for (int i = e.length; i-- > 0;) {
                        Assertions.assertArrayEquals(p1.gradients[i], p2.gradients1[i]);
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
  void functionComputesBackgroundGradient1With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    Assumptions.assumeTrue(f2.evaluatesBackground());
    functionComputesTargetGradient1With2Peaks(CubicSplineFunction.BACKGROUND);
    functionComputesTargetGradient1With2Peaks(
        CubicSplineFunction.BACKGROUND + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesSignalGradient1With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    Assumptions.assumeTrue(f2.evaluatesSignal());
    functionComputesTargetGradient1With2Peaks(CubicSplineFunction.SIGNAL);
    functionComputesTargetGradient1With2Peaks(
        CubicSplineFunction.SIGNAL + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesXGradient1With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesTargetGradient1With2Peaks(CubicSplineFunction.X_POSITION);
    functionComputesTargetGradient1With2Peaks(
        CubicSplineFunction.X_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesYGradient1With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesTargetGradient1With2Peaks(CubicSplineFunction.Y_POSITION);
    functionComputesTargetGradient1With2Peaks(
        CubicSplineFunction.Y_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesZGradient1With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    Assumptions.assumeTrue(f2.evaluatesZ());
    functionComputesTargetGradient1With2Peaks(CubicSplineFunction.Z_POSITION);
    functionComputesTargetGradient1With2Peaks(
        CubicSplineFunction.Z_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  private void functionComputesTargetGradient1With2Peaks(int targetParameter) {
    final int gradientIndex = findGradientIndex(f2, targetParameter);

    final Statistics s = new Statistics();

    final StandardValueProcedure p1a = new StandardValueProcedure();
    final StandardValueProcedure p1b = new StandardValueProcedure();
    final StandardGradient1Procedure p2 = new StandardGradient1Procedure();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              // Peak 2
              for (final double signal2 : testsignal2) {
                for (final double cx2 : testcx2) {
                  for (final double cy2 : testcy2) {
                    for (final double cz2 : testcz2) {
                      final double[] a = createParameters(background, signal1, cx1, cy1, cz1,
                          signal2, cx2, cy2, cz2);

                      // System.out.println(java.util.Arrays.toString(a));

                      // Evaluate all gradients
                      p2.getValues(f2, a);

                      // Numerically solve gradient.
                      // Calculate the step size h to be an exact numerical representation
                      final double xx = a[targetParameter];

                      // Get h to minimise roundoff error
                      final double h = Precision.representableDelta(xx, stepH);

                      // Evaluate at (x+h) and (x-h)
                      a[targetParameter] = xx + h;
                      p1a.getValues(f2, a);

                      a[targetParameter] = xx - h;
                      p1b.getValues(f2, a);

                      // Only test close to the XY centre
                      for (final int x : testx) {
                        for (final int y : testy) {
                          final int i = y * maxx + x;
                          final double high = p1a.values[i];
                          final double low = p1b.values[i];

                          final double gradient = (high - low) / (2 * h);
                          final double dyda = p2.gradients[i][gradientIndex];
                          final double error = DoubleEquality.relativeError(gradient, dyda);
                          s.add(error);

                          Assertions.assertTrue((gradient * dyda) >= 0,
                              () -> String.format("%s sign != %s", gradient, dyda));
                          // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f == [%d] %f? (%g)", x,
                          // y, gradient, gradientIndex, dyda, error);
                          Assertions.assertTrue(eq.almostEqualRelativeOrAbsolute(gradient, dyda),
                              () -> String.format("%s != %s", gradient, dyda));
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
      return String.format("functionComputesTargetGradient1With2Peaks %s %s (error %s +/- %s)",
          f1.getClass().getSimpleName(), CubicSplineFunction.getName(targetParameter),
          MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  @Test
  void functionComputesBackgroundGradient2With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    Assumptions.assumeTrue(f2.evaluatesBackground());
    functionComputesTargetGradient2With2Peaks(CubicSplineFunction.BACKGROUND);
    functionComputesTargetGradient2With2Peaks(
        CubicSplineFunction.BACKGROUND + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesSignalGradient2With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    Assumptions.assumeTrue(f2.evaluatesSignal());
    functionComputesTargetGradient2With2Peaks(CubicSplineFunction.SIGNAL);
    functionComputesTargetGradient2With2Peaks(
        CubicSplineFunction.SIGNAL + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesXGradient2With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesTargetGradient2With2Peaks(CubicSplineFunction.X_POSITION);
    functionComputesTargetGradient2With2Peaks(
        CubicSplineFunction.X_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesYGradient2With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    functionComputesTargetGradient2With2Peaks(CubicSplineFunction.Y_POSITION);
    functionComputesTargetGradient2With2Peaks(
        CubicSplineFunction.Y_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  @Test
  void functionComputesZGradient2With2Peaks() {
    Assumptions.assumeTrue(null != f2);
    Assumptions.assumeTrue(f2.evaluatesZ());
    functionComputesTargetGradient1With2Peaks(CubicSplineFunction.Z_POSITION);
    functionComputesTargetGradient1With2Peaks(
        CubicSplineFunction.Z_POSITION + CubicSplineFunction.PARAMETERS_PER_PEAK);
  }

  private void functionComputesTargetGradient2With2Peaks(int targetParameter) {
    final int gradientIndex = findGradientIndex(f2, targetParameter);

    final Statistics s = new Statistics();

    final StandardGradient1Procedure p1a = new StandardGradient1Procedure();
    final StandardGradient1Procedure p1b = new StandardGradient1Procedure();
    final StandardGradient2Procedure p2 = new StandardGradient2Procedure();

    for (final double background : testbackground) {
      // Peak 1
      for (final double signal1 : testsignal1) {
        for (final double cx1 : testcx1) {
          for (final double cy1 : testcy1) {
            for (final double cz1 : testcz1) {
              // Peak 2
              for (final double signal2 : testsignal2) {
                for (final double cx2 : testcx2) {
                  for (final double cy2 : testcy2) {
                    for (final double cz2 : testcz2) {
                      final double[] a = createParameters(background, signal1, cx1, cy1, cz1,
                          signal2, cx2, cy2, cz2);

                      // System.out.println(java.util.Arrays.toString(a));

                      f2.initialise2(a);
                      final boolean test = !f2.isNodeBoundary(gradientIndex);
                      // Comment out when printing errors
                      if (!test) {
                        continue;
                      }

                      // Evaluate all gradients
                      p2.getValues(f2, a);

                      // Numerically solve gradient.
                      // Calculate the step size h to be an exact numerical representation
                      final double xx = a[targetParameter];

                      // Get h to minimise roundoff error
                      final double h = Precision.representableDelta(xx, stepH);

                      // Evaluate at (x+h) and (x-h)
                      a[targetParameter] = xx + h;
                      p1a.getValues(f2, a);

                      a[targetParameter] = xx - h;
                      p1b.getValues(f2, a);

                      // Only test close to the XY centre
                      for (final int x : testx) {
                        for (final int y : testy) {
                          final int i = y * maxx + x;
                          final double high = p1a.gradients[i][gradientIndex];
                          final double low = p1b.gradients[i][gradientIndex];

                          final double gradient = (high - low) / (2 * h);
                          final double d2yda2 = p2.gradients2[i][gradientIndex];
                          final double error = DoubleEquality.relativeError(gradient, d2yda2);
                          // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f == [%d] %f? (%g)", x,
                          // y, gradient, gradientIndex, d2yda2, error);
                          if (test) {
                            s.add(error);
                            if ((gradient * d2yda2) < 0) {
                              Assertions.fail(String.format("%s sign != %s", gradient, d2yda2));
                            }
                            // logger.fine(FunctionUtils.getSupplier("[%d,%d] %f == [%d] %f? (%g)",
                            // x, y, gradient, gradientIndex, d2yda2, error);
                            if (!eq.almostEqualRelativeOrAbsolute(gradient, d2yda2)) {
                              Assertions.fail(String.format("%s != %s", gradient, d2yda2));
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
      return String.format("functionComputesTargetGradient2With2Peaks %s %s (error %s +/- %s)",
          f1.getClass().getSimpleName(), CubicSplineFunction.getName(targetParameter),
          MathUtils.rounded(s.getMean()), MathUtils.rounded(s.getStandardDeviation()));
    });
  }

  protected double[] createParameters(double... args) {
    return args;
  }
}
