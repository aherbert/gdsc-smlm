/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.smlm.function.IcsiFastLog.DataType;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.function.Supplier;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"unused", "javadoc"})
public class FastLogTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(FastLogTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  IcsiFastLog icsiLog = IcsiFastLog.create(DataType.BOTH);
  FFastLog ffastLog = new FFastLog();
  DFastLog dfastLog = new DFastLog();
  TurboLog turboLog = new TurboLog();

  //@formatter:off
  private static class MathLog extends FastLog {
    @Override
    public double getBase()  { return Math.E;}
    @Override
    public double getScale() { return LN2; }
    @Override
    public int getN() { return 52; }

    @Override
    public float log(float x) {  return (float) Math.log(x);  }
    @Override
    public float log2(float x) { return (float) (log(x) / LN2);  }
    @Override
    public float log(double x) { return (float) Math.log(x); }
    @Override
    public float log2(double x) { return (float) (log(x) / LN2); }
    @Override
    public double logD(double x) { return Math.log(x); }
    @Override
    public double  log2D(double x) { return (log(x) / LN2); }

    @Override
    public float fastLog(float x) { return log(x); }
    @Override
    public float fastLog2(float x) { return log2(x); }
    @Override
    public float fastLog(double x) { return log(x);  }
    @Override
    public float fastLog2(double x) { return log2(x); }
    @Override
    public double fastLogD(double x) { return log(x);  }
    @Override
    public double fastLog2D(double x) { return log2(x); }
  }
  private static class FastMathLog extends MathLog {
    @Override
    public float log(float x) {  return (float) FastMath.log(x);  }
    @Override
    public float log(double x) { return (float) FastMath.log(x); }
    @Override
    public double logD(double x) { return FastMath.log(x); }
  }
  private abstract static class BaseTestLog {
    FastLog fl;
    String name;
    BaseTestLog(FastLog fl) {
      this.name = this.getClass().getSimpleName() + " " + fl.getClass().getSimpleName();
      this.fl=fl; }
    abstract float log(float x);
    abstract double log(double x);
    int getN() { return fl.getN(); }
  }
  private static class TestLog extends BaseTestLog {
    TestLog(FastLog fl) { super(fl); }
    @Override
    float log(float x) { return fl.log(x); }
    @Override
    double log(double x) { return fl.logD(x); }
  }
  private static class TestFastLog extends BaseTestLog {
    TestFastLog(FastLog fl) { super(fl); }
    @Override
    float log(float x) { return fl.fastLog(x); }
    @Override
    double log(double x) { return fl.fastLogD(x); }
  }
  // To test Math.log(1+x).
  // This is what is used in the MLE LVM gradient calculator
  private static class Test1PLog extends BaseTestLog {
    Test1PLog(FastLog fl) { super(fl); }
    @Override
    float log(float x) { return (float) Math.log(1+x); }
    @Override
    double log(double x) { return Math.log(1+x); }
  }
  private static class TestLog1P extends BaseTestLog {
    TestLog1P(FastLog fl) { super(fl); }
    @Override
    float log(float x) { return (float) Math.log1p(x); }
    @Override
    double log(double x) { return Math.log1p(x); }
  }
  private static class TestLog1PApache extends BaseTestLog {
    TestLog1PApache(FastLog fl) { super(fl); }
    @Override
    float log(float x) { return (float) FastMath.log1p(x); }
    @Override
    double log(double x) { return FastMath.log1p(x); }
  }
  //@formatter:on

  @Test
  public void canComputeFFastLog_fastLog() {
    canComputeLog(new TestFastLog(ffastLog), false);
  }

  @Test
  public void canComputeFFastLog_log() {
    canComputeLog(new TestLog(ffastLog), true);
  }

  @Test
  public void canComputeDFastLog_fastLog() {
    canComputeLog(new TestFastLog(dfastLog), false);
  }

  @Test
  public void canComputeDFastLog_log() {
    canComputeLog(new TestLog(dfastLog), true);
  }

  @Test
  public void canComputeICSCFastLog_fastLog() {
    canComputeLog(new TestFastLog(icsiLog), false);
  }

  @Test
  public void canComputeICSCFastLog_log() {
    canComputeLog(new TestLog(icsiLog), true);
  }

  @Test
  public void canComputeTurboLog_fastLog() {
    canComputeLog(new TestFastLog(turboLog), false);
  }

  @Test
  public void canComputeTurboLog_log() {
    canComputeLog(new TestLog(turboLog), true);
  }

  private static void canComputeLog(BaseTestLog fl, boolean edgeCases) {
    testLog(fl, Float.NaN, edgeCases);
    testLog(fl, Float.NEGATIVE_INFINITY, edgeCases);
    testLog(fl, -Float.MAX_VALUE, edgeCases);
    testLog(fl, -Float.MIN_VALUE, edgeCases);
    testLog(fl, -2f, edgeCases);
    testLog(fl, -1f, edgeCases);
    testLog(fl, -0f, edgeCases);
    testLog(fl, 0f, edgeCases);
    testLog(fl, Float.MIN_VALUE, false); // Not enough precision to test this
    testLog(fl, 1e-10f, true);
    testLog(fl, 1f, true);
    testLog(fl, 2f, true);
    testLog(fl, 2048f, true);
    testLog(fl, Float.MAX_VALUE, true);
    testLog(fl, Float.POSITIVE_INFINITY, edgeCases);
  }

  private static void testLog(BaseTestLog fl, float value, boolean test) {
    final float e = (float) Math.log(value);
    final float o = fl.log(value);
    final float error = FloatEquality.relativeError(e, o);
    logger.log(uk.ac.sussex.gdsc.test.utils.TestLogUtils.getRecord(Level.INFO,
        "%s v=%g : %fl vs %s (%g)", fl.name, value, e, o, error));
    if (test) {
      if (Double.isNaN(e) && Double.isNaN(o)) {
        return;
      }
      if (e == o) {
        return;
      }
      Assertions.assertTrue(error < 1e-4f);
    }
  }

  @Test
  public void canComputeDoubleFFast_fastLog() {
    canComputeDoubleLog(new TestFastLog(ffastLog), false);
  }

  @Test
  public void canComputeDoubleFFastLog_log() {
    canComputeDoubleLog(new TestLog(ffastLog), true);
  }

  @Test
  public void canComputeDoubleDFast_fastLog() {
    canComputeDoubleLog(new TestFastLog(dfastLog), false);
  }

  @Test
  public void canComputeDoubleDFastLog_log() {
    canComputeDoubleLog(new TestLog(dfastLog), true);
  }

  @Test
  public void canComputeDoubleICSCFast_fastLog() {
    canComputeDoubleLog(new TestFastLog(icsiLog), false);
  }

  @Test
  public void canComputeDoubleICSCFastLog_log() {
    canComputeDoubleLog(new TestLog(icsiLog), true);
  }

  @Test
  public void canComputeDoubleTurbo_fastLog() {
    canComputeDoubleLog(new TestFastLog(turboLog), false);
  }

  @Test
  public void canComputeDoubleTurboLog_log() {
    canComputeDoubleLog(new TestLog(turboLog), true);
  }

  private static void canComputeDoubleLog(BaseTestLog fl, boolean edgeCases) {
    testDoubleLog(fl, Double.NaN, edgeCases);
    testDoubleLog(fl, Double.NEGATIVE_INFINITY, edgeCases);
    testDoubleLog(fl, -Double.MAX_VALUE, edgeCases);
    testDoubleLog(fl, -Double.MIN_VALUE, edgeCases);
    testDoubleLog(fl, -2d, edgeCases);
    testDoubleLog(fl, -1d, edgeCases);
    testDoubleLog(fl, -0d, edgeCases);
    testDoubleLog(fl, 0d, edgeCases);
    testDoubleLog(fl, Double.MIN_VALUE, false); // Not enough precision to test this
    testDoubleLog(fl, 1e-10, true);
    testDoubleLog(fl, 1d, true);
    testDoubleLog(fl, 2d, true);
    testDoubleLog(fl, 2048d, true);
    testDoubleLog(fl, Double.MAX_VALUE / 2, true);
    testDoubleLog(fl, Double.MAX_VALUE, true);
    testDoubleLog(fl, Double.POSITIVE_INFINITY, edgeCases);
  }

  private static void testDoubleLog(BaseTestLog fl, double value, boolean test) {
    final double e = Math.log(value);
    final double o = fl.log(value);
    final double error = DoubleEquality.relativeError(e, o);
    logger.log(uk.ac.sussex.gdsc.test.utils.TestLogUtils.getRecord(Level.INFO,
        "%s v=%g : %fl vs %s (%g)", fl.name, value, e, o, error));
    if (test) {
      if (Double.isNaN(e) && Double.isNaN(o)) {
        return;
      }
      if (e == o) {
        return;
      }
      Assertions.assertTrue(error < 1e-4);
    }
  }

  // A robust error test using a uniform random number, report min,av,max,sd of the error.
  // The error 'should' always be negative as the truncation rounds down. However the table
  // pre-computes using an exponent offset which can lead to rounding up.

  @SeededTest
  public void canTestFloatError(RandomSeed seed) {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    // All float values is a lot so we do a representative set
    final float[] d = generateRandomFloats(seed, 1000000);
    final float[] logD = new float[d.length];
    for (int i = 0; i < d.length; i++) {
      logD[i] = (float) Math.log(d[i]);
    }

    final int min = 0;
    final int max = 23;
    // int min = 13, max = 13;

    // for (int n = min; n <= max; n++)
    // {
    // canTestFloatError(new TestFastLog(ICSIFastLog.create(n, DataType.FLOAT)), d, logD);
    // }
    // for (int n = min; n <= max; n++)
    // {
    // canTestFloatError(new TestFastLog(new FFastLog(n)), d, logD);
    // }
    for (int n = min; n <= max; n++) {
      runCanTestFloatError(new TestFastLog(new DFastLog(n)), d, logD);
    }
    for (int n = min; n <= max; n++) {
      runCanTestFloatError(new TestFastLog(new TurboLog(n)), d, logD);
    }
    for (int n = min; n <= max; n++) {
      runCanTestFloatError(new TestFastLog(new TurboLog2(n)), d, logD);
    }
  }

  private static float[] generateRandomFloats(RandomSeed seed, int n) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());
    final float[] d = new float[n];
    for (int i = 0; i < d.length; i++) {
      d[i] = nextUniformFloat(rng);
    }
    return d;
  }

  private static float nextUniformFloat(UniformRandomProvider rng) {
    int uniform = rng.nextInt();
    // Mask out sign and the last bit of the exponent (avoid infinity and NaN)
    uniform = BitFlagUtils.unset(uniform, 0x80000000 | 0x00800000);
    // assert ((u >> 23) & 0xff) < 255;
    return Float.intBitsToFloat(uniform);
  }

  private static double nextUniformDouble(UniformRandomProvider rng) {
    long uniform = rng.nextLong();
    // Mask out sign and the last bit of the exponent (avoid infinity and NaN)
    uniform &= ~(0x8000000000000000L | 0x0010000000000000L);
    // assert ((u >> 52) & 0x7ffL) < 2047;
    return Double.longBitsToDouble(uniform);
  }

  @Test
  public void canTestFloatErrorRange() {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    final TurboList<TestFastLog> test = new TurboList<>();
    final int n = 13;
    test.add(new TestFastLog(IcsiFastLog.create(n, DataType.FLOAT)));
    test.add(new TestFastLog(new FFastLog(n)));
    test.add(new TestFastLog(new DFastLog(n)));
    test.add(new TestFastLog(new TurboLog(n)));
    test.add(new TestFastLog(new TurboLog2(n)));

    // Full range in blocks.
    // Only when the number is around 1 or min value are there significant errors
    final float[] d = null;
    final float[] logD = null;

    // All
    // testFloatErrorRange(test, n, d, logD, 0, 255, 0);

    // Only a problem around min value and x==1
    // testFloatErrorRange(test, n, d, logD, 0, 2, 0);
    testFloatErrorRange(test, n, d, logD, 125, 130, 0);
    // testFloatErrorRange(test, n, d, logD, 253, 255, 0);
  }

  private static void testFloatErrorRange(TurboList<TestFastLog> test, int n, float[] data,
      float[] logD, int mine, int maxe, int ee) {
    for (int e = mine; e < maxe; e += ee + 1) {
      data = generateFloats(e, e + ee, data);
      if (logD == null || logD.length < data.length) {
        logD = new float[data.length];
      }
      for (int i = 0; i < data.length; i++) {
        logD[i] = (float) Math.log(data[i]);
      }
      logger.info(getSupplier("e=%d-%d", e, e + ee));
      for (final TestFastLog fl : test) {
        runCanTestFloatError(fl, data, logD);
      }
    }
  }

  private static float[] generateFloats(int mine, int maxe, float[] data) {
    // Mantissa = 23-bit, Exponent = 8-bit
    final int mbits = 23;
    mine = MathUtils.clip(0, 255, mine);
    maxe = MathUtils.clip(0, 255, maxe);
    if (mine > maxe) {
      throw new IllegalStateException();
    }
    final int mn = (1 << mbits);
    final int n = mn * (maxe - mine + 1);
    if (data == null || data.length < n) {
      data = new float[n];
    }
    int index = 0;
    for (int m = 0; m < mn; m++) {
      for (int e = mine; e <= maxe; e++) {
        final int bits = m | (e << 23);
        final float v = Float.intBitsToFloat(bits);
        // logger.fine(FunctionUtils.getSupplier("%g = %s", v, Integer.toBinaryString(bits));
        data[index++] = v;
      }
    }
    return data;
  }

  private static class FPair {
    int index;
    float value;
  }

  private static class DPair {
    int index;
    double value;
  }

  private static class Stats {
    double min;
    double max;
    double sum;
    double ss;
    double minv;
    double maxv;
    int count = 1;

    Stats(double exp, double value) {
      min = max = sum = exp;
      minv = maxv = value;
      ss = exp * exp;
    }

    void add(double exp, double value) {
      if (min > exp) {
        min = exp;
        minv = value;
      } else if (max < exp) {
        max = exp;
        maxv = value;
      }
      sum += exp;
      ss += exp * exp;
      count++;
    }

    double getMean() {
      return sum / count;
    }

    double getSD() {
      final double sd = ss - (sum * sum) / count;
      if (sd > 0.0) {
        return Math.sqrt(sd / (count - 1));
      }
      return 0.0;
    }

    String summary() {
      return String.format("min=%s (%s), max=%s (%s), mean=%s, sd=%g", min, minv, max, maxv,
          getMean(), getSD());
    }
  }

  private static void runCanTestFloatError(BaseTestLog fl, float[] data, float[] logD) {
    final FPair pair = new FPair();
    if (!next(fl, pair, data)) {
      return;
    }

    float value = logD[pair.index - 1];
    double delta = value - pair.value;
    delta = Math.abs(delta);
    // if (delta > 1)
    // {
    // //logger.fine(FunctionUtils.getSupplier("Big error: %fl %fl", v, d[pair.i-1]);
    // }
    final Stats s1 = new Stats(delta, data[pair.index - 1]);
    final Stats s2 = (value != 0) ? new Stats(Math.abs(delta / value), data[pair.index - 1])
        : new Stats(0, data[pair.index - 1]);
    while (next(fl, pair, data)) {
      value = logD[pair.index - 1];
      delta = value - pair.value;
      delta = Math.abs(delta);
      // if (delta > 5)
      // {
      // //logger.fine(FunctionUtils.getSupplier("Big error: [%g] %fl %fl %fl", d[pair.i - 1], v,
      // pair.fl, v));
      // }
      s1.add(delta, data[pair.index - 1]);
      if (value != 0) {
        s2.add(Math.abs(delta / value), data[pair.index - 1]);
      }
    }
    logger.info(getSupplier("%s, n=%d, c=%d : %s : relative %s", fl.name, fl.getN(), s1.count,
        s1.summary(), s2.summary()));
  }

  private static boolean next(BaseTestLog fl, FPair pair, float[] data) {
    while (pair.index < data.length) {
      final float x = data[pair.index++];
      if (x == 0) {
        continue;
      }
      pair.value = fl.log(x);
      if (pair.value != Float.NEGATIVE_INFINITY) {
        return true;
        // logger.fine(FunctionUtils.getSupplier("%g", x);
      }
    }
    return false;
  }

  private static boolean next(BaseTestLog fl, DPair pair, double[] data) {
    while (pair.index < data.length) {
      final double x = data[pair.index++];
      if (x == 0) {
        continue;
      }
      pair.value = fl.log(x);
      if (pair.value != Double.NEGATIVE_INFINITY) {
        return true;
        // logger.fine(FunctionUtils.getSupplier("%g", x);
      }
    }
    return false;
  }

  @SeededTest
  public void canTestDoubleError(RandomSeed seed) {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    // All float values is a lot so we do a representative set
    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());
    final double lower = Double.MIN_VALUE;
    final double upper = Double.MAX_VALUE;
    final double[] d = new double[10000000];
    final double[] logD = new double[d.length];
    for (int i = 0; i < d.length; i++) {
      final double v = nextUniformDouble(rng);
      d[i] = v;
      logD[i] = Math.log(v);
    }

    // int min = 0, max = 23;
    final int min = 4;
    final int max = 13;

    // for (int n = min; n <= max; n++)
    // {
    // canTestDoubleError(new TestFastLog(ICSIFastLog.create(n, DataType.DOUBLE)), d, logD);
    // }
    // for (int n = min; n <= max; n++)
    // {
    // canTestDoubleError(new TestFastLog(new DFastLog(n)), d, logD);
    // }
    // for (int n = min; n <= max; n++)
    // {
    // canTestDoubleError(new TestFastLog(new FFastLog(n)), d, logD);
    // }
    for (int n = min; n <= max; n++) {
      runCanTestDoubleError(new TestFastLog(new TurboLog(n)), d, logD);
      runCanTestDoubleError(new TestFastLog(new TurboLog2(n)), d, logD);
    }
  }

  @SeededTest
  public void canTestDoubleErrorLog1P(RandomSeed seed) {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    // All float values is a lot so we do a representative set
    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());
    final double lower = Double.MIN_VALUE;
    final double upper = Double.MAX_VALUE;
    final double[] d = new double[100000];
    final double[] logD = new double[d.length];
    for (int i = 0; i < d.length; i++) {
      final double v = nextUniformDouble(rng);
      d[i] = v;
      logD[i] = Math.log1p(v);
    }

    runCanTestDoubleError(new Test1PLog(new MathLog()), d, logD);
    runCanTestDoubleError(new TestLog1P(new MathLog()), d, logD);
  }

  @SeededTest
  public void canTestDoubleErrorRange(RandomSeed seed) {
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());

    final TurboList<TestFastLog> test = new TurboList<>();
    final int n = 13;
    test.add(new TestFastLog(IcsiFastLog.create(n, DataType.DOUBLE)));
    test.add(new TestFastLog(new FFastLog(n)));
    test.add(new TestFastLog(new DFastLog(n)));
    test.add(new TestFastLog(new TurboLog(n)));

    // Full range in blocks.
    // Only when the number is around 1 or min value are there significant errors
    final double[] d = new double[10000000];
    final double[] logD = null;

    // All
    // testDoubleErrorRange(test, n, d, logD, 0, 255, 0);

    // Only a problem around min value and x==1
    // testDoubleErrorRange(rng, test, n, d, logD, 0, 2, 0);
    testDoubleErrorRange(rng, test, n, d, logD, 1021, 1026, 0);
    testDoubleErrorRange(rng, test, n, d, logD, 2045, 2047, 0);
  }

  private static void testDoubleErrorRange(UniformRandomProvider rng, TurboList<TestFastLog> test,
      int n, double[] data, double[] logD, int mine, int maxe, int ee) {
    for (int e = mine; e < maxe; e += ee + 1) {
      data = generateDoubles(rng, e, e + ee, data);
      if (logD == null || logD.length < data.length) {
        logD = new double[data.length];
      }
      for (int i = 0; i < data.length; i++) {
        logD[i] = Math.log(data[i]);
      }
      logger.info(getSupplier("e=%d-%d", e, e + ee));
      for (final TestFastLog fl : test) {
        runCanTestDoubleError(fl, data, logD);
      }
    }
  }

  private static double[] generateDoubles(UniformRandomProvider rng, int mine, int maxe,
      double[] data) {
    // Mantissa = 52-bit, Exponent = 11-bit
    mine = MathUtils.clip(0, 2047, mine);
    maxe = MathUtils.clip(0, 2047, maxe);
    if (mine > maxe) {
      throw new IllegalStateException();
    }
    int index = 0;
    while (index < data.length) {
      // Only generate the mantissa
      final long m = rng.nextLong() & 0xfffffffffffffL;

      for (long e = mine; e <= maxe && index < data.length; e++) {
        final long bits = m | (e << 52);
        final double v = Double.longBitsToDouble(bits);
        // logger.fine(FunctionUtils.getSupplier("%g = %s", v, Long.toBinaryString(bits));
        data[index++] = v;
      }
    }
    return data;
  }

  private static void runCanTestDoubleError(BaseTestLog fl, double[] data, double[] logD) {
    final DPair pair = new DPair();
    if (!next(fl, pair, data)) {
      return;
    }

    double value = logD[pair.index - 1];
    double delta = value - pair.value;
    delta = Math.abs(delta);
    final Stats s1 = new Stats(delta, data[pair.index - 1]);
    final Stats s2 = (value != 0) ? new Stats(Math.abs(delta / value), data[pair.index - 1])
        : new Stats(0, data[pair.index - 1]);
    while (next(fl, pair, data)) {
      value = logD[pair.index - 1];
      // logger.fine(FunctionUtils.getSupplier("%g vs %g", v, pair.fl);
      delta = value - pair.value;
      delta = Math.abs(delta);
      s1.add(delta, data[pair.index - 1]);
      if (value != 0) {
        s2.add(Math.abs(delta / value), data[pair.index - 1]);
      }
    }
    logger.info(getSupplier("%s, n=%d, c=%d : %s : relative %s", fl.name, fl.getN(), s1.count,
        s1.summary(), s2.summary()));
  }

  // Speed test of float/double version verses the Math.log.
  private abstract static class DummyTimingTask extends BaseTimingTask {
    BaseTestLog log;
    int qvalue;

    public DummyTimingTask(String name, BaseTestLog log, int qvalue) {
      super(name + " " + log.name + " q=" + qvalue);
      this.log = log;
      this.qvalue = qvalue;
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int index) {
      return null;
    }
  }

  private static class FloatTimingTask extends DummyTimingTask {
    float[] values;

    public FloatTimingTask(BaseTestLog log, int qvalue, float[] values) {
      super("log(float)", log, qvalue);
      this.values = values;
    }

    @Override
    public Object run(Object data) {
      final float[] rng = new float[values.length];
      for (int i = 0; i < values.length; i++) {
        rng[i] = log.log(values[i]);
      }
      return rng;
    }
  }

  @SpeedTag
  @SeededTest
  public void canTestFloatSpeed(RandomSeed seed) {
    // No assertions, this is just a report
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());
    final float[] x = new float[1000000];
    for (int i = 0; i < x.length; i++) {
      x[i] = nextUniformFloat(rng);
    }

    final TimingService ts = new TimingService(5);
    ts.execute(new FloatTimingTask(new TestLog(new MathLog()), 0, x));
    ts.execute(new FloatTimingTask(new TestLog(new FastMathLog()), 0, x));
    for (final int q : new int[] {11}) {
      final int n = 23 - q;
      final IcsiFastLog fl = IcsiFastLog.create(n, DataType.FLOAT);
      ts.execute(new FloatTimingTask(new TestLog(fl), q, x));
      ts.execute(new FloatTimingTask(new TestFastLog(fl), q, x));
      final FFastLog ff = new FFastLog(n);
      ts.execute(new FloatTimingTask(new TestLog(ff), q, x));
      ts.execute(new FloatTimingTask(new TestFastLog(ff), q, x));
      final DFastLog df = new DFastLog(n);
      ts.execute(new FloatTimingTask(new TestLog(df), q, x));
      ts.execute(new FloatTimingTask(new TestFastLog(df), q, x));
      final TurboLog tf = new TurboLog(n);
      ts.execute(new FloatTimingTask(new TestLog(tf), q, x));
      ts.execute(new FloatTimingTask(new TestFastLog(tf), q, x));
      // TurboLog2 tf2 = new TurboLog2(n);
      // ts.execute(new FloatTimingTask(new TestLog(tf2), q, x));
      // ts.execute(new FloatTimingTask(new TestFastLog(tf2), q, x));
      // For the same precision we can reduce n
      final TurboLog2 tf3 = new TurboLog2(n - 1);
      ts.execute(new FloatTimingTask(new TestLog(tf3), q + 1, x));
      ts.execute(new FloatTimingTask(new TestFastLog(tf3), q + 1, x));
    }

    final int size = ts.getSize();
    ts.repeat(size);
    logger.info(ts.getReport(size));
  }

  private static class DoubleTimingTask extends DummyTimingTask {
    double[] values;

    public DoubleTimingTask(BaseTestLog log, int qvalue, double[] values) {
      super("log(double)", log, qvalue);
      this.values = values;
    }

    @Override
    public Object run(Object data) {
      final double[] rng = new double[values.length];
      for (int i = 0; i < values.length; i++) {
        rng[i] = log.log(values[i]);
      }
      return rng;
    }
  }

  private static class DoubleToFloatTimingTask extends DummyTimingTask {
    double[] values;
    float[] xf;

    public DoubleToFloatTimingTask(BaseTestLog log, int qvalue, double[] values, float[] xf) {
      super("log((float)double)", log, qvalue);
      this.values = values;
      this.xf = xf;
    }

    @Override
    public Object run(Object data) {
      final double[] rng = new double[values.length];
      for (int i = 0; i < values.length; i++) {
        rng[i] = log.log((float) values[i]);
      }
      // rng[i] = log.log(xf[i]);
      return rng;
    }
  }

  @SpeedTag
  @SeededTest
  public void canTestDoubleSpeed(RandomSeed seed) {
    // No assertions, this is just a report
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());
    final double[] values = new double[1000000];
    for (int i = 0; i < values.length; i++) {
      values[i] = nextUniformDouble(rng);
    }

    final TimingService ts = new TimingService(5);
    ts.execute(new DoubleTimingTask(new TestLog(new MathLog()), 0, values));
    ts.execute(new DoubleTimingTask(new TestLog(new FastMathLog()), 0, values));
    //// Test min acceptable precision
    // TurboLog2 tf3 = new TurboLog2(8);
    // ts.execute(new DoubleTimingTask(new TestLog(tf3), 15, x));
    // ts.execute(new DoubleTimingTask(new TestFastLog(tf3), 15, x));
    for (final int q : new int[] {11}) {
      final int n = 23 - q;
      final IcsiFastLog fl = IcsiFastLog.create(n, DataType.DOUBLE);
      ts.execute(new DoubleTimingTask(new TestLog(fl), q, values));
      ts.execute(new DoubleTimingTask(new TestFastLog(fl), q, values));
      final DFastLog df = new DFastLog(n);
      // ts.execute(new DoubleTimingTask(new DFastLog_log2(fl), q, x));
      // ts.execute(new DoubleTimingTask(new DTestFastLog2(fl), q, x));
      ts.execute(new DoubleTimingTask(new TestLog(df), q, values));
      ts.execute(new DoubleTimingTask(new TestFastLog(df), q, values));
      final TurboLog tf = new TurboLog(n);
      ts.execute(new DoubleTimingTask(new TestLog(tf), q, values));
      ts.execute(new DoubleTimingTask(new TestFastLog(tf), q, values));
      // Test same precision
      final TurboLog2 tf2 = new TurboLog2(n - 1);
      ts.execute(new DoubleTimingTask(new TestLog(tf2), q + 1, values));
      ts.execute(new DoubleTimingTask(new TestFastLog(tf2), q + 1, values));

      //// Min acceptable precision. This is usually the same speed
      //// showing the precomputed table is optimally used for moderate n.
      // ts.execute(new DoubleTimingTask(new TestLog(tf3), 15, x));
      // ts.execute(new DoubleTimingTask(new TestFastLog(tf3), 15, x));
    }

    final int size = ts.getSize();
    ts.repeat(size);
    logger.info(ts.getReport(size));
  }

  @SpeedTag
  @SeededTest
  public void canTestDoubleSpeedLog1P(RandomSeed seed) {
    // No assertions, this is just a report
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());
    final double[] x = new double[1000000];
    for (int i = 0; i < x.length; i++) {
      x[i] = nextUniformDouble(rng);
    }

    final MathLog fl = new MathLog();

    final TimingService ts = new TimingService(5);
    // ts.execute(new DoubleTimingTask(new TestLog(fl), 0, x));
    ts.execute(new DoubleTimingTask(new Test1PLog(fl), 0, x));
    ts.execute(new DoubleTimingTask(new TestLog1P(fl), 0, x));
    ts.execute(new DoubleTimingTask(new TestLog1PApache(fl), 0, x));
    // ts.execute(new DoubleTimingTask(new TestLog(fl), 0, x));
    ts.execute(new DoubleTimingTask(new Test1PLog(fl), 0, x));
    ts.execute(new DoubleTimingTask(new TestLog1P(fl), 0, x));
    ts.execute(new DoubleTimingTask(new TestLog1PApache(fl), 0, x));

    final int size = ts.getSize();
    ts.repeat(size);
    logger.info(ts.getReport(size));
  }

  @SpeedTag
  @SeededTest
  public void canTestFloatVsDoubleSpeed(RandomSeed seed) {
    // No assertions, this is just a report
    Assumptions.assumeTrue(logger.isLoggable(Level.INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final UniformRandomProvider rng = RngUtils.create(seed.getSeedAsLong());
    final double[] x = new double[1000000];
    final float[] xf = new float[x.length];
    for (int i = 0; i < x.length; i++) {
      x[i] = nextUniformFloat(rng);
      xf[i] = (float) x[i];
    }

    final TimingService ts = new TimingService(5);
    ts.execute(new DoubleTimingTask(new TestLog(new MathLog()), 0, x));
    ts.execute(new DoubleTimingTask(new TestLog(new FastMathLog()), 0, x));
    for (final int q : new int[] {11}) {
      final int n = 23 - q;
      final IcsiFastLog ff = IcsiFastLog.create(n, DataType.FLOAT);
      final IcsiFastLog fd = IcsiFastLog.create(n, DataType.DOUBLE);
      final IcsiFastLog ff2 = IcsiFastLog.create(n, DataType.FLOAT);
      ts.execute(new DoubleToFloatTimingTask(new TestLog(ff), q, x, xf));
      ts.execute(new DoubleToFloatTimingTask(new TestFastLog(ff), q, x, xf));
      ts.execute(new FloatTimingTask(new TestLog(ff2), q, xf));
      ts.execute(new FloatTimingTask(new TestFastLog(ff2), q, xf));
      ts.execute(new DoubleTimingTask(new TestLog(fd), q, x));
      ts.execute(new DoubleTimingTask(new TestFastLog(fd), q, x));

      // ts.execute(new DoubleToFloatTimingTask(new TestLog(ff), q, x, xf));
      // ts.execute(new DoubleToFloatTimingTask(new TestFastLog(ff), q, x, xf));
      // ts.execute(new FloatTimingTask(new TestLog(ff2), q, xf));
      // ts.execute(new FloatTimingTask(new TestFastLog(ff2), q, xf));
      // ts.execute(new DoubleTimingTask(new TestLog(fd), q, x));
      // ts.execute(new DoubleTimingTask(new TestFastLog(fd), q, x));

      final TurboLog tf = new TurboLog(n);
      ts.execute(new DoubleToFloatTimingTask(new TestLog(tf), q, x, xf));
      ts.execute(new DoubleToFloatTimingTask(new TestFastLog(tf), q, x, xf));
      ts.execute(new FloatTimingTask(new TestLog(tf), q, xf));
      ts.execute(new FloatTimingTask(new TestFastLog(tf), q, xf));
      ts.execute(new DoubleTimingTask(new TestLog(tf), q, x));
      ts.execute(new DoubleTimingTask(new TestFastLog(tf), q, x));

      final TurboLog2 tf2 = new TurboLog2(n);
      ts.execute(new DoubleToFloatTimingTask(new TestLog(tf2), q, x, xf));
      ts.execute(new DoubleToFloatTimingTask(new TestFastLog(tf2), q, x, xf));
      ts.execute(new FloatTimingTask(new TestLog(tf2), q, xf));
      ts.execute(new FloatTimingTask(new TestFastLog(tf2), q, xf));
      ts.execute(new DoubleTimingTask(new TestLog(tf2), q, x));
      ts.execute(new DoubleTimingTask(new TestFastLog(tf2), q, x));

      // Slower as the look-up table is bigger
      final FFastLog f1 = new FFastLog(n);
      final DFastLog f2 = new DFastLog(n);
      ts.execute(new FloatTimingTask(new TestLog(f1), q, xf));
      ts.execute(new FloatTimingTask(new TestFastLog(f1), q, xf));
      ts.execute(new DoubleTimingTask(new TestLog(f2), q, x));
      ts.execute(new DoubleTimingTask(new TestFastLog(f2), q, x));
    }

    final int size = ts.getSize();
    ts.repeat(size);
    logger.info(ts.getReport(size));
  }

  /**
   * Get a supplier for the string using the format and arguments.
   *
   * <p>This can be used where it is not convenient to create a lambda function directly because the
   * arguments are not effectively final.
   *
   * <p>Returns a function of:
   *
   * <pre>
   * {@code
   * () -> String.format(format, args);
   * }
   * </pre>
   *
   * @param format the format
   * @param args the arguments
   * @return the supplier
   */
  static final Supplier<String> getSupplier(final String format, final Object... args) {
    return () -> String.format(format, args);
  }
}
