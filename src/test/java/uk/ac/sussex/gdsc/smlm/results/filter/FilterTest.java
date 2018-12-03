package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.core.utils.XmlUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingResult;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class FilterTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(FilterTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @SeededTest
  public void canCompareMultiFilter(RandomSeed seed) {
    final UniformRandomProvider UniformRandomProvider = RngUtils.create(seed.getSeedAsLong());
    final MultiFilter f = new MultiFilter(0, 0, 0, 0, 0, 0, 0, 0, 0);
    for (int i = 1000; i-- > 0;) {
      final MultiFilter f1 =
          (MultiFilter) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final MultiFilter f2 =
          (MultiFilter) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final int e = f1.weakest((Filter) f2);
      final int o = f1.weakest(f2);
      Assertions.assertEquals(e, o);
    }
  }

  @SeededTest
  public void canCompareMultiFilter2(RandomSeed seed) {
    final UniformRandomProvider UniformRandomProvider = RngUtils.create(seed.getSeedAsLong());
    final MultiFilter2 f = new MultiFilter2(0, 0, 0, 0, 0, 0, 0, 0, 0);
    for (int i = 1000; i-- > 0;) {
      final MultiFilter2 f1 =
          (MultiFilter2) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final MultiFilter2 f2 =
          (MultiFilter2) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final int e = f1.weakest((Filter) f2);
      final int o = f1.weakest(f2);
      Assertions.assertEquals(e, o);
    }
  }

  @SeededTest
  public void directCompareMultiFilterIsFaster(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final UniformRandomProvider UniformRandomProvider = RngUtils.create(seed.getSeedAsLong());
    final MultiFilter f1 = new MultiFilter(0, 0, 0, 0, 0, 0, 0, 0, 0);
    final MultiFilter2 f2 = new MultiFilter2(0, 0, 0, 0, 0, 0, 0, 0, 0);

    final double[][][] data = new double[1000][][];
    for (int i = data.length; i-- > 0;) {
      data[i] = new double[][] {random(f1.getNumberOfParameters(), UniformRandomProvider),
          random(f1.getNumberOfParameters(), UniformRandomProvider)};
    }

    final TimingService ts = new TimingService();

    ts.execute(new BaseTimingTask("MultiFilter") {
      @Override
      public Object getData(int i) {
        return new MultiFilter[] {(MultiFilter) f1.create(data[i][0]),
            (MultiFilter) f1.create(data[i][1])};
      }

      @Override
      public Object run(Object data) {
        final MultiFilter f1 = ((MultiFilter[]) data)[0];
        final MultiFilter f2 = ((MultiFilter[]) data)[1];
        f1.weakest((Filter) f2);
        return null;
      }

      @Override
      public int getSize() {
        return data.length;
      }
    });

    ts.execute(new BaseTimingTask("MultiFilter direct") {
      @Override
      public Object getData(int i) {
        return new MultiFilter[] {(MultiFilter) f1.create(data[i][0]),
            (MultiFilter) f1.create(data[i][1])};
      }

      @Override
      public Object run(Object data) {
        final MultiFilter f1 = ((MultiFilter[]) data)[0];
        final MultiFilter f2 = ((MultiFilter[]) data)[1];
        f1.weakest(f2);
        return null;
      }

      @Override
      public int getSize() {
        return data.length;
      }
    });

    ts.execute(new BaseTimingTask("MultiFilter2") {
      @Override
      public Object getData(int i) {
        return new MultiFilter2[] {(MultiFilter2) f2.create(data[i][0]),
            (MultiFilter2) f2.create(data[i][1])};
      }

      @Override
      public Object run(Object data) {
        final MultiFilter2 f1 = ((MultiFilter2[]) data)[0];
        final MultiFilter2 f2 = ((MultiFilter2[]) data)[1];
        f1.weakest((Filter) f2);
        return null;
      }

      @Override
      public int getSize() {
        return data.length;
      }
    });

    ts.execute(new BaseTimingTask("MultiFilter2 direct") {
      @Override
      public Object getData(int i) {
        return new MultiFilter2[] {(MultiFilter2) f2.create(data[i][0]),
            (MultiFilter2) f2.create(data[i][1])};
      }

      @Override
      public Object run(Object data) {
        final MultiFilter2 f1 = ((MultiFilter2[]) data)[0];
        final MultiFilter2 f2 = ((MultiFilter2[]) data)[1];
        f1.weakest(f2);
        return null;
      }

      @Override
      public int getSize() {
        return data.length;
      }
    });

    ts.check();

    final int size = ts.repeat();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    for (int i = 0; i < size; i += 2) {
      final TimingResult slow = ts.get(-(i + 2));
      final TimingResult fast = ts.get(-(i + 1));
      Assertions.assertTrue(slow.getMin() > fast.getMin());
    }
  }

  private static double[] random(int n, UniformRandomProvider r) {
    final double[] p = new double[n];
    while (n-- > 0) {
      p[n] = r.nextInt(3);
    }
    return p;
  }

  @SeededTest
  public void canSerialiseMultiFilter(RandomSeed seed) {
    // Check the XStream serialisation supports inheritance
    final UniformRandomProvider UniformRandomProvider = RngUtils.create(seed.getSeedAsLong());
    testSerialisation(new MultiFilter(0, 0, 0, 0, 0, 0, 0, 0, 0), UniformRandomProvider);
    testSerialisation(new MultiFilter2(0, 0, 0, 0, 0, 0, 0, 0, 0), UniformRandomProvider);
    testSerialisation(new MultiFilterCRLB(0, 0, 0, 0, 0, 0, 0, 0, 0), UniformRandomProvider);
  }

  private static void testSerialisation(MultiFilter f,
      UniformRandomProvider UniformRandomProvider) {
    for (int i = 10; i-- > 0;) {
      final MultiFilter f1 =
          (MultiFilter) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final String xml = f1.toXML();
      logger.log(TestLogUtils.getRecord(Level.FINE, XmlUtils.prettyPrintXml(xml)));
      final MultiFilter f2 = (MultiFilter) Filter.fromXML(xml);
      Assertions.assertTrue(f1.getClass().equals(f2.getClass()));
      Assertions.assertEquals(f1, f2);
    }
  }
}
