package uk.ac.sussex.gdsc.smlm.math3.distribution;

import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({ "javadoc" })
public class CustomGammaDistributionTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(CustomGammaDistributionTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    private abstract class MyTimingTask extends BaseTimingTask
    {
        RandomSeed seed;
        UniformRandomProvider r;
        double shape = 0.5;
        double scale = 300;
        int n = 1000, m = 10;

        public MyTimingTask(String name, RandomSeed seed)
        {
            super(name);
            this.seed = seed;
        }

        @Override
        public int getSize()
        {
            return 1;
        }

        @Override
        public Object getData(int i)
        {
            r = RngUtils.create(seed.getSeedAsLong());
            shape = 0.5;
            return null;
        }
    }

    private class StaticTimingTask extends MyTimingTask
    {
        public StaticTimingTask(RandomSeed seed)
        {
            super("RandomDataGenerator", seed);
        }

        @Override
        public Object run(Object data)
        {
            final RandomDataGenerator rdg = new RandomDataGenerator(new RandomGeneratorAdapter(r));
            final double[] e = new double[n * m];
            for (int i = 0, k = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++, k++)
                    e[k] = rdg.nextGamma(shape, scale);
                shape += 1;
            }
            return e;
        }
    }

    private class InstanceTimingTask extends MyTimingTask
    {
        public InstanceTimingTask(RandomSeed seed)
        {
            super("Instance", seed);
        }

        @Override
        public Object run(Object data)
        {
            final CustomGammaDistribution dist = new CustomGammaDistribution(new RandomGeneratorAdapter(r), 1, scale);
            final double[] e = new double[n * m];
            for (int i = 0, k = 0; i < n; i++)
            {
                dist.setShape(shape);
                for (int j = 0; j < m; j++, k++)
                    e[k] = dist.sample();
                shape += 1;
            }
            return e;
        }
    }

    @SeededTest
    public void canCreateSamples(RandomSeed seed)
    {
        final StaticTimingTask t1 = new StaticTimingTask(seed);
        t1.getData(0);
        final double[] e = (double[]) t1.run(null);

        final InstanceTimingTask t2 = new InstanceTimingTask(seed);
        t2.getData(0);
        final double[] o = (double[]) t2.run(null);

        Assertions.assertArrayEquals(e, o);
    }

    @SpeedTag
    @SeededTest
    public void customDistributionIsFaster(RandomSeed seed)
    {
        Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

        final TimingService ts = new TimingService(5);
        ts.execute(new StaticTimingTask(seed));
        ts.execute(new InstanceTimingTask(seed));

        final int size = ts.getSize();
        ts.repeat(size);
        if (logger.isLoggable(Level.INFO))
            logger.info(ts.getReport(size));

        Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
    }
}
