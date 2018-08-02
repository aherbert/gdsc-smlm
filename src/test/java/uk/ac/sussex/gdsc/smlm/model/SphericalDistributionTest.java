package uk.ac.sussex.gdsc.smlm.model;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.logging.Logger;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.smlm.ij.results.IJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class SphericalDistributionTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(SphericalDistributionTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    @SeededTest
    public void canSampleUsingTransformationMethod(RandomSeed seed)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final double radius = 10 + rg.nextDouble() * 10;
        final SphericalDistribution dist = new SphericalDistribution(radius, new RandomGeneratorAdapter(rg));
        dist.setUseRejectionMethod(false);
        for (int i = 100; i-- > 0;)
            dist.next();
    }

    @SeededTest
    public void canSampleUsingRejectionMethod(RandomSeed seed)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final double radius = 10 + rg.nextDouble() * 10;
        final SphericalDistribution dist = new SphericalDistribution(radius, new RandomGeneratorAdapter(rg));
        dist.setUseRejectionMethod(true);
        for (int i = 100; i-- > 0;)
            dist.next();
    }

    @SeededTest
    public void rejectionMethodIsFasterThanTransformationMethod(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity();

        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final double radius = 10 + rg.nextDouble() * 10;
        final SphericalDistribution dist = new SphericalDistribution(radius, new RandomGeneratorAdapter(rg));
        dist.setUseRejectionMethod(false);
        for (int i = 100; i-- > 0;)
            dist.next();
        dist.setUseRejectionMethod(true);
        for (int i = 100; i-- > 0;)
            dist.next();

        dist.setUseRejectionMethod(false);
        final long time1 = getRunTime(dist);
        dist.setUseRejectionMethod(true);
        final long time2 = getRunTime(dist);
        ExtraAssertions.assertTrue(time1 > time2, "Rejection = %d, Transformation = %d", time2, time1);
        TestLog.info(logger, "Rejection = %d, Transformation = %d", time2, time1);
    }

    private static long getRunTime(SphericalDistribution dist)
    {
        final long start = System.nanoTime();
        for (int i = 1000000; i-- > 0;)
            dist.next();
        return System.nanoTime() - start;
    }

    // These are not tests. They draw an image and use classes outside the package.
    // Comment out for production code.

    //@SeededTest
    public void rejectionMethodSamplesEvenly(RandomSeed seed)
    {
        drawImage(seed, true);
    }

    //@SeededTest
    public void transformationMethodSamplesEvenly(RandomSeed seed)
    {
        drawImage(seed, false);
    }

    private static void drawImage(RandomSeed seed, boolean useRejctionMethod)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults results = new MemoryPeakResults();
        results.setSortAfterEnd(true);
        final int radius = 10;
        final Rectangle bounds = new Rectangle(0, 0, radius * 2, radius * 2);
        final SphericalDistribution dist = new SphericalDistribution(radius, new RandomGeneratorAdapter(rg));
        dist.setUseRejectionMethod(useRejctionMethod);
        final float scale = 10;
        results.begin();
        final float intensity = 1;
        for (int i = 100000; i-- > 0;)
        {
            final double[] xyz = dist.next();
            final int frame = (int) (1 + scale * radius + Math.round(scale * xyz[2]));
            final float x = radius + (float) xyz[0];
            final float y = radius + (float) xyz[1];
            results.add(new PeakResult(frame, x, y, intensity));
        }
        results.end();
        final IJImagePeakResults image = new IJImagePeakResults(
                (useRejctionMethod) ? "Rejection Method" : "Transformation Method", bounds, scale);
        image.setRollingWindowSize(1);
        image.begin();
        image.addAll(Arrays.asList(results.toArray()));
        // Place breakpoint here in debug mode to view the image.
        // It should have an even colour through the stack.
        image.end();
    }
}
