package uk.ac.sussex.gdsc.smlm.ga;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import uk.ac.sussex.gdsc.test.BaseTimingTask;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingResult;
import uk.ac.sussex.gdsc.test.TimingService;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

@SuppressWarnings({ "javadoc" })
public class RampedSelectionStrategyTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(RampedSelectionStrategyTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

	@Test
	public void canSearchUsingActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length - 1; i++)
		{
			final long key = sum[i];
			final int j = RampedSelectionStrategy.search(sum, key);
			Assertions.assertEquals(i + 1, j);
		}
	}

	@Test
	public void canBinarySearchUsingActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length - 1; i++)
		{
			final long key = sum[i];
			final int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assertions.assertEquals(i + 1, j);
		}
	}

	@Test
	public void canSearchUsingNotActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length; i++)
		{
			final long key = sum[i] - 1;
			final int j = RampedSelectionStrategy.search(sum, key);
			Assertions.assertEquals(i, j);
		}
	}

	@Test
	public void canBinarySearchUsingNotActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length; i++)
		{
			final long key = sum[i] - 1;
			final int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assertions.assertEquals(i, j);
		}
	}

	@Test
	public void binarySearchEqualsSearch()
	{
		final long[] sum = RampedSelectionStrategy.createSum(100);
		for (int key = (int) sum[sum.length - 1]; key-- > 0;)
		{
			final int i = RampedSelectionStrategy.search(sum, key);
			final int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assertions.assertEquals(i, j);
		}
	}

	@SpeedTag
	@Test
	public void speedTest50()
	{
		ExtraAssumptions.assumeLowComplexity();
		speedTest(50, false, 10);
	}

	@SpeedTag
	@Test
	public void speedTest200()
	{
		ExtraAssumptions.assumeMediumComplexity();
		speedTest(200, true, 5);
	}

	@SpeedTag
	@Test
	public void speedTest1000()
	{
		ExtraAssumptions.assumeMediumComplexity();
		speedTest(1000, true, 2);
	}

	// Too slow for common use
	@SpeedTag
	@Test
	public void speedTest5000()
	{
		ExtraAssumptions.assumeHighComplexity();
		speedTest(5000, true, 1);
	}

	private static void speedTest(final int size, boolean faster, int runs)
	{
		final long[] sum = RampedSelectionStrategy.createSum(size);

		final TimingService ts = new TimingService(runs);

		ts.execute(new BaseTimingTask("search" + size)
		{
			@Override
			public Object getData(int i)
			{
				return sum;
			}

			@Override
			public Object run(Object data)
			{
				for (int key = (int) sum[sum.length - 1]; key-- > 0;)
					RampedSelectionStrategy.search(sum, key);
				return null;
			}
			@Override
			public int getSize()
			{
				return 1;
			}
		});

		ts.execute(new BaseTimingTask("binarySearch" + size)
		{
			@Override
			public Object getData(int i)
			{
				return sum[i];
			}

			@Override
			public Object run(Object data)
			{
				for (int key = (int) sum[sum.length - 1]; key-- > 0;)
					RampedSelectionStrategy.binarySearch(sum, key);
				return null;
			}

			@Override
			public int getSize()
			{
				return 1;
			}
		});

		final int n = ts.repeat();
		ts.repeat(n);
		if (logger.isLoggable(Level.INFO))
			ts.report();

		final TimingResult slow = ts.get((faster) ? ts.getSize() - 2 : ts.getSize() - 1);
		final TimingResult fast = ts.get((faster) ? ts.getSize() - 1 : ts.getSize() - 2);
		TestLog.logSpeedTestResult(logger,slow, fast);
	}
}
