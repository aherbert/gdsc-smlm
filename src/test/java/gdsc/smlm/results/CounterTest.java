package gdsc.smlm.results;

import org.junit.Assert;
import org.junit.Test;

public class CounterTest
{
	@Test
	public void canIncrementAndGet()
	{
		Counter c = new Counter();
		Assert.assertEquals(1, c.incrementAndGet());
		Assert.assertEquals(2, c.incrementAndGet());
		Assert.assertEquals(3, c.incrementAndGet());
		
		c = new Counter();
		Assert.assertEquals(10, c.incrementAndGet(10));
		Assert.assertEquals(20, c.incrementAndGet(10));
		Assert.assertEquals(30, c.incrementAndGet(10));
	}

	@Test
	public void canGetAndIncrement()
	{
		Counter c = new Counter();
		Assert.assertEquals(0, c.getAndIncrement());
		Assert.assertEquals(1, c.getAndIncrement());
		Assert.assertEquals(2, c.getAndIncrement());
		
		c = new Counter();
		Assert.assertEquals(0, c.getAndIncrement(10));
		Assert.assertEquals(10, c.getAndIncrement(10));
		Assert.assertEquals(20, c.getAndIncrement(10));
	}
}
