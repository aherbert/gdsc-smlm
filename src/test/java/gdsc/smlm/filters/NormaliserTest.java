package gdsc.smlm.filters;

import org.junit.Assert;
import org.junit.Test;

public class NormaliserTest
{
	private gdsc.core.utils.Random rand;

	int[] primes = new int[] { 113, 97, 53 };
	int[] boxSizes = new int[] { 15, 9, 5, 3, 2, 1 };

	private float[] floatCreateData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			//data[i] = i;
			data[i] = rand.next();
		//rand.shuffle(data);

		return data;
	}

	@Test
	public void nonNormaliserCanCopyToOutDataWithBorder()
	{
		rand = new gdsc.core.utils.Random(-300519);

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = floatCreateData(width, height);
				for (int boxSize : boxSizes)
				{
					// Assume fixed normaliser works
					FixedNormaliser n = new FixedNormaliser(1);
					NonNormaliser nn = new NonNormaliser();
					float[] e = new float[data.length];
					float[] o = new float[data.length];
					n.normalise(data, e, width, height, boxSize);
					nn.normalise(data, o, width, height, boxSize);
					Assert.assertArrayEquals(o, e, 0);
					return;
				}
			}
	}
}