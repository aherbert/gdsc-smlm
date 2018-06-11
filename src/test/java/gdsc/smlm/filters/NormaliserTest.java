/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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
				}
			}
	}
}
