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

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.test.TestSettings;

public class NormaliserTest extends AbstractFilterTest
{
	@Test
	public void nonNormaliserCanCopyToOutDataWithBorder()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(rg,width, height);
				for (int boxSize : boxSizes)
				{
					String msg = null; //String.format("%dx%d : border=%d", width, height, boxSize);
					//System.out.println(msg);
					
					// Assume fixed normaliser works
					FixedNormaliser n = new FixedNormaliser(1);
					NonNormaliser nn = new NonNormaliser();
					float[] e = new float[data.length];
					float[] o = new float[data.length];
					n.normalise(data, e, width, height, boxSize);
					nn.normalise(data, o, width, height, boxSize);
					Assert.assertArrayEquals(msg, o, e, 0);
				}
			}
	}
}
