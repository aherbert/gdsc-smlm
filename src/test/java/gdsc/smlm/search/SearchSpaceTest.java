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
package gdsc.smlm.search;

import org.junit.Assert;
import org.junit.Test;

@SuppressWarnings({ "javadoc" })
public class SearchSpaceTest
{
	@Test
	public void canEnumerateSearchSpace()
	{
		SearchDimension d1 = new SearchDimension(0, 10, 1, 1);
		SearchDimension d2 = new SearchDimension(0, 10, 0.5, 2, 2.5, 7.5);
		double[][] ss = SearchSpace.createSearchSpace(createDimensions(d1, d2));
		Assert.assertEquals(d1.getMaxLength() * d2.getMaxLength(), ss.length);
		//for (double[] p : ss)
		//	System.out.println(java.util.Arrays.toString(p));
	}

	@Test
	public void canEnumerateRefineSpace()
	{
		SearchDimension d1 = new SearchDimension(0, 10, 2, 10);
		SearchDimension d2 = new SearchDimension(0, 10, 1, 10);

		double[][] ss = SearchSpace.createRefineSpace(createDimensions(d1, d2), new double[] { 0, 5 });
		//for (double[] p : ss)
		//	System.out.println(java.util.Arrays.toString(p));
		Assert.assertEquals(2 + 3, ss.length);

		ss = SearchSpace.createRefineSpace(createDimensions(d1, d2), new double[] { 4, 5 });
		//for (double[] p : ss)
		//	System.out.println(java.util.Arrays.toString(p));
		Assert.assertEquals(3 + 3, ss.length);

		ss = SearchSpace.createRefineSpace(createDimensions(d1, d2), new double[] { 10, 5 });
		//for (double[] p : ss)
		//	System.out.println(java.util.Arrays.toString(p));
		Assert.assertEquals(2 + 3, ss.length);
	}

	@Test
	public void canMoveCentre()
	{
		SearchDimension d1 = new SearchDimension(0, 10, 0, 1, 2.5, 7.5);
		double[] v1 = d1.values();
		Assert.assertTrue(d1.isAtBounds(0));
		Assert.assertTrue(d1.isAtBounds(v1[0]));
		Assert.assertTrue(d1.isAtBounds(v1[v1.length - 1]));
		Assert.assertFalse(d1.isAtBounds(5));

		d1.setCentre(0);
		Assert.assertTrue(d1.isAtBounds(0));
		Assert.assertTrue(d1.isAtBounds(5));

		double[] v2 = d1.values();
		//System.out.println(java.util.Arrays.toString(v1));
		//System.out.println(java.util.Arrays.toString(v2));
		Assert.assertTrue(v1.length > v2.length);
	}

	@Test
	public void canReduceSearchSpace()
	{
		SearchDimension d1 = new SearchDimension(0, 10, 0, 1);
		d1.setCentre(0);
		double[] v1 = d1.values();
		d1.reduce();
		double[] v2 = d1.values();
		//System.out.println(java.util.Arrays.toString(v1));
		//System.out.println(java.util.Arrays.toString(v2));
		for (int i = 0; i < v1.length; i++)
			Assert.assertEquals(v1[i] * d1.getReduceFactor(), v2[i], 0);
	}

	@Test
	public void canPadRange()
	{
		SearchDimension d1 = new SearchDimension(0, 10, 0, 1, 2.5, 7.5);
		d1.setPad(false);

		double[] v1 = d1.values();
		Assert.assertTrue(d1.isAtBounds(0));

		d1.setCentre(0);
		Assert.assertTrue(d1.isAtBounds(0));

		double[] v2 = d1.values();
		//System.out.println(java.util.Arrays.toString(v1));
		//System.out.println(java.util.Arrays.toString(v2));
		Assert.assertTrue(v1.length > v2.length);

		d1.setCentre(10);
		v2 = d1.values();
		//System.out.println(java.util.Arrays.toString(v2));
		Assert.assertTrue(v1.length > v2.length);

		d1.setPad(true);

		double[] v3 = d1.values();
		//System.out.println(java.util.Arrays.toString(v3));
		Assert.assertTrue(v1.length == v3.length);

		d1.setCentre(0);

		v3 = d1.values();
		//System.out.println(java.util.Arrays.toString(v3));
		Assert.assertTrue(v1.length == v3.length);
	}

	private static SearchDimension[] createDimensions(SearchDimension... d)
	{
		return d;
	}
}
