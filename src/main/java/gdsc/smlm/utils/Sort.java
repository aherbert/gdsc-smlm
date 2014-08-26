package gdsc.smlm.utils;

import java.util.Arrays;
import java.util.Comparator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Provides sorting functionality
 */
public class Sort
{
	/**
	 * Sorts the indices in descending order of their values
	 * 
	 * @param indices
	 * @param values
	 * @return The indices
	 */
	public static int[] sort(int[] indices, final int[] values)
	{
		// Convert data for sorting
		int[][] data = new int[indices.length][2];
		for (int i = indices.length; i-- > 0;)
		{
			data[i][0] = values[indices[i]];
			data[i][1] = indices[i];
		}

		Arrays.sort(data, new Comparator<int[]>()
		{
			public int compare(int[] o1, int[] o2)
			{
				// Largest first
				return o2[0] - o1[0];
			}
		});

		// Copy back
		for (int i = indices.length; i-- > 0;)
		{
			indices[i] = data[i][1];
		}

		return indices;
	}
	
	/**
	 * Sorts the indices in descending order of their values
	 * 
	 * @param indices
	 * @param values
	 * @return The indices
	 */
	public static int[] sort(int[] indices, final float[] values)
	{
		// Convert data for sorting
		float[][] data = new float[indices.length][2];
		for (int i = indices.length; i-- > 0;)
		{
			data[i][0] = values[indices[i]];
			data[i][1] = indices[i];
		}

		Arrays.sort(data, new Comparator<float[]>()
		{
			public int compare(float[] o1, float[] o2)
			{
				// Largest first
				if (o1[0] > o2[0])
					return -1;
				if (o1[0] < o2[0])
					return 1;
				return 0;
			}
		});

		// Copy back
		for (int i = indices.length; i-- > 0;)
		{
			indices[i] = (int) data[i][1];
		}

		return indices;
	}

	/**
	 * Sorts the indices in descending order of their values
	 * 
	 * @param indices
	 * @param values
	 * @return The indices
	 */
	public static int[] sort(int[] indices, final double[] values)
	{
		// Convert data for sorting
		double[][] data = new double[indices.length][2];
		for (int i = indices.length; i-- > 0;)
		{
			data[i][0] = values[indices[i]];
			data[i][1] = indices[i];
		}

		Arrays.sort(data, new Comparator<double[]>()
		{
			public int compare(double[] o1, double[] o2)
			{
				// Largest first
				if (o1[0] > o2[0])
					return -1;
				if (o1[0] < o2[0])
					return 1;
				return 0;
			}
		});

		// Copy back
		for (int i = indices.length; i-- > 0;)
		{
			indices[i] = (int) data[i][1];
		}

		return indices;
	}

	/**
	 * Sorts the indices in descending order of their values. Does not evaluate
	 * equivalence (returns only 1 or -1 to the sort routine)
	 * <p>
	 * Warning: This method can cause sorting violations within the sort implementation if there are elements with equal
	 * values - use with caution.
	 * 
	 * @param indices
	 * @param values
	 * @return The indices
	 */
	public static int[] sortIgnoreEqual(int[] indices, final float[] values)
	{
		// Convert data for sorting
		float[][] data = new float[indices.length][2];
		for (int i = indices.length; i-- > 0;)
		{
			data[i][0] = values[indices[i]];
			data[i][1] = indices[i];
		}

		Arrays.sort(data, new Comparator<float[]>()
		{
			public int compare(float[] o1, float[] o2)
			{
				if (o1[0] > o2[0])
					return -1;
				return 1;
			}
		});

		// Copy back
		for (int i = indices.length; i-- > 0;)
		{
			indices[i] = (int) data[i][1];
		}

		return indices;
	}
}
