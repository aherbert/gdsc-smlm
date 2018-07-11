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
package gdsc.smlm.fitting.linear;

/**
 * Solves (one) linear equation, a x = b
 */
public class GaussJordan
{
	private int max_row, max_col;

	/**
	 * @return True if OK
	 */
	private boolean find_pivot(float[][] a, int[] piv)
	{
		float max = 0;

		for (int i = 0; i < piv.length; i++)
		{
			if (piv[i] != 1)
			{
				for (int j = 0; j < piv.length; j++)
				{
					if (piv[j] == 0)
					{
						if (Math.abs(a[i][j]) >= max)
						{
							max = Math.abs(a[i][j]);
							max_row = i;
							max_col = j;
						}
					}
					else if (piv[j] > 1)
					{
						// This should not happen, i.e. a second pivot around a column
						return false;
					}
				}
			}
		}

		piv[max_col]++;

		return true;
	}

	private void interchange_rows_vector(float[][] a, float[] b)
	{
		for (int j = a[max_row].length; j-- > 0;)
		{
			float tmp = a[max_row][j];
			a[max_row][j] = a[max_col][j];
			a[max_col][j] = tmp;
		}

		float tmp = b[max_row];
		b[max_row] = b[max_col];
		b[max_col] = tmp;
	}

	private boolean pivot_vector(float[][] a, float[] b)
	{
		if (a[max_col][max_col] == 0)
			return false;

		float piv_inv = 1 / a[max_col][max_col];

		a[max_col][max_col] = 1;

		for (int i = 0; i < a[max_col].length; i++)
			a[max_col][i] *= piv_inv;
		b[max_col] *= piv_inv;

		for (int i = 0; i < a[max_col].length; i++)
		{
			if (i != max_col)
			{
				float x = a[i][max_col];
				a[i][max_col] = 0;

				for (int j = 0; j < a[max_col].length; j++)
					a[i][j] -= x * a[max_col][j];

				b[i] -= x * b[max_col];
			}
		}

		return true;
	}

	private void unscramble_vector(float[][] a, int[] row, int[] col)
	{
		for (int j = row.length; j-- > 0;)
		{
			if (row[j] != col[j])
			{
				for (int i = row.length; i-- > 0;)
				{
					float tmp = a[i][row[j]];
					a[i][row[j]] = a[i][col[j]];
					a[i][col[j]] = tmp;
				}
			}
		}
	}

	/**
	 * Solves (one) linear equation, a x = b, for x[n].
	 * <p>
	 * On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(float[][] a, float[] b)
	{
		int[] piv = new int[b.length];
		int[] row = new int[b.length];
		int[] col = new int[b.length];
		return solve(a, b, piv, row, col);
	}

	/**
	 * Solves (one) linear equation, a x = b, for x[n].
	 * <p>
	 * On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
	 * <p>
	 * piv[n], row[n], col[n] (all ints) are used for storage
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(float[][] a, float[] b, int[] piv, int[] row, int[] col)
	{
		max_row = 0;
		max_col = 0;

		for (int i = 0; i < piv.length; i++)
			piv[i] = 0;

		for (int i = 0; i < piv.length; i++)
		{
			if (!find_pivot(a, piv))
			{
				return false;
			}

			if (max_row != max_col)
				interchange_rows_vector(a, b);

			row[i] = max_row;
			col[i] = max_col;

			if (!pivot_vector(a, b))
			{
				return false;
			}
		}

		unscramble_vector(a, row, col);
		return true;
	}

	/**
	 * @return True if OK
	 */
	private boolean find_pivot(double[][] a, int[] piv)
	{
		double max = 0;

		for (int i = 0; i < piv.length; i++)
		{
			if (piv[i] != 1)
			{
				for (int j = 0; j < piv.length; j++)
				{
					if (piv[j] == 0)
					{
						if (Math.abs(a[i][j]) >= max)
						{
							max = Math.abs(a[i][j]);
							max_row = i;
							max_col = j;
						}
					}
					else if (piv[j] > 1)
					{
						// This should not happen, i.e. a second pivot around a column
						return false;
					}
				}
			}
		}

		piv[max_col]++;

		return true;
	}

	private void interchange_rows_vector(double[][] a, double[] b)
	{
		for (int j = a[max_row].length; j-- > 0;)
		{
			double tmp = a[max_row][j];
			a[max_row][j] = a[max_col][j];
			a[max_col][j] = tmp;
		}

		double tmp = b[max_row];
		b[max_row] = b[max_col];
		b[max_col] = tmp;
	}

	private boolean pivot_vector(double[][] a, double[] b)
	{
		if (a[max_col][max_col] == 0)
			return false;

		double piv_inv = 1 / a[max_col][max_col];

		a[max_col][max_col] = 1;

		for (int i = 0; i < a[max_col].length; i++)
			a[max_col][i] *= piv_inv;
		b[max_col] *= piv_inv;

		for (int i = 0; i < a[max_col].length; i++)
		{
			if (i != max_col)
			{
				double x = a[i][max_col];
				a[i][max_col] = 0;

				for (int j = 0; j < a[max_col].length; j++)
					a[i][j] -= x * a[max_col][j];

				b[i] -= x * b[max_col];
			}
		}

		return true;
	}

	private void unscramble_vector(double[][] a, int[] row, int[] col)
	{
		for (int j = row.length; j-- > 0;)
		{
			if (row[j] != col[j])
			{
				for (int i = row.length; i-- > 0;)
				{
					double tmp = a[i][row[j]];
					a[i][row[j]] = a[i][col[j]];
					a[i][col[j]] = tmp;
				}
			}
		}
	}

	/**
	 * Solves (one) linear equation, a x = b, for x[n].
	 * <p>
	 * On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(double[][] a, double[] b)
	{
		int[] piv = new int[b.length];
		int[] row = new int[b.length];
		int[] col = new int[b.length];
		return solve(a, b, piv, row, col);
	}

	/**
	 * Solves (one) linear equation, a x = b, for x[n].
	 * <p>
	 * On input have a[n][n], b[n]. On output these replaced by a_inverse[n][n], x[n].
	 * <p>
	 * piv[n], row[n], col[n] (all ints) are used for storage
	 *
	 * @return False if the equation is singular (no solution)
	 */
	public boolean solve(double[][] a, double[] b, int[] piv, int[] row, int[] col)
	{
		max_row = 0;
		max_col = 0;

		for (int i = 0; i < piv.length; i++)
			piv[i] = 0;

		for (int i = 0; i < piv.length; i++)
		{
			if (!find_pivot(a, piv))
			{
				return false;
			}

			if (max_row != max_col)
				interchange_rows_vector(a, b);

			row[i] = max_row;
			col[i] = max_col;

			if (!pivot_vector(a, b))
			{
				return false;
			}
		}

		unscramble_vector(a, row, col);
		return true;
	}
}
