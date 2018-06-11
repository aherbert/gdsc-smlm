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
package gdsc.smlm.utils;

import org.ejml.data.DenseMatrix64F;

import gdsc.core.utils.Maths;

/**
 * Compute the inertia tensor for a 3D object
 * 
 * @author Alex Herbert
 */
public class Tensor2D
{
	private final double[] com;

	private final DenseMatrix64F tensor;
	private final double[] eigenValues;
	private final double[][] eigenVectors;

	/**
	 * Instantiates a new tensor 3D.
	 *
	 * @param data
	 *            the data (packed in YX order)
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 */
	public Tensor2D(float[] data, int w, int h)
	{
		// Compute centre-of-mass
		double cx = 0;
		double cy = 0;
		double sumXY = 0;
		for (int y = 0, j = 0; y < h; y++)
		{
			double sumX = 0;
			for (int x = 0; x < w; x++)
			{
				float f = data[j++];
				sumX += f;
				cx += f * x;
			}
			sumXY += sumX;
			cy += sumX * y;
		}
		cx = cx / sumXY;
		cy = cy / sumXY;
		com = new double[] { cx, cy };

		// Compute tensor
		// https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
		tensor = new DenseMatrix64F(2, 2);
		for (int y = 0, j = 0; y < h; y++)
		{
			final double dy = y - cy;
			final double dy2 = dy * dy;
			double summ = 0;
			double summdx = 0;
			for (int x = 0; x < w; x++)
			{
				final double m = data[j++];
				final double dx = x - cx;
				final double dx2 = dx * dx;

				//tensor[0][0] += m * dy2;
				//tensor[0][1] -= m * dx * dy;
				//tensor[1][1] += m * dx2;

				summ += m;
				summdx += m * dx;
				//tensor.data[0] += m * dy2;
				//tensor.data[1] -= m * dx * dy;
				tensor.data[3] += m * dx2;
			}
			tensor.data[0] += summ * dy2;
			tensor.data[1] -= summdx * dy;
		}

		// Inertia tensor is symmetric    
		tensor.data[2] = tensor.data[1];

		// Eigen decompose

		// Algebra solution
		//http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html

		final double a = tensor.data[0];
		final double b = tensor.data[1];
		final double c = tensor.data[2];
		final double d = tensor.data[3];
		final double T_2 = (a + d) / 2.0; // Trace / 2 
		final double D = a * d - b * c; // Determinant
		final double root_T2_4_D = Math.sqrt(T_2 * T_2 - D);
		final double L1 = T_2 + root_T2_4_D;
		final double L2 = T_2 - root_T2_4_D;

		eigenValues = new double[] { L1, L2 };
		eigenVectors = new double[2][2];

		if (c != 0)
		{
			eigenVectors[0][0] = L1 - d;
			eigenVectors[0][1] = c;
			eigenVectors[1][0] = L2 - d;
			eigenVectors[1][1] = c;
			normalise(eigenVectors);
		}
		else if (b != 0)
		{
			eigenVectors[0][0] = b;
			eigenVectors[0][1] = L1 - a;
			eigenVectors[1][0] = b;
			eigenVectors[1][1] = L2 - a;
			normalise(eigenVectors);
		}
		else
		{
			// b==0, c==0
			eigenVectors[0][0] = 1;
			//eigenVectors[0][1] = 0;
			//eigenVectors[1][0] = 0;
			eigenVectors[1][1] = 1;
		}

		//// EJML
		//EigenDecomposition<DenseMatrix64F> decomp = DecompositionFactory.eig(2, true, true);
		//
		//if (!decomp.decompose(tensor))
		//{
		//	eigenValues = null;
		//	eigenVectors = null;
		//	return;
		//}
		//
		//eigenValues = new double[2];
		//eigenVectors = new double[2][];
		//for (int i = 0; i < 2; i++)
		//{
		//	eigenValues[i] = decomp.getEigenvalue(i).real;
		//	eigenVectors[i] = decomp.getEigenVector(i).data;
		//}

		// Sort
		if (eigenValues[1] > eigenValues[0])
		{
			double tmp = eigenValues[1];
			eigenValues[1] = eigenValues[0];
			eigenValues[0] = tmp;
			double[] tmp2 = eigenVectors[1];
			eigenVectors[1] = eigenVectors[0];
			eigenVectors[0] = tmp2;
		}
	}

	private static void normalise(double[][] eigenVectors)
	{
		for (int i = 0; i < 2; i++)
		{
			double l = Math.sqrt(Maths.pow2(eigenVectors[i][0]) + Maths.pow2(eigenVectors[i][1]));
			if (l > 0)
			{
				eigenVectors[i][0] /= l;
				eigenVectors[i][1] /= l;
			}
		}
	}

	/**
	 * Gets the centre of mass.
	 *
	 * @return the centre of mass
	 */
	public double[] getCentreOfMass()
	{
		return com;
	}

	/**
	 * Gets the tensor.
	 *
	 * @return the tensor
	 */
	public DenseMatrix64F getTensor()
	{
		return tensor;
	}

	/**
	 * Checks for a succesfull Eigen decomposition.
	 *
	 * @return true, if successful
	 */
	public boolean hasDecomposition()
	{
		return eigenValues != null;
	}

	/**
	 * Gets the eigen values. These are sorted from large to small.
	 * <p>
	 * Note: The minor moment of inertia will be around the longest axis of the object
	 *
	 * @return the eigen values
	 */
	public double[] getEigenValues()
	{
		return eigenValues;
	}

	/**
	 * Gets the eigen vectors.
	 *
	 * @return the eigen vectors
	 */
	public double[][] getEigenVectors()
	{
		return eigenVectors;
	}
}
