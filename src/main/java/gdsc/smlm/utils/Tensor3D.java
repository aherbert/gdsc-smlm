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
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.EigenDecomposition;


/**
 * Compute the inertia tensor for a 3D object
 * 
 * @author Alex Herbert
 */
public class Tensor3D
{
	private final double[] com;

	private final DenseMatrix64F tensor;
	private final double[] eigenValues;
	private final double[][] eigenVectors;

	/**
	 * Instantiates a new tensor 3D.
	 *
	 * @param data
	 *            the data (series of z slices packed in YX order)
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 */
	public Tensor3D(float[][] data, int w, int h)
	{
		// Compute centre-of-mass
		double cx = 0;
		double cy = 0;
		double cz = 0;
		double sumXYZ = 0;
		for (int z = 0; z < data.length; z++)
		{
			final float[] d = data[z];
			double sumXY = 0;
			for (int y = 0, j = 0; y < h; y++)
			{
				double sumX = 0;
				for (int x = 0; x < w; x++)
				{
					float f = d[j++];
					sumX += f;
					cx += f * x;
				}
				sumXY += sumX;
				cy += sumX * y;
			}
			cz += sumXY * z;
			sumXYZ += sumXY;
		}
		cx = cx / sumXYZ;
		cy = cy / sumXYZ;
		cz = cz / sumXYZ;
		com = new double[] { cx, cy, cz };

		// Compute tensor
		// https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
		tensor = new DenseMatrix64F(3, 3);
		for (int z = 0; z < data.length; z++)
		{
			final double dz = z - cz;
			final double dz2 = dz * dz;
			final float[] d = data[z];
			for (int y = 0, j = 0; y < h; y++)
			{
				final double dy = y - cy;
				final double dy2 = dy * dy;
				double summ = 0;
				double summdx = 0;
				double summdx2 = 0;
				for (int x = 0; x < w; x++)
				{
					final double m = d[j++];
					final double dx = x - cx;
					final double dx2 = dx * dx;

					//tensor[0][0] += m * (dy2 + dz2);
					//tensor[0][1] -= m * dx * dy;
					//tensor[0][2] -= m * dx * dz;
					//tensor[1][1] += m * (dx2 + dz2);
					//tensor[1][2] -= m * dy * dz;
					//tensor[2][2] += m * (dx2 + dy2);

					summ += m;
					summdx += m * dx;
					summdx2 += m * dx2;
					//tensor.data[0] += m * (dy2 + dz2);
					//tensor.data[1] -= m * dx * dy;
					//tensor.data[2] -= m * dx * dz;
					//tensor.data[4] += m * (dx2 + dz2);
					//tensor.data[5] -= m * dy * dz;
					//tensor.data[8] += m * (dx2 + dy2);
				}
				tensor.data[0] += summ * (dy2 + dz2);
				tensor.data[1] -= summdx * dy;
				tensor.data[2] -= summdx * dz;
				tensor.data[4] += summdx2 + summ * dz2;
				tensor.data[5] -= summ * dy * dz;
				tensor.data[8] += summdx2 + summ * dy2;
			}
		}

		// Inertia tensor is symmetric    
		//tensor[1][0] = tensor[0][1];
		//tensor[2][0] = tensor[0][2];
		//tensor[2][1] = tensor[1][2];
		tensor.data[3] = tensor.data[1];
		tensor.data[6] = tensor.data[2];
		tensor.data[7] = tensor.data[5];

		// Eigen decompose
		EigenDecomposition<DenseMatrix64F> decomp = DecompositionFactory.eig(3, true, true);

		if (!decomp.decompose(tensor))
		{
			eigenValues = null;
			eigenVectors = null;
			return;
		}

		eigenValues = new double[3];
		eigenVectors = new double[3][];
		for (int i = 0; i < 3; i++)
		{
			eigenValues[i] = decomp.getEigenvalue(i).real;
			eigenVectors[i] = decomp.getEigenVector(i).data;
		}

		sort3xN(eigenValues, eigenVectors);
	}

	/**
	 * Vector sorting routine for 3xn set of vectors. On output the weights (and corresponding vectors) will be in
	 * descending order.
	 * 
	 * @param w
	 *            Vector weights
	 * @param v
	 *            Vectors
	 */
	private static void sort3xN(double[] w, double[][] v)
	{
		int k, j, i;
		double p;

		for (i = 3; i-- > 0;)
		{
			p = w[k = i];
			for (j = i; j-- > 0;)
				if (w[j] <= p)
					p = w[k = j];
			if (k != i)
			{
				w[k] = w[i];
				w[i] = p;
				double[] vv = v[k];
				v[k] = v[i];
				v[i] = vv;
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
