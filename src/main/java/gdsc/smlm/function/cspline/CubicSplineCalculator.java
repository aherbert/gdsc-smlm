package gdsc.smlm.function.cspline;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;

import gdsc.core.data.TrivalueProvider;
import gdsc.core.math.interpolation.CubicSplinePosition;
import gdsc.core.math.interpolation.CustomTricubicFunction;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes the cubic spline coefficients for a 3D cubic spline from interpolated points
 */
public class CubicSplineCalculator
{
	private static final DenseMatrix64F A;
	static
	{
		A = new DenseMatrix64F(64, 64);
		CubicSplinePosition[] s = new CubicSplinePosition[4];
		for (int i = 0; i < 4; i++)
			s[i] = new CubicSplinePosition((double) i / 3);
		int c = 0;
		//int c2 = 0;
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
				{
					double[] t = CustomTricubicFunction.computePowerTable(s[i], s[j], s[k]);
					System.arraycopy(t, 0, A.data, c, 64);
					c += 64;
					//for (int x = 0; x < 64; x++)
					//	A.set(c2, x, t[x]);
					//c2++;
				}
	}

	private final LinearSolver<DenseMatrix64F> solver;

	/**
	 * Instantiates a new cubic spline calculator.
	 */
	public CubicSplineCalculator()
	{
		solver = LinearSolverFactory.linear(64);
		// Note: Linear solver should not modify A or B for this to work!
		if (!solver.setA(A) || solver.modifiesA() || solver.modifiesB())
			throw new IllegalStateException("Unable to create linear solver");
	}

	/**
	 * Compute the coefficients given the spline node value at interpolated points. The value should be interpolated at
	 * [0,1/3,2/3,1] in each dimension.
	 *
	 * @param value
	 *            the value
	 * @return the coefficients (or null if computation failed)
	 */
	public double[] compute(double[][][] value)
	{
		DenseMatrix64F B = new DenseMatrix64F(64, 1);
		int c = 0;
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					B.data[c++] = value[i][j][k];
		solver.solve(B, B);
		return B.data;
	}

	/**
	 * Compute the coefficients given the spline node value at interpolated points. The value should be interpolated at
	 * [0,1/3,2/3,1] in each dimension.
	 *
	 * @param value
	 *            the value
	 * @return the coefficients (or null if computation failed)
	 */
	public double[] compute(TrivalueProvider value)
	{
		DenseMatrix64F B = new DenseMatrix64F(64, 1);
		int c = 0;
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					B.data[c++] = value.get(i, j, k);
		solver.solve(B, B);
		return B.data;
	}
}