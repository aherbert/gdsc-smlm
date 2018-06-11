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
package gdsc.smlm.fitting.nonlinear.gradient;


/**
 * Helper functions for the gradient procedures
 */
class GradientProcedureHelper
{
	static void initialiseWorkingMatrix4(double[] data)
	{
		data[0] = 0;
		data[1] = 0;
		data[2] = 0;
		data[3] = 0;
		data[4] = 0;
		data[5] = 0;
		data[6] = 0;
		data[7] = 0;
		data[8] = 0;
		data[9] = 0;
	}

	static void initialiseWorkingMatrix5(double[] data)
	{
		data[0] = 0;
		data[1] = 0;
		data[2] = 0;
		data[3] = 0;
		data[4] = 0;
		data[5] = 0;
		data[6] = 0;
		data[7] = 0;
		data[8] = 0;
		data[9] = 0;
		data[10] = 0;
		data[11] = 0;
		data[12] = 0;
		data[13] = 0;
		data[14] = 0;
	}

	static void initialiseWorkingMatrix6(double[] data)
	{
		data[0] = 0;
		data[1] = 0;
		data[2] = 0;
		data[3] = 0;
		data[4] = 0;
		data[5] = 0;
		data[6] = 0;
		data[7] = 0;
		data[8] = 0;
		data[9] = 0;
		data[10] = 0;
		data[11] = 0;
		data[12] = 0;
		data[13] = 0;
		data[14] = 0;
		data[15] = 0;
		data[16] = 0;
		data[17] = 0;
		data[18] = 0;
		data[19] = 0;
		data[20] = 0;
	}

	static void getMatrix(double[] data, double[][] matrix, int n)
	{
		// Generate symmetric matrix
		for (int j = 0, i = 0; j < n; j++)
			for (int k = 0; k <= j; k++)
			{
				//System.out.printf("matrix[%d][%d] = data[%d];\n", j, k, i);
				//if (j != k)
				//	System.out.printf("matrix[%d][%d] = data[%d];\n", k, j, i);
				matrix[j][k] = matrix[k][j] = data[i++];
			}
		//throw new RuntimeException();
	}

	static void getMatrix4(double[] data, double[][] matrix)
	{
		matrix[0][0] = data[0];
		matrix[1][0] = data[1];
		matrix[0][1] = data[1];
		matrix[1][1] = data[2];
		matrix[2][0] = data[3];
		matrix[0][2] = data[3];
		matrix[2][1] = data[4];
		matrix[1][2] = data[4];
		matrix[2][2] = data[5];
		matrix[3][0] = data[6];
		matrix[0][3] = data[6];
		matrix[3][1] = data[7];
		matrix[1][3] = data[7];
		matrix[3][2] = data[8];
		matrix[2][3] = data[8];
		matrix[3][3] = data[9];
	}

	static void getMatrix5(double[] data, double[][] matrix)
	{
		matrix[0][0] = data[0];
		matrix[1][0] = data[1];
		matrix[0][1] = data[1];
		matrix[1][1] = data[2];
		matrix[2][0] = data[3];
		matrix[0][2] = data[3];
		matrix[2][1] = data[4];
		matrix[1][2] = data[4];
		matrix[2][2] = data[5];
		matrix[3][0] = data[6];
		matrix[0][3] = data[6];
		matrix[3][1] = data[7];
		matrix[1][3] = data[7];
		matrix[3][2] = data[8];
		matrix[2][3] = data[8];
		matrix[3][3] = data[9];
		matrix[4][0] = data[10];
		matrix[0][4] = data[10];
		matrix[4][1] = data[11];
		matrix[1][4] = data[11];
		matrix[4][2] = data[12];
		matrix[2][4] = data[12];
		matrix[4][3] = data[13];
		matrix[3][4] = data[13];
		matrix[4][4] = data[14];
	}

	static void getMatrix6(double[] data, double[][] matrix)
	{
		// Generate symmetric matrix
		matrix[0][0] = data[0];
		matrix[1][0] = data[1];
		matrix[0][1] = data[1];
		matrix[1][1] = data[2];
		matrix[2][0] = data[3];
		matrix[0][2] = data[3];
		matrix[2][1] = data[4];
		matrix[1][2] = data[4];
		matrix[2][2] = data[5];
		matrix[3][0] = data[6];
		matrix[0][3] = data[6];
		matrix[3][1] = data[7];
		matrix[1][3] = data[7];
		matrix[3][2] = data[8];
		matrix[2][3] = data[8];
		matrix[3][3] = data[9];
		matrix[4][0] = data[10];
		matrix[0][4] = data[10];
		matrix[4][1] = data[11];
		matrix[1][4] = data[11];
		matrix[4][2] = data[12];
		matrix[2][4] = data[12];
		matrix[4][3] = data[13];
		matrix[3][4] = data[13];
		matrix[4][4] = data[14];
		matrix[5][0] = data[15];
		matrix[0][5] = data[15];
		matrix[5][1] = data[16];
		matrix[1][5] = data[16];
		matrix[5][2] = data[17];
		matrix[2][5] = data[17];
		matrix[5][3] = data[18];
		matrix[3][5] = data[18];
		matrix[5][4] = data[19];
		matrix[4][5] = data[19];
		matrix[5][5] = data[20];
	}

	static void getMatrix(double[] data, double[] matrix, int n)
	{
		// Generate symmetric matrix
		for (int j = 0, i = 0; j < n; j++)
			for (int k = 0; k <= j; k++)
			{
				//System.out.printf("matrix[%d] = data[%d];\n", j * n + k, i);
				//if (j != k)
				//	System.out.printf("matrix[%d] = data[%d];\n", k * n + j, i);
				matrix[j * n + k] = matrix[k * n + j] = data[i++];
			}
		//throw new RuntimeException();
	}

	static void getMatrix4(double[] data, double[] matrix)
	{
		matrix[0] = data[0];
		matrix[4] = data[1];
		matrix[1] = data[1];
		matrix[5] = data[2];
		matrix[8] = data[3];
		matrix[2] = data[3];
		matrix[9] = data[4];
		matrix[6] = data[4];
		matrix[10] = data[5];
		matrix[12] = data[6];
		matrix[3] = data[6];
		matrix[13] = data[7];
		matrix[7] = data[7];
		matrix[14] = data[8];
		matrix[11] = data[8];
		matrix[15] = data[9];
	}

	static void getMatrix5(double[] data, double[] matrix)
	{
		matrix[0] = data[0];
		matrix[5] = data[1];
		matrix[1] = data[1];
		matrix[6] = data[2];
		matrix[10] = data[3];
		matrix[2] = data[3];
		matrix[11] = data[4];
		matrix[7] = data[4];
		matrix[12] = data[5];
		matrix[15] = data[6];
		matrix[3] = data[6];
		matrix[16] = data[7];
		matrix[8] = data[7];
		matrix[17] = data[8];
		matrix[13] = data[8];
		matrix[18] = data[9];
		matrix[20] = data[10];
		matrix[4] = data[10];
		matrix[21] = data[11];
		matrix[9] = data[11];
		matrix[22] = data[12];
		matrix[14] = data[12];
		matrix[23] = data[13];
		matrix[19] = data[13];
		matrix[24] = data[14];
	}

	static void getMatrix6(double[] data, double[] matrix)
	{
		// Generate symmetric matrix
		matrix[0] = data[0];
		matrix[6] = data[1];
		matrix[1] = data[1];
		matrix[7] = data[2];
		matrix[12] = data[3];
		matrix[2] = data[3];
		matrix[13] = data[4];
		matrix[8] = data[4];
		matrix[14] = data[5];
		matrix[18] = data[6];
		matrix[3] = data[6];
		matrix[19] = data[7];
		matrix[9] = data[7];
		matrix[20] = data[8];
		matrix[15] = data[8];
		matrix[21] = data[9];
		matrix[24] = data[10];
		matrix[4] = data[10];
		matrix[25] = data[11];
		matrix[10] = data[11];
		matrix[26] = data[12];
		matrix[16] = data[12];
		matrix[27] = data[13];
		matrix[22] = data[13];
		matrix[28] = data[14];
		matrix[30] = data[15];
		matrix[5] = data[15];
		matrix[31] = data[16];
		matrix[11] = data[16];
		matrix[32] = data[17];
		matrix[17] = data[17];
		matrix[33] = data[18];
		matrix[23] = data[18];
		matrix[34] = data[19];
		matrix[29] = data[19];
		matrix[35] = data[20];
	}
}
