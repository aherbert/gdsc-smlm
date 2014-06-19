package gdsc.smlm.fitting.linear;

import java.util.Arrays;

import gdsc.smlm.fitting.linear.GaussJordan;

import org.junit.Assert;
import org.junit.Test;

public class GaussJordanTest
{
	@Test
	public void canSolveLinearEquation()
	{
		GaussJordan solver = new GaussJordan();
		
		// Solves (one) linear equation, a x = b, for x[n]
		
		// Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
		float[][] a = new float[][] { 
				new float[] { 1, 2, 3 }, 
				new float[] { 4, 5, 6 }, 
				new float[] { 7, 8, 10 } 
				};
		float[] b = new float[] { 3, 3, 4 };
		float[] expecteds = new float[] { -2, 1, 1 };
		
		boolean result = solver.solve(a, b);
		
		Assert.assertTrue(result);
		Assert.assertArrayEquals(expecteds, b, 1e-4f);
		
		log("x = %s\n", Arrays.toString(b));
		for (int i=0; i<b.length; i++)
			log("a[%d] = %s\n", i, Arrays.toString(a[i]));
	}
	
	void log(String format, Object... args)
	{
		System.out.printf(format, args);
	}
}
