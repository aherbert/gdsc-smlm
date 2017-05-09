package gdsc.smlm.fitting.nonlinear.gradient;

//import java.util.Arrays;

import org.apache.commons.math3.random.Well19937c;

import gdsc.core.ij.Utils;
import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.ValueProcedure;

class FakeGradientFunction implements Gradient2Function, Gradient1Function, NonLinearFunction
{
	private final int maxx, n, nparams;
	private final PseudoRandomGenerator r;
	private final double[] dy_da;

	FakeGradientFunction(int maxx, int nparams)
	{
		this(maxx, nparams, 1000, 30051977);
	}

	FakeGradientFunction(int maxx, int nparams, int randomSize, int randomSeed)
	{
		this.maxx = maxx;
		this.n = maxx * maxx;
		this.nparams = nparams;
		this.r = new PseudoRandomGenerator(randomSize, new Well19937c(randomSeed));
		this.dy_da = new double[nparams];
	}

	public int size()
	{
		return 0;
	}

	public void initialise(double[] a)
	{
		int seed = 0;
		for (int i = a.length; i-- > 0;)
			seed += hashCode(a[i]);
		//System.out.printf("seed = %d\n", seed);
		r.setSeed(seed);
	}

	public void initialise0(double[] a)
	{
		initialise(a);
	}

	public void initialise1(double[] a)
	{
		initialise(a);
	}

	public void initialise2(double[] a)
	{
		initialise(a);
	}
	
	private int hashCode(double d)
	{
		long bits = Double.doubleToLongBits(d);
		return (int) (bits ^ (bits >>> 32));
	}

	public int[] gradientIndices()
	{
		return Utils.newArray(nparams, 0, 1);
	}

	public int getNumberOfGradients()
	{
		return nparams;
	}

	public void forEach(ValueProcedure procedure)
	{
		for (int i = 0; i < n; i++)
			procedure.execute(r.nextDouble());
	}

	public void forEach(Gradient1Procedure procedure)
	{
		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				for (int j = nparams; j-- > 0;)
					dy_da[j] = r.nextDouble() * x + y;
				//System.out.println(Arrays.toString(dy_da));
				procedure.execute(r.nextDouble(), dy_da);
			}
		}
	}
	
	public void forEach(Gradient2Procedure procedure)
	{
		final double[] d2y_da2 = new double[nparams];
		
		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				for (int j = nparams; j-- > 0;)
				{
					dy_da[j] = r.nextDouble() * x + y;
					d2y_da2[j] = r.nextDouble() * x + y;
				}
				//System.out.println(Arrays.toString(dy_da));
				procedure.execute(r.nextDouble(), dy_da, d2y_da2);
			}
		}
	}

	public double eval(int i, double[] dy_da)
	{
		// Unpack the predictor to the 2D coordinates
		int y = i / maxx;
		int x = i % maxx;
		for (int j = nparams; j-- > 0;)
			dy_da[j] = r.nextDouble() * x + y;
		//System.out.println(Arrays.toString(dy_da));
		return r.nextDouble();
	}

	public double eval(int x)
	{
		return r.nextDouble();
	}

	public double eval(int x, double[] dyda, double[] w)
	{
		throw new NotImplementedException();
	}

	public double evalw(int x, double[] w)
	{
		throw new NotImplementedException();
	}

	public boolean canComputeWeights()
	{
		return false;
	}
}