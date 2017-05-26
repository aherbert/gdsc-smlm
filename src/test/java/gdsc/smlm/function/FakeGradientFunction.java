package gdsc.smlm.function;

//import java.util.Arrays;

import org.apache.commons.math3.random.Well19937c;

import gdsc.core.ij.Utils;
import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.PseudoRandomSequence;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.ValueProcedure;

public class FakeGradientFunction
		implements ExtendedGradient2Function, Gradient2Function, Gradient1Function, NonLinearFunction
{
	private final int maxx, n, nparams;
	private final PseudoRandomSequence r;
	private final double[] dy_da;

	public FakeGradientFunction(int maxx, int nparams)
	{
		this(maxx, nparams, 1000, 30051977, 10.0);
	}

	public FakeGradientFunction(int maxx, int nparams, double scale)
	{
		this(maxx, nparams, 1000, 30051977, scale);
	}

	public FakeGradientFunction(int maxx, int nparams, int randomSize, int randomSeed, double scale)
	{
		this.maxx = maxx;
		this.n = maxx * maxx;
		this.nparams = nparams;
		this.r = new PseudoRandomSequence(randomSize, new Well19937c(randomSeed), scale);
		this.dy_da = new double[nparams];
	}

	public int size()
	{
		return n;
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

	public void initialiseExtended2(double[] a)
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
		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
			for (int x = 0; x < maxx; x++)
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

	public void forEach(ExtendedGradient2Procedure procedure)
	{
		final double[] d2y_dadb = new double[nparams * nparams];

		// Simulate a 2D forEach
		for (int y = 0; y < maxx; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				for (int j = nparams; j-- > 0;)
					dy_da[j] = r.nextDouble() * x + y;
				for (int j = d2y_dadb.length; j-- > 0;)
					d2y_dadb[j] = r.nextDouble() * x + y;

				//System.out.println(Arrays.toString(dy_da));
				procedure.executeExtended(r.nextDouble(), dy_da, d2y_dadb);
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