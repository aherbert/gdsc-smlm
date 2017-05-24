package gdsc.smlm.function;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.ValueProcedure;

public class DummyGradientFunction implements Gradient1Function, Gradient2Function
{
	int n;

	public DummyGradientFunction(int n)
	{
		this.n = n;
	}

	public int size()
	{
		return 0;
	}

	public void initialise0(double[] a)
	{
	}

	public void initialise1(double[] a)
	{
	}

	public void initialise2(double[] a)
	{
	}

	public int[] gradientIndices()
	{
		return null;
	}

	public int getNumberOfGradients()
	{
		return n;
	}

	public void forEach(ValueProcedure procedure)
	{
	}

	public void forEach(Gradient1Procedure procedure)
	{
	}

	public void forEach(Gradient2Procedure procedure)
	{
	}

	public void initialise(double[] a)
	{
	}
}