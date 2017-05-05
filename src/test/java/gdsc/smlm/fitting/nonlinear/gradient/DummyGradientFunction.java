package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.ValueProcedure;

class DummyGradientFunction implements Gradient1Function
{
	int n;

	DummyGradientFunction(int n)
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
}