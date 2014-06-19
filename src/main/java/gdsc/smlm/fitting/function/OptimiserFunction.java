package gdsc.smlm.fitting.function;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Allow optimisation using Apache Commons Math 3 Optimiser
 */
public abstract class OptimiserFunction
{
	protected List<Double> x = null;
	protected List<Double> y = null;

	public void addPoint(double x, double y)
	{
		if (this.x == null)
		{
			this.x = new ArrayList<Double>();
			this.y = new ArrayList<Double>();
		}
		this.x.add(x);
		this.y.add(y);
	}

	public void addData(float[] x, float[] y)
	{
		this.x = new ArrayList<Double>();
		this.y = new ArrayList<Double>();
		for (int i = 0; i < x.length; i++)
		{
			this.x.add((double) x[i]);
			this.y.add((double) y[i]);
		}
	}

	public void addData(double[] x, double[] y)
	{
		this.x = new ArrayList<Double>();
		this.y = new ArrayList<Double>();
		for (int i = 0; i < x.length; i++)
		{
			this.x.add(x[i]);
			this.y.add(y[i]);
		}
	}

	private double[] getValues(List<Double> list)
	{
		double[] values = new double[list.size()];
		for (int i = 0; i < list.size(); i++)
		{
			values[i] = list.get(i).doubleValue();
		}
		return values;
	}

	public double[] getX()
	{
		return getValues(x);
	}

	public double[] getY()
	{
		return getValues(y);
	}

	public double[] getWeights()
	{
		double[] w = new double[y.size()];
		Arrays.fill(w, 1);
		return w;
	}

	public int size()
	{
		return x.size();
	}
}