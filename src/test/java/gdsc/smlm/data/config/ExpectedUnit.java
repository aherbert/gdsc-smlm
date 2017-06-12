package gdsc.smlm.data.config;

public class ExpectedUnit<T>
{
	final T u;
	final double value;

	public ExpectedUnit(T u, double value)
	{
		this.u = u;
		this.value = value;
	}
}
