package gdsc.smlm.units;

public class ExpectedUnit<T extends Unit>
{
	final T u;
	final double value;

	public ExpectedUnit(T u, double value)
	{
		this.u = u;
		this.value = value;
	}
}
