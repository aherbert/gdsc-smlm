package gdsc.smlm.data.units;

import gdsc.smlm.data.units.Unit;

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
