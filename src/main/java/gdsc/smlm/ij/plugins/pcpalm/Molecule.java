package gdsc.smlm.ij.plugins.pcpalm;

/**
 * Used to store all the information required for the PC-PALM analysis
 */
public class Molecule
{
	public double x, y, precision, photons;
	
	// Used to construct a single linked list of molecules
	public Molecule next = null;

	public Molecule(double x, double y, double precision, double photons)
	{
		this.x = x;
		this.y = y;
		this.precision = precision;
		this.photons = photons;
	}

	public double distance(Molecule other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return Math.sqrt(dx * dx + dy * dy);
	}

	public double distance2(Molecule other)
	{
		final double dx = x - other.x;
		final double dy = y - other.y;
		return dx * dx + dy * dy;
	}
}