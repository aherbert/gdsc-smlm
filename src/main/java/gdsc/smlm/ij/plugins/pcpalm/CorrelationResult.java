package gdsc.smlm.ij.plugins.pcpalm;

import gdsc.smlm.results.ImageSource;

/**
 * Used to store the correlation (g(r)) result for the PC-PALM analysis
 */
public class CorrelationResult
{
	public int id;
	public ImageSource source;
	public double minx, miny, maxx, maxy, nmPerPixel, peakDensity;
	public double n;
	public boolean binaryImage;
	public double[][] gr;
	public boolean spatialDomain;

	public CorrelationResult(int id, ImageSource source, double minx, double miny, double maxx, double maxy, double uniquePoints,
			double nmPerPixel, double peakDensity, boolean binaryImage, double[][] gr, boolean spatialDomain)
	{
		this.id = id;
		this.source = source;
		this.minx = minx;
		this.miny = miny;
		this.maxx = maxx;
		this.maxy = maxy;
		this.n = uniquePoints;
		this.nmPerPixel = nmPerPixel;
		this.peakDensity = peakDensity;
		this.binaryImage = binaryImage;
		this.gr = gr;
		this.spatialDomain = spatialDomain;
	}
}