/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.ij.plugins.pcpalm;

import gdsc.smlm.results.ImageSource;

/**
 * Used to store the correlation (g(r)) result for the PC-PALM analysis
 */
public class CorrelationResult implements Comparable<CorrelationResult>
{
	public int id;
	public ImageSource source;
	public double minx, miny, maxx, maxy, nmPerPixel, peakDensity;
	public double n;
	public boolean binaryImage;
	public double[][] gr;
	public boolean spatialDomain;

	public CorrelationResult(int id, ImageSource source, double minx, double miny, double maxx, double maxy,
			double uniquePoints, double nmPerPixel, double peakDensity, boolean binaryImage, double[][] gr,
			boolean spatialDomain)
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(CorrelationResult o)
	{
		return id - o.id;
	}
}
