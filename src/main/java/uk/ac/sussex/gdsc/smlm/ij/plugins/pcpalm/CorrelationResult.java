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
package uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm;

import uk.ac.sussex.gdsc.smlm.results.ImageSource;

/**
 * Used to store the correlation (g(r)) result for the PC-PALM analysis.
 */
public class CorrelationResult implements Comparable<CorrelationResult>
{
    /** The id. */
    public int id;

    /** The source. */
    public ImageSource source;

    /** The minx. */
    public double minx;

    /** The miny. */
    public double miny;

    /** The maxx. */
    public double maxx;

    /** The maxy. */
    public double maxy;

    /** The nm per pixel. */
    public double nmPerPixel;

    /** The peak density. */
    public double peakDensity;

    /** The number of unique points. */
    public double n;

    /**
     * Set to true if pixels in the image have a 1/0 value. Otherwise it is assumed
     * the image has continuous real values, i.e. pixels can have 1 or more localisations.
     */
    public boolean binaryImage;

    /** The correlation curve. */
    public double[][] gr;

    /** The spatial domain. */
    public boolean spatialDomain;

    /**
     * Instantiates a new correlation result.
     *
     * @param id
     *            the id
     * @param source
     *            the source
     * @param minx
     *            the minx
     * @param miny
     *            the miny
     * @param maxx
     *            the maxx
     * @param maxy
     *            the maxy
     * @param uniquePoints
     *            the unique points
     * @param nmPerPixel
     *            the nm per pixel
     * @param peakDensity
     *            the peak density
     * @param binaryImage
     *            Set to true if pixels in the image have a 1/0 value. Otherwise it is assumed
     *            the image has continuous real values, i.e. pixels can have 1 or more localisations.
     * @param gr
     *            the correlation curve
     * @param spatialDomain
     *            the spatial domain
     */
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
    @Override
    public int compareTo(CorrelationResult o)
    {
        return Integer.compare(id, o.id);
    }
}
