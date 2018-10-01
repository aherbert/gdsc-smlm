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
package uk.ac.sussex.gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a coordinate range
 */
public class CoordinateFilter extends DirectFilter
{
    /** The default increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface. */
    public static final double DEFAULT_INCREMENT = 0.01;
    /** The default range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface. */
    public static final double DEFAULT_RANGE = 1;

    @XStreamAsAttribute
    private final float minX;
    @XStreamAsAttribute
    private final float maxX;
    @XStreamAsAttribute
    private final float minY;
    @XStreamAsAttribute
    private final float maxY;

    /**
     * Instantiates a new coordinate filter.
     *
     * @param minX
     *            the min X
     * @param maxX
     *            the max X
     * @param minY
     *            the min Y
     * @param maxY
     *            the max Y
     */
    public CoordinateFilter(float minX, float maxX, float minY, float maxY)
    {
        if (maxX < minX)
        {
            final float f = maxX;
            maxX = minX;
            minX = f;
        }
        this.minX = minX;
        this.maxX = maxX;
        if (maxY < minY)
        {
            final float f = maxY;
            maxY = minY;
            minY = f;
        }
        this.minY = minY;
        this.maxY = maxY;
    }

    @Override
    protected String generateName()
    {
        return "X " + minX + "-" + maxX + ", Y " + minY + "-" + maxY;
    }

    @Override
    public void setup(MemoryPeakResults peakResults)
    {
        // Ignore
    }

    @Override
    public boolean accept(PeakResult peak)
    {
        return peak.getXPosition() >= minX && peak.getXPosition() <= maxX && peak.getYPosition() >= minY &&
                peak.getYPosition() <= maxY;
    }

    @Override
    public int getValidationFlags()
    {
        return V_X | V_Y;
    }

    @Override
    public int validate(final PreprocessedPeakResult peak)
    {
        if (peak.getX() < minX || peak.getX() > maxX)
            return V_X;
        if (peak.getY() < minY || peak.getY() > maxY)
            return V_Y;
        return 0;
    }

    /** {@inheritDoc} */
    @Override
    public String getDescription()
    {
        return "Filter results using a coordinate range.";
    }

    /** {@inheritDoc} */
    @Override
    public int getNumberOfParameters()
    {
        return 4;
    }

    /** {@inheritDoc} */
    @Override
    protected double getParameterValueInternal(int index)
    {
        switch (index)
        {
            case 0:
                return minX;
            case 1:
                return maxX;
            case 2:
                return minY;
            default:
                return maxY;
        }
    }

    /** {@inheritDoc} */
    @Override
    public double getParameterIncrement(int index)
    {
        checkIndex(index);
        return DEFAULT_INCREMENT;
    }

    /** {@inheritDoc} */
    @Override
    public double getDisabledParameterValue(int index)
    {
        checkIndex(index);
        switch (index)
        {
            case 0:
                return Double.NEGATIVE_INFINITY;
            case 1:
                return Double.POSITIVE_INFINITY;
            case 2:
                return Double.NEGATIVE_INFINITY;
            default:
                return Double.POSITIVE_INFINITY;
        }
    }

    /** {@inheritDoc} */
    @Override
    public ParameterType getParameterType(int index)
    {
        checkIndex(index);
        switch (index)
        {
            case 0:
                return ParameterType.MIN_X;
            case 1:
                return ParameterType.MAX_X;
            case 2:
                return ParameterType.MIN_Y;
            default:
                return ParameterType.MAX_Y;
        }
    }

    /** {@inheritDoc} */
    @Override
    public Filter adjustParameter(int index, double delta)
    {
        checkIndex(index);
        switch (index)
        {
            case 0:
                return new CoordinateFilter(updateParameter(minX, delta, DEFAULT_RANGE), maxX, minY, maxY);
            case 1:
                return new CoordinateFilter(minX, updateParameter(maxX, delta, DEFAULT_RANGE), minY, maxY);
            case 2:
                return new CoordinateFilter(minX, maxX, updateParameter(minY, delta, DEFAULT_RANGE), maxY);
            default:
                return new CoordinateFilter(minX, maxX, minY, updateParameter(maxY, delta, DEFAULT_RANGE));
        }
    }

    /** {@inheritDoc} */
    @Override
    public Filter create(double... parameters)
    {
        return new CoordinateFilter((float) parameters[0], (float) parameters[1], (float) parameters[2],
                (float) parameters[3]);
    }

    /** {@inheritDoc} */
    @Override
    public void weakestParameters(double[] parameters)
    {
        setMin(parameters, 0, minX);
        setMax(parameters, 1, maxX);
        setMin(parameters, 2, minY);
        setMax(parameters, 3, maxY);
    }

    /** {@inheritDoc} */
    @Override
    public int lowerBoundOrientation(int index)
    {
        return (index == 1 || index == 3) ? 1 : -1;
    }

    /** {@inheritDoc} */
    @Override
    public int length()
    {
        return 4;
    }

    /** {@inheritDoc} */
    @Override
    public double[] sequence()
    {
        // Ignore the mode parameters
        return new double[] { minX, maxX, minY, maxY };
    }

    /** {@inheritDoc} */
    @Override
    public double[] mutationStepRange()
    {
        return new double[] { DEFAULT_RANGE, DEFAULT_RANGE, DEFAULT_RANGE, DEFAULT_RANGE };
    }
}
