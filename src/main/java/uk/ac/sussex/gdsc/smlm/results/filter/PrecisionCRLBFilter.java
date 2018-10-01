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
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a precision threshold. Calculates the precision using the Cramér-Rao lower bound (CRLB) of the
 * variance of estimators of the fit parameter. The variance for the fitted X and Y position is averaged to produce a
 * localisation precision.
 */
public class PrecisionCRLBFilter extends DirectFilter implements IMultiFilter
{
    @XStreamAsAttribute
    private final double precision;
    @XStreamOmitField
    private float variance;

    /**
     * Instantiates a new precision CRLB filter.
     *
     * @param precision
     *            the precision
     */
    public PrecisionCRLBFilter(double precision)
    {
        this.precision = Math.max(0, precision);
    }

    @Override
    public void setup(MemoryPeakResults peakResults)
    {
        // Add the 2-fold scale factor here:
        // (varX + varY)/2 < precision^2
        // (varX + varY) < precision^2 * 2
        variance = Filter.getUpperSquaredLimit(precision) * 2f;
    }

    @Override
    public boolean accept(PeakResult peak)
    {
        // Use the estimated parameter deviations for the peak
        if (peak.hasParameterDeviations())
        {
            final float vx = peak.getParameterDeviation(PeakResult.X);
            final float vy = peak.getParameterDeviation(PeakResult.Y);
            return (vx * vx + vy * vy) <= variance;
        }
        return true;
    }

    @Override
    public int getValidationFlags()
    {
        return V_LOCATION_VARIANCE_CRLB;
    }

    @Override
    public int validate(final PreprocessedPeakResult peak)
    {
        if (peak.getLocationVarianceCRLB() > variance)
            return V_LOCATION_VARIANCE_CRLB;
        return 0;
    }

    /** {@inheritDoc} */
    @Override
    public String getDescription()
    {
        return "Filter results using an upper precision threshold (uses fitted parameter variance).";
    }

    /** {@inheritDoc} */
    @Override
    public boolean requiresParameterDeviations()
    {
        return true;
    }

    /** {@inheritDoc} */
    @Override
    public int getNumberOfParameters()
    {
        return 1;
    }

    /** {@inheritDoc} */
    @Override
    protected double getParameterValueInternal(int index)
    {
        return precision;
    }

    /** {@inheritDoc} */
    @Override
    public double getParameterIncrement(int index)
    {
        checkIndex(index);
        return PrecisionFilter.DEFAULT_INCREMENT;
    }

    /** {@inheritDoc} */
    @Override
    public ParameterType getParameterType(int index)
    {
        checkIndex(index);
        return ParameterType.PRECISION_CRLB;
    }

    /** {@inheritDoc} */
    @Override
    public Filter adjustParameter(int index, double delta)
    {
        checkIndex(index);
        return new PrecisionCRLBFilter(updateParameter(precision, delta, PrecisionFilter.DEFAULT_RANGE));
    }

    /** {@inheritDoc} */
    @Override
    public Filter create(double... parameters)
    {
        return new PrecisionCRLBFilter(parameters[0]);
    }

    /** {@inheritDoc} */
    @Override
    public void weakestParameters(double[] parameters)
    {
        setMax(parameters, 0, precision);
    }

    /** {@inheritDoc} */
    @Override
    public int lowerBoundOrientation(int index)
    {
        return 1;
    }

    /** {@inheritDoc} */
    @Override
    public double[] upperLimit()
    {
        return new double[] { PrecisionFilter.UPPER_LIMIT };
    }

    /** {@inheritDoc} */
    @Override
    public double[] mutationStepRange()
    {
        return new double[] { PrecisionFilter.DEFAULT_RANGE };
    }

    @Override
    public double getSignal()
    {
        return 0;
    }

    @Override
    public double getSNR()
    {
        return 0;
    }

    @Override
    public double getMinWidth()
    {
        return 0;
    }

    @Override
    public double getMaxWidth()
    {
        return 0;
    }

    @Override
    public double getShift()
    {
        return 0;
    }

    @Override
    public double getEShift()
    {
        return 0;
    }

    @Override
    public double getPrecision()
    {
        return precision;
    }

    @Override
    public PrecisionType getPrecisionType()
    {
        return PrecisionType.CRLB;
    }

    @Override
    public double getMinZ()
    {
        return 0;
    }

    @Override
    public double getMaxZ()
    {
        return 0;
    }
}
