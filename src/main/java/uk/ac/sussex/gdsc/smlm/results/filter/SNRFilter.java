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
 * Filter results using a signal-to-noise ratio (SNR) threshold
 */
public class SNRFilter extends DirectFilter implements IMultiFilter
{
    /** The default increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface. */
    public static final double DEFAULT_INCREMENT = 1;
    /** The default range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface. */
    public static final double DEFAULT_RANGE = 10;

    @XStreamAsAttribute
    private final float snr;

    /**
     * Instantiates a new signal-to-noise ratio (SNR) filter.
     *
     * @param snr
     *            the signal-to-noise ratio (SNR)
     */
    public SNRFilter(float snr)
    {
        this.snr = Math.max(0, snr);
    }

    @Override
    public void setup(MemoryPeakResults peakResults)
    {
        // Do nothing
    }

    @Override
    public boolean accept(PeakResult peak)
    {
        return peak.getSNR() >= this.snr;
    }

    @Override
    public int getValidationFlags()
    {
        return V_SNR;
    }

    @Override
    public int validate(final PreprocessedPeakResult peak)
    {
        if (peak.getSNR() < this.snr)
            return V_SNR;
        return 0;
    }

    /** {@inheritDoc} */
    @Override
    public String getDescription()
    {
        return "Filter results using a lower SNR threshold.";
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
        return snr;
    }

    /** {@inheritDoc} */
    @Override
    public double getParameterIncrement(int index)
    {
        checkIndex(index);
        return SNRFilter.DEFAULT_INCREMENT;
    }

    /** {@inheritDoc} */
    @Override
    public ParameterType getParameterType(int index)
    {
        checkIndex(index);
        return ParameterType.SNR;
    }

    /** {@inheritDoc} */
    @Override
    public Filter adjustParameter(int index, double delta)
    {
        checkIndex(index);
        return new SNRFilter(updateParameter(snr, delta, DEFAULT_RANGE));
    }

    /** {@inheritDoc} */
    @Override
    public Filter create(double... parameters)
    {
        return new SNRFilter((float) parameters[0]);
    }

    /** {@inheritDoc} */
    @Override
    public void weakestParameters(double[] parameters)
    {
        setMin(parameters, 0, snr);
    }

    /** {@inheritDoc} */
    @Override
    public double[] mutationStepRange()
    {
        return new double[] { DEFAULT_RANGE };
    }

    @Override
    public double getSignal()
    {
        return 0;
    }

    @Override
    public double getSNR()
    {
        return snr;
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
        return 0;
    }

    @Override
    public PrecisionType getPrecisionType()
    {
        return PrecisionType.NONE;
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
