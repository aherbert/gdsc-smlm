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
package gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Filter results using a signal-to-noise ratio (SNR) threshold
 */
public class SNRFilter extends DirectFilter implements IMultiFilter
{
	public static final double DEFAULT_INCREMENT = 1;
	public static final double DEFAULT_RANGE = 10;

	@XStreamAsAttribute
	final float snr;

	public SNRFilter(float snr)
	{
		this.snr = Math.max(0, snr);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a lower SNR threshold.";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
		return snr;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterIncrement(int)
	 */
	@Override
	public double getParameterIncrement(int index)
	{
		checkIndex(index);
		return SNRFilter.DEFAULT_INCREMENT;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterType(int)
	 */
	@Override
	public ParameterType getParameterType(int index)
	{
		checkIndex(index);
		return ParameterType.SNR;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#adjustParameter(int, double)
	 */
	@Override
	public Filter adjustParameter(int index, double delta)
	{
		checkIndex(index);
		return new SNRFilter(updateParameter(snr, delta, DEFAULT_RANGE));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new SNRFilter((float) parameters[0]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, snr);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
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
