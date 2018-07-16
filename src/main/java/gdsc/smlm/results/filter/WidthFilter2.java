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
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Filter results using a width range. Assumes width is identical on the X and Y axis.
 */
public class WidthFilter2 extends DirectFilter implements IMultiFilter
{
	/** The default increment for the min width. Used for {@link gdsc.smlm.ga.Chromosome} interface. */
	public static final double DEFAULT_MIN_INCREMENT = 0.02;
	/** The default range for the min width. Used for {@link gdsc.smlm.ga.Chromosome} interface. */
	public static final double DEFAULT_MIN_RANGE = 1;

	/** The min width. */
	@XStreamAsAttribute
	protected final double minWidth;

	/** The max width. */
	@XStreamAsAttribute
	protected final double maxWidth;

	/** The lower sigma threshold. */
	@XStreamOmitField
	protected float lowerSigmaThreshold;

	/** The upper sigma threshold. */
	@XStreamOmitField
	protected float upperSigmaThreshold;

	/** The width enabled. */
	@XStreamOmitField
	protected boolean widthEnabled;

	/** The calculator. */
	@XStreamOmitField
	protected Gaussian2DPeakResultCalculator calculator;

	/**
	 * Instantiates a new width filter 2.
	 *
	 * @param minWidth
	 *            the min width
	 * @param maxWidth
	 *            the max width
	 */
	public WidthFilter2(double minWidth, double maxWidth)
	{
		// Only swap if max width is enabled
		if (maxWidth != 0 && maxWidth < minWidth)
		{
			final double f = maxWidth;
			maxWidth = minWidth;
			minWidth = f;
		}
		this.minWidth = Math.max(0, minWidth);
		this.maxWidth = Math.max(0, maxWidth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#generateName()
	 */
	@Override
	protected String generateName()
	{
		return "Width " + minWidth + "-" + maxWidth;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#setup(gdsc.smlm.results.MemoryPeakResults)
	 */
	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(), 0);

		// Set the width limit
		lowerSigmaThreshold = 0;
		upperSigmaThreshold = Float.POSITIVE_INFINITY;
		double s = PSFHelper.getGaussian2DWx(peakResults.getPSF());
		lowerSigmaThreshold = (float) (s * minWidth);
		upperSigmaThreshold = Filter.getUpperLimit(s * maxWidth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#setup()
	 */
	@Override
	public void setup()
	{
		setup(minWidth, maxWidth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#setup(int)
	 */
	@Override
	public void setup(int flags)
	{
		if (areSet(flags, IDirectFilter.NO_WIDTH))
			widthEnabled = false;
		else
			setup(minWidth, maxWidth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#setup(int, gdsc.smlm.results.filter.FilterSetupData[])
	 */
	@Override
	public void setup(int flags, FilterSetupData... filterSetupData)
	{
		setup(flags);
	}

	/**
	 * Setup the filter.
	 *
	 * @param minWidth
	 *            the min width
	 * @param maxWidth
	 *            the max width
	 */
	protected void setup(final double minWidth, double maxWidth)
	{
		widthEnabled = false;
		if (maxWidth > 1 && maxWidth != Double.POSITIVE_INFINITY)
		{
			upperSigmaThreshold = Filter.getUpperLimit(maxWidth);
			widthEnabled = upperSigmaThreshold != Float.POSITIVE_INFINITY;
		}
		else
			upperSigmaThreshold = Float.POSITIVE_INFINITY;
		if (minWidth < 1)
		{
			widthEnabled = true;
			lowerSigmaThreshold = (float) minWidth;
		}
		else
			lowerSigmaThreshold = 0f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#getFilterSetupFlags()
	 */
	@Override
	public int getFilterSetupFlags() throws IllegalStateException
	{
		return (widthEnabled) ? 0 : IDirectFilter.NO_WIDTH;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#accept(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean accept(PeakResult peak)
	{
		final float sd = calculator.getStandardDeviation(peak.getParameters());
		return sd <= upperSigmaThreshold && sd >= lowerSigmaThreshold;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#getValidationFlags()
	 */
	@Override
	public int getValidationFlags()
	{
		return V_X_SD_FACTOR;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#validate(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (widthEnabled)
		{
			if (peak.getXSDFactor() > upperSigmaThreshold || peak.getXSDFactor() < lowerSigmaThreshold)
				return V_X_SD_FACTOR;
		}
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
		return "Filter results using a width range. (Width is relative to initial peak width.)";
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 2;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
		switch (index)
		{
			case 0:
				return minWidth;
			default:
				return maxWidth;
		}
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
		switch (index)
		{
			case 0:
				return WidthFilter2.DEFAULT_MIN_INCREMENT;
			default:
				return WidthFilter.DEFAULT_INCREMENT;
		}
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
		switch (index)
		{
			case 0:
				return ParameterType.MIN_WIDTH;
			default:
				return ParameterType.MAX_WIDTH;
		}
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
		switch (index)
		{
			case 0:
				return new WidthFilter2(updateParameter(minWidth, delta, DEFAULT_MIN_RANGE), maxWidth);
			default:
				return new WidthFilter2(minWidth, updateParameter(maxWidth, delta, WidthFilter.DEFAULT_RANGE));
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new WidthFilter2(parameters[0], parameters[1]);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, minWidth);
		setMax(parameters, 1, maxWidth);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return (index == 1) ? 1 : -1;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT };
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	@Override
	public double[] mutationStepRange()
	{
		return new double[] { WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getSignal()
	 */
	@Override
	public double getSignal()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getSNR()
	 */
	@Override
	public double getSNR()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getMinWidth()
	 */
	@Override
	public double getMinWidth()
	{
		return minWidth;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getMaxWidth()
	 */
	@Override
	public double getMaxWidth()
	{
		return maxWidth;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getShift()
	 */
	@Override
	public double getShift()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getEShift()
	 */
	@Override
	public double getEShift()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getPrecision()
	 */
	@Override
	public double getPrecision()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getPrecisionType()
	 */
	@Override
	public PrecisionType getPrecisionType()
	{
		return PrecisionType.NONE;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getMinZ()
	 */
	@Override
	public double getMinZ()
	{
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IMultiFilter#getMaxZ()
	 */
	@Override
	public double getMaxZ()
	{
		return 0;
	}
}
