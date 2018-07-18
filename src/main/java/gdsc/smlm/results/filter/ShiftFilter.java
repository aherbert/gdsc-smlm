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
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Filter results using a X/Y coordinate shift. This filter requires that the result X and Y coordinates are reported
 * relative to their initial positions.
 */
public class ShiftFilter extends DirectFilter implements IMultiFilter
{
	/** The default increment. Used for {@link gdsc.smlm.ga.Chromosome} interface. */
	public static final double DEFAULT_INCREMENT = 0.05;
	/** The default range. Used for {@link gdsc.smlm.ga.Chromosome} interface. */
	public static final double DEFAULT_RANGE = 10;
	/** The default limit. Used for {@link gdsc.smlm.ga.Chromosome} interface. */
	public static final double UPPER_LIMIT = 5;

	@XStreamAsAttribute
	private final double shift;
	@XStreamOmitField
	private float offsetx;
	@XStreamOmitField
	private float offsety;
	@XStreamOmitField
	private float shift2;
	@XStreamOmitField
	private boolean shiftEnabled;

	/**
	 * Instantiates a new shift filter.
	 *
	 * @param shift
	 *            the shift
	 */
	public ShiftFilter(double shift)
	{
		this.shift = Math.max(0, shift);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		// Set the shift limit
		final double[] s = PSFHelper.getGaussian2DWxWy(peakResults.getPSF());
		offsetx = getUpperLimit(s[0] * shift);
		offsety = getUpperLimit(s[1] * shift);
	}

	@Override
	public void setup()
	{
		setup(shift);
	}

	@Override
	public void setup(int flags)
	{
		if (areSet(flags, IDirectFilter.NO_SHIFT))
			shiftEnabled = false;
		else
			setup(shift);
	}

	@Override
	public void setup(int flags, FilterSetupData... filterSetupData)
	{
		if (areSet(flags, IDirectFilter.NO_SHIFT))
		{
			shiftEnabled = false;
			return;
		}

		for (int i = filterSetupData.length; i-- > 0;)
			if (filterSetupData[i] instanceof ShiftFilterSetupData)
			{
				setup(((ShiftFilterSetupData) filterSetupData[i]).shift);
				return;
			}
		// Default
		setup(shift);
	}

	private void setup(final double shift)
	{
		shift2 = getUpperSquaredLimit(shift);
		shiftEnabled = (shift2 != Float.POSITIVE_INFINITY);
	}

	@Override
	public int getFilterSetupFlags() throws IllegalStateException
	{
		return (shiftEnabled) ? 0 : IDirectFilter.NO_SHIFT;
	}

	@Override
	public FilterSetupData[] getFilterSetupData() throws IllegalStateException
	{
		if (shiftEnabled && shift2 != Float.POSITIVE_INFINITY)
		{
			if (shift2 == getUpperSquaredLimit(shift))
				// This is the default so ignore
				return null;
			return getFilterSetupData(new ShiftFilterSetupData(Math.sqrt(shift2)));
		}
		return null;
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		return Math.abs(peak.getXShift()) <= offsetx && Math.abs(peak.getYShift()) <= offsety;
	}

	@Override
	public int getValidationFlags()
	{
		return V_X_RELATIVE_SHIFT | V_Y_RELATIVE_SHIFT;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (shiftEnabled)
		{
			if (peak.getXRelativeShift2() > shift2)
				return V_X_RELATIVE_SHIFT;
			if (peak.getYRelativeShift2() > shift2)
				return V_Y_RELATIVE_SHIFT;
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
		return "Filter results using a shift factor. (X/Y shift is relative to initial peak width.)";
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
		return shift;
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
		return DEFAULT_INCREMENT;
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
		return ParameterType.SHIFT;
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
		return new ShiftFilter(updateParameter(shift, delta, DEFAULT_RANGE));
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new ShiftFilter(parameters[0]);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMax(parameters, 0, shift);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { UPPER_LIMIT };
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
		return shift;
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
