package gdsc.smlm.results.filter;

import java.util.Arrays;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift, precision and z-depth.
 */
public class MultiFilter extends DirectFilter implements IMultiFilter
{
	@XStreamAsAttribute
	final double signal;
	@XStreamAsAttribute
	final float snr;
	@XStreamAsAttribute
	final double minWidth;
	@XStreamAsAttribute
	final double maxWidth;
	@XStreamAsAttribute
	final double shift;
	@XStreamAsAttribute
	final double eshift;
	@XStreamAsAttribute
	final double precision;
	@XStreamAsAttribute
	final float minZ;
	@XStreamAsAttribute
	final float maxZ;

	@XStreamOmitField
	float signalThreshold;
	@XStreamOmitField
	float lowerSigmaThreshold;
	@XStreamOmitField
	float upperSigmaThreshold;
	@XStreamOmitField
	float offsetx;
	@XStreamOmitField
	float offsety;
	@XStreamOmitField
	float eoffset;
	@XStreamOmitField
	double variance;
	@XStreamOmitField
	Gaussian2DPeakResultCalculator calculator;
	@XStreamOmitField
	boolean widthEnabled;
	@XStreamOmitField
	int flags;
	@XStreamOmitField
	boolean zEnabled;
	@XStreamOmitField
	MultiFilterComponentSet components = null;
	@XStreamOmitField
	MultiFilterComponentSet components_NoWidth_Shift = null;
	@XStreamOmitField
	MultiFilterComponentSet components_Width_Shift = null;
	@XStreamOmitField
	MultiFilterComponentSet components_NoWidth_NoShift = null;
	@XStreamOmitField
	MultiFilterComponentSet components_Width_NoShift = null;
	@XStreamOmitField
	MultiFilterComponentSet components_Shift0 = null;
	@XStreamOmitField
	private int filterSetupFlags;
	@XStreamOmitField
	private FilterSetupData[] filterSetupData;

	public MultiFilter(double signal, float snr, double minWidth, double maxWidth, double shift, double eshift,
			double precision, float minZ, float maxZ)
	{
		this.signal = Math.max(0, signal);
		this.snr = Math.max(0, snr);
		// Only swap if max width is enabled
		if (maxWidth != 0 && maxWidth < minWidth)
		{
			final double f = maxWidth;
			maxWidth = minWidth;
			minWidth = f;
		}
		this.minWidth = Math.max(0, minWidth);
		this.maxWidth = Math.max(0, maxWidth);
		this.shift = Math.max(0, shift);
		this.eshift = Math.max(0, eshift);
		this.precision = Math.max(0, precision);
		this.minZ = minZ;
		this.maxZ = maxZ;
	}

	@Override
	protected String generateName()
	{
		return String.format(
				"Multi: Signal=%.1f, SNR=%.1f, Width=%.2f-%.2f, Shift=%.2f, EShift=%.2f, Precision=%.1f, Width=%.2f-%.2f",
				signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		setupCalculator(peakResults);

		signalThreshold = (float) (signal);

		// Set the width limit
		lowerSigmaThreshold = 0;
		upperSigmaThreshold = Float.POSITIVE_INFINITY;
		// Set the shift limit. The calculator can support both 1/2 axis widths
		// when extracting the Standard Deviation from the parameters.
		double[] s = PSFHelper.getGaussian2DWxWy(peakResults.getPSF());
		double s0 = (s[0] == s[1]) ? s[0] : Gaussian2DPeakResultHelper.getStandardDeviation(s[0], s[1]);
		lowerSigmaThreshold = (float) (s0 * minWidth);
		upperSigmaThreshold = Filter.getUpperLimit(s0 * maxWidth);
		offsetx = getUpperLimit(s[0] * shift);
		offsety = getUpperLimit(s[1] * shift);
		// Convert to squared distance
		eoffset = Filter.getUpperSquaredLimit(s0 * eshift);

		// Configure the precision limit
		variance = Filter.getDUpperSquaredLimit(precision);

		zEnabled = isZEnabled();
	}

	private boolean isZEnabled()
	{
		return (minZ != 0 || maxZ != 0) && minZ <= maxZ;
	}

	protected void setupCalculator(MemoryPeakResults peakResults)
	{
		calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(),
				Gaussian2DPeakResultHelper.LSE_PRECISION);
	}

	@Override
	public void setup()
	{
		filterSetupFlags = 0;
		this.filterSetupData = null;
		setup(true, true, 0);
	}

	@Override
	public void setup(int flags)
	{
		filterSetupFlags = flags;
		this.filterSetupData = null;
		setup(!areSet(flags, DirectFilter.NO_WIDTH), !areSet(flags, DirectFilter.NO_SHIFT),
				// Pass through the flags that are recognised
				flags & (DirectFilter.XY_WIDTH | DirectFilter.NO_Z));
	}

	private void setup(final boolean widthEnabled, final boolean shiftEnabled, final int flags)
	{
		// Note: The filter caches the combinations that are likely to be turned on/off:
		// width filtering and shift filtering
		// Other filters related to the PSF (XY widths, z-depth) are assumed to be constant.
		
		if (components_Width_Shift == null || this.flags != flags)
		{
			// Store this in case the filter is setup with different flags 
			this.flags = flags;

			// Create the components we require
			final MultiFilterComponent[] components1 = new MultiFilterComponent[7];
			int s1 = 0;
			Class<? extends MultiFilterComponent> widthComponentClass = null;
			Class<? extends MultiFilterComponent> shiftComponentClass = null;

			// Current order of filter power obtained from BenchmarkFilterAnalysis:
			// SNR, Max Width, Precision, Shift, Min width
			if (isFiniteStrictlyPositive(snr))
			{
				components1[s1++] = new MultiFilterSNRComponent(snr);
			}
			if ((maxWidth > 1 && maxWidth != Double.POSITIVE_INFINITY) || (minWidth > 0 && minWidth < 1))
			{
				// Handle the width being 1/2 axis variable.
				if (areSet(flags, DirectFilter.XY_WIDTH))
					components1[s1++] = new MultiFilterXYWidthComponent(minWidth, maxWidth);
				else
					components1[s1++] = new MultiFilterWidthComponent(minWidth, maxWidth);
				widthComponentClass = components1[s1 - 1].getClass();
			}
			if (isFiniteStrictlyPositive(precision))
			{
				components1[s1++] = createPrecisionComponent();
			}
			if (isFiniteStrictlyPositive(shift))
			{
				components1[s1++] = new MultiFilterShiftComponent(shift);
				shiftComponentClass = components1[s1 - 1].getClass();
			}
			if (isFiniteStrictlyPositive(signal))
			{
				components1[s1++] = new MultiFilterSignalComponent(signal);
			}
			if (isFiniteStrictlyPositive(eshift))
			{
				components1[s1++] = new MultiFilterEShiftComponent(eshift);
			}
			if (isZEnabled() && !areSet(flags, DirectFilter.NO_Z))
			{
				components1[s1++] = new MultiFilterZComponent(minZ, maxZ);
			}

			final MultiFilterComponent[] components2 = MultiFilter.remove(components1, s1, widthComponentClass);
			final MultiFilterComponent[] components3 = MultiFilter.remove(components1, s1, shiftComponentClass);

			final MultiFilterComponent[] components4 = MultiFilter.remove(components2, components2.length,
					shiftComponentClass);

			components_Width_Shift = MultiFilterComponentSetFactory.create(components1, s1);
			components_NoWidth_Shift = MultiFilterComponentSetFactory.create(components2, components2.length);
			components_Width_NoShift = MultiFilterComponentSetFactory.create(components3, components3.length);
			components_NoWidth_NoShift = MultiFilterComponentSetFactory.create(components4, components4.length);

			MultiFilterComponent[] data = new MultiFilterComponent[components3.length + 1];
			System.arraycopy(components3, 0, data, 1, components3.length);
			components_Shift0 = MultiFilterComponentSetFactory.create(data, data.length);
		}

		if (widthEnabled)
		{
			components = (shiftEnabled) ? components_Width_Shift : components_NoWidth_Shift;
		}
		else
		{
			components = (shiftEnabled) ? components_NoWidth_Shift : components_NoWidth_NoShift;
		}

		//		// This is the legacy support for all components together
		//		this.widthEnabled = widthEnabled;
		//		signalThreshold = (float) signal;
		//		if (widthEnabled)
		//		{
		//			lowerSigmaThreshold = (float) minWidth;
		//			upperSigmaThreshold = Filter.getUpperLimit(maxWidth);
		//		}
		//		offset = Filter.getUpperSquaredLimit(shift);
		//		eoffset = Filter.getUpperSquaredLimit(eshift);
		//		variance = Filter.getDUpperSquaredLimit(precision);
	}

	protected MultiFilterComponent createPrecisionComponent()
	{
		return new MultiFilterVarianceComponent(precision);
	}

	static MultiFilterComponent[] remove(MultiFilterComponent[] in, int size,
			Class<? extends MultiFilterComponent> clazz)
	{
		if (clazz == null)
			return in;
		MultiFilterComponent[] out = new MultiFilterComponent[size];
		int length = 0;
		for (int i = 0; i < size; i++)
		{
			if (in[i].getClass() != clazz)
				out[length++] = in[i];
		}
		return Arrays.copyOf(out, length);
	}

	@Override
	public void setup(int flags, FilterSetupData... filterSetupData)
	{
		setup(flags);
		for (int i = filterSetupData.length; i-- > 0;)
		{
			if (filterSetupData[i] instanceof ShiftFilterSetupData)
			{
				this.filterSetupData = getFilterSetupData(filterSetupData[i]);
				double shift = ((ShiftFilterSetupData) filterSetupData[i]).shift;
				if (shift > 0)
				{
					components_Shift0.replace0(new MultiFilterShiftComponent(shift));
					components = components_Shift0;
				}
				else
				{
					components = components_Width_NoShift;
				}
				return;
			}
		}
	}

	@Override
	public int getFilterSetupFlags() throws IllegalStateException
	{
		return filterSetupFlags;
	}

	@Override
	public FilterSetupData[] getFilterSetupData() throws IllegalStateException
	{
		return filterSetupData;
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		// Current order of filter power obtained from BenchmarkFilterAnalysis:
		// SNR, Max Width, Precision, Shift, Min width

		if (peak.getSNR() < this.snr)
			return false;

		// Width
		final float sd = calculator.getStandardDeviation(peak.getParameters());
		if (sd > upperSigmaThreshold || sd < lowerSigmaThreshold)
			return false;

		// Precision
		if (getVariance(peak) > variance)
			return false;

		// Shift
		if (Math.abs(peak.getXShift()) > offsetx || Math.abs(peak.getYShift()) > offsety)
			return false;

		if (peak.getSignal() < signalThreshold)
			return false;

		if (zEnabled && (peak.getZPosition() < minZ || peak.getZPosition() > maxZ))
			return false;

		// Euclidian shift
		final float dx = peak.getXPosition();
		final float dy = peak.getYPosition();
		if (dx * dx + dy * dy > eoffset)
			return false;

		return true;
	}

	protected double getVariance(PeakResult peak)
	{
		return calculator.getLSEVariance(peak.getParameters(), peak.getNoise());
	}

	public int getValidationFlags()
	{
		if (components == null)
			setup();
		return components.getValidationFlags();
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		return components.validate(peak);

		//		// This is the legacy support for all components together
		//		
		//		// Current order of filter power obtained from BenchmarkFilterAnalysis:
		//		// Precision, Max Width, SNR, Shift, Min width
		//
		//		if (peak.getLocationVariance() > variance)
		//			return V_LOCATION_VARIANCE;
		//
		//		if (widthEnabled)
		//		{
		//			final float xsdf = peak.getXSDFactor();
		//			if (xsdf > upperSigmaThreshold || xsdf < lowerSigmaThreshold)
		//				return V_X_SD_FACTOR;
		//		}
		//
		//		if (peak.getSNR() < this.snr)
		//			return V_SNR;
		//
		//		// Shift
		//		final float xs2 = peak.getXRelativeShift2();
		//		if (xs2 > offset)
		//			return V_X_RELATIVE_SHIFT;
		//		final float ys2 = peak.getYRelativeShift2();
		//		if (ys2 > offset)
		//			return V_Y_RELATIVE_SHIFT;
		//
		//		if (peak.getPhotons() < signal)
		//			return V_PHOTONS;
		//
		//		// Euclidian shift
		//		if (xs2 + ys2 > eoffset)
		//			return V_X_RELATIVE_SHIFT | V_Y_RELATIVE_SHIFT;
		//
		//		return 0;
	}

	@Override
	public double getNumericalValue()
	{
		// This is not the first parameter so override
		return snr;
	}

	@Override
	public String getNumericalValueName()
	{
		// This is not the first parameter so override
		return ParameterType.SNR.toString();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using multiple thresholds: Signal, SNR, width, shift, Euclidian shift, precision and Z-depth";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 9;
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
				return signal;
			case 1:
				return snr;
			case 2:
				return minWidth;
			case 3:
				return maxWidth;
			case 4:
				return shift;
			case 5:
				return eshift;
			case 6:
				return precision;
			case 7:
				return minZ;
			default:
				return maxZ;
		}
	}

	@Override
	public double[] getParameters()
	{
		return new double[] { signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ };
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
				return SignalFilter.DEFAULT_INCREMENT;
			case 1:
				return SNRFilter.DEFAULT_INCREMENT;
			case 2:
				return WidthFilter2.DEFAULT_MIN_INCREMENT;
			case 3:
				return WidthFilter.DEFAULT_INCREMENT;
			case 4:
				return ShiftFilter.DEFAULT_INCREMENT;
			case 5:
				return EShiftFilter.DEFAULT_INCREMENT;
			case 6:
				return PrecisionFilter.DEFAULT_INCREMENT;
			case 7:
				return ZCoordinateFilter.DEFAULT_INCREMENT;
			default:
				return ZCoordinateFilter.DEFAULT_INCREMENT;
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
				return ParameterType.SIGNAL;
			case 1:
				return ParameterType.SNR;
			case 2:
				return ParameterType.MIN_WIDTH;
			case 3:
				return ParameterType.MAX_WIDTH;
			case 4:
				return ParameterType.SHIFT;
			case 5:
				return ParameterType.ESHIFT;
			case 6:
				return getPrecisionParamaterType();
			case 7:
				return ParameterType.MIN_Z;
			default:
				return ParameterType.MAX_Z;
		}
	}

	protected ParameterType getPrecisionParamaterType()
	{
		return ParameterType.PRECISION;
	}

	static double[] defaultRange = new double[] { SignalFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE,
			WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE,
			EShiftFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE, ZCoordinateFilter.DEFAULT_RANGE,
			ZCoordinateFilter.DEFAULT_RANGE };

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#adjustParameter(int, double)
	 */
	@Override
	public Filter adjustParameter(int index, double delta)
	{
		checkIndex(index);
		double[] params = new double[] { signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ };
		params[index] = updateParameter(params[index], delta, defaultRange[index]);
		return new MultiFilter(params[0], (float) params[1], params[2], params[3], params[4], params[5], params[6],
				(float) params[7], (float) params[8]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new MultiFilter(parameters[0], (float) parameters[1], parameters[2], parameters[3], parameters[4],
				parameters[5], parameters[6], (float) parameters[7], (float) parameters[8]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, signal);
		setMin(parameters, 1, snr);
		setMin(parameters, 2, minWidth);
		setMax(parameters, 3, maxWidth);
		setMax(parameters, 4, shift);
		setMax(parameters, 5, eshift);
		setMax(parameters, 6, precision);
		setMin(parameters, 7, minZ);
		setMax(parameters, 8, maxZ);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return (index < 3 || index == 7) ? -1 : 1;
	}

	/**
	 * Compare to the other filter, count the number of weakest parameters. If negative then this filter has more weak
	 * parameters. If positive then this filter has less weak parameters. If the same or the number of parameters do not
	 * match then return 0. If the other filter is null return -1.
	 * 
	 * @param o
	 *            The other filter
	 * @return the count difference
	 */
	public int weakest(MultiFilter o)
	{
		if (o == null)
			return -1;

		// Count the number of weakest
		//@formatter:off
		return 
			compareMin(signal, o.signal) +
    		compareMin(snr, o.snr) +
    		compareMin(minWidth, o.minWidth) +
    		compareMax(maxWidth, o.maxWidth) +
    		compareMax(shift, o.shift) +
    		compareMax(eshift, o.eshift) + 
    		compareMax(precision, o.precision) +
    		compareMin(minZ, o.minZ) +
    		compareMax(maxZ, o.maxZ);
		//@formatter:on
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, WidthFilter.UPPER_LIMIT,
				WidthFilter.UPPER_LIMIT, ShiftFilter.UPPER_LIMIT, EShiftFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT,
				ZCoordinateFilter.UPPER_LIMIT, ZCoordinateFilter.UPPER_LIMIT };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return defaultRange;
	}

	public double getSignal()
	{
		return signal;
	}

	public double getSNR()
	{
		return snr;
	}

	public double getMinWidth()
	{
		return minWidth;
	}

	public double getMaxWidth()
	{
		return maxWidth;
	}

	public double getShift()
	{
		return shift;
	}

	public double getEShift()
	{
		return eshift;
	}

	public double getPrecision()
	{
		return precision;
	}

	public PrecisionType getPrecisionType()
	{
		return PrecisionType.ESTIMATE;
	}

	public double getMinZ()
	{
		return minZ;
	}

	public double getMaxZ()
	{
		return maxZ;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#initialiseState()
	 */
	@Override
	protected void initialiseState()
	{
		// This is run after a clone() occurs.
		// Q. Can the setup state be maintained?

		//components_Width_Shift = null;

		// Replace any object that is manipulated by the instance
		if (components_Shift0 != null)
		{
			boolean update = components == components_Shift0;
			components_Shift0 = components_Shift0.clone();
			if (update)
				components = components_Shift0;
		}
	}
}