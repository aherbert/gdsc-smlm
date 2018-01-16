package gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift and precision. Calculates the
 * precision using the true fitted background if a bias is provided.
 */
public class MultiFilter2 extends DirectFilter implements IMultiFilter
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

	@XStreamOmitField
	float signalThreshold;
	@XStreamOmitField
	float lowerSigmaThreshold;
	@XStreamOmitField
	float upperSigmaThreshold;
	@XStreamOmitField
	float offset;
	@XStreamOmitField
	float eoffset;
	@XStreamOmitField
	double variance;
	@XStreamOmitField
	boolean useBackground = false;
	@XStreamOmitField
	private Gaussian2DPeakResultCalculator calculator;
	@XStreamOmitField
	boolean widthEnabled;
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

	public MultiFilter2(double signal, float snr, double minWidth, double maxWidth, double shift, double eshift,
			double precision)
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
	}

	@Override
	protected String generateName()
	{
		return String.format("Multi2: Signal=%.1f, SNR=%.1f, Width=%.2f-%.2f, Shift=%.2f, EShift=%.2f, Precision2=%.1f",
				signal, snr, minWidth, maxWidth, shift, eshift, precision);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		try
		{
			calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(),
					Gaussian2DPeakResultHelper.LSE_PRECISION_X);
			useBackground = true;
		}
		catch (ConfigurationException e)
		{
			calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(),
					Gaussian2DPeakResultHelper.LSE_PRECISION);
			useBackground = false;
		}

		signalThreshold = (float) (signal);

		// Set the width limit
		lowerSigmaThreshold = 0;
		upperSigmaThreshold = Float.POSITIVE_INFINITY;
		// Set the shift limit
		double s = PSFHelper.getGaussian2DWx(peakResults.getPSF());
		lowerSigmaThreshold = (float) (s * minWidth);
		upperSigmaThreshold = Filter.getUpperLimit(s * maxWidth);
		offset = Filter.getUpperLimit(s * shift);
		// Convert to squared distance
		eoffset = Filter.getUpperSquaredLimit(s * eshift);

		// Configure the precision limit
		variance = Filter.getDUpperSquaredLimit(precision);
	}

	@Override
	public void setup()
	{
		setup(true, true);
	}

	@Override
	public void setup(int flags)
	{
		setup(!areSet(flags, DirectFilter.NO_WIDTH), !areSet(flags, DirectFilter.NO_SHIFT));
	}

	private void setup(final boolean widthEnabled, final boolean shiftEnabled)
	{
		if (components_Width_Shift == null)
		{
			// Create the components we require
			final MultiFilterComponent[] components1 = new MultiFilterComponent[6];
			int s1 = 0;

			// Current order of filter power obtained from BenchmarkFilterAnalysis:
			// SNR, Max Width, Precision, Shift, Min width
			if (snr != 0)
			{
				components1[s1++] = new MultiFilterSNRComponent(snr);
			}
			if (maxWidth != 0 || minWidth != 0)
			{
				components1[s1++] = new MultiFilterWidthComponent(minWidth, maxWidth);
			}
			if (precision != 0)
			{
				components1[s1++] = new MultiFilterVariance2Component(precision);
			}
			if (shift != 0)
			{
				components1[s1++] = new MultiFilterShiftComponent(shift);
			}
			if (signal != 0)
			{
				components1[s1++] = new MultiFilterSignalComponent(signal);
			}
			if (eshift != 0)
			{
				components1[s1++] = new MultiFilterEShiftComponent(eshift);
			}

			final MultiFilterComponent[] components2 = MultiFilter.remove(components1, s1,
					MultiFilterWidthComponent.class);
			final MultiFilterComponent[] components3 = MultiFilter.remove(components1, s1,
					MultiFilterShiftComponent.class);

			final MultiFilterComponent[] components4 = MultiFilter.remove(components2, components2.length,
					MultiFilterShiftComponent.class);

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

	@Override
	public void setup(FilterSetupData... filterSetupData)
	{
		setup(true, true);
		for (int i = filterSetupData.length; i-- > 0;)
		{
			if (filterSetupData[i] instanceof ShiftFilterSetupData)
			{
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
	public boolean accept(PeakResult peak)
	{
		// Current order of filter power obtained from BenchmarkFilterAnalysis:
		// SNR, Max Width, Precision, Shift, Min width

		if (SNRFilter.getSNR(peak) < this.snr)
			return false;

		// Width
		final float sd = calculator.getStandardDeviation(peak.getParameters());
		if (sd > upperSigmaThreshold || sd < lowerSigmaThreshold)
			return false;

		// Precision
		if (useBackground)
		{
			if (calculator.getLSEVariance(peak.getParameters()) > variance)
				return false;
		}
		else
		{
			if (calculator.getLSEVariance(peak.getParameters(), peak.getNoise()) > variance)
				return false;
		}

		// Shift
		if (Math.abs(peak.getXPosition()) > offset || Math.abs(peak.getYPosition()) > offset)
			return false;

		if (peak.getSignal() < signalThreshold)
			return false;

		// Euclidian shift
		final float dx = peak.getXPosition();
		final float dy = peak.getYPosition();
		if (dx * dx + dy * dy > eoffset)
			return false;

		return true;
	}

	public int getValidationFlags()
	{
		return components.getValidationFlags();
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		return components.validate(peak);

		//		// Current order of filter power obtained from BenchmarkFilterAnalysis:
		//		// Precision, Max Width, SNR, Shift, Min width
		//
		//		if (peak.getLocationVariance2() > variance)
		//			return V_LOCATION_VARIANCE2;
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
		return "Filter results using multiple thresholds: Signal, SNR, width, shift, Euclidian shift and precision (uses fitted background to set noise)";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 7;
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
			default:
				return precision;
		}
	}

	@Override
	public double[] getParameters()
	{
		return new double[] { signal, snr, minWidth, maxWidth, shift, eshift, precision };
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
			default:
				return PrecisionFilter.DEFAULT_INCREMENT;
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
			default:
				return ParameterType.PRECISION2;
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
		double[] params = new double[] { signal, snr, minWidth, maxWidth, shift, eshift, precision };
		params[index] = updateParameter(params[index], delta, MultiFilter.defaultRange[index]);
		return new MultiFilter2(params[0], (float) params[1], params[2], params[3], params[4], params[5], params[6]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new MultiFilter2(parameters[0], (float) parameters[1], parameters[2], parameters[3], parameters[4],
				parameters[5], parameters[6]);
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
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return (index >= 3) ? 1 : -1;
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
    		compareMax(precision, o.precision);
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
				WidthFilter.UPPER_LIMIT, ShiftFilter.UPPER_LIMIT, EShiftFilter.UPPER_LIMIT,
				PrecisionFilter.UPPER_LIMIT };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return MultiFilter.defaultRange;
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
		return PrecisionType.ESTIMATE_USING_LOCAL_BACKGROUND;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#initialiseState()
	 */
	@Override
	protected void initialiseState()
	{
		components_Width_Shift = null;
	}
}