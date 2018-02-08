package gdsc.smlm.results.filter;

import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift, precision and z-depth.
 * Calculates the precision using the fitted X/Y parameter deviations.
 */
public class MultiFilterCRLB extends MultiFilter implements IMultiFilter
{
	public MultiFilterCRLB(double signal, float snr, double minWidth, double maxWidth, double shift, double eshift,
			double precision, float minZ, float maxZ)
	{
		super(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
	}

	@Override
	protected String generateName()
	{
		return String.format(
				"MultiCRLB: Signal=%.1f, SNR=%.1f, Width=%.2f-%.2f, Shift=%.2f, EShift=%.2f, Precision=%.1f, Width=%.2f-%.2f",
				signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
	}

	@Override
	protected void setupCalculator(MemoryPeakResults peakResults)
	{
		calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(), 0);
	}

	@Override
	protected MultiFilterComponent createPrecisionComponent()
	{
		return new MultiFilterVarianceCRLBComponent(precision);
	}
	

	@Override
	protected double getVariance(PeakResult peak)
	{
		if (peak.hasParameterDeviations())
		{
			float vx = peak.getParameterDeviation(PeakResult.X);
			float vy = peak.getParameterDeviation(PeakResult.Y);
			return (vx * vx + vy * vy);
		}
		return variance; // Return the current limit
	}	

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using multiple thresholds: Signal, SNR, width, shift, Euclidian shift, precision (uses fitted parameter variance) and Z-depth";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#requiresParameterDeviations()
	 */
	@Override
	public boolean requiresParameterDeviations()
	{
		return (precision != 0);
	}

	@Override
	protected ParameterType getPrecisionParamaterType()
	{
		return ParameterType.PRECISION_CRLB;
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
		return new MultiFilterCRLB(params[0], (float) params[1], params[2], params[3], params[4], params[5], params[6],
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
		return new MultiFilterCRLB(parameters[0], (float) parameters[1], parameters[2], parameters[3], parameters[4],
				parameters[5], parameters[6], (float) parameters[7], (float) parameters[8]);
	}

	@Override
	public PrecisionType getPrecisionType()
	{
		return PrecisionType.CRLB;
	}
}