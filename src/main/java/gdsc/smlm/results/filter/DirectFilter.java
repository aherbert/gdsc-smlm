package gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Support direct filtering of PreprocessedPeakResult objects.
 * <p>
 * The decision to support for filtering as both a DirectFilter and Filter at the same time is left to the implementing
 * class. It is not a requirement.
 */
public abstract class DirectFilter extends Filter implements IDirectFilter
{
	@XStreamOmitField
	private int result = 0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#setup()
	 */
	public void setup()
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#setup(int)
	 */
	public void setup(final int flags)
	{
	}

	/**
	 * Check if all of the given bits are set in the flags
	 * 
	 * @param flags
	 * @param bits
	 * @return True if all are set
	 */
	public static boolean areSet(final int flags, final int bits)
	{
		return (flags & bits) == bits;
	}

	/**
	 * Check if any of the given bits are set in the flags
	 * 
	 * @param flags
	 * @param bits
	 * @return True if any are set
	 */
	public static boolean anySet(final int flags, final int bits)
	{
		return (flags & bits) != 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#accept(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	final public boolean accept(final PreprocessedPeakResult peak)
	{
		return (result = validate(peak)) == 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#validate(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	public abstract int validate(final PreprocessedPeakResult peak);

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getFilterType()
	 */
	public FilterType getFilterType()
	{
		return FilterType.DIRECT;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#getResult()
	 */
	public int getResult()
	{
		return result;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#copy()
	 */
	public IDirectFilter copy()
	{
		return (IDirectFilter) clone();
	}

	/**
	 * Generate a status message using all the properties of the peak that failed validation
	 * 
	 * @param peak
	 *            The peak
	 * @param flags
	 *            The validation flags
	 * @return The message
	 */
	public static String getStatusMessage(PreprocessedPeakResult peak, int flags)
	{
		if (flags == 0)
			return "";
		StringBuilder sb = new StringBuilder();
		//@formatter:off
		if (areSet(flags, V_AMPLITUDE))	         append(sb, "Amplitude",        peak.getAmplitude());
		if (areSet(flags, V_PHOTONS))            append(sb, "Photon",           peak.getPhotons());
		if (areSet(flags, V_SNR))                append(sb, "SNR",              peak.getSNR());
		if (areSet(flags, V_NOISE))              append(sb, "Noise",            peak.getNoise());
		if (areSet(flags, V_LOCATION_VARIANCE))  append(sb, "Precision",        peak.getLocationVariance());
		if (areSet(flags, V_LOCATION_VARIANCE2)) append(sb, "Precision2",       peak.getLocationVariance2());
		if (areSet(flags, V_SD))                 append(sb, "SD",               peak.getSD());
		if (areSet(flags, V_BACKGROUND))         append(sb, "Background",       peak.getBackground());
		if (areSet(flags, V_AMPLITUDE))          append(sb, "Amplitude",        peak.getAmplitude());
		if (areSet(flags, V_ANGLE))              append(sb, "Angle",            peak.getAngle());
		if (areSet(flags, V_X))                  append(sb, "X",                peak.getX());
		if (areSet(flags, V_Y))                  append(sb, "Y",                peak.getY());
		if (areSet(flags, V_X_RELATIVE_SHIFT))   append(sb, "X Relative Shift", peak.getXRelativeShift2());
		if (areSet(flags, V_Y_RELATIVE_SHIFT))   append(sb, "Y Relative Shift", peak.getYRelativeShift2());
		if (areSet(flags, V_X_SD))               append(sb, "X SD",             peak.getXSD());
		if (areSet(flags, V_Y_SD))               append(sb, "Y SD",             peak.getYSD());
		if (areSet(flags, V_X_SD_FACTOR))        append(sb, "X SD Factor",      peak.getXSDFactor());
		if (areSet(flags, V_Y_SD_FACTOR))        append(sb, "Y SD Factor",      peak.getYSDFactor());
		//@formatter:on
		return sb.toString();
	}

	private static void append(StringBuilder sb, String name, double value)
	{
		if (sb.length() != 0)
			sb.append("; ");
		sb.append(name);
		sb.append('=');
		sb.append(value);
	}

	/**
	 * Generate a message using all flags that are set
	 * 
	 * @param flags
	 *            The validation flags
	 * @return The message
	 */
	public static String getFlagMessage(int flags)
	{
		if (flags == 0)
			return "";
		StringBuilder sb = new StringBuilder();
		//@formatter:off
		if (areSet(flags, V_AMPLITUDE))	         append(sb, "Amplitude");
		if (areSet(flags, V_PHOTONS))            append(sb, "Photon");
		if (areSet(flags, V_SNR))                append(sb, "SNR");
		if (areSet(flags, V_NOISE))              append(sb, "Noise");
		if (areSet(flags, V_LOCATION_VARIANCE))  append(sb, "Precision");
		if (areSet(flags, V_LOCATION_VARIANCE2)) append(sb, "Precision2");
		if (areSet(flags, V_SD))                 append(sb, "SD");
		if (areSet(flags, V_BACKGROUND))         append(sb, "Background");
		if (areSet(flags, V_AMPLITUDE))          append(sb, "Amplitude");
		if (areSet(flags, V_ANGLE))              append(sb, "Angle");
		if (areSet(flags, V_X))                  append(sb, "X");
		if (areSet(flags, V_Y))                  append(sb, "Y");
		if (areSet(flags, V_X_RELATIVE_SHIFT))   append(sb, "X Relative Shift");
		if (areSet(flags, V_Y_RELATIVE_SHIFT))   append(sb, "Y Relative Shift");
		if (areSet(flags, V_X_SD))               append(sb, "X SD");
		if (areSet(flags, V_Y_SD))               append(sb, "Y SD");
		if (areSet(flags, V_X_SD_FACTOR))        append(sb, "X SD Factor");
		if (areSet(flags, V_Y_SD_FACTOR))        append(sb, "Y SD Factor");
		//@formatter:on
		return sb.toString();
	}

	private static void append(StringBuilder sb, String name)
	{
		if (sb.length() != 0)
			sb.append("; ");
		sb.append(name);
	}
}