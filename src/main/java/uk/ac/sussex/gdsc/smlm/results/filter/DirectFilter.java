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

import com.thoughtworks.xstream.annotations.XStreamOmitField;

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
    @XStreamOmitField
    private float strength = Float.NaN;

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#setup()
     */
    @Override
    public void setup()
    {
        // Do nothing
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#setup(int)
     */
    @Override
    public void setup(final int flags)
    {
        // Do nothing
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#setup(int,
     * uk.ac.sussex.gdsc.smlm.results.filter.FilterSetupData[])
     */
    @Override
    public void setup(int flags, FilterSetupData... filterSetupData)
    {
        // Do nothing
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#getFlags()
     */
    @Override
    public int getFilterSetupFlags() throws IllegalStateException
    {
        return 0;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#getFilterSetupData()
     */
    @Override
    public FilterSetupData[] getFilterSetupData() throws IllegalStateException
    {
        return null;
    }

    /**
     * Convenience method to convert variable input data to ann array.
     *
     * @param data
     *            the data
     * @return the filter setup data
     */
    protected static FilterSetupData[] getFilterSetupData(FilterSetupData... data)
    {
        return data;
    }

    /**
     * Check if all of the given bits are set in the flags.
     *
     * @param flags
     *            the flags
     * @param bits
     *            the bits
     * @return True if all are set
     */
    public static boolean areSet(final int flags, final int bits)
    {
        return (flags & bits) == bits;
    }

    /**
     * Check if any of the given bits are set in the flags.
     *
     * @param flags
     *            the flags
     * @param bits
     *            the bits
     * @return True if any are set
     */
    public static boolean anySet(final int flags, final int bits)
    {
        return (flags & bits) != 0;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#accept(uk.ac.sussex.gdsc.smlm.results.filter.
     * PreprocessedPeakResult)
     */
    @Override
    final public boolean accept(final PreprocessedPeakResult peak)
    {
        return (result = validate(peak)) == 0;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#validate(uk.ac.sussex.gdsc.smlm.results.filter.
     * PreprocessedPeakResult)
     */
    @Override
    public abstract int validate(final PreprocessedPeakResult peak);

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.Filter#getFilterType()
     */
    @Override
    public FilterType getFilterType()
    {
        return FilterType.DIRECT;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#getResult()
     */
    @Override
    public int getResult()
    {
        return result;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.filter.IDirectFilter#copy()
     */
    @Override
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
        final StringBuilder sb = new StringBuilder();
        //@formatter:off
		if (areSet(flags, V_AMPLITUDE))	             append(sb, "Amplitude",        peak.getAmplitude());
		if (areSet(flags, V_PHOTONS))                append(sb, "Signal",           peak.getSignal());
		if (areSet(flags, V_SNR))                    append(sb, "SNR",              peak.getSNR());
		if (areSet(flags, V_NOISE))                  append(sb, "Noise",            peak.getNoise());
		if (areSet(flags, V_LOCATION_VARIANCE))      append(sb, "Precision",        Math.sqrt(peak.getLocationVariance()));
		if (areSet(flags, V_LOCATION_VARIANCE2))     append(sb, "Precision2",       Math.sqrt(peak.getLocationVariance2()));
		if (areSet(flags, V_LOCATION_VARIANCE_CRLB)) append(sb, "Precision CRLB",   Math.sqrt(peak.getLocationVarianceCRLB()));
		if (areSet(flags, V_SD))                     append(sb, "SD",               peak.getSD());
		if (areSet(flags, V_BACKGROUND))             append(sb, "Background",       peak.getBackground());
		if (areSet(flags, V_AMPLITUDE))              append(sb, "Amplitude",        peak.getAmplitude());
		if (areSet(flags, V_ANGLE))                  append(sb, "Angle",            peak.getAngle());
		if (areSet(flags, V_X))                      append(sb, "X",                peak.getX());
		if (areSet(flags, V_Y))                      append(sb, "Y",                peak.getY());
		if (areSet(flags, V_Z))                      append(sb, "Z",                peak.getZ());
		if (areSet(flags, V_X_RELATIVE_SHIFT))       append(sb, "X Relative Shift", Math.sqrt(peak.getXRelativeShift2()));
		if (areSet(flags, V_Y_RELATIVE_SHIFT))       append(sb, "Y Relative Shift", Math.sqrt(peak.getYRelativeShift2()));
		if (areSet(flags, V_X_SD))                   append(sb, "X SD",             peak.getXSD());
		if (areSet(flags, V_Y_SD))                   append(sb, "Y SD",             peak.getYSD());
		if (areSet(flags, V_X_SD_FACTOR))            append(sb, "X SD Factor",      peak.getXSDFactor());
		if (areSet(flags, V_Y_SD_FACTOR))            append(sb, "Y SD Factor",      peak.getYSDFactor());
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
        final StringBuilder sb = new StringBuilder();
        //@formatter:off
		if (areSet(flags, V_AMPLITUDE))	             append(sb, "Amplitude");
		if (areSet(flags, V_PHOTONS))                append(sb, "Signal");
		if (areSet(flags, V_SNR))                    append(sb, "SNR");
		if (areSet(flags, V_NOISE))                  append(sb, "Noise");
		if (areSet(flags, V_LOCATION_VARIANCE))      append(sb, "Precision");
		if (areSet(flags, V_LOCATION_VARIANCE2))     append(sb, "Precision2");
		if (areSet(flags, V_LOCATION_VARIANCE_CRLB)) append(sb, "Precision CRLB");
		if (areSet(flags, V_SD))                     append(sb, "SD");
		if (areSet(flags, V_BACKGROUND))             append(sb, "Background");
		if (areSet(flags, V_AMPLITUDE))              append(sb, "Amplitude");
		if (areSet(flags, V_ANGLE))                  append(sb, "Angle");
		if (areSet(flags, V_X))                      append(sb, "X");
		if (areSet(flags, V_Y))                      append(sb, "Y");
		if (areSet(flags, V_Z))                      append(sb, "Z");
		if (areSet(flags, V_X_RELATIVE_SHIFT))       append(sb, "X Relative Shift");
		if (areSet(flags, V_Y_RELATIVE_SHIFT))       append(sb, "Y Relative Shift");
		if (areSet(flags, V_X_SD))                   append(sb, "X SD");
		if (areSet(flags, V_Y_SD))                   append(sb, "Y SD");
		if (areSet(flags, V_X_SD_FACTOR))            append(sb, "X SD Factor");
		if (areSet(flags, V_Y_SD_FACTOR))            append(sb, "Y SD Factor");
		//@formatter:on
        return sb.toString();
    }

    private static void append(StringBuilder sb, String name)
    {
        if (sb.length() != 0)
            sb.append("; ");
        sb.append(name);
    }

    /**
     * Gets the strength. This is used in the {@link #weakestUnsafe(Filter)} method before the parameters are compared.
     *
     * @return the strength
     */
    public float getStrength()
    {
        return strength;
    }

    /**
     * Sets the strength. This is used in the {@link #weakestUnsafe(DirectFilter)} method before the parameters are
     * compared.
     *
     * @param strength
     *            the new strength
     */
    public void setStrength(float strength)
    {
        this.strength = strength;
    }

    /**
     * Compute strength using the limits on the parameters. Strength is computed using the distance from the lower/upper
     * bounds inside the range for each parameter. The bounds to choose is determined by the method
     * {@link #lowerBoundOrientation(int)}.
     * <p>
     * Warning: No checks are made for the input arrays to be null or incorrect length.
     * <p>
     * If lower is equal or above upper then this index is ignored.
     *
     * @param lower
     *            the lower limit
     * @param upper
     *            the upper limit
     * @return the strength
     */
    public float computeStrength(double[] lower, double[] upper)
    {
        double s = 0;
        final double[] p = getParameters();
        for (int i = 0; i < p.length; i++)
        {
            final double range = upper[i] - lower[i];
            if (range <= 0)
                // Ignore this as it will produce bad strength data
                continue;
            final int o = lowerBoundOrientation(i);
            if (o < 0)
                s += (p[i] - lower[i]) / range;
            else if (o > 0)
                s += (upper[i] - p[i]) / range;
        }
        return (float) s;
    }

    /**
     * Get the lower bound orientation for the given parameter index. This is the orientation of the bound for the
     * weakest parameter. E.g. if a parameter must be low to be weak, the orientation is -1.
     *
     * @param index
     *            the index
     * @return the lower bound orientation.
     */
    public int lowerBoundOrientation(int index)
    {
        return -1;
    }

    /**
     * Compare to the other filter using the strength property and return the weakest. If equal (or the strength is not
     * set) then default to the {@link #weakest(Filter)} method.
     * <p>
     * This method does not check for null or if the other filter has a different number of parameters.
     *
     * @param o
     *            The other filter
     * @return the count difference
     */
    public int weakestUnsafe(DirectFilter o)
    {
        //		// This should not happen if used correctly
        //		if (Float.isNaN(strength) || Float.isNaN(o.strength))
        //		{
        //			System.out.println("No strength (nan)");
        //			return super.weakestUnsafe(o);
        //		}
        //		if (Float.isInfinite(strength) || Float.isInfinite(o.strength))
        //		{
        //			System.out.println("No strength (inf)");
        //			return super.weakestUnsafe(o);
        //		}

        //System.out.println("weakestUnsafe");

        //int result = super.weakestUnsafe(o);

        if (this.strength < o.strength)
            //if (result >= 0)
            //	System.out.println("Strength not same as weakest");
            return -1;
        if (this.strength > o.strength)
            //if (result <= 0)
            //	System.out.println("Strength not same as weakest");
            return 1;
        return super.weakestUnsafe(o);
    }
}
