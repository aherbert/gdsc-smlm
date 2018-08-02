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
package uk.ac.sussex.gdsc.smlm.data.config;

import uk.ac.sussex.gdsc.core.data.utils.AddMultiplyTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.IdentityTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.MultiplyAddTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.MultiplyTypeConverter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;

/**
 * Factory for creating unit converters
 */
public class UnitConverterFactory
{
    /**
     * Creates the converter.
     *
     * @param from
     *            the from
     * @param to
     *            the to
     * @param nmPerPixel
     *            the nm per pixel
     * @return the unit converter
     * @throws ConversionException
     *             if a converter cannot be created
     */
    public static TypeConverter<DistanceUnit> createConverter(DistanceUnit from, DistanceUnit to, double nmPerPixel)
            throws ConversionException
    {
        if (from == to)
            return new IdentityTypeConverter<>(from);
        if (from == null)
            throw new ConversionException("from must not be null");
        if (to == null)
            throw new ConversionException("to must not be null");

        switch (from)
        {
            case NM:
                switch (to)
                {
                    case PIXEL:
                        return new MultiplyTypeConverter<>(from, to, 1.0 / checkNmPerPixel(nmPerPixel));
                    case UM:
                        return new MultiplyTypeConverter<>(from, to, 1e-3);
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            case PIXEL:
                switch (to)
                {
                    case NM:
                        return new MultiplyTypeConverter<>(from, to, checkNmPerPixel(nmPerPixel));
                    case UM:
                        return new MultiplyTypeConverter<>(from, to, checkNmPerPixel(nmPerPixel) / 1e3);
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            case UM:
                switch (to)
                {
                    case NM:
                        return new MultiplyTypeConverter<>(from, to, 1e3);
                    case PIXEL:
                        return new MultiplyTypeConverter<>(from, to, 1e3 / checkNmPerPixel(nmPerPixel));
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            default:
                throw new ConversionException(from + " to " + to);
        }
    }

    private static double checkNmPerPixel(double nmPerPixel)
    {
        if (!(nmPerPixel > 0 && nmPerPixel <= java.lang.Double.MAX_VALUE))
            throw new ConversionException("nm/pixel must be positive");
        return nmPerPixel;
    }

    /**
     * Creates the converter.
     *
     * @param from
     *            the from
     * @param to
     *            the to
     * @param msPerFrame
     *            the ms per frame
     * @return the unit converter
     * @throws ConversionException
     *             if a converter cannot be created
     */
    public static TypeConverter<TimeUnit> createConverter(TimeUnit from, TimeUnit to, double msPerFrame)
            throws ConversionException
    {
        if (from == to)
            return new IdentityTypeConverter<>(from);
        if (from == null)
            throw new ConversionException("from must not be null");
        if (to == null)
            throw new ConversionException("to must not be null");

        switch (from)
        {
            case MILLISECOND:
                switch (to)
                {
                    case FRAME:
                        return new MultiplyTypeConverter<>(from, to, 1.0 / checkMsPerFrame(msPerFrame));
                    case SECOND:
                        return new MultiplyTypeConverter<>(from, to, 1e-3);
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            case FRAME:
                switch (to)
                {
                    case MILLISECOND:
                        return new MultiplyTypeConverter<>(from, to, checkMsPerFrame(msPerFrame));
                    case SECOND:
                        return new MultiplyTypeConverter<>(from, to, checkMsPerFrame(msPerFrame) / 1e3);
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            case SECOND:
                switch (to)
                {
                    case MILLISECOND:
                        return new MultiplyTypeConverter<>(from, to, 1e3);
                    case FRAME:
                        return new MultiplyTypeConverter<>(from, to, 1e3 / checkMsPerFrame(msPerFrame));
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            default:
                throw new ConversionException(from + " to " + to);
        }
    }

    private static double checkMsPerFrame(double msPerFrame)
    {
        if (!(msPerFrame > 0 && msPerFrame <= java.lang.Double.MAX_VALUE))
            throw new ConversionException("ms/frame must be positive");
        return msPerFrame;
    }

    /**
     * Creates the converter.
     *
     * @param from
     *            the from
     * @param to
     *            the to
     * @return the unit converter
     * @throws ConversionException
     *             if a converter cannot be created
     */
    public static TypeConverter<AngleUnit> createConverter(AngleUnit from, AngleUnit to) throws ConversionException
    {
        if (from == to)
            return new IdentityTypeConverter<>(from);
        if (from == null)
            throw new ConversionException("from must not be null");
        if (to == null)
            throw new ConversionException("to must not be null");

        switch (from)
        {
            case RADIAN:
                switch (to)
                {
                    case DEGREE:
                        return new MultiplyTypeConverter<>(from, to, 180.0 / Math.PI);
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            case DEGREE:
                switch (to)
                {
                    case RADIAN:
                        return new MultiplyTypeConverter<>(from, to, Math.PI / 180.0);
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            default:
                throw new ConversionException(from + " to " + to);
        }
    }

    /**
     * Creates the converter.
     *
     * @param from
     *            the from
     * @param to
     *            the to
     * @param offset
     *            the offset
     * @param countPerPhoton
     *            the count per photon
     * @return the unit converter
     * @throws ConversionException
     *             if a converter cannot be created
     */
    public static TypeConverter<IntensityUnit> createConverter(IntensityUnit from, IntensityUnit to, double offset,
            double countPerPhoton) throws ConversionException
    {
        if (from == to)
            return new IdentityTypeConverter<>(from);
        if (from == null)
            throw new ConversionException("from must not be null");
        if (to == null)
            throw new ConversionException("to must not be null");

        switch (from)
        {
            case COUNT:
                switch (to)
                {
                    case PHOTON:
                        return new AddMultiplyTypeConverter<>(from, to, checkOffset(offset, -1),
                                1.0 / checkCountPerPhoton(countPerPhoton));
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            case PHOTON:
                switch (to)
                {
                    case COUNT:
                        return new MultiplyAddTypeConverter<>(from, to, checkCountPerPhoton(countPerPhoton),
                                checkOffset(offset, 1));
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            default:
                throw new ConversionException(from + " to " + to);
        }
    }

    private static double checkOffset(double offset, int sign)
    {
        if (!(offset >= 0 && offset <= java.lang.Double.MAX_VALUE))
            throw new ConversionException("count/photon must be positive");
        return offset * sign;
    }

    private static double checkCountPerPhoton(double countPerPhoton)
    {
        if (!(countPerPhoton > 0 && countPerPhoton <= java.lang.Double.MAX_VALUE))
            throw new ConversionException("count/photon must be positive");
        return countPerPhoton;
    }

    /**
     * Creates the converter.
     *
     * @param from
     *            the from
     * @param to
     *            the to
     * @param countPerPhoton
     *            the count per photon
     * @return the unit converter
     * @throws ConversionException
     *             if a converter cannot be created
     */
    public static TypeConverter<IntensityUnit> createConverter(IntensityUnit from, IntensityUnit to,
            double countPerPhoton) throws ConversionException
    {
        if (from == to)
            return new IdentityTypeConverter<>(from);
        if (from == null)
            throw new ConversionException("from must not be null");
        if (to == null)
            throw new ConversionException("to must not be null");

        switch (from)
        {
            case COUNT:
                switch (to)
                {
                    case PHOTON:
                        return new MultiplyTypeConverter<>(from, to, 1.0 / checkCountPerPhoton(countPerPhoton));
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            case PHOTON:
                switch (to)
                {
                    case COUNT:
                        return new MultiplyTypeConverter<>(from, to, checkCountPerPhoton(countPerPhoton));
                    default:
                        throw new ConversionException(from + " to " + to);
                }

            default:
                throw new ConversionException(from + " to " + to);
        }
    }
}
