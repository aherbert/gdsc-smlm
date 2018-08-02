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
package uk.ac.sussex.gdsc.smlm.results;

import java.awt.Rectangle;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;

import com.google.protobuf.InvalidProtocolBufferException;
import com.google.protobuf.util.JsonFormat;
import com.google.protobuf.util.JsonFormat.Parser;

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.FitMode;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.FluorophoreType;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.IntensityUnits;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.LocationUnits;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.ROI;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.Spot;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.SpotList;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.ThetaUnits;

/**
 * Reads the fit results from file using the Tagged Spot File (TSF) format.
 * <p>
 * Has only limited support for TSF in that only 1 channel, position, slice and fluorophore type can be read into a
 * dataset.
 *
 * @author Alex Herbert
 */
public class TSFPeakResultsReader
{
    private String filename = null;
    private SpotList spotList = null;
    private boolean readHeader = false;
    private boolean isGDSC = false;
    private boolean isMulti = false;

    private int channel = 1;
    private int slice = 0;
    private int position = 0;
    private int fluorophoreType = 1;

    /**
     * Instantiates a new TSF peak results reader.
     *
     * @param filename
     *            the filename
     */
    public TSFPeakResultsReader(String filename)
    {
        this.filename = filename;
    }

    /**
     * Read the TSF header.
     *
     * @return the spot list header
     */
    public SpotList readHeader()
    {
        if (readHeader)
            return spotList;
        readHeader = true;

        try (FileInputStream fi = new FileInputStream(filename))
        {
            try (DataInputStream di = new DataInputStream(fi))
            {
                // the file has an initial 0, then the offset (as long)
                // to the position of spotList
                final int magic = di.readInt();
                if (magic != 0)
                    throw new RuntimeException("Magic number is not 0 (required for a TSF file)");
                if (fi.available() == 0)
                    throw new RuntimeException("Cannot read offset");
                final long offset = di.readLong();
                if (offset == 0)
                    throw new RuntimeException("Offset is 0, cannot find header data in this file");
                fi.skip(offset);
                spotList = SpotList.parseDelimitedFrom(fi);

                // We can do special processing for a TSF file we created
                isGDSC = (spotList.getApplicationId() == TSFPeakResultsWriter.APPLICATION_ID);
                isMulti = isMulti(spotList);
            }
        }
        catch (final Exception e)
        {
            System.err.println("Failed to read SpotList message");
            e.printStackTrace();
        }

        return spotList;
    }

    /**
     * Checks if the header data in the SpotList contains multiple channels, slices, positions, or fluorophores. These
     * are categories that can be use to filter the data when reading.
     *
     * @param spotList
     *            the spot list
     * @return true, if is multi
     */
    public static boolean isMulti(SpotList spotList)
    {
        if (spotList.getNrChannels() > 1)
            return true;
        if (spotList.getNrSlices() > 1)
            return true;
        if (spotList.getNrPos() > 1)
            return true;
        if (spotList.getFluorophoreTypesCount() > 1)
            return true;
        return false;
    }

    /**
     * Checks if is a binary TSF file by attempting to read the SpotList header.
     *
     * @param filename
     *            the filename
     * @return true, if is a TSF file
     */
    public static boolean isTSF(String filename)
    {
        try (FileInputStream fi = new FileInputStream(filename))
        {
            try (DataInputStream di = new DataInputStream(fi))
            {
                // the file has an initial 0, then the offset (as long)
                // to the position of spotList
                final int magic = di.readInt();
                if (magic != 0)
                    // Magic number should be zero
                    return false;
                if (fi.available() == 0)
                    // No more contents
                    return false;
                final long offset = di.readLong();
                if (offset == 0)
                    // No offset record
                    return false;
                fi.skip(offset);
                final SpotList spotList = SpotList.parseDelimitedFrom(fi);
                if (spotList != null)
                    return true;
            }
        }
        catch (final Exception e)
        {
            // Fail
        }

        return false;
    }

    /**
     * Read the results from the TSF file into memory
     *
     * @return The results set (or null if an error occurred)
     */
    public MemoryPeakResults read()
    {
        readHeader();
        if (spotList == null)
            return null;

        final MemoryPeakResults results = createResults();

        // Used in the exception handler to check the correct number of spots were read
        long expectedSpots = -1;

        // Read the messages that contain the spot data
        try (FileInputStream fi = new FileInputStream(filename))
        {
            fi.skip(12); // size of int + size of long

            final LocationUnits locationUnits = spotList.getLocationUnits();
            boolean locationUnitsWarning = false;

            final IntensityUnits intensityUnits = spotList.getIntensityUnits();
            boolean intensityUnitsWarning = false;

            FitMode fitMode = FitMode.ONEAXIS;
            if (spotList.hasFitMode())
                fitMode = spotList.getFitMode();

            final boolean filterPosition = position > 0;
            final boolean filterSlice = slice > 0;

            // Set up to read two-axis and theta data
            final PSF psf = PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
            final int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
            final int isx = indices[0];
            final int isy = indices[1];
            final int ia = PSFHelper.getGaussian2DAngleIndex(psf);
            int nParams = PeakResult.STANDARD_PARAMETERS;
            switch (fitMode)
            {
                case TWOAXISANDTHETA:
                    nParams++;
                case TWOAXIS:
                    nParams++;
                case ONEAXIS:
                default:
                    nParams++;
                    break;
            }

            expectedSpots = getExpectedSpots();
            long totalSpots = 0;
            while (fi.available() > 0 && (totalSpots < expectedSpots || expectedSpots == 0))
            {
                totalSpots++;

                final Spot spot = Spot.parseDelimitedFrom(fi);

                // Only read the specified channel, position, slice and fluorophore type
                if (spot.getChannel() != channel)
                    continue;
                if (filterPosition && spot.hasPos() && spot.getPos() != position)
                    continue;
                if (filterSlice && spot.hasSlice() && spot.getSlice() != slice)
                    continue;
                if (spot.hasFluorophoreType() && spot.getFluorophoreType() != fluorophoreType)
                    continue;

                if (spot.hasLocationUnits() && !locationUnitsWarning && spot.getLocationUnits() != locationUnits)
                {
                    System.err.println("Spot message has different location units, the units will be ignored: " +
                            spot.getLocationUnits());
                    locationUnitsWarning = true;
                }
                if (spot.hasIntensityUnits() && !intensityUnitsWarning && spot.getIntensityUnits() != intensityUnits)
                {
                    System.err.println("Spot message has different intensity units, the units will be ignored: " +
                            spot.getIntensityUnits());
                    intensityUnitsWarning = true;
                }

                // Required fields
                final int frame = spot.getFrame();

                final float[] params = new float[nParams];
                params[PeakResult.X] = spot.getX();
                params[PeakResult.Y] = spot.getY();
                params[PeakResult.INTENSITY] = spot.getIntensity();
                if (spot.hasBackground())
                    params[PeakResult.BACKGROUND] = spot.getBackground();
                if (spot.hasZ())
                    params[PeakResult.Z] = spot.getZ();

                // Support different Gaussian shapes
                if (fitMode == FitMode.ONEAXIS)
                    params[isx] = (float) (spot.getWidth() / Gaussian2DFunction.SD_TO_FWHM_FACTOR);
                else
                {
                    if (!spot.hasA())
                        params[isx] = params[isy] = (float) (spot.getWidth() / Gaussian2DFunction.SD_TO_FWHM_FACTOR);
                    else
                    {
                        final double a = Math.sqrt(spot.getA());
                        final double sd = spot.getWidth() / Gaussian2DFunction.SD_TO_FWHM_FACTOR;
                        params[isx] = (float) (sd * a);
                        params[isy] = (float) (sd / a);
                    }

                    if (fitMode == FitMode.TWOAXISANDTHETA && spot.hasTheta())
                        params[ia] = spot.getTheta();
                }

                // We can use the original position in pixels used for fitting
                final int origX = (spot.hasXPosition()) ? spot.getXPosition() : 0;
                final int origY = (spot.hasYPosition()) ? spot.getYPosition() : 0;

                // Q. Should we use the required field 'molecule'?

                float origValue = 0;
                double error = 0;
                float noise = 0;
                float meanIntensity = 0;
                float[] paramsStdDev = null;
                int id = 0;
                int endFrame = frame;

                if (isGDSC)
                {
                    error = spot.getError();
                    noise = spot.getNoise();
                    meanIntensity = spot.getMeanIntensity();
                    id = spot.getCluster();
                    origValue = spot.getOriginalValue();
                    endFrame = spot.getEndFrame();
                    if (spot.getParamStdDevsCount() != 0)
                    {
                        paramsStdDev = new float[spot.getParamStdDevsCount()];
                        for (int i = 0; i < paramsStdDev.length; i++)
                            paramsStdDev[i] = spot.getParamStdDevs(i);
                    }
                }
                // Use the standard cluster field for the ID
                else if (spot.hasCluster())
                    id = spot.getCluster();

                // Allow storing any of the optional attributes
                final AttributePeakResult peakResult = new AttributePeakResult(frame, origX, origY, origValue, error,
                        noise, meanIntensity, params, paramsStdDev);

                peakResult.setEndFrame(endFrame);
                peakResult.setId(id);
                if (spot.hasXPrecision() || spot.hasYPrecision())
                {
                    // Use the average. Note this is not the Euclidean distance since we divide by n
                    double p = 0;
                    int n = 0;
                    if (spot.hasXPrecision())
                    {
                        p = spot.getXPrecision() * spot.getXPrecision();
                        n++;
                    }
                    if (spot.hasYPrecision())
                    {
                        p += spot.getYPrecision() * spot.getYPrecision();
                        n++;
                    }
                    peakResult.setPrecision(Math.sqrt(p / n));
                }
                results.add(peakResult);
            }
        }
        catch (final IOException e)
        {
            System.err.println("Failed to read TSF file: " + filename);
            e.printStackTrace();

            if (expectedSpots == -1)
                // No attempt to read the spots was made.
                // The exception was created during set-up.
                return null;

            // If expectedSpots==0 then the number of spots was unknown and the file was
            // read until the EOF.
            // Only fail if there is a number of expected spots.
            if (expectedSpots != 0)
            {
                System.err.println("Unexpected error in reading Spot messages, no results will be returned");
                return null;
            }
        }

        // Do log a warning if the expected spots does not match the size.
        // The spots may be from multiple channels. etc.

        return results;
    }

    private long getExpectedSpots()
    {
        if (spotList.hasNrSpots())
            return spotList.getNrSpots();
        return 0;
    }

    private static uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType[] cameraTypeMap;
    private static uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit[] thetaUnitsMap;
    private static uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit[] locationUnitsMap;
    private static uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit[] intensityUnitsMap;
    static
    {
        // These should have 1:1 mapping. We can extends the TSF proto if necessary.
        cameraTypeMap = new uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType[CameraType.values().length];
        cameraTypeMap[CameraType.CCD.ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType.CCD;
        cameraTypeMap[CameraType.EMCCD
                .ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType.EMCCD;
        cameraTypeMap[CameraType.SCMOS
                .ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType.SCMOS;
        thetaUnitsMap = new uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit[ThetaUnits.values().length];
        thetaUnitsMap[ThetaUnits.RADIANS.ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit.RADIAN;
        thetaUnitsMap[ThetaUnits.DEGREES.ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit.DEGREE;
        locationUnitsMap = new uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit[LocationUnits
                .values().length];
        locationUnitsMap[LocationUnits.NM.ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit.NM;
        locationUnitsMap[LocationUnits.UM.ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit.UM;
        locationUnitsMap[LocationUnits.PIXELS
                .ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit.PIXEL;
        intensityUnitsMap = new uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit[IntensityUnits
                .values().length];
        intensityUnitsMap[IntensityUnits.COUNTS
                .ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit.COUNT;
        intensityUnitsMap[IntensityUnits.PHOTONS
                .ordinal()] = uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit.PHOTON;
    }

    private MemoryPeakResults createResults()
    {
        // Limit the capacity since we may not need all the spots
        int capacity = 1000;
        if (spotList.hasNrSpots())
            capacity = (int) Math.min(100000, spotList.getNrSpots());
        final MemoryPeakResults results = new MemoryPeakResults(capacity);

        // Create the type of Gaussian PSF
        if (spotList.hasFitMode())
            switch (spotList.getFitMode())
            {
                case ONEAXIS:
                    results.setPSF(PSFHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D));
                    break;
                case TWOAXIS:
                    results.setPSF(PSFHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D));
                    break;
                case TWOAXISANDTHETA:
                    results.setPSF(PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D));
                    break;
                default:
                    break;
            }

        // Generic reconstruction
        String name;
        if (spotList.hasName())
            name = spotList.getName();
        else
            name = new File(filename).getName();
        // Append these if not using the defaults
        if (channel != 1 || slice != 0 || position != 0 || fluorophoreType != 1)
            name = String.format("%s c=%d, s=%d, p=%d, ft=%d", name, channel, slice, position, fluorophoreType);
        results.setName(name);

        //		if (spotList.hasNrPixelsX() && spotList.hasNrPixelsY())
        //		{
        //			// Do not do this. The size of the camera may not map to the data bounds due
        //			// to the support for position offsets.
        //			results.setBounds(new Rectangle(0, 0, spotList.getNrPixelsX(), spotList.getNrPixelsY()));
        //		}

        final CalibrationWriter cal = new CalibrationWriter();

        if (spotList.hasPixelSize())
            cal.setNmPerPixel(spotList.getPixelSize());
        if (spotList.getEcfCount() >= channel)
        {
            // ECF is per channel
            final double ecf = spotList.getEcf(channel - 1);
            // QE is per fluorophore type
            final double qe = (spotList.getQeCount() >= fluorophoreType) ? spotList.getQe(fluorophoreType - 1) : 1;
            // e-/photon / e-/count => count/photon
            cal.setCountPerPhoton(qe / ecf);
            cal.setQuantumEfficiency(qe);
        }

        if (isGDSC)
        {
            // Special processing for results we created to allow
            // perfect dataset reconstruction

            if (spotList.hasSource())
                // Deserialise
                results.setSource(ImageSource.fromXML(spotList.getSource()));

            if (spotList.hasRoi())
            {
                final ROI roi = spotList.getRoi();
                if (roi.hasX() && roi.hasY() && roi.hasXWidth() && roi.hasYWidth())
                    results.setBounds(new Rectangle(roi.getX(), roi.getY(), roi.getXWidth(), roi.getYWidth()));
            }

            if (spotList.hasGain())
                cal.setCountPerPhoton(spotList.getGain());
            if (spotList.hasExposureTime())
                cal.setExposureTime(spotList.getExposureTime());
            if (spotList.hasReadNoise())
                cal.setReadNoise(spotList.getReadNoise());
            if (spotList.hasBias())
                cal.setBias(spotList.getBias());
            if (spotList.hasCameraType())
                cal.setCameraType(cameraTypeMap[spotList.getCameraType().ordinal()]);
            else
                cal.setCameraType(null);

            if (spotList.hasConfiguration())
                results.setConfiguration(spotList.getConfiguration());
            // Allow restoring the GDSC PSF exactly
            if (spotList.hasPSF())
                try
                {
                    final Parser parser = JsonFormat.parser();
                    final PSF.Builder psfBuilder = PSF.newBuilder();
                    parser.merge(spotList.getPSF(), psfBuilder);
                    results.setPSF(psfBuilder.build());
                }
                catch (final InvalidProtocolBufferException e)
                {
                    // This should be OK
                    System.err.println("Unable to deserialise the PSF settings");
                }
        }

        if (spotList.hasLocationUnits())
        {
            cal.setDistanceUnit(locationUnitsMap[spotList.getLocationUnits().ordinal()]);
            if (!spotList.hasPixelSize() && spotList.getLocationUnits() != LocationUnits.PIXELS)
                System.err.println(
                        "TSF location units are not pixels and no pixel size calibration is available. The dataset will be constructed in the native units: " +
                                spotList.getLocationUnits());
        }
        else
            cal.setDistanceUnit(null);

        if (spotList.hasIntensityUnits())
        {
            cal.setIntensityUnit(intensityUnitsMap[spotList.getIntensityUnits().ordinal()]);
            if (!spotList.hasGain() && spotList.getIntensityUnits() != IntensityUnits.COUNTS)
                System.err.println(
                        "TSF intensity units are not counts and no gain calibration is available. The dataset will be constructed in the native units: " +
                                spotList.getIntensityUnits());
        }
        else
            cal.setIntensityUnit(null);

        if (spotList.hasThetaUnits())
            cal.setAngleUnit(thetaUnitsMap[spotList.getThetaUnits().ordinal()]);
        else
            cal.setAngleUnit(null);

        results.setCalibration(cal.getCalibration());

        return results;
    }

    /**
     * Gets the channel to read.
     *
     * @return the channel
     */
    public int getChannel()
    {
        return channel;
    }

    /**
     * Sets the channel to read.
     *
     * @param channel
     *            the new channel
     */
    public void setChannel(int channel)
    {
        this.channel = get1Based(channel);
    }

    private static int get1Based(int value)
    {
        return (value < 1) ? 1 : value;
    }

    /**
     * Gets the slice to read.
     *
     * @return the slice
     */
    public int getSlice()
    {
        return slice;
    }

    /**
     * Sets the slice to read.
     *
     * @param slice
     *            the new slice
     */
    public void setSlice(int slice)
    {
        this.slice = get1Based(slice);
    }

    /**
     * Gets the position to read.
     *
     * @return the position
     */
    public int getPosition()
    {
        return position;
    }

    /**
     * Sets the position to read.
     *
     * @param position
     *            the new position
     */
    public void setPosition(int position)
    {
        this.position = get1Based(position);
    }

    /**
     * Gets the fluorophore type to read.
     *
     * @return the fluorophore type
     */
    public int getFluorophoreType()
    {
        return fluorophoreType;
    }

    /**
     * Sets the fluorophore type to read.
     *
     * @param fluorophoreType
     *            the new fluorophore type
     */
    public void setFluorophoreType(int fluorophoreType)
    {
        this.fluorophoreType = get1Based(fluorophoreType);
    }

    /**
     * Checks if is a GDSC TSF file.
     *
     * @return true, if is GDSC TSF file
     */
    public boolean isGDSC()
    {
        readHeader();
        return isGDSC;
    }

    /**
     * Checks if the header data in the SpotList contains multiple channels, slices, positions, or fluorophores. These
     * are categories that can be use to filter the data when reading.
     *
     * @return true, if the data contains multiple categories of localisations
     */
    public boolean isMulti()
    {
        readHeader();
        return isMulti;
    }

    /**
     * Gets the options for reading the results.
     *
     * @return the options
     */
    public ResultOption[] getOptions()
    {
        if (!isMulti())
            return null;

        final ResultOption[] options = new ResultOption[4];
        int count = 0;

        if (spotList.getNrChannels() > 1)
            options[count++] = createOption(1, "Channel", spotList.getNrChannels(), 1, false);
        if (spotList.getNrSlices() > 1)
            options[count++] = createOption(2, "Slice", spotList.getNrSlices(), 0, true);
        if (spotList.getNrPos() > 1)
            options[count++] = createOption(3, "Position", spotList.getNrPos(), 0, true);
        if (spotList.getFluorophoreTypesCount() > 1)
        {
            // Build a string for the allowed value to provide space for the description
            final String[] values = new String[spotList.getFluorophoreTypesCount()];
            for (int i = 0; i < spotList.getFluorophoreTypesCount(); i++)
            {
                final FluorophoreType type = spotList.getFluorophoreTypes(i);
                String value = Integer.toString(type.getId());
                if (type.hasDescription())
                    value += ":" + type.getDescription();
                if (type.hasIsFiducial())
                    value += ":fiducial=" + type.getIsFiducial();
                values[i] = value;
            }
            options[count++] = new ResultOption(4, "Fluorophore type", values[0], values);
        }

        return Arrays.copyOf(options, count);
    }

    private static ResultOption createOption(int id, String name, int total, int value, boolean allowZero)
    {
        final Integer[] values = new Integer[total + ((allowZero) ? 1 : 0)];
        for (int v = (allowZero) ? 0 : 1, i = 0; v <= total; v++)
            values[i++] = v;
        return new ResultOption(id, name, new Integer(value), values);
    }

    /**
     * Sets the options for reading the results.
     *
     * @param options
     *            the new options for reading the results
     */
    public void setOptions(ResultOption[] options)
    {
        if (options == null)
            return;
        for (final ResultOption option : options)
            switch (option.id)
            {
                case 1:
                    setChannel((Integer) option.getValue());
                    break;
                case 2:
                    setSlice((Integer) option.getValue());
                    break;
                case 3:
                    setPosition((Integer) option.getValue());
                    break;
                case 4:
                    String value = (String) option.getValue();
                    // Remove the appended description
                    final int index = value.indexOf(':');
                    if (index != -1)
                        value = value.substring(0, index);
                    setFluorophoreType(Integer.parseInt(value));
                    break;
            }
    }
}
