package uk.ac.sussex.gdsc.smlm.results;

import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.internal.ArrayComparisonFailure;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;

import uk.ac.sussex.gdsc.core.utils.NotImplementedException;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;

@SuppressWarnings({ "javadoc" })
public class PeakResultsReaderTest
{
    private static Logger logger;

    @BeforeAll
    public static void beforeAll()
    {
        logger = Logger.getLogger(PeakResultsReaderTest.class.getName());
    }

    @AfterAll
    public static void afterAll()
    {
        logger = null;
    }

    static final boolean[] onOff = new boolean[] { true, false };

    // TODO - Add tests to compare writing to a IJTablePeakResults, saving the TextPanel contents to file and then reading.

    // -=-=-=-=-

    @SeededTest
    public void writeTextMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, false, false);
    }

    @SeededTest
    public void writeSequentialTextMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, true, ResultsFileFormat.TEXT, false, false, false, false, false);
    }

    @SeededTest
    public void writeTextWithDeviationsMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, true, false, false, false, false);
    }

    @SeededTest
    public void writeTextWithEndFrameMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, true, false, false, false);
    }

    @SeededTest
    public void writeTextWithIdMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, true, false, false);
    }

    @SeededTest
    public void writeTextWithPrecisionMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, true, false);
    }

    @SeededTest
    public void writeTextWithCombinationsMatchesRead(RandomSeed seed)
    {
        writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.TEXT, false);
    }

    // -=-=-=-=-

    @SeededTest
    public void writeBinaryMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, false, false);
    }

    @SeededTest
    public void writeSequentialBinaryMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, true, ResultsFileFormat.BINARY, false, false, false, false, false);
    }

    @SeededTest
    public void writeBinaryWithDeviationsMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, true, false, false, false, false);
    }

    @SeededTest
    public void writeBinaryWithEndFrameMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, true, false, false, false);
    }

    @SeededTest
    public void writeBinaryWithIdMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, true, false, false);
    }

    @SeededTest
    public void writeBinaryWithPrecisionMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, true, false);
    }

    @SeededTest
    public void writeBinaryWithCombinationsMatchesRead(RandomSeed seed)
    {
        writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.BINARY, false);
    }

    // -=-=-=-=-

    @SeededTest
    public void writeTextWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, false, true);
    }

    @SeededTest
    public void writeSequentialTextWithSortMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, true, ResultsFileFormat.TEXT, false, false, false, false, true);
    }

    @SeededTest
    public void writeTextWithDeviationsWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, true, false, false, false, true);
    }

    @SeededTest
    public void writeTextWithEndFrameWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, true, false, false, true);
    }

    @SeededTest
    public void writeTextWithIdWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, true, false, true);
    }

    @SeededTest
    public void writeTextWithPrecisionWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, true, true);
    }

    @SeededTest
    public void writeTextWithCombinationsWithSortMatchesRead(RandomSeed seed)
    {
        writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.TEXT, true);
    }

    // -=-=-=-=-

    @SeededTest
    public void writeBinaryWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, false, true);
    }

    @SeededTest
    public void writeSequentialBinaryWithSortMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, true, ResultsFileFormat.BINARY, false, false, false, false, true);
    }

    @SeededTest
    public void writeBinaryWithDeviationsWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, true, false, false, false, true);
    }

    @SeededTest
    public void writeBinaryWithEndFrameWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, true, false, false, true);
    }

    @SeededTest
    public void writeBinaryWithIdWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, true, false, true);
    }

    @SeededTest
    public void writeBinaryWithPrecisionWithSortMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, true, true);
    }

    @SeededTest
    public void writeBinaryWithCombinationsWithSortMatchesRead(RandomSeed seed)
    {
        writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.BINARY, true);
    }

    // -=-=-=-=-

    // Note: For MALK we cannot do all the tests as the format only contains X,Y,T,I

    @SeededTest
    public void writeMALKMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, false, ResultsFileFormat.MALK, false, false, false, false, false);
    }

    @SeededTest
    public void writeSequentialMALKMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, true, ResultsFileFormat.MALK, false, false, false, false, false);
    }

    @SeededTest
    public void writeMALKWithSortMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, false, ResultsFileFormat.MALK, false, false, false, true, false);
    }

    @SeededTest
    public void writeSequentialMALKWithSortMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, true, ResultsFileFormat.MALK, false, false, false, true, false);
    }

    // -=-=-=-=-

    // Note: For TSF we cannot specify as binary because the widths are converted into a
    // different format and then back again.

    @SeededTest
    public void writeTSFMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, false, false, false, false);
    }

    @SeededTest
    public void writeSequentialTSFMatchesRead(RandomSeed seed)
    {
        writeMatchesRead(seed, true, ResultsFileFormat.TSF, false, false, false, false, false);
    }

    @SeededTest
    public void writeTSFWithDeviationsMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TSF, true, false, false, false, false);
    }

    @SeededTest
    public void writeTSFWithEndFrameMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, true, false, false, false);
    }

    @SeededTest
    public void writeTSFWithIdMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, false, true, false, false);
    }

    @SeededTest
    public void writeTSFWithPrecisionMatchesRead(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, false, false, true, false);
    }

    @SeededTest
    public void writeTSFWithCombinationsMatchesRead(RandomSeed seed)
    {
        writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.TSF, false);
    }

    // -=-=-=-=-

    @SeededTest
    public void readWithScannerMatchesNonScanner(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, false, false, false, false);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithDeviations(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, true, false, false, false, false);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithEndFrame(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, true, false, false, false);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithId(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, false, true, false, false);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithPrecision(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, false, false, true, false);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithCombinations(RandomSeed seed)
    {
        ExtraAssumptions.assumeLowComplexity(); // Scanner is not a default so do not always test
        readWithScannerMatchesNonScannerWithCombinations(seed, false);
    }

    // -=-=-=-=-

    @SeededTest
    public void readWithScannerMatchesNonScannerWithSort(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, false, false, false, true);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithDeviationsWithSort(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, true, false, false, false, true);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithEndFrameWithSort(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, true, false, false, true);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithIdWithSort(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, false, true, false, true);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithPrecisionWithSort(RandomSeed seed)
    {
        ExtraAssumptions.assumeMediumComplexity(); // The combinations test overlaps this
        readWithScannerMatchesNonScanner(seed, false, false, false, true, true);
    }

    @SeededTest
    public void readWithScannerMatchesNonScannerWithCombinationsWithSort(RandomSeed seed)
    {
        ExtraAssumptions.assumeLowComplexity(); // Scanner is not a default so do not always test
        readWithScannerMatchesNonScannerWithCombinations(seed, true);
    }

    // -=-=-=-=-

    @SeededTest
    public void readTextWithNonScannerIsFasterThanScanner(RandomSeed seed)
    {
        readWith2IsFasterThan1(seed, false, false, false, false, ResultsFileFormat.TEXT, true, ResultsFileFormat.TEXT,
                false, 1);
    }

    @SeededTest
    public void readTextWithNonScannerIsFasterThanScannerWithDeviationsWithEndFrameWithIdWithPrecision(RandomSeed seed)
    {
        readWith2IsFasterThan1(seed, true, true, true, true, ResultsFileFormat.TEXT, true, ResultsFileFormat.TEXT,
                false, 1);
    }

    @SeededTest
    public void readWithMALKIsFasterThanText(RandomSeed seed)
    {
        readWith2IsFasterThan1(seed, false, false, false, false, ResultsFileFormat.TEXT, false, ResultsFileFormat.MALK,
                false, 2);
    }

    @SeededTest
    public void readWithBinaryIsFasterThanText(RandomSeed seed)
    {
        readWith2IsFasterThan1(seed, false, false, false, false, ResultsFileFormat.TEXT, false,
                ResultsFileFormat.BINARY, false, 2);
    }

    @SeededTest
    public void readWithBinaryIsFasterThanTSF(RandomSeed seed)
    {
        readWith2IsFasterThan1(seed, false, false, false, false, ResultsFileFormat.TSF, false, ResultsFileFormat.BINARY,
                false, 20);
    }

    @SeededTest
    public void canConvertMalkToNMAndPhotons(RandomSeed seed)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults out = createResults(rg, 200, false, false, false, false);

        // Output in pixel and count
        final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
        cal.setDistanceUnit(DistanceUnit.PIXEL);
        cal.setIntensityUnit(IntensityUnit.COUNT);
        out.setCalibration(cal.getCalibration());
        out.setPSF(PSFHelper.create(PSFType.CUSTOM));

        final String filename = createFile();

        writeFile(false, ResultsFileFormat.MALK, false, false, false, false, false, out, filename);

        final MemoryPeakResults in = readFile(filename, false);

        // Change to nm and photon for the validation
        out.convertToUnits(DistanceUnit.NM, IntensityUnit.PHOTON, null);

        checkEqual(ResultsFileFormat.MALK, false, false, false, false, false, out, in);
    }

    @SeededTest
    public void writeTextWithComputedPrecisionMatchesRead(RandomSeed seed)
    {
        // Create without precision
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults results = createResults(rg, 200, false, false, false, false);
        // Ensure units are OK for computing precision
        final CalibrationWriter cw = results.getCalibrationWriter();
        cw.setIntensityUnit(IntensityUnit.PHOTON);
        cw.setDistanceUnit(DistanceUnit.PIXEL);
        results.setCalibration(cw.getCalibration());

        final String filename = createFile();

        final TextFilePeakResults out = new TextFilePeakResults(filename, false, false, false, true);
        // Compute precision
        out.setComputePrecision(true);
        out.copySettings(results);
        out.begin();
        out.addAll(Arrays.asList(results.toArray()));
        out.end();

        final MemoryPeakResults in = readFile(filename, false);

        checkEqual(ResultsFileFormat.TEXT, false, false, false, false, false, results, in);
    }

    @SeededTest
    public void canReadTextIntoPreferredUnits(RandomSeed seed)
    {
        canReadIntoPreferredUnits(seed, ResultsFileFormat.TEXT);
    }

    @SeededTest
    public void canReadBinaryIntoPreferredUnits(RandomSeed seed)
    {
        canReadIntoPreferredUnits(seed, ResultsFileFormat.BINARY);
    }

    @SeededTest
    public void canReadTSFIntoPreferredUnits(RandomSeed seed)
    {
        canReadIntoPreferredUnits(seed, ResultsFileFormat.TSF);
    }

    private static void canReadIntoPreferredUnits(RandomSeed seed, ResultsFileFormat fileFormat)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults out = createResults(rg, 200, false, false, false, false);

        // Output in nm and count
        final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
        cal.setDistanceUnit(DistanceUnit.NM);
        cal.setIntensityUnit(IntensityUnit.COUNT);
        if (fileFormat == ResultsFileFormat.TSF)
            // For now just support using the native float TSF value
            cal.setNmPerPixel((float) cal.getNmPerPixel());
        out.setCalibration(cal.getCalibration());

        final String filename = createFile();

        writeFile(false, fileFormat, false, false, false, false, false, out, filename);

        final MemoryPeakResults in = readFile(filename, false, false);

        // Change to preferred units
        out.convertToUnits(MemoryPeakResults.PREFERRED_DISTANCE_UNIT, MemoryPeakResults.PREFERRED_INTENSITY_UNIT,
                MemoryPeakResults.PREFERRED_ANGLE_UNIT);

        checkEqual(fileFormat, false, false, false, false, false, out, in);
    }

    @SeededTest
    public void canReadTextAndSimplifyGaussian2DPSF(RandomSeed seed)
    {
        canReadAndSimplifyGaussian2DPSF(seed, ResultsFileFormat.TEXT);
    }

    @SeededTest
    public void canReadBinaryAndSimplifyGaussian2DPSF(RandomSeed seed)
    {
        canReadAndSimplifyGaussian2DPSF(seed, ResultsFileFormat.BINARY);
    }

    @SeededTest
    public void canReadTSFAndSimplifyGaussian2DPSF(RandomSeed seed)
    {
        canReadAndSimplifyGaussian2DPSF(seed, ResultsFileFormat.TSF);
    }

    private static void canReadAndSimplifyGaussian2DPSF(RandomSeed seed, ResultsFileFormat fileFormat)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults out = createResults(rg, 1, false, false, false, false);

        final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
        cal.setDistanceUnit(MemoryPeakResults.PREFERRED_DISTANCE_UNIT);
        cal.setIntensityUnit(MemoryPeakResults.PREFERRED_INTENSITY_UNIT);
        cal.setAngleUnit(MemoryPeakResults.PREFERRED_ANGLE_UNIT);
        if (fileFormat == ResultsFileFormat.TSF)
            // For now just support using the native float TSF value
            cal.setNmPerPixel((float) cal.getNmPerPixel());
        out.setCalibration(cal.getCalibration());

        // Remove angle
        final int ia = PSFHelper.getGaussian2DAngleIndex(out.getPSF());
        out.forEach(new PeakResultProcedure()
        {
            @Override
            public void execute(PeakResult peakResult)
            {
                peakResult.getParameters()[ia] = 0;
            }
        });

        final String filename = createFile();

        writeFile(false, fileFormat, false, false, false, false, false, out, filename);

        MemoryPeakResults in = readFile(filename, false, false);

        // Change to two-axis PSF
        out.setPSF(PSFHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D));
        final int twoAxisLength = PSFHelper.getParameterCount(out.getPSF()) + PeakResult.STANDARD_PARAMETERS;
        out.forEach(new PeakResultProcedure()
        {
            @Override
            public void execute(PeakResult peakResult)
            {
                peakResult.resizeParameters(twoAxisLength);
            }
        });

        checkEqual(fileFormat, false, false, false, false, false, out, in);

        // Remove sy
        final int[] indices = PSFHelper.getGaussian2DWxWyIndices(out.getPSF());
        final int isx = indices[0];
        final int isy = indices[1];
        out.forEach(new PeakResultProcedure()
        {
            @Override
            public void execute(PeakResult peakResult)
            {
                final float[] p = peakResult.getParameters();
                p[isy] = p[isx];
            }
        });

        writeFile(false, fileFormat, false, false, false, false, false, out, filename);

        in = readFile(filename, false, false);

        // Change to one-axis PSF
        out.setPSF(PSFHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D));
        final int oneAxisLength = PSFHelper.getParameterCount(out.getPSF()) + PeakResult.STANDARD_PARAMETERS;
        out.forEach(new PeakResultProcedure()
        {
            @Override
            public void execute(PeakResult peakResult)
            {
                peakResult.resizeParameters(oneAxisLength);
            }
        });

        checkEqual(fileFormat, false, false, false, false, false, out, in);
    }

    private static void readWith2IsFasterThan1(RandomSeed seed, boolean showDeviations, boolean showEndFrame,
            boolean showId, boolean showPrecision, ResultsFileFormat f1, boolean useScanner1, ResultsFileFormat f2,
            boolean useScanner2, int loops)
    {
        ExtraAssumptions.assume(logger, Level.INFO);
        ExtraAssumptions.assume(TestComplexity.HIGH);

        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults out = createResults(rg, 20000, showDeviations, showEndFrame, showId, showPrecision);
        final String filename = createFile();

        writeFile(false, f1, showDeviations, showEndFrame, showId, showPrecision, false, out, filename);
        final long time1 = getReadTime(filename, useScanner1, loops);

        writeFile(false, f2, showDeviations, showEndFrame, showId, showPrecision, false, out, filename);
        final long time2 = getReadTime(filename, useScanner2, loops);

        if (useScanner1 != useScanner2)
            logger.info(TestLog.getSupplier("%s (scan=%b) is %.2fx faster than %s (scan=%b)", f2, useScanner2,
                    (double) time1 / time2, f1, useScanner1));
        else
            logger.info(TestLog.getSupplier("%s is %.2fx faster than %s", f2, (double) time1 / time2, f1));
        ExtraAssertions.assertTrue(time2 < time1, "%s (%d) is not faster than %s (%d)", f2, time2, f1, time1);
    }

    // -=-=-=-=-

    private static void writeMatchesRead(RandomSeed seed, boolean sequential, ResultsFileFormat fileFormat,
            boolean showDeviations, boolean showEndFrame, boolean showId, boolean showPrecision, boolean sort)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults out = createResults(rg, 200, showDeviations, showEndFrame, showId, showPrecision);
        if (fileFormat == ResultsFileFormat.MALK)
        {
            final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
            cal.setDistanceUnit(DistanceUnit.NM);
            cal.setIntensityUnit(IntensityUnit.PHOTON);
            out.setCalibration(cal.getCalibration());
            out.setPSF(PSFHelper.create(PSFType.CUSTOM));
        }
        if (fileFormat == ResultsFileFormat.TSF)
        {
            final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
            // For now just support using the native float TSF datatype
            cal.setNmPerPixel((float) cal.getNmPerPixel());
            out.setCalibration(cal.getCalibration());
            // TSF converts the width parameters so make sure they are not zero
            out.forEach(new PeakResultProcedure()
            {

                @Override
                public void execute(PeakResult peakResult)
                {
                    check(peakResult.getParameters());
                    if (showDeviations)
                        check(peakResult.getParameterDeviations());
                }

                private void check(float[] parameters)
                {
                    for (int i = 0; i < parameters.length; i++)
                        if (parameters[i] == 0)
                            parameters[i] = 0.1f;
                }
            });
        }

        //		System.out.println(out.getCalibration());
        //		System.out.println(out.getPSF().toString());
        //
        //		System.out.println(TextFormat.shortDebugString(out.getCalibration()));
        //
        //		try
        //		{
        //			Printer printer = JsonFormat.printer()
        //					.omittingInsignificantWhitespace()
        //					//.includingDefaultValueFields()
        //					;
        //			System.out.println(printer.print(out.getCalibration()));
        //			System.out.println(printer.print(out.getPSF()));
        //		}
        //		catch (InvalidProtocolBufferException e)
        //		{
        //			// This shouldn't happen so throw it
        //		}

        final String filename = createFile();

        writeFile(sequential, fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, out, filename);

        final MemoryPeakResults in = readFile(filename, false);

        checkEqual(fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, out, in);
    }

    private static void writeWithCombinationsMatchesRead(RandomSeed seed, boolean sequential,
            ResultsFileFormat fileFormat, boolean sort)
    {
        for (final boolean showDeviations : onOff)
            for (final boolean showEndFrame : onOff)
                for (final boolean showId : onOff)
                    for (final boolean showPrecision : onOff)
                    {
                        if (count(showDeviations, showEndFrame, showId, showPrecision) < 2)
                            continue;
                        writeMatchesRead(seed, sequential, fileFormat, showDeviations, showEndFrame, showId,
                                showPrecision, sort);
                    }
    }

    private static int count(boolean... flags)
    {
        int c = 0;
        for (final boolean flag : flags)
            if (flag)
                c++;
        return c;
    }

    private static void readWithScannerMatchesNonScanner(RandomSeed seed, boolean showDeviations, boolean showEndFrame,
            boolean showId, boolean showPrecision, boolean sort)
    {
        final UniformRandomProvider rg = TestSettings.getRandomGenerator(seed.getSeed());
        final MemoryPeakResults out = createResults(rg, 1000, showDeviations, showEndFrame, showId, showPrecision);
        final String filename = createFile();

        final ResultsFileFormat fileFormat = ResultsFileFormat.TEXT;
        writeFile(false, fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, out, filename);

        final MemoryPeakResults in = readFile(filename, false);
        final MemoryPeakResults in2 = readFile(filename, true);

        checkEqual(fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, in, in2);
    }

    private static void readWithScannerMatchesNonScannerWithCombinations(RandomSeed seed, boolean sort)
    {
        for (final boolean showDeviations : onOff)
            for (final boolean showEndFrame : onOff)
                for (final boolean showId : onOff)
                    for (final boolean showPrecision : onOff)
                    {
                        if (count(showDeviations, showEndFrame, showId, showPrecision) < 2)
                            continue;
                        readWithScannerMatchesNonScanner(seed, showDeviations, showEndFrame, showId, showPrecision,
                                sort);
                    }
    }

    private static void checkEqual(ResultsFileFormat fileFormat, boolean showDeviations, boolean showEndFrame,
            boolean showId, boolean showPrecision, boolean sort, MemoryPeakResults expectedResults,
            MemoryPeakResults actualResults) throws ArrayComparisonFailure
    {
        Assertions.assertNotNull(actualResults, "Input results are null");
        Assertions.assertEquals(expectedResults.size(), actualResults.size(), "Size differ");

        final PeakResult[] expected = expectedResults.toArray();
        final PeakResult[] actual = actualResults.toArray();
        if (sort)
            // Results should be sorted by time
            Arrays.sort(expected, new Comparator<PeakResult>()
            {
                @Override
                public int compare(PeakResult o1, PeakResult o2)
                {
                    return o1.getFrame() - o2.getFrame();
                }
            });

        // TSF requires the bias be subtracted
        //		double bias = expectedResults.getCalibration().getBias();

        for (int i = 0; i < actualResults.size(); i++)
        {
            final PeakResult p1 = expected[i];
            final PeakResult p2 = actual[i];

            ExtraAssertions.assertEquals(p1.getFrame(), p2.getFrame(), "Peak mismatch @ [%d]", i);

            if (fileFormat == ResultsFileFormat.MALK)
            {
                final double delta = 1e-5f;
                ExtraAssertions.assertEqualsRelative(p1.getXPosition(), p2.getXPosition(), delta, "X @ [%d]", i);
                ExtraAssertions.assertEqualsRelative(p1.getYPosition(), p2.getYPosition(), delta, "Y @ [%d]", i);
                ExtraAssertions.assertEqualsRelative(p1.getIntensity(), p2.getIntensity(), delta, "Signal @ " + i);
                continue;
            }

            ExtraAssertions.assertEquals(p1.getOrigX(), p2.getOrigX(), "Orig X mismatch @ [%d]", i);
            ExtraAssertions.assertEquals(p1.getOrigY(), p2.getOrigY(), "Orig Y mismatch @ [%d]", i);
            ExtraAssertions.assertNotNull(p2.getParameters(), "Params is null @ [%d]", i);
            if (showEndFrame)
                ExtraAssertions.assertEquals(p1.getEndFrame(), p2.getEndFrame(), "End frame mismatch @ [%d]", i);
            if (showId)
                ExtraAssertions.assertEquals(p1.getId(), p2.getId(), "ID mismatch @ [%d]", i);
            if (showDeviations)
                ExtraAssertions.assertNotNull(p2.getParameterDeviations(), "Deviations @ [%d]", i);

            // Binary should be exact for float numbers
            if (fileFormat == ResultsFileFormat.BINARY)
            {
                ExtraAssertions.assertEquals(p1.getOrigValue(), p2.getOrigValue(), "Orig value mismatch @ [%d]", i);
                ExtraAssertions.assertEquals(p1.getError(), p2.getError(), "Error mismatch @ [%d]", i);
                ExtraAssertions.assertEquals(p1.getNoise(), p2.getNoise(), "Noise mismatch @ [%d]", i);
                ExtraAssertions.assertEquals(p1.getMeanIntensity(), p2.getMeanIntensity(),
                        "Mean intensity mismatch @ [%d]", i);
                ExtraAssertions.assertArrayEquals(p1.getParameters(), p2.getParameters(), "Params mismatch @ [%d]", i);
                if (showDeviations)
                    ExtraAssertions.assertArrayEquals(p1.getParameterDeviations(), p2.getParameterDeviations(),
                            "Params StdDev mismatch @ [%d]", i);
                if (showPrecision)
                    ExtraAssertions.assertEquals(p1.getPrecision(), p2.getPrecision(), "Precision mismatch @ [%d]", i);
                continue;
            }

            // Otherwise have an error
            final double delta = 1e-5f;
            ExtraAssertions.assertEqualsRelative(p1.getOrigValue(), p2.getOrigValue(), delta,
                    "Orig value mismatch @ [%d]", i);
            ExtraAssertions.assertEqualsRelative(p1.getError(), p2.getError(), delta, "Error mismatch @ [%d]", i);
            ExtraAssertions.assertEqualsRelative(p1.getNoise(), p2.getNoise(), delta, "Noise mismatch @ [%d]", i);
            ExtraAssertions.assertEqualsRelative(p1.getMeanIntensity(), p2.getMeanIntensity(), delta,
                    "Mean intensity mismatch @ [%d]", i);
            ExtraAssertions.assertArrayEqualsRelative(p1.getParameters(), p2.getParameters(), delta,
                    "Params mismatch @ [%d]", i);
            if (showDeviations)
                ExtraAssertions.assertArrayEqualsRelative(p1.getParameterDeviations(), p2.getParameterDeviations(),
                        delta, "Params StdDev mismatch @ [%d]", i);
            if (showPrecision)
                ExtraAssertions.assertEqualsRelative(p1.getPrecision(), p2.getPrecision(), delta,
                        "Precision mismatch @ [%d]", i);
        }

        // Check the header information
        Assertions.assertEquals(expectedResults.getName(), actualResults.getName(), "Name");
        Assertions.assertEquals(expectedResults.getConfiguration(), actualResults.getConfiguration(), "Configuration");

        final Rectangle r1 = expectedResults.getBounds();
        final Rectangle r2 = actualResults.getBounds();
        if (r1 != null)
        {
            Assertions.assertNotNull(r2, "Bounds");
            Assertions.assertEquals(r1.x, r2.x, "Bounds x");
            Assertions.assertEquals(r1.y, r2.y, "Bounds y");
            Assertions.assertEquals(r1.width, r2.width, "Bounds width");
            Assertions.assertEquals(r1.height, r2.height, "Bounds height");
        }
        else
            Assertions.assertNull(r2, "Bounds");

        final Calibration c1 = expectedResults.getCalibration();
        final Calibration c2 = actualResults.getCalibration();
        //try
        //{
        //	Printer printer = JsonFormat.printer().omittingInsignificantWhitespace()
        //	//.includingDefaultValueFields()
        //	;
        //	System.out.println(printer.print(c1));
        //	System.out.println(printer.print(c2));
        //}
        //catch (InvalidProtocolBufferException e)
        //{
        //	// This shouldn't happen so throw it
        //}

        if (c1 != null)
        {
            Assertions.assertNotNull(c2, "Calibration");
            Assertions.assertTrue(c1.equals(c2), "Calibration");
        }
        else
            Assertions.assertNull(c2, "Calibration");

        final PSF p1 = expectedResults.getPSF();
        final PSF p2 = actualResults.getPSF();
        if (p1 != null)
        {
            Assertions.assertNotNull(p2, "PSF");
            Assertions.assertTrue(p1.equals(p2), "PSF");
        }
        else
            Assertions.assertNull(p2, "PSF");
    }

    private static MemoryPeakResults createResults(UniformRandomProvider rg, int i, boolean showDeviations,
            boolean showEndFrame, boolean showId, boolean showPrecision)
    {
        final double bias = rg.nextDouble();

        final boolean extended = showEndFrame || showId || showPrecision;

        final MemoryPeakResults results = new MemoryPeakResults(
                PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D));
        while (i-- > 0)
        {
            final int startFrame = rg.nextInt(i + 1);
            final int origX = rg.nextInt(256);
            final int origY = rg.nextInt(256);
            final float origValue = rg.nextFloat();
            final double error = rg.nextDouble();
            final float noise = rg.nextFloat();
            final float meanIntensity = rg.nextFloat();
            final float[] params = createData(rg);
            final float[] paramsStdDev = (showDeviations) ? createData(rg) : null;
            if (extended)
            {
                final AttributePeakResult r = new AttributePeakResult(startFrame, origX, origY, origValue, error, noise,
                        meanIntensity, params, paramsStdDev);
                if (showEndFrame)
                    r.setEndFrame(startFrame + rg.nextInt(10));
                if (showId)
                    r.setId(i + 1);
                if (showPrecision)
                    r.setPrecision(rg.nextDouble());
                results.add(r);
            }
            else
                results.add(startFrame, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev);
        }
        results.setName(Long.toString(rg.nextLong()));
        results.setConfiguration(Long.toString(rg.nextLong()));
        results.setBounds(new Rectangle(rg.nextInt(10), rg.nextInt(10), rg.nextInt(100), rg.nextInt(100)));
        final CalibrationWriter cal = new CalibrationWriter();
        cal.setNmPerPixel(rg.nextDouble());
        cal.setCountPerPhoton(rg.nextDouble());
        cal.setExposureTime(rg.nextDouble());
        cal.setReadNoise(rg.nextDouble());
        cal.setBias(bias);
        cal.setQuantumEfficiency(rg.nextDouble());
        // Subtract 1 to avoid the additional UNRECOGNISED enum value
        cal.setCameraType(CameraType.values()[rg.nextInt(CameraType.values().length - 1)]);
        cal.setDistanceUnit(DistanceUnit.values()[rg.nextInt(DistanceUnit.values().length - 1)]);
        cal.setIntensityUnit(IntensityUnit.values()[rg.nextInt(IntensityUnit.values().length - 1)]);
        cal.setAngleUnit(AngleUnit.values()[rg.nextInt(AngleUnit.values().length - 1)]);
        results.setCalibration(cal.getCalibration());
        return results;
    }

    private static float[] createData(UniformRandomProvider rg)
    {
        return Gaussian2DPeakResultHelper.createTwoAxisAndAngleParams(rg.nextFloat(), rg.nextFloat(), rg.nextFloat(),
                rg.nextFloat(), rg.nextFloat(), rg.nextFloat(), rg.nextFloat(), rg.nextFloat());
    }

    private static String createFile()
    {
        File file;
        try
        {
            file = File.createTempFile("test", null);
            file.deleteOnExit();
            final String filename = file.getPath();
            return filename;
        }
        catch (final IOException e)
        {
            Assertions.fail("Cannot create temp files for IO testing");
        }
        return null; // Allow compilation but the assert will stop the code
    }

    private static void writeFile(boolean sequential, ResultsFileFormat fileFormat, boolean showDeviations,
            boolean showEndFrame, boolean showId, boolean showPrecision, boolean sort, MemoryPeakResults results,
            String filename)
    {
        final PeakResults out;
        switch (fileFormat)
        {
            case BINARY:
                out = new BinaryFilePeakResults(filename, showDeviations, showEndFrame, showId, showPrecision);
                break;
            case TEXT:
                out = new TextFilePeakResults(filename, showDeviations, showEndFrame, showId, showPrecision);
                break;
            case TSF:
                out = new TSFPeakResultsWriter(filename);
                break;
            case MALK:
                out = new MALKFilePeakResults(filename);
                break;
            default:
                throw new NotImplementedException("Unsupported file format: " + fileFormat);
        }
        out.copySettings(results);
        if (sort && out instanceof FilePeakResults)
            ((FilePeakResults) out).setSortAfterEnd(sort);
        out.begin();

        // TODO - option to test adding using:
        // add(peak);
        // addAll(PeakResult[])

        if (sequential)
            results.forEach(new PeakResultProcedure()
            {
                @Override
                public void execute(PeakResult peak)
                {
                    out.add(peak.getFrame(), peak.getOrigX(), peak.getOrigY(), peak.getOrigValue(), peak.getError(),
                            peak.getNoise(), peak.getMeanIntensity(), peak.getParameters(),
                            peak.getParameterDeviations());
                }
            });
        else
            out.addAll(Arrays.asList(results.toArray()));
        out.end();
    }

    private static MemoryPeakResults readFile(String filename, boolean useScanner)
    {
        return readFile(filename, useScanner, true);
    }

    private static MemoryPeakResults readFile(String filename, boolean useScanner, boolean rawResults)
    {
        final PeakResultsReader reader = new PeakResultsReader(filename);
        reader.setUseScanner(useScanner);
        reader.setRawResults(rawResults);
        final MemoryPeakResults in = reader.getResults();
        return in;
    }

    private static long getReadTime(String filename, final boolean useScanner, final int loops)
    {
        // Initialise reading code
        readFile(filename, useScanner);

        final long start = System.nanoTime();
        for (int i = loops; i-- > 0;)
            readFile(filename, useScanner);
        final long time = System.nanoTime() - start;
        return time;
    }
}
