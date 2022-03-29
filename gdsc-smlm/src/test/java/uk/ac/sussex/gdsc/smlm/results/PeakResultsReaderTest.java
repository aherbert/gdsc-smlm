/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.utils.rng.RadixStringSampler;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.api.function.FloatFloatBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogging.TestLevel;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.FormatSupplier;
import uk.ac.sussex.gdsc.test.utils.functions.ObjectArrayFormatSupplier;

@SuppressWarnings({"javadoc"})
class PeakResultsReaderTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(PeakResultsReaderTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  private static final boolean[] TRUE_FALSE = new boolean[] {true, false};

  // -=-=-=-=-

  @SeededTest
  void writeTextMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, false, false, false);
  }

  @SeededTest
  void writeSequentialTextMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, true, ResultsFileFormat.TEXT, false, false, false, false, false, false);
  }

  @SeededTest
  void writeTextWithDeviationsMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, true, false, false, false, false, false);
  }

  @SeededTest
  void writeTextWithEndFrameMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, true, false, false, false, false);
  }

  @SeededTest
  void writeTextWithIdMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, true, false, false, false);
  }

  @SeededTest
  void writeTextWithPrecisionMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, true, false, false);
  }

  @SeededTest
  void writeTextWithCategoryMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, false, true, false);
  }

  @SeededTest
  void writeTextWithCombinationsMatchesRead(RandomSeed seed) {
    writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.TEXT, false);
  }

  // -=-=-=-=-

  @SeededTest
  void writeBinaryMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, false, false,
        false);
  }

  @SeededTest
  void writeSequentialBinaryMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, true, ResultsFileFormat.BINARY, false, false, false, false, false,
        false);
  }

  @SeededTest
  void writeBinaryWithDeviationsMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, true, false, false, false, false,
        false);
  }

  @SeededTest
  void writeBinaryWithEndFrameMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, true, false, false, false,
        false);
  }

  @SeededTest
  void writeBinaryWithIdMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, true, false, false,
        false);
  }

  @SeededTest
  void writeBinaryWithPrecisionMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, true, false,
        false);
  }

  @SeededTest
  void writeBinaryWithCategoryPrecisionMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, false, true,
        false);
  }

  @SeededTest
  void writeBinaryWithCombinationsMatchesRead(RandomSeed seed) {
    writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.BINARY, false);
  }

  // -=-=-=-=-

  @SeededTest
  void writeTextWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, false, false, true);
  }

  @SeededTest
  void writeSequentialTextWithSortMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, true, ResultsFileFormat.TEXT, false, false, false, false, false, true);
  }

  @SeededTest
  void writeTextWithDeviationsWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, true, false, false, false, false, true);
  }

  @SeededTest
  void writeTextWithEndFrameWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, true, false, false, false, true);
  }

  @SeededTest
  void writeTextWithIdWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, true, false, false, true);
  }

  @SeededTest
  void writeTextWithPrecisionWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, true, false, true);
  }

  @SeededTest
  void writeTextWithCategoryWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TEXT, false, false, false, false, true, true);
  }

  @SeededTest
  void writeTextWithCombinationsWithSortMatchesRead(RandomSeed seed) {
    writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.TEXT, true);
  }

  // -=-=-=-=-

  @SeededTest
  void writeBinaryWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, false, false,
        true);
  }

  @SeededTest
  void writeSequentialBinaryWithSortMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, true, ResultsFileFormat.BINARY, false, false, false, false, false, true);
  }

  @SeededTest
  void writeBinaryWithDeviationsWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, true, false, false, false, false, true);
  }

  @SeededTest
  void writeBinaryWithEndFrameWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, true, false, false, false, true);
  }

  @SeededTest
  void writeBinaryWithIdWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, true, false, false, true);
  }

  @SeededTest
  void writeBinaryWithPrecisionWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, true, false, true);
  }

  @SeededTest
  void writeBinaryWithCategoryWithSortMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.BINARY, false, false, false, false, true, true);
  }

  @SeededTest
  void writeBinaryWithCombinationsWithSortMatchesRead(RandomSeed seed) {
    writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.BINARY, true);
  }

  // -=-=-=-=-

  // Note: For Malk we cannot do all the tests as the format only contains X,Y,T,I

  @SeededTest
  void writeMalkMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, false, ResultsFileFormat.MALK, false, false, false, false, false, false);
  }

  @SeededTest
  void writeSequentialMalkMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, true, ResultsFileFormat.MALK, false, false, false, false, false, false);
  }

  @SeededTest
  void writeMalkWithSortMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, false, ResultsFileFormat.MALK, false, false, false, false, false, true);
  }

  @SeededTest
  void writeSequentialMalkWithSortMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, true, ResultsFileFormat.MALK, false, false, false, false, false, true);
  }

  // -=-=-=-=-

  // Note: For Tsf we cannot specify as binary because the widths are converted into a
  // different format and then back again.

  @SeededTest
  void writeTsfMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, false, false, false, false, false);
  }

  @SeededTest
  void writeSequentialTsfMatchesRead(RandomSeed seed) {
    writeMatchesRead(seed, true, ResultsFileFormat.TSF, false, false, false, false, false, false);
  }

  @SeededTest
  void writeTsfWithDeviationsMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TSF, true, false, false, false, false, false);
  }

  @SeededTest
  void writeTsfWithEndFrameMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, true, false, false, false, false);
  }

  @SeededTest
  void writeTsfWithIdMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, false, true, false, false, false);
  }

  @SeededTest
  void writeTsfWithPrecisionMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, false, false, true, false, false);
  }

  @SeededTest
  void writeTsfWithCategoryMatchesRead(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeMatchesRead(seed, false, ResultsFileFormat.TSF, false, false, false, false, true, false);
  }

  @SeededTest
  void writeTsfWithCombinationsMatchesRead(RandomSeed seed) {
    writeWithCombinationsMatchesRead(seed, false, ResultsFileFormat.TSF, false);
  }

  // -=-=-=-=-

  @SeededTest
  void readWithScannerMatchesNonScanner(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, false, false, false, false);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithDeviations(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, true, false, false, false, false, false);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithEndFrame(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, true, false, false, false, false);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithId(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, true, false, false, false);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithPrecision(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, false, true, false, false);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithCategory(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, false, false, true, false);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithCombinations(RandomSeed seed) {
    // Scanner is not a default so do not always test
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    checkScannerMatchesNonScannerWithCombinations(seed, false);
  }

  // -=-=-=-=-

  @SeededTest
  void readWithScannerMatchesNonScannerWithSort(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, false, false, false, true);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithDeviationsWithSort(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, true, false, false, false, false, true);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithEndFrameWithSort(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, true, false, false, false, true);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithIdWithSort(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, true, false, false, true);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithPrecisionWithSort(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, false, true, false, true);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithCategoryWithSort(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    checkScannerMatchesNonScanner(seed, false, false, false, false, true, true);
  }

  @SeededTest
  void readWithScannerMatchesNonScannerWithCombinationsWithSort(RandomSeed seed) {
    // Scanner is not a default so do not always test
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    checkScannerMatchesNonScannerWithCombinations(seed, true);
  }

  // -=-=-=-=-

  @SeededTest
  void readTextWithNonScannerIsFasterThanScanner(RandomSeed seed) {
    readWith2IsFasterThan1(seed, false, false, false, false, false, ResultsFileFormat.TEXT, true,
        ResultsFileFormat.TEXT, false, 1);
  }

  @SeededTest
  public void readTextWithNonScannerIsFasterThanScannerWithExtended(RandomSeed seed) {
    readWith2IsFasterThan1(seed, true, true, true, true, true, ResultsFileFormat.TEXT, true,
        ResultsFileFormat.TEXT, false, 1);
  }

  @SeededTest
  void readWithMalkIsFasterThanText(RandomSeed seed) {
    readWith2IsFasterThan1(seed, false, false, false, false, false, ResultsFileFormat.TEXT, false,
        ResultsFileFormat.MALK, false, 2);
  }

  @SeededTest
  void readWithBinaryIsFasterThanText(RandomSeed seed) {
    readWith2IsFasterThan1(seed, false, false, false, false, false, ResultsFileFormat.TEXT, false,
        ResultsFileFormat.BINARY, false, 2);
  }

  @SeededTest
  void readWithBinaryIsFasterThanTsf(RandomSeed seed) {
    readWith2IsFasterThan1(seed, false, false, false, false, false, ResultsFileFormat.TSF, false,
        ResultsFileFormat.BINARY, false, 20);
  }

  @SeededTest
  void canConvertMalkToNmAndPhotons(RandomSeed seed) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final MemoryPeakResults out = createResults(rg, 200, false, false, false, false, false);

    // Output in pixel and count
    final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
    cal.setDistanceUnit(DistanceUnit.PIXEL);
    cal.setIntensityUnit(IntensityUnit.COUNT);
    out.setCalibration(cal.getCalibration());
    out.setPsf(PsfHelper.create(PSFType.CUSTOM));

    final String filename = createFile();

    writeFile(false, ResultsFileFormat.MALK, false, false, false, false, false, false, out,
        filename);

    final MemoryPeakResults in = readFile(filename, false);

    // Change to nm and photon for the validation
    out.convertToUnits(DistanceUnit.NM, IntensityUnit.PHOTON, null);

    checkEqual(ResultsFileFormat.MALK, false, false, false, false, false, false, out, in);
  }

  @SeededTest
  void writeTextWithComputedPrecisionMatchesRead(RandomSeed seed) {
    // Create without precision
    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final MemoryPeakResults results = createResults(rg, 200, false, false, false, false, false);
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

    checkEqual(ResultsFileFormat.TEXT, false, false, false, true, false, false, results, in);
  }

  @SeededTest
  void canReadTextIntoPreferredUnits(RandomSeed seed) {
    canReadIntoPreferredUnits(seed, ResultsFileFormat.TEXT);
  }

  @SeededTest
  void canReadBinaryIntoPreferredUnits(RandomSeed seed) {
    canReadIntoPreferredUnits(seed, ResultsFileFormat.BINARY);
  }

  @SeededTest
  void canReadTsfIntoPreferredUnits(RandomSeed seed) {
    canReadIntoPreferredUnits(seed, ResultsFileFormat.TSF);
  }

  private static void canReadIntoPreferredUnits(RandomSeed seed, ResultsFileFormat fileFormat) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final MemoryPeakResults out = createResults(rg, 200, false, false, false, false, false);

    // Output in nm and count
    final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
    cal.setDistanceUnit(DistanceUnit.NM);
    cal.setIntensityUnit(IntensityUnit.COUNT);
    if (fileFormat == ResultsFileFormat.TSF) {
      // For now just support using the native float TSF value
      cal.setNmPerPixel((float) cal.getNmPerPixel());
    }
    out.setCalibration(cal.getCalibration());

    final String filename = createFile();

    writeFile(false, fileFormat, false, false, false, false, false, false, out, filename);

    final MemoryPeakResults in = readFile(filename, false, false);

    // Change to preferred units
    out.convertToUnits(MemoryPeakResults.PREFERRED_DISTANCE_UNIT,
        MemoryPeakResults.PREFERRED_INTENSITY_UNIT, MemoryPeakResults.PREFERRED_ANGLE_UNIT);

    checkEqual(fileFormat, false, false, false, false, false, false, out, in);
  }

  @SeededTest
  void canReadTextAndSimplifyGaussian2DPsf(RandomSeed seed) {
    canReadAndSimplifyGaussian2DPsf(seed, ResultsFileFormat.TEXT);
  }

  @SeededTest
  void canReadBinaryAndSimplifyGaussian2DPsf(RandomSeed seed) {
    canReadAndSimplifyGaussian2DPsf(seed, ResultsFileFormat.BINARY);
  }

  @SeededTest
  void canReadTsfAndSimplifyGaussian2DPsf(RandomSeed seed) {
    canReadAndSimplifyGaussian2DPsf(seed, ResultsFileFormat.TSF);
  }

  private static void canReadAndSimplifyGaussian2DPsf(RandomSeed seed,
      ResultsFileFormat fileFormat) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final MemoryPeakResults out = createResults(rg, 1, false, false, false, false, false);

    final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
    cal.setDistanceUnit(MemoryPeakResults.PREFERRED_DISTANCE_UNIT);
    cal.setIntensityUnit(MemoryPeakResults.PREFERRED_INTENSITY_UNIT);
    cal.setAngleUnit(MemoryPeakResults.PREFERRED_ANGLE_UNIT);
    if (fileFormat == ResultsFileFormat.TSF) {
      // For now just support using the native float TSF value
      cal.setNmPerPixel((float) cal.getNmPerPixel());
    }
    out.setCalibration(cal.getCalibration());

    // Remove angle
    final int ia = PsfHelper.getGaussian2DAngleIndex(out.getPsf());
    out.forEach(new PeakResultProcedure() {
      @Override
      public void execute(PeakResult peakResult) {
        peakResult.getParameters()[ia] = 0;
      }
    });

    final String filename = createFile();

    writeFile(false, fileFormat, false, false, false, false, false, false, out, filename);

    MemoryPeakResults in = readFile(filename, false, false);

    // Change to two-axis PSF
    out.setPsf(PsfHelper.create(PSFType.TWO_AXIS_GAUSSIAN_2D));
    final int twoAxisLength =
        PsfHelper.getParameterCount(out.getPsf()) + PeakResult.STANDARD_PARAMETERS;
    out.forEach(new PeakResultProcedure() {
      @Override
      public void execute(PeakResult peakResult) {
        peakResult.resizeParameters(twoAxisLength);
      }
    });

    checkEqual(fileFormat, false, false, false, false, false, false, out, in);

    // Remove sy
    final int[] indices = PsfHelper.getGaussian2DWxWyIndices(out.getPsf());
    final int isx = indices[0];
    final int isy = indices[1];
    out.forEach(new PeakResultProcedure() {
      @Override
      public void execute(PeakResult peakResult) {
        final float[] p = peakResult.getParameters();
        p[isy] = p[isx];
      }
    });

    writeFile(false, fileFormat, false, false, false, false, false, false, out, filename);

    in = readFile(filename, false, false);

    // Change to one-axis PSF
    out.setPsf(PsfHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D));
    final int oneAxisLength =
        PsfHelper.getParameterCount(out.getPsf()) + PeakResult.STANDARD_PARAMETERS;
    out.forEach(new PeakResultProcedure() {
      @Override
      public void execute(PeakResult peakResult) {
        peakResult.resizeParameters(oneAxisLength);
      }
    });

    checkEqual(fileFormat, false, false, false, false, false, false, out, in);
  }

  private static void readWith2IsFasterThan1(RandomSeed seed, boolean showDeviations,
      boolean showEndFrame, boolean showId, boolean showPrecision, boolean showCategory,
      ResultsFileFormat f1, boolean useScanner1, ResultsFileFormat f2, boolean useScanner2,
      int loops) {
    Assumptions.assumeTrue(logger.isLoggable(TestLevel.TEST_INFO));
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));

    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final MemoryPeakResults out =
        createResults(rg, 20000, showDeviations, showEndFrame, showId, showPrecision, showCategory);
    final String filename = createFile();

    writeFile(false, f1, showDeviations, showEndFrame, showId, showPrecision, showCategory, false,
        out, filename);
    final long time1 = getReadTime(filename, useScanner1, loops);

    writeFile(false, f2, showDeviations, showEndFrame, showId, showPrecision, showCategory, false,
        out, filename);
    final long time2 = getReadTime(filename, useScanner2, loops);

    if (useScanner1 != useScanner2) {
      logger.log(TestLevel.TEST_INFO,
          FormatSupplier.getSupplier("%s (scan=%b) is %.2fx faster than %s (scan=%b)", f2,
              useScanner2, (double) time1 / time2, f1, useScanner1));
    } else {
      logger.log(TestLevel.TEST_INFO,
          FormatSupplier.getSupplier("%s is %.2fx faster than %s", f2, (double) time1 / time2, f1));
    }
    Assertions.assertTrue(time2 < time1,
        () -> String.format("%s (%d) is not faster than %s (%d)", f2, time2, f1, time1));
  }

  // -=-=-=-=-

  private static void writeMatchesRead(RandomSeed seed, boolean sequential,
      ResultsFileFormat fileFormat, boolean showDeviations, boolean showEndFrame, boolean showId,
      boolean showPrecision, boolean showCategory, boolean sort) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final MemoryPeakResults out =
        createResults(rg, 200, showDeviations, showEndFrame, showId, showPrecision, showCategory);
    if (fileFormat == ResultsFileFormat.MALK) {
      final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
      cal.setDistanceUnit(DistanceUnit.NM);
      cal.setIntensityUnit(IntensityUnit.PHOTON);
      out.setCalibration(cal.getCalibration());
      out.setPsf(PsfHelper.create(PSFType.CUSTOM));
    }
    if (fileFormat == ResultsFileFormat.TSF) {
      final CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
      // For now just support using the native float TSF datatype
      cal.setNmPerPixel((float) cal.getNmPerPixel());
      out.setCalibration(cal.getCalibration());
      // TSF converts the width parameters so make sure they are not zero
      out.forEach(new PeakResultProcedure() {

        @Override
        public void execute(PeakResult peakResult) {
          check(peakResult.getParameters());
          if (showDeviations) {
            check(peakResult.getParameterDeviations());
          }
        }

        private void check(float[] parameters) {
          for (int i = 0; i < parameters.length; i++) {
            if (parameters[i] == 0) {
              parameters[i] = 0.1f;
            }
          }
        }
      });
    }

    // System.out.println(out.getCalibration());
    // System.out.println(out.getPSF().toString());
    //
    // System.out.println(TextFormat.shortDebugString(out.getCalibration()));
    //
    // try
    // {
    // Printer printer = JsonFormat.printer()
    // .omittingInsignificantWhitespace()
    // //.includingDefaultValueFields()
    // ;
    // System.out.println(printer.print(out.getCalibration()));
    // System.out.println(printer.print(out.getPSF()));
    // }
    // catch (InvalidProtocolBufferException e)
    // {
    // // This shouldn't happen so throw it
    // }

    final String filename = createFile();

    writeFile(sequential, fileFormat, showDeviations, showEndFrame, showId, showPrecision,
        showCategory, sort, out, filename);

    final MemoryPeakResults in = readFile(filename, false);

    checkEqual(fileFormat, showDeviations, showEndFrame, showId, showPrecision, showCategory, sort,
        out, in);
  }

  private static void writeWithCombinationsMatchesRead(RandomSeed seed, boolean sequential,
      ResultsFileFormat fileFormat, boolean sort) {
    for (final boolean showDeviations : TRUE_FALSE) {
      for (final boolean showEndFrame : TRUE_FALSE) {
        for (final boolean showId : TRUE_FALSE) {
          for (final boolean showPrecision : TRUE_FALSE) {
            for (final boolean showCategory : TRUE_FALSE) {
              if (count(showDeviations, showEndFrame, showId, showPrecision, showCategory) < 2) {
                continue;
              }
              writeMatchesRead(seed, sequential, fileFormat, showDeviations, showEndFrame, showId,
                  showPrecision, showCategory, sort);
            }
          }
        }
      }
    }
  }

  private static int count(boolean... flags) {
    int count = 0;
    for (final boolean flag : flags) {
      if (flag) {
        count++;
      }
    }
    return count;
  }

  private static void checkScannerMatchesNonScanner(RandomSeed seed, boolean showDeviations,
      boolean showEndFrame, boolean showId, boolean showPrecision, boolean showCategory,
      boolean sort) {
    final UniformRandomProvider rg = RngFactory.create(seed.get());
    final MemoryPeakResults out =
        createResults(rg, 1000, showDeviations, showEndFrame, showId, showPrecision, showCategory);
    final String filename = createFile();

    final ResultsFileFormat fileFormat = ResultsFileFormat.TEXT;
    writeFile(false, fileFormat, showDeviations, showEndFrame, showId, showPrecision, showCategory,
        sort, out, filename);

    final MemoryPeakResults in = readFile(filename, false);
    final MemoryPeakResults in2 = readFile(filename, true);

    checkEqual(fileFormat, showDeviations, showEndFrame, showId, showPrecision, showCategory, sort,
        in, in2);
  }

  private static void checkScannerMatchesNonScannerWithCombinations(RandomSeed seed, boolean sort) {
    for (final boolean showDeviations : TRUE_FALSE) {
      for (final boolean showEndFrame : TRUE_FALSE) {
        for (final boolean showId : TRUE_FALSE) {
          for (final boolean showPrecision : TRUE_FALSE) {
            for (final boolean showCategory : TRUE_FALSE) {
              if (count(showDeviations, showEndFrame, showId, showPrecision, showCategory) < 2) {
                continue;
              }
              checkScannerMatchesNonScanner(seed, showDeviations, showEndFrame, showId,
                  showPrecision, showCategory, sort);
            }
          }
        }
      }
    }
  }

  private static void checkEqual(ResultsFileFormat fileFormat, boolean showDeviations,
      boolean showEndFrame, boolean showId, boolean showPrecision, boolean showCategory,
      boolean sort, MemoryPeakResults expectedResults, MemoryPeakResults actualResults) {
    Assertions.assertNotNull(actualResults, "Input results are null");
    Assertions.assertEquals(expectedResults.size(), actualResults.size(), "Size differ");

    final PeakResult[] expected = expectedResults.toArray();
    final PeakResult[] actual = actualResults.toArray();
    if (sort) {
      // Results should be sorted by time
      Arrays.sort(expected, (o1, o2) -> o1.getFrame() - o2.getFrame());
    }

    // TSF requires the bias be subtracted
    // double bias = expectedResults.getCalibration().getBias();

    final DoubleDoubleBiPredicate deltaD = Predicates.doublesIsRelativelyCloseTo(1e-5);
    final FloatFloatBiPredicate deltaF = Predicates.floatsIsRelativelyCloseTo(1e-5);

    for (int i = 0; i < actualResults.size(); i++) {
      final PeakResult p1 = expected[i];
      final PeakResult p2 = actual[i];
      final ObjectArrayFormatSupplier msg = new ObjectArrayFormatSupplier("%s @ [" + i + "]", 1);

      Assertions.assertEquals(p1.getFrame(), p2.getFrame(), msg.set(0, "Peak"));

      if (fileFormat == ResultsFileFormat.MALK) {
        TestAssertions.assertTest(p1.getXPosition(), p2.getXPosition(), deltaF, msg.set(0, "X"));
        TestAssertions.assertTest(p1.getYPosition(), p2.getYPosition(), deltaF, msg.set(0, "Y"));
        TestAssertions.assertTest(p1.getIntensity(), p2.getIntensity(), deltaF,
            msg.set(0, "Intensity"));
        continue;
      }

      Assertions.assertEquals(p1.getOrigX(), p2.getOrigX(), msg.set(0, "Orig X"));
      Assertions.assertEquals(p1.getOrigY(), p2.getOrigY(), msg.set(0, "Orig Y"));
      Assertions.assertNotNull(p2.getParameters(), msg.set(0, "Params is null"));
      if (showEndFrame) {
        Assertions.assertEquals(p1.getEndFrame(), p2.getEndFrame(), msg.set(0, "End frame"));
      }
      if (showId) {
        Assertions.assertEquals(p1.getId(), p2.getId(), msg.set(0, "ID"));
      }
      if (showDeviations) {
        Assertions.assertNotNull(p2.getParameterDeviations(), msg.set(0, "Deviations"));
      }
      if (showCategory) {
        Assertions.assertEquals(p1.getCategory(), p2.getCategory(), msg.set(0, "Category"));
      }

      // Binary should be exact for float numbers
      if (fileFormat == ResultsFileFormat.BINARY) {
        Assertions.assertEquals(p1.getOrigValue(), p2.getOrigValue(), msg.set(0, "Orig value"));
        Assertions.assertEquals(p1.getError(), p2.getError(), msg.set(0, "Error"));
        Assertions.assertEquals(p1.getNoise(), p2.getNoise(), msg.set(0, "Noise"));
        Assertions.assertEquals(p1.getMeanIntensity(), p2.getMeanIntensity(),
            msg.set(0, "Mean intensity"));
        Assertions.assertArrayEquals(p1.getParameters(), p2.getParameters(), msg.set(0, "Params"));
        if (showDeviations) {
          Assertions.assertArrayEquals(p1.getParameterDeviations(), p2.getParameterDeviations(),
              msg.set(0, "Params StdDev"));
        }
        if (showPrecision) {
          Assertions.assertEquals(p1.getPrecision(), p2.getPrecision(), msg.set(0, "Precision"));
        }
        continue;
      }

      // Otherwise have an error
      TestAssertions.assertTest(p1.getOrigValue(), p2.getOrigValue(), deltaF,
          msg.set(0, "Orig value"));
      TestAssertions.assertTest(p1.getError(), p2.getError(), deltaD, msg.set(0, "Error"));
      TestAssertions.assertTest(p1.getNoise(), p2.getNoise(), deltaF, msg.set(0, "Noise"));
      TestAssertions.assertTest(p1.getMeanIntensity(), p2.getMeanIntensity(), deltaF,
          msg.set(0, "Mean intensity"));
      TestAssertions.assertArrayTest(p1.getParameters(), p2.getParameters(), deltaF,
          msg.set(0, "Params"));
      if (showDeviations) {
        TestAssertions.assertArrayTest(p1.getParameterDeviations(), p2.getParameterDeviations(),
            deltaF, msg.set(0, "Params StdDev"));
      }
      if (showPrecision) {
        // Handle NaN precisions
        final double pa = p1.getPrecision();
        final double pb = p2.getPrecision();
        if (!Double.isNaN(pa) || !Double.isNaN(pb)) {
          TestAssertions.assertTest(p1.getPrecision(), p2.getPrecision(), deltaD,
              msg.set(0, "Precision"));
        }
      }
    }

    // Check the header information
    Assertions.assertEquals(expectedResults.getName(), actualResults.getName(), "Name");
    Assertions.assertEquals(expectedResults.getConfiguration(), actualResults.getConfiguration(),
        "Configuration");

    final Rectangle r1 = expectedResults.getBounds();
    final Rectangle r2 = actualResults.getBounds();
    if (r1 != null) {
      Assertions.assertNotNull(r2, "Bounds");
      Assertions.assertEquals(r1.x, r2.x, "Bounds x");
      Assertions.assertEquals(r1.y, r2.y, "Bounds y");
      Assertions.assertEquals(r1.width, r2.width, "Bounds width");
      Assertions.assertEquals(r1.height, r2.height, "Bounds height");
    } else {
      Assertions.assertNull(r2, "Bounds");
    }

    final Calibration c1 = expectedResults.getCalibration();
    final Calibration c2 = actualResults.getCalibration();
    // try
    // {
    // Printer printer = JsonFormat.printer().omittingInsignificantWhitespace()
    // //.includingDefaultValueFields()
    // ;
    // System.out.println(printer.print(c1));
    // System.out.println(printer.print(c2));
    // }
    // catch (InvalidProtocolBufferException e)
    // {
    // // This shouldn't happen so throw it
    // }

    if (c1 != null) {
      Assertions.assertNotNull(c2, "Calibration");
      // Be lenient and allow no TimeUnit to match TimeUnit.FRAME
      boolean ok = c1.equals(c2);
      if (!ok && new CalibrationReader(c1).getTimeUnitValue() == TimeUnit.TIME_UNIT_NA_VALUE) {
        switch (fileFormat) {
          case BINARY:
          case MALK:
          case TEXT:
          case TSF:
            final CalibrationWriter writer = new CalibrationWriter(c1);
            writer.setTimeUnit(TimeUnit.FRAME);
            ok = writer.getCalibration().equals(c2);
            break;
          default:
            // Do not assume frames for other file formats
            break;
        }
      }
      Assertions.assertTrue(ok, "Calibration");
    } else {
      Assertions.assertNull(c2, "Calibration");
    }

    final PSF p1 = expectedResults.getPsf();
    final PSF p2 = actualResults.getPsf();
    if (p1 != null) {
      Assertions.assertNotNull(p2, "PSF");
      Assertions.assertTrue(p1.equals(p2), "PSF");
    } else {
      Assertions.assertNull(p2, "PSF");
    }
  }

  private static MemoryPeakResults createResults(UniformRandomProvider rg, int size,
      boolean showDeviations, boolean showEndFrame, boolean showId, boolean showPrecision,
      boolean showCategory) {
    final double bias = rg.nextDouble();

    final boolean extended = showEndFrame || showId || showPrecision || showCategory;

    final MemoryPeakResults results =
        new MemoryPeakResults(PsfHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D));
    while (size-- > 0) {
      final int startFrame = rg.nextInt(size + 1);
      final int origX = rg.nextInt(256);
      final int origY = rg.nextInt(256);
      final float origValue = rg.nextFloat();
      final double error = rg.nextDouble();
      final float noise = rg.nextFloat();
      final float meanIntensity = rg.nextFloat();
      final float[] params = createData(rg);
      final float[] paramsStdDev = (showDeviations) ? createData(rg) : null;
      if (extended) {
        final AttributePeakResult r = new AttributePeakResult(startFrame, origX, origY, origValue,
            error, noise, meanIntensity, params, paramsStdDev);
        if (showEndFrame) {
          r.setEndFrame(startFrame + rg.nextInt(10));
        }
        if (showId) {
          r.setId(size + 1);
        }
        if (showPrecision) {
          r.setPrecision(rg.nextDouble());
        }
        if (showCategory) {
          r.setCategory(size & 64);
        }
        results.add(r);
      } else {
        results.add(startFrame, origX, origY, origValue, error, noise, meanIntensity, params,
            paramsStdDev);
      }
    }
    results.setName(RadixStringSampler.nextBase64String(rg, 16));
    results.setConfiguration(RadixStringSampler.nextBase64String(rg, 16));
    results
        .setBounds(new Rectangle(rg.nextInt(10), rg.nextInt(10), rg.nextInt(100), rg.nextInt(100)));
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

  private static float[] createData(UniformRandomProvider rg) {
    return Gaussian2DPeakResultHelper.createTwoAxisAndAngleParams(rg.nextFloat(), rg.nextFloat(),
        rg.nextFloat(), rg.nextFloat(), rg.nextFloat(), rg.nextFloat(), rg.nextFloat(),
        rg.nextFloat());
  }

  private static String createFile() {
    File file;
    try {
      file = File.createTempFile("test", null);
      file.deleteOnExit();
      final String filename = file.getPath();
      return filename;
    } catch (final IOException ex) {
      Assertions.fail("Cannot create temp files for IO testing");
    }
    return null; // Allow compilation but the assert will stop the code
  }

  private static void writeFile(boolean sequential, ResultsFileFormat fileFormat,
      boolean showDeviations, boolean showEndFrame, boolean showId, boolean showPrecision,
      boolean showCategory, boolean sort, MemoryPeakResults results, String filename) {
    final PeakResults out;
    switch (fileFormat) {
      case BINARY:
        out = new BinaryFilePeakResults(filename, showDeviations, showEndFrame, showId,
            showPrecision, showCategory);
        break;
      case TEXT:
        out = new TextFilePeakResults(filename, showDeviations, showEndFrame, showId, showPrecision,
            showCategory);
        break;
      case TSF:
        out = new TsfPeakResultsWriter(filename);
        break;
      case MALK:
        out = new MalkFilePeakResults(filename);
        break;
      default:
        throw new NotImplementedException("Unsupported file format: " + fileFormat);
    }
    out.copySettings(results);
    if (sort && out instanceof FilePeakResults) {
      ((FilePeakResults) out).setSortAfterEnd(sort);
    }
    out.begin();

    // TODO - option to test adding using:
    // add(peak);
    // addAll(PeakResult[])

    if (sequential) {
      results.forEach(new PeakResultProcedure() {
        @Override
        public void execute(PeakResult peak) {
          out.add(peak.getFrame(), peak.getOrigX(), peak.getOrigY(), peak.getOrigValue(),
              peak.getError(), peak.getNoise(), peak.getMeanIntensity(), peak.getParameters(),
              peak.getParameterDeviations());
        }
      });
    } else {
      out.addAll(Arrays.asList(results.toArray()));
    }
    out.end();
  }

  private static MemoryPeakResults readFile(String filename, boolean useScanner) {
    return readFile(filename, useScanner, true);
  }

  private static MemoryPeakResults readFile(String filename, boolean useScanner,
      boolean rawResults) {
    final PeakResultsReader reader = new PeakResultsReader(filename);
    reader.setUseScanner(useScanner);
    reader.setRawResults(rawResults);
    return reader.getResults();
  }

  private static long getReadTime(String filename, final boolean useScanner, final int loops) {
    // Initialise reading code
    readFile(filename, useScanner);

    final long start = System.nanoTime();
    for (int i = loops; i-- > 0;) {
      readFile(filename, useScanner);
    }
    final long time = System.nanoTime() - start;
    return time;
  }
}
