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
package gdsc.smlm.results;

import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.core.utils.NotImplementedException;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.test.TestAssert;
import gdsc.test.TestSettings;
import gdsc.test.TestSettings.LogLevel;
import gdsc.test.TestSettings.TestComplexity;

@SuppressWarnings({ "javadoc" })
public class PeakResultsReaderTest
{
	static final boolean[] onOff = new boolean[] { true, false };

	// TODO - Add tests to compare writing to a IJTablePeakResults, saving the TextPanel contents to file and then reading.

	// -=-=-=-=-

	@Test
	public void writeTextMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, false, false, false);
	}

	@Test
	public void writeSequentialTextMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.TEXT, false, false, false, false, false);
	}

	@Test
	public void writeTextWithDeviationsMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, true, false, false, false, false);
	}

	@Test
	public void writeTextWithEndFrameMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, true, false, false, false);
	}

	@Test
	public void writeTextWithIdMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, true, false, false);
	}

	@Test
	public void writeTextWithPrecisionMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, false, true, false);
	}

	@Test
	public void writeTextWithCombinationsMatchesRead()
	{
		writeWithCombinationsMatchesRead(false, ResultsFileFormat.TEXT, false);
	}

	// -=-=-=-=-

	@Test
	public void writeBinaryMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, false, false, false);
	}

	@Test
	public void writeSequentialBinaryMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.BINARY, false, false, false, false, false);
	}

	@Test
	public void writeBinaryWithDeviationsMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, true, false, false, false, false);
	}

	@Test
	public void writeBinaryWithEndFrameMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, true, false, false, false);
	}

	@Test
	public void writeBinaryWithIdMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, true, false, false);
	}

	@Test
	public void writeBinaryWithPrecisionMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, false, true, false);
	}

	@Test
	public void writeBinaryWithCombinationsMatchesRead()
	{
		writeWithCombinationsMatchesRead(false, ResultsFileFormat.BINARY, false);
	}

	// -=-=-=-=-

	@Test
	public void writeTextWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, false, false, true);
	}

	@Test
	public void writeSequentialTextWithSortMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.TEXT, false, false, false, false, true);
	}

	@Test
	public void writeTextWithDeviationsWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, true, false, false, false, true);
	}

	@Test
	public void writeTextWithEndFrameWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, true, false, false, true);
	}

	@Test
	public void writeTextWithIdWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, true, false, true);
	}

	@Test
	public void writeTextWithPrecisionWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, false, true, true);
	}

	@Test
	public void writeTextWithCombinationsWithSortMatchesRead()
	{
		writeWithCombinationsMatchesRead(false, ResultsFileFormat.TEXT, true);
	}

	// -=-=-=-=-

	@Test
	public void writeBinaryWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, false, false, true);
	}

	@Test
	public void writeSequentialBinaryWithSortMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.BINARY, false, false, false, false, true);
	}

	@Test
	public void writeBinaryWithDeviationsWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, true, false, false, false, true);
	}

	@Test
	public void writeBinaryWithEndFrameWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, true, false, false, true);
	}

	@Test
	public void writeBinaryWithIdWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, true, false, true);
	}

	@Test
	public void writeBinaryWithPrecisionWithSortMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, false, true, true);
	}

	@Test
	public void writeBinaryWithCombinationsWithSortMatchesRead()
	{
		writeWithCombinationsMatchesRead(false, ResultsFileFormat.BINARY, true);
	}

	// -=-=-=-=-

	// Note: For MALK we cannot do all the tests as the format only contains X,Y,T,I

	@Test
	public void writeMALKMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.MALK, false, false, false, false, false);
	}

	@Test
	public void writeSequentialMALKMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.MALK, false, false, false, false, false);
	}

	@Test
	public void writeMALKWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.MALK, false, false, false, true, false);
	}

	@Test
	public void writeSequentialMALKWithSortMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.MALK, false, false, false, true, false);
	}

	// -=-=-=-=-

	// Note: For TSF we cannot specify as binary because the widths are converted into a
	// different format and then back again.

	@Test
	public void writeTSFMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TSF, false, false, false, false, false);
	}

	@Test
	public void writeSequentialTSFMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.TSF, false, false, false, false, false);
	}

	@Test
	public void writeTSFWithDeviationsMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TSF, true, false, false, false, false);
	}

	@Test
	public void writeTSFWithEndFrameMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TSF, false, true, false, false, false);
	}

	@Test
	public void writeTSFWithIdMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TSF, false, false, true, false, false);
	}

	@Test
	public void writeTSFWithPrecisionMatchesRead()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		writeMatchesRead(false, ResultsFileFormat.TSF, false, false, false, true, false);
	}

	@Test
	public void writeTSFWithCombinationsMatchesRead()
	{
		writeWithCombinationsMatchesRead(false, ResultsFileFormat.TSF, false);
	}

	// -=-=-=-=-

	@Test
	public void readWithScannerMatchesNonScanner()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, false, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviations()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(true, false, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithEndFrame()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, true, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithId()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, false, true, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithPrecision()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, false, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithCombinations()
	{
		TestSettings.assumeLowComplexity(); // Scanner is not a default so do not always test
		readWithScannerMatchesNonScannerWithCombinations(false);
	}

	// -=-=-=-=-

	@Test
	public void readWithScannerMatchesNonScannerWithSort()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, false, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithSort()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(true, false, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithEndFrameWithSort()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, true, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithIdWithSort()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, false, true, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithPrecisionWithSort()
	{
		TestSettings.assumeMediumComplexity(); // The combinations test overlaps this
		readWithScannerMatchesNonScanner(false, false, false, true, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithCombinationsWithSort()
	{
		TestSettings.assumeLowComplexity(); // Scanner is not a default so do not always test
		readWithScannerMatchesNonScannerWithCombinations(true);
	}

	// -=-=-=-=-

	@Test
	public void readTextWithNonScannerIsFasterThanScanner()
	{
		readWith2IsFasterThan1(false, false, false, false, ResultsFileFormat.TEXT, true, ResultsFileFormat.TEXT, false,
				1);
	}

	@Test
	public void readTextWithNonScannerIsFasterThanScannerWithDeviationsWithEndFrameWithIdWithPrecision()
	{
		readWith2IsFasterThan1(true, true, true, true, ResultsFileFormat.TEXT, true, ResultsFileFormat.TEXT, false, 1);
	}

	@Test
	public void readWithMALKIsFasterThanText()
	{
		readWith2IsFasterThan1(false, false, false, false, ResultsFileFormat.TEXT, false, ResultsFileFormat.MALK, false,
				2);
	}

	@Test
	public void readWithBinaryIsFasterThanText()
	{
		readWith2IsFasterThan1(false, false, false, false, ResultsFileFormat.TEXT, false, ResultsFileFormat.BINARY,
				false, 2);
	}

	@Test
	public void readWithBinaryIsFasterThanTSF()
	{
		readWith2IsFasterThan1(false, false, false, false, ResultsFileFormat.TSF, false, ResultsFileFormat.BINARY,
				false, 20);
	}

	@Test
	public void canConvertMalkToNMAndPhotons()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MemoryPeakResults out = createResults(rg, 200, false, false, false, false);

		// Output in pixel and count
		CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
		cal.setDistanceUnit(DistanceUnit.PIXEL);
		cal.setIntensityUnit(IntensityUnit.COUNT);
		out.setCalibration(cal.getCalibration());
		out.setPSF(PSFHelper.create(PSFType.CUSTOM));

		String filename = createFile();

		writeFile(false, ResultsFileFormat.MALK, false, false, false, false, false, out, filename);

		MemoryPeakResults in = readFile(filename, false);

		// Change to nm and photon for the validation
		out.convertToUnits(DistanceUnit.NM, IntensityUnit.PHOTON, null);

		checkEqual(ResultsFileFormat.MALK, false, false, false, false, false, out, in);
	}

	@Test
	public void writeTextWithComputedPrecisionMatchesRead()
	{
		// Create without precision
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MemoryPeakResults results = createResults(rg, 200, false, false, false, false);
		// Ensure units are OK for computing precision
		CalibrationWriter cw = results.getCalibrationWriter();
		cw.setIntensityUnit(IntensityUnit.PHOTON);
		cw.setDistanceUnit(DistanceUnit.PIXEL);
		results.setCalibration(cw.getCalibration());

		String filename = createFile();

		TextFilePeakResults out = new TextFilePeakResults(filename, false, false, false, true);
		// Compute precision
		out.setComputePrecision(true);
		out.copySettings(results);
		out.begin();
		out.addAll(Arrays.asList(results.toArray()));
		out.end();

		MemoryPeakResults in = readFile(filename, false);

		checkEqual(ResultsFileFormat.TEXT, false, false, false, false, false, results, in);
	}

	@Test
	public void canReadTextIntoPreferredUnits()
	{
		canReadIntoPreferredUnits(ResultsFileFormat.TEXT);
	}

	@Test
	public void canReadBinaryIntoPreferredUnits()
	{
		canReadIntoPreferredUnits(ResultsFileFormat.BINARY);
	}

	@Test
	public void canReadTSFIntoPreferredUnits()
	{
		canReadIntoPreferredUnits(ResultsFileFormat.TSF);
	}

	private void canReadIntoPreferredUnits(ResultsFileFormat fileFormat)
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MemoryPeakResults out = createResults(rg, 200, false, false, false, false);

		// Output in nm and count
		CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
		cal.setDistanceUnit(DistanceUnit.NM);
		cal.setIntensityUnit(IntensityUnit.COUNT);
		if (fileFormat == ResultsFileFormat.TSF)
		{
			// For now just support using the native float TSF value
			cal.setNmPerPixel((float) cal.getNmPerPixel());
		}
		out.setCalibration(cal.getCalibration());

		String filename = createFile();

		writeFile(false, fileFormat, false, false, false, false, false, out, filename);

		MemoryPeakResults in = readFile(filename, false, false);

		// Change to preferred units
		out.convertToUnits(MemoryPeakResults.PREFERRED_DISTANCE_UNIT, MemoryPeakResults.PREFERRED_INTENSITY_UNIT,
				MemoryPeakResults.PREFERRED_ANGLE_UNIT);

		checkEqual(fileFormat, false, false, false, false, false, out, in);
	}

	@Test
	public void canReadTextAndSimplifyGaussian2DPSF()
	{
		canReadAndSimplifyGaussian2DPSF(ResultsFileFormat.TEXT);
	}

	@Test
	public void canReadBinaryAndSimplifyGaussian2DPSF()
	{
		canReadAndSimplifyGaussian2DPSF(ResultsFileFormat.BINARY);
	}

	@Test
	public void canReadTSFAndSimplifyGaussian2DPSF()
	{
		canReadAndSimplifyGaussian2DPSF(ResultsFileFormat.TSF);
	}

	private void canReadAndSimplifyGaussian2DPSF(ResultsFileFormat fileFormat)
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MemoryPeakResults out = createResults(rg, 1, false, false, false, false);

		CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
		cal.setDistanceUnit(MemoryPeakResults.PREFERRED_DISTANCE_UNIT);
		cal.setIntensityUnit(MemoryPeakResults.PREFERRED_INTENSITY_UNIT);
		cal.setAngleUnit(MemoryPeakResults.PREFERRED_ANGLE_UNIT);
		if (fileFormat == ResultsFileFormat.TSF)
		{
			// For now just support using the native float TSF value
			cal.setNmPerPixel((float) cal.getNmPerPixel());
		}
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

		String filename = createFile();

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
				float[] p = peakResult.getParameters();
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

	private void readWith2IsFasterThan1(boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean showPrecision, ResultsFileFormat f1, boolean useScanner1, ResultsFileFormat f2, boolean useScanner2,
			int loops)
	{
		TestSettings.assume(LogLevel.INFO, TestComplexity.HIGH);

		RandomGenerator rg = TestSettings.getRandomGenerator();
		MemoryPeakResults out = createResults(rg, 20000, showDeviations, showEndFrame, showId, showPrecision);
		String filename = createFile();

		writeFile(false, f1, showDeviations, showEndFrame, showId, showPrecision, false, out, filename);
		long time1 = getReadTime(filename, useScanner1, loops);

		writeFile(false, f2, showDeviations, showEndFrame, showId, showPrecision, false, out, filename);
		long time2 = getReadTime(filename, useScanner2, loops);

		if (useScanner1 != useScanner2)
			System.out.printf("%s (scan=%b) is %.2fx faster than %s (scan=%b)\n", f2, useScanner2,
					(double) time1 / time2, f1, useScanner1);
		else
			System.out.printf("%s is %.2fx faster than %s\n", f2, (double) time1 / time2, f1);
		TestAssert.assertTrue(time2 < time1, "%s (%d) is not faster than %s (%d)", f2, time2, f1, time1);
	}

	// -=-=-=-=-

	private void writeMatchesRead(boolean sequential, ResultsFileFormat fileFormat, boolean showDeviations,
			boolean showEndFrame, boolean showId, boolean showPrecision, boolean sort)
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MemoryPeakResults out = createResults(rg, 200, showDeviations, showEndFrame, showId, showPrecision);
		if (fileFormat == ResultsFileFormat.MALK)
		{
			CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
			cal.setDistanceUnit(DistanceUnit.NM);
			cal.setIntensityUnit(IntensityUnit.PHOTON);
			out.setCalibration(cal.getCalibration());
			out.setPSF(PSFHelper.create(PSFType.CUSTOM));
		}
		if (fileFormat == ResultsFileFormat.TSF)
		{
			CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
			// For now just support using the native float TSF datatype
			cal.setNmPerPixel((float) cal.getNmPerPixel());
			out.setCalibration(cal.getCalibration());
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

		String filename = createFile();

		writeFile(sequential, fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, out, filename);

		MemoryPeakResults in = readFile(filename, false);

		checkEqual(fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, out, in);
	}

	private void writeWithCombinationsMatchesRead(boolean sequential, ResultsFileFormat fileFormat, boolean sort)
	{
		for (boolean showDeviations : onOff)
		{
			for (boolean showEndFrame : onOff)
			{
				for (boolean showId : onOff)
				{
					for (boolean showPrecision : onOff)
					{
						if (count(showDeviations, showEndFrame, showId, showPrecision) < 2)
							continue;
						writeMatchesRead(sequential, fileFormat, showDeviations, showEndFrame, showId, showPrecision,
								sort);
					}
				}
			}
		}
	}

	private int count(boolean... flags)
	{
		int c = 0;
		for (boolean flag : flags)
			if (flag)
				c++;
		return c;
	}

	private void readWithScannerMatchesNonScanner(boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean showPrecision, boolean sort)
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		readWithScannerMatchesNonScanner(rg, showDeviations, showEndFrame, showId, showPrecision, sort);
	}

	private void readWithScannerMatchesNonScanner(RandomGenerator rg, boolean showDeviations, boolean showEndFrame,
			boolean showId, boolean showPrecision, boolean sort)
	{
		MemoryPeakResults out = createResults(rg, 1000, showDeviations, showEndFrame, showId, showPrecision);
		String filename = createFile();

		ResultsFileFormat fileFormat = ResultsFileFormat.TEXT;
		writeFile(false, fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, out, filename);

		MemoryPeakResults in = readFile(filename, false);
		MemoryPeakResults in2 = readFile(filename, true);

		checkEqual(fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, in, in2);
	}

	private void readWithScannerMatchesNonScannerWithCombinations(boolean sort)
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		for (boolean showDeviations : onOff)
		{
			for (boolean showEndFrame : onOff)
			{
				for (boolean showId : onOff)
				{
					for (boolean showPrecision : onOff)
					{
						if (count(showDeviations, showEndFrame, showId, showPrecision) < 2)
							continue;
						readWithScannerMatchesNonScanner(rg, showDeviations, showEndFrame, showId, showPrecision, sort);
					}
				}
			}
		}
	}

	private void checkEqual(ResultsFileFormat fileFormat, boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean showPrecision, boolean sort, MemoryPeakResults expectedResults, MemoryPeakResults actualResults)
			throws ArrayComparisonFailure
	{
		Assert.assertNotNull("Input results are null", actualResults);
		Assert.assertEquals("Size differ", expectedResults.size(), actualResults.size());

		final double delta = (fileFormat == ResultsFileFormat.BINARY) ? 0 : 1e-5f;

		PeakResult[] expected = expectedResults.toArray();
		PeakResult[] actual = actualResults.toArray();
		if (sort)
		{
			// Results should be sorted by time
			Arrays.sort(expected, new Comparator<PeakResult>()
			{
				@Override
				public int compare(PeakResult o1, PeakResult o2)
				{
					return o1.getFrame() - o2.getFrame();
				}
			});
		}

		// TSF requires the bias be subtracted
		//		double bias = expectedResults.getCalibration().getBias();

		for (int i = 0; i < actualResults.size(); i++)
		{
			PeakResult p1 = expected[i];
			PeakResult p2 = actual[i];

			TestAssert.assertEquals(p1.getFrame(), p2.getFrame(), "Peak mismatch @ [%d]", i);

			if (fileFormat == ResultsFileFormat.MALK)
			{
				TestAssert.assertEqualsRelative(p1.getXPosition(), p2.getXPosition(), delta, "X @ [%d]", i);
				TestAssert.assertEqualsRelative(p1.getYPosition(), p2.getYPosition(), delta, "Y @ [%d]", i);
				TestAssert.assertEqualsRelative(p1.getIntensity(), p2.getIntensity(), delta, "Signal @ " + i);
				continue;
			}

			TestAssert.assertEquals(p1.getOrigX(), p2.getOrigX(), "Orig X mismatch @ [%d]", i);
			TestAssert.assertEquals(p1.getOrigY(), p2.getOrigY(), "Orig Y mismatch @ [%d]", i);
			TestAssert.assertEqualsRelative(p1.getOrigValue(), p2.getOrigValue(), delta, "Orig value mismatch @ [%d]",
					i);
			TestAssert.assertEqualsRelative(p1.getError(), p2.getError(), 1e-6, "Error mismatch @ [%d]", i);
			TestAssert.assertEqualsRelative(p1.getNoise(), p2.getNoise(), delta, "Noise mismatch @ [%d]", i);
			TestAssert.assertEqualsRelative(p1.getMeanIntensity(), p2.getMeanIntensity(), delta,
					"Mean intensity mismatch @ [%d]", i);
			TestAssert.assertNotNull(p2.getParameters(), "Params is null @ [%d]", i);
			TestAssert.assertArrayEqualsRelative(p1.getParameters(), p2.getParameters(), delta,
					"Params mismatch @ [%d]", i);
			if (showDeviations)
			{
				Assert.assertNotNull(p2.getParameterDeviations());
				TestAssert.assertArrayEqualsRelative(p1.getParameterDeviations(), p2.getParameterDeviations(), delta,
						"Params StdDev mismatch @ [%d]", i);
			}
			if (showEndFrame)
			{
				TestAssert.assertEquals(p1.getEndFrame(), p2.getEndFrame(), "End frame mismatch @ [%d]", i);
			}
			if (showId)
			{
				TestAssert.assertEquals(p1.getId(), p2.getId(), "ID mismatch @ [%d]", i);
			}
			if (showPrecision)
			{
				TestAssert.assertEqualsRelative(p1.getPrecision(), p2.getPrecision(), delta,
						"Precision mismatch @ [%d]", i);
			}
		}

		// Check the header information
		Assert.assertEquals("Name", expectedResults.getName(), actualResults.getName());
		Assert.assertEquals("Configuration", expectedResults.getConfiguration(), actualResults.getConfiguration());

		Rectangle r1 = expectedResults.getBounds();
		Rectangle r2 = actualResults.getBounds();
		if (r1 != null)
		{
			Assert.assertNotNull("Bounds", r2);
			Assert.assertEquals("Bounds x", r1.x, r2.x);
			Assert.assertEquals("Bounds y", r1.y, r2.y);
			Assert.assertEquals("Bounds width", r1.width, r2.width);
			Assert.assertEquals("Bounds height", r1.height, r2.height);
		}
		else
		{
			Assert.assertNull("Bounds", r2);
		}

		Calibration c1 = expectedResults.getCalibration();
		Calibration c2 = actualResults.getCalibration();
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
			Assert.assertNotNull("Calibration", c2);
			Assert.assertTrue("Calibration", c1.equals(c2));
		}
		else
		{
			Assert.assertNull("Calibration", c2);
		}

		PSF p1 = expectedResults.getPSF();
		PSF p2 = actualResults.getPSF();
		if (p1 != null)
		{
			Assert.assertNotNull("PSF", p2);
			Assert.assertTrue("PSF", p1.equals(p2));
		}
		else
		{
			Assert.assertNull("PSF", p2);
		}
	}

	private MemoryPeakResults createResults(RandomGenerator rg, int i, boolean showDeviations, boolean showEndFrame,
			boolean showId, boolean showPrecision)
	{
		double bias = rg.nextDouble();

		boolean extended = showEndFrame || showId || showPrecision;

		MemoryPeakResults results = new MemoryPeakResults(PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D));
		while (i-- > 0)
		{
			int startFrame = rg.nextInt(i + 1);
			int origX = rg.nextInt(256);
			int origY = rg.nextInt(256);
			float origValue = rg.nextFloat();
			double error = rg.nextDouble();
			float noise = rg.nextFloat();
			float meanIntensity = rg.nextFloat();
			float[] params = createData(rg);
			float[] paramsStdDev = (showDeviations) ? createData(rg) : null;
			if (extended)
			{
				AttributePeakResult r = new AttributePeakResult(startFrame, origX, origY, origValue, error, noise,
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
		CalibrationWriter cal = new CalibrationWriter();
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

	private float[] createData(RandomGenerator rg)
	{
		return Gaussian2DPeakResultHelper.createTwoAxisAndAngleParams(rg.nextFloat(), rg.nextFloat(), rg.nextFloat(),
				rg.nextFloat(), rg.nextFloat(), rg.nextFloat(), rg.nextFloat(), rg.nextFloat());
	}

	private String createFile()
	{
		File file;
		try
		{
			file = File.createTempFile("test", null);
			file.deleteOnExit();
			String filename = file.getPath();
			return filename;
		}
		catch (IOException e)
		{
			Assert.fail("Cannot create temp files for IO testing");
		}
		return null; // Allow compilation but the assert will stop the code
	}

	private void writeFile(boolean sequential, ResultsFileFormat fileFormat, boolean showDeviations,
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
		{
			((FilePeakResults) out).setSortAfterEnd(sort);
		}
		out.begin();

		// TODO - option to test adding using:
		// add(peak);
		// addAll(PeakResult[])

		if (sequential)
		{
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
		}
		else
		{
			out.addAll(Arrays.asList(results.toArray()));
		}
		out.end();
	}

	private MemoryPeakResults readFile(String filename, boolean useScanner)
	{
		return readFile(filename, useScanner, true);
	}

	private MemoryPeakResults readFile(String filename, boolean useScanner, boolean rawResults)
	{
		PeakResultsReader reader = new PeakResultsReader(filename);
		reader.setUseScanner(useScanner);
		reader.setRawResults(rawResults);
		MemoryPeakResults in = reader.getResults();
		return in;
	}

	private long getReadTime(String filename, final boolean useScanner, final int loops)
	{
		// Initialise reading code
		readFile(filename, useScanner);

		long start = System.nanoTime();
		for (int i = loops; i-- > 0;)
			readFile(filename, useScanner);
		long time = System.nanoTime() - start;
		return time;
	}
}
