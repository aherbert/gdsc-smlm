package gdsc.smlm.results;

import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.Random;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import gdsc.smlm.results.procedures.PeakResultProcedure;

public class PeakResultsReaderTest
{
	private gdsc.core.utils.Random rand = new Random();
	static final boolean[] onOff = new boolean[] { true, false };

	// TODO - Add tests to compare writing to a IJTablePeakResults, saving the TextPanel contents to file and then reading.

	// -=-=-=-=-

	@Test
	public void writeTextMatchesRead()
	{
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
		writeMatchesRead(false, ResultsFileFormat.TEXT, true, false, false, false, false);
	}

	@Test
	public void writeTextWithEndFrameMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, true, false, false, false);
	}

	@Test
	public void writeTextWithIdMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, true, false, false);
	}

	@Test
	public void writeTextWithPrecisionMatchesRead()
	{
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
		writeMatchesRead(false, ResultsFileFormat.BINARY, true, false, false, false, false);
	}

	@Test
	public void writeBinaryWithEndFrameMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, true, false, false, false);
	}

	@Test
	public void writeBinaryWithIdMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, true, false, false);
	}

	@Test
	public void writeBinaryWithPrecisionMatchesRead()
	{
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
		writeMatchesRead(false, ResultsFileFormat.TEXT, true, false, false, false, true);
	}

	@Test
	public void writeTextWithEndFrameWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, true, false, false, true);
	}

	@Test
	public void writeTextWithIdWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.TEXT, false, false, true, false, true);
	}

	@Test
	public void writeTextWithPrecisionWithSortMatchesRead()
	{
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
		writeMatchesRead(false, ResultsFileFormat.BINARY, true, false, false, false, true);
	}

	@Test
	public void writeBinaryWithEndFrameWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, true, false, false, true);
	}

	@Test
	public void writeBinaryWithIdWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.BINARY, false, false, true, false, true);
	}

	@Test
	public void writeBinaryWithPrecisionWithSortMatchesRead()
	{
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
		writeMatchesRead(false, ResultsFileFormat.TSF, true, false, false, false, false);
	}

	@Test
	public void writeTSFWithEndFrameMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.TSF, false, true, false, false, false);
	}

	@Test
	public void writeTSFWithIdMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.TSF, false, false, true, false, false);
	}

	@Test
	public void writeTSFWithPrecisionMatchesRead()
	{
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
		readWithScannerMatchesNonScanner(false, false, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviations()

	{
		readWithScannerMatchesNonScanner(true, false, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithEndFrame()

	{
		readWithScannerMatchesNonScanner(false, true, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithId()

	{
		readWithScannerMatchesNonScanner(false, false, true, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithPrecision()

	{
		readWithScannerMatchesNonScanner(false, false, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithCombinations()
	{
		readWithScannerMatchesNonScannerWithCombinations(false);
	}

	// -=-=-=-=-

	@Test
	public void readWithScannerMatchesNonScannerWithSort()
	{
		readWithScannerMatchesNonScanner(false, false, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithSort()

	{
		readWithScannerMatchesNonScanner(true, false, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithEndFrameWithSort()

	{
		readWithScannerMatchesNonScanner(false, true, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithIdWithSort()

	{
		readWithScannerMatchesNonScanner(false, false, true, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithPrecisionWithSort()

	{
		readWithScannerMatchesNonScanner(false, false, false, true, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithCombinationsWithSort()
	{
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
		MemoryPeakResults out = createResults(200, false, false, false, false);

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
		MemoryPeakResults out = createResults(200, false, false, false, false);

		// Output in nm and count
		CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
		cal.setDistanceUnit(DistanceUnit.NM);
		cal.setIntensityUnit(IntensityUnit.COUNT);
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
		MemoryPeakResults out = createResults(200, false, false, false, false);

		CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
		cal.setDistanceUnit(MemoryPeakResults.PREFERRED_DISTANCE_UNIT);
		cal.setIntensityUnit(MemoryPeakResults.PREFERRED_INTENSITY_UNIT);
		cal.setAngleUnit(MemoryPeakResults.PREFERRED_ANGLE_UNIT);
		out.setCalibration(cal.getCalibration());

		// Remove angle
		final int ia = PSFHelper.getGaussian2DAngleIndex(out.getPSF());
		out.forEach(new PeakResultProcedure()
		{
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
		MemoryPeakResults out = createResults(20000, showDeviations, showEndFrame, showId, showPrecision);
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
		Assert.assertTrue(String.format(f1 + " is slower (%d > %d) than " + f2, time2, time1), time2 < time1);
	}

	// -=-=-=-=-

	private void writeMatchesRead(boolean sequential, ResultsFileFormat fileFormat, boolean showDeviations,
			boolean showEndFrame, boolean showId, boolean showPrecision, boolean sort)
	{
		MemoryPeakResults out = createResults(200, showDeviations, showEndFrame, showId, showPrecision);
		if (fileFormat == ResultsFileFormat.MALK)
		{
			CalibrationWriter cal = new CalibrationWriter(out.getCalibration());
			cal.setDistanceUnit(DistanceUnit.NM);
			cal.setIntensityUnit(IntensityUnit.PHOTON);
			out.setCalibration(cal.getCalibration());
			out.setPSF(PSFHelper.create(PSFType.CUSTOM));
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
		MemoryPeakResults out = createResults(1000, showDeviations, showEndFrame, showId, showPrecision);
		String filename = createFile();

		ResultsFileFormat fileFormat = ResultsFileFormat.TEXT;
		writeFile(false, fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, out, filename);

		MemoryPeakResults in = readFile(filename, false);
		MemoryPeakResults in2 = readFile(filename, true);

		checkEqual(fileFormat, showDeviations, showEndFrame, showId, showPrecision, sort, in, in2);
	}

	private void readWithScannerMatchesNonScannerWithCombinations(boolean sort)
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
						readWithScannerMatchesNonScanner(showDeviations, showEndFrame, showId, showPrecision, sort);
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

		final float delta = (fileFormat == ResultsFileFormat.BINARY) ? 0 : 1e-6f;

		PeakResult[] expected = expectedResults.toArray();
		PeakResult[] actual = actualResults.toArray();
		if (sort)
		{
			// Results should be sorted by time
			Arrays.sort(expected, new Comparator<PeakResult>()
			{
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

			Assert.assertEquals("Peak mismatch @ " + i, p1.getFrame(), p2.getFrame());

			if (fileFormat == ResultsFileFormat.MALK)
			{
				Assert.assertEquals("X @ " + i, p1.getXPosition(), p2.getXPosition(), delta);
				Assert.assertEquals("Y @ " + i, p1.getYPosition(), p2.getYPosition(), delta);
				Assert.assertEquals("Signal @ " + i, p1.getSignal(), p2.getSignal(), delta);
				continue;
			}

			Assert.assertEquals("Orig X mismatch @ " + i, p1.getOrigX(), p2.getOrigX());
			Assert.assertEquals("Orig Y mismatch @ " + i, p1.getOrigY(), p2.getOrigY());
			Assert.assertEquals("Orig value mismatch @ " + i, p1.getOrigValue(), p2.getOrigValue(), delta);
			Assert.assertEquals("Error mismatch @ " + i, p1.getError(), p2.getError(), 1e-6);
			Assert.assertEquals("Noise mismatch @ " + i, p1.getNoise(), p2.getNoise(), delta);
			Assert.assertNotNull("Params is null @ " + i, p2.getParameters());
			Assert.assertArrayEquals("Params mismatch @ " + i, p1.getParameters(), p2.getParameters(), delta);
			if (showDeviations)
			{
				Assert.assertNotNull(p2.getParameterDeviations());
				Assert.assertArrayEquals("Params StdDev mismatch @ " + i, p1.getParameterDeviations(),
						p2.getParameterDeviations(), delta);
			}
			if (showEndFrame)
			{
				Assert.assertEquals("End frame mismatch @ " + i, p1.getEndFrame(), p2.getEndFrame());
			}
			if (showId)
			{
				Assert.assertEquals("ID mismatch @ " + i, p1.getId(), p2.getId());
			}
			if (showPrecision)
			{
				Assert.assertEquals("Precision mismatch @ " + i, p1.getPrecision(), p2.getPrecision(), delta);
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
		//		try
		//		{
		//			Printer printer = JsonFormat.printer().omittingInsignificantWhitespace()
		//			//.includingDefaultValueFields()
		//			;
		//			System.out.println(printer.print(c1));
		//			System.out.println(printer.print(c2));
		//		}
		//		catch (InvalidProtocolBufferException e)
		//		{
		//			// This shouldn't happen so throw it
		//		}

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

	private MemoryPeakResults createResults(int i, boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean showPrecision)
	{
		double bias = rand.next();

		boolean extended = showEndFrame || showId || showPrecision;

		MemoryPeakResults results = new MemoryPeakResults(PSFHelper.create(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D));
		while (i-- > 0)
		{
			int startFrame = (int) (i * rand.next());
			int origX = (int) (i * rand.next());
			int origY = (int) (i * rand.next());
			float origValue = rand.next();
			double error = rand.next();
			float noise = rand.next();
			float[] params = createData();
			float[] paramsStdDev = (showDeviations) ? createData() : null;
			if (extended)
			{
				AttributePeakResult r = new AttributePeakResult(startFrame, origX, origY, origValue, error, noise,
						params, paramsStdDev);
				if (showEndFrame)
					r.setEndFrame(startFrame + (int) (10 * rand.next()));
				if (showId)
					r.setId(i + 1);
				if (showPrecision)
					r.setPrecision(rand.next());
				results.add(r);
			}
			else
				results.add(startFrame, origX, origY, origValue, error, noise, params, paramsStdDev);
		}
		results.setName(Float.toString(rand.next()) + Float.toString(rand.next()));
		results.setConfiguration(Float.toString(rand.next()) + Float.toString(rand.next()));
		results.setBounds(new Rectangle((int) (10 * rand.next()), (int) (10 * rand.next()), (int) (100 * rand.next()),
				(int) (100 * rand.next())));
		CalibrationWriter cal = new CalibrationWriter();
		cal.setNmPerPixel(rand.next());
		cal.setCountPerPhoton(rand.next());
		cal.setExposureTime(rand.next());
		cal.setReadNoise(rand.next());
		cal.setBias(bias);
		cal.setQuantumEfficiency(rand.next());
		// Subtract 1 to avoid the additional UNRECOGNISED enum value
		cal.setCameraType(CameraType.values()[rand.nextInt(CameraType.values().length - 1)]);
		cal.setDistanceUnit(DistanceUnit.values()[rand.nextInt(DistanceUnit.values().length - 1)]);
		cal.setIntensityUnit(IntensityUnit.values()[rand.nextInt(IntensityUnit.values().length - 1)]);
		cal.setAngleUnit(AngleUnit.values()[rand.nextInt(AngleUnit.values().length - 1)]);
		results.setCalibration(cal.getCalibration());
		return results;
	}

	private float[] createData()
	{
		return Gaussian2DPeakResultHelper.createTwoAxisAndAngleParams(rand.next(), rand.next(), rand.next(),
				rand.next(), rand.next(), rand.next(), rand.next(), rand.next());
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
				public void execute(PeakResult peak)
				{
					out.add(peak.getFrame(), peak.getOrigX(), peak.getOrigY(), peak.getOrigValue(), peak.getError(),
							peak.getNoise(), peak.getParameters(), peak.getParameterDeviations());
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
