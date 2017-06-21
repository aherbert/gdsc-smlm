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
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.SMLMSettings.*;
import gdsc.smlm.data.config.SMLMSettings.PSFType;
import gdsc.smlm.ij.results.ResultsFileFormat;
import gdsc.smlm.results.procedures.PeakResultProcedure;

public class PeakResultsReaderTest
{
	private gdsc.core.utils.Random rand = new Random();

	// TODO - Add tests to compare writing to a IJTablePeakResults, saving the TextPanel contents to file and then reading.

	// -=-=-=-=-

	@Test
	public void writeTextMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, false, false, false, false, false);
	}

	@Test
	public void writeSequentialTextMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.GDSC_TEXT, false, false, false, false, false);
	}

	@Test
	public void writeTextWithDeviationsMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, false, false, false, false);
	}

	@Test
	public void writeTextWithEndFrameMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, false, true, false, false, false);
	}

	@Test
	public void writeTextWithIdMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, false, false, true, false, false);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, true, false, false, false);
	}

	@Test
	public void writeTextWithDeviationsAndIdMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, false, true, false, false);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameAndIdMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, true, true, false, false);
	}

	// -=-=-=-=-

	@Test
	public void writeBinaryMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, false, false, false, false, true);
	}

	@Test
	public void writeSequentialBinaryMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.GDSC_BINARY, false, false, false, false, true);
	}

	@Test
	public void writeBinaryWithDeviationsMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, false, false, false, true);
	}

	@Test
	public void writeBinaryWithEndFrameMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, false, true, false, false, true);
	}

	@Test
	public void writeBinaryWithIdMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, false, false, true, false, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, true, false, false, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndIdMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, false, true, false, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameAndIdMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, true, true, false, true);
	}

	// -=-=-=-=-

	@Test
	public void writeTextWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, false, false, false, true, false);
	}

	@Test
	public void writeSequentialTextWithSortMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.GDSC_TEXT, false, false, false, true, false);
	}

	@Test
	public void writeTextWithDeviationsWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, false, false, true, false);
	}

	@Test
	public void writeTextWithEndFrameWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, false, true, false, true, false);
	}

	@Test
	public void writeTextWithIdWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, false, false, true, true, false);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameWithSortMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, true, false, true, false);
	}

	@Test
	public void writeTextWithDeviationsAndIdWithSortMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, false, true, true, false);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameAndIdWithSortMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_TEXT, true, true, true, true, false);
	}

	// -=-=-=-=-

	@Test
	public void writeBinaryWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, false, false, false, true, true);
	}

	@Test
	public void writeSequentialBinaryWithSortMatchesRead()
	{
		writeMatchesRead(true, ResultsFileFormat.GDSC_BINARY, false, false, false, true, true);
	}

	@Test
	public void writeBinaryWithDeviationsWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, false, false, true, true);
	}

	@Test
	public void writeBinaryWithEndFrameWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, false, true, false, true, true);
	}

	@Test
	public void writeBinaryWithIdWithSortMatchesRead()
	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, false, false, true, true, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameWithSortMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, true, false, true, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndIdWithSortMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, false, true, true, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameAndIdWithSortMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.GDSC_BINARY, true, true, true, true, true);
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
	public void writeTSFWithDeviationsAndEndFrameMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.TSF, true, true, false, false, false);
	}

	@Test
	public void writeTSFWithDeviationsAndIdMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.TSF, true, false, true, false, false);
	}

	@Test
	public void writeTSFWithDeviationsAndEndFrameAndIdMatchesRead()

	{
		writeMatchesRead(false, ResultsFileFormat.TSF, true, true, true, false, false);
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
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrame()

	{
		readWithScannerMatchesNonScanner(true, true, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithId()

	{
		readWithScannerMatchesNonScanner(true, false, true, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrameWithId()

	{
		readWithScannerMatchesNonScanner(true, true, true, false, false);
	}

	// -=-=-=-=-

	@Test
	public void readWithScannerMatchesNonScannerWithSort()
	{
		readWithScannerMatchesNonScanner(false, false, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithSort()

	{
		readWithScannerMatchesNonScanner(true, false, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithEndFrameWithSort()

	{
		readWithScannerMatchesNonScanner(false, true, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithIdWithSort()

	{
		readWithScannerMatchesNonScanner(false, false, true, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrameWithSort()

	{
		readWithScannerMatchesNonScanner(true, true, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithIdWithSort()

	{
		readWithScannerMatchesNonScanner(true, false, true, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrameWithIdWithSort()

	{
		readWithScannerMatchesNonScanner(true, true, true, true, false);
	}

	// -=-=-=-=-

	@Test
	public void readTextWithNonScannerIsFasterThanScanner()
	{
		readWith2IsFasterThan1(false, false, false, ResultsFileFormat.GDSC_TEXT, true, ResultsFileFormat.GDSC_TEXT,
				false, 1);
	}

	@Test
	public void readTextWithNonScannerIsFasterThanScannerWithDeviationsWithEndFrameWithId()
	{
		readWith2IsFasterThan1(true, true, true, ResultsFileFormat.GDSC_TEXT, true, ResultsFileFormat.GDSC_TEXT, false,
				1);
	}

	@Test
	public void readWithMALKIsFasterThanText()
	{
		readWith2IsFasterThan1(false, false, false, ResultsFileFormat.GDSC_TEXT, false, ResultsFileFormat.MALK, false,
				2);
	}

	@Test
	public void readWithBinaryIsFasterThanText()
	{
		readWith2IsFasterThan1(false, false, false, ResultsFileFormat.GDSC_TEXT, false, ResultsFileFormat.GDSC_BINARY,
				false, 2);
	}

	@Test
	public void readWithBinaryIsFasterThanTSF()
	{
		readWith2IsFasterThan1(false, false, false, ResultsFileFormat.TSF, false, ResultsFileFormat.GDSC_BINARY, false,
				20);
	}

	private void readWith2IsFasterThan1(boolean showDeviations, boolean showEndFrame, boolean showId,
			ResultsFileFormat f1, boolean useScanner1, ResultsFileFormat f2, boolean useScanner2, int loops)
	{
		MemoryPeakResults out = createResults(20000, showDeviations, showEndFrame, showId);
		String filename = createFile();

		writeFile(false, f1, showDeviations, showEndFrame, showId, false, out, filename);
		long time1 = getReadTime(filename, useScanner1, loops);

		writeFile(false, f2, showDeviations, showEndFrame, showId, false, out, filename);
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
			boolean showEndFrame, boolean showId, boolean sort, boolean binary)
	{
		MemoryPeakResults out = createResults(200, showDeviations, showEndFrame, showId);
		if (fileFormat == ResultsFileFormat.MALK)
		{
			Calibration cal = out.getCalibration();
			cal.setDistanceUnit(DistanceUnit.NM);
			cal.setIntensityUnit(IntensityUnit.PHOTON);
			out.setPSF(PSFHelper.create(PSFType.CUSTOM));			
		}
		
		String filename = createFile();

		writeFile(sequential, fileFormat, showDeviations, showEndFrame, showId, sort, out, filename);

		MemoryPeakResults in = readFile(filename, false);

		checkEqual(fileFormat, showDeviations, showEndFrame, showId, sort, binary, out, in);
	}

	private void readWithScannerMatchesNonScanner(boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean sort, boolean binary)
	{
		MemoryPeakResults out = createResults(1000, showDeviations, showEndFrame, showId);
		String filename = createFile();

		ResultsFileFormat fileFormat = ResultsFileFormat.GDSC_TEXT;
		writeFile(false, fileFormat, showDeviations, showEndFrame, showId, sort, out, filename);

		MemoryPeakResults in = readFile(filename, false);
		MemoryPeakResults in2 = readFile(filename, true);

		checkEqual(fileFormat, showDeviations, showEndFrame, showId, sort, binary, in, in2);
	}

	private void checkEqual(ResultsFileFormat fileFormat, boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean sort, boolean binary, MemoryPeakResults expectedResults, MemoryPeakResults actualResults)
			throws ArrayComparisonFailure
	{
		Assert.assertNotNull("Input results are null", actualResults);
		Assert.assertEquals("Size differ", expectedResults.size(), actualResults.size());

		final float delta = (binary) ? 0 : 1e-6f;

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

			Assert.assertEquals("Orig X mismatch @ " + i, p1.origX, p2.origX);
			Assert.assertEquals("Orig Y mismatch @ " + i, p1.origY, p2.origY);
			Assert.assertEquals("Orig value mismatch @ " + i, p1.origValue, p2.origValue, delta);
			Assert.assertEquals("Error mismatch @ " + i, p1.error, p2.error, 1e-6);
			Assert.assertEquals("Noise mismatch @ " + i, p1.noise, p2.noise, delta);
			Assert.assertNotNull("Params is null @ " + i, p2.params);
			Assert.assertArrayEquals("Params mismatch @ " + i, p1.params, p2.params, delta);
			if (showDeviations)
			{
				Assert.assertNotNull(p2.paramsStdDev);
				Assert.assertArrayEquals("Params StdDev mismatch @ " + i, p1.paramsStdDev, p2.paramsStdDev, delta);
			}
			if (showEndFrame)
			{
				Assert.assertEquals("End frame mismatch @ " + i, p1.getEndFrame(), p2.getEndFrame());
			}
			if (showId)
			{
				Assert.assertEquals("ID mismatch @ " + i, p1.getId(), p2.getId());
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
		if (c1 != null)
		{
			Assert.assertNotNull("Calibration", c2);
			Assert.assertEquals("Calibration nmPerPixel", c1.getNmPerPixel(), c2.getNmPerPixel(), 1e-6);
			Assert.assertEquals("Calibration gain", c1.getGain(), c2.getGain(), 1e-6);
			Assert.assertEquals("Calibration exposureTime", c1.getExposureTime(), c2.getExposureTime(), 1e-6);
			Assert.assertEquals("Calibration readNoise", c1.getReadNoise(), c2.getReadNoise(), 1e-6);
			Assert.assertEquals("Calibration bias", c1.getBias(), c2.getBias(), 1e-6);
			Assert.assertEquals("Calibration amplification", c1.getAmplification(), c2.getAmplification(), 1e-6);
			Assert.assertEquals("Calibration CameraType", c1.getCameraType(), c2.getCameraType());
			Assert.assertEquals("Calibration DistanceUnit", c1.getDistanceUnit(), c2.getDistanceUnit());
			Assert.assertEquals("Calibration IntensityUnit", c1.getIntensityUnit(), c2.getIntensityUnit());
			Assert.assertEquals("Calibration AngleUnit", c1.getAngleUnit(), c2.getAngleUnit());
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

	private MemoryPeakResults createResults(int i, boolean showDeviations, boolean showEndFrame, boolean showId)
	{
		double bias = rand.next();

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
			params[PeakResult.BACKGROUND] += bias;
			float[] paramsStdDev = (showDeviations) ? createData() : null;
			if (showEndFrame || showId)
				results.add(new ExtendedPeakResult(startFrame, origX, origY, origValue, error, noise, params,
						paramsStdDev, startFrame + (int) (10 * rand.next()), i + 1));
			else
				results.add(startFrame, origX, origY, origValue, error, noise, params, paramsStdDev);
		}
		results.setName(Float.toString(rand.next()) + Float.toString(rand.next()));
		results.setConfiguration(Float.toString(rand.next()) + Float.toString(rand.next()));
		results.setBounds(new Rectangle((int) (10 * rand.next()), (int) (10 * rand.next()), (int) (100 * rand.next()),
				(int) (100 * rand.next())));
		Calibration cal = new Calibration();
		cal.setNmPerPixel(rand.next());
		cal.setGain(rand.next());
		cal.setExposureTime(rand.next());
		cal.setReadNoise(rand.next());
		cal.setBias(bias);
		cal.setAmplification(rand.next());
		// Subtract 1 to avoid the additional UNRECOGNISED enum value
		cal.setCameraType(CameraType.values()[rand.nextInt(CameraType.values().length-1)]);
		cal.setDistanceUnit(DistanceUnit.values()[rand.nextInt(DistanceUnit.values().length-1)]);
		cal.setIntensityUnit(IntensityUnit.values()[rand.nextInt(IntensityUnit.values().length-1)]);
		cal.setAngleUnit(AngleUnit.values()[rand.nextInt(AngleUnit.values().length-1)]);
		results.setCalibration(cal);
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
			boolean showEndFrame, boolean showId, boolean sort, MemoryPeakResults results, String filename)
	{
		final PeakResults out;
		switch (fileFormat)
		{
			case GDSC_BINARY:
				out = new BinaryFilePeakResults(filename, showDeviations, showEndFrame, showId);
				break;
			case GDSC_TEXT:
				out = new TextFilePeakResults(filename, showDeviations, showEndFrame, showId);
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
					out.add(peak.getFrame(), peak.origX, peak.origY, peak.origValue, peak.error, peak.noise,
							peak.params, peak.paramsStdDev);
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
