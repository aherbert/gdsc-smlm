package gdsc.smlm.results;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.Random;
import gdsc.smlm.ij.results.ResultsFileFormat;

public class PeakResultsReaderTest
{
	private gdsc.core.utils.Random rand = new Random();

	// TODO - Add tests to compare writing to a IJTablePeakResults, saving the TextPanel contents to file and then reading.

	// -=-=-=-=-

	@Test
	public void writeTextMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, false, false, false, false);
	}

	@Test
	public void writeTextWithDeviationsMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, false, false, false);
	}

	@Test
	public void writeTextWithEndFrameMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, false, true, false, false);
	}

	@Test
	public void writeTextWithIdMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, false, false, true, false);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, true, false, false);
	}

	@Test
	public void writeTextWithDeviationsAndIdMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, false, true, false);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameAndIdMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, true, true, false);
	}

	// -=-=-=-=-

	@Test
	public void writeBinaryMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, false, false, false, false);
	}

	@Test
	public void writeBinaryWithDeviationsMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, false, false, false);
	}

	@Test
	public void writeBinaryWithEndFrameMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, false, true, false, false);
	}

	@Test
	public void writeBinaryWithIdMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, false, false, true, false);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, true, false, false);
	}

	@Test
	public void writeBinaryWithDeviationsAndIdMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, false, true, false);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameAndIdMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, true, true, false);
	}

	// -=-=-=-=-

	@Test
	public void writeTextWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, false, false, false, true);
	}

	@Test
	public void writeTextWithDeviationsWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, false, false, true);
	}

	@Test
	public void writeTextWithEndFrameWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, false, true, false, true);
	}

	@Test
	public void writeTextWithIdWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, false, false, true, true);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameWithSortMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, true, false, true);
	}

	@Test
	public void writeTextWithDeviationsAndIdWithSortMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, false, true, true);
	}

	@Test
	public void writeTextWithDeviationsAndEndFrameAndIdWithSortMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_TEXT, true, true, true, true);
	}

	// -=-=-=-=-

	@Test
	public void writeBinaryWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, false, false, false, true);
	}

	@Test
	public void writeBinaryWithDeviationsWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, false, false, true);
	}

	@Test
	public void writeBinaryWithEndFrameWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, false, true, false, true);
	}

	@Test
	public void writeBinaryWithIdWithSortMatchesRead()
	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, false, false, true, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameWithSortMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, true, false, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndIdWithSortMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, false, true, true);
	}

	@Test
	public void writeBinaryWithDeviationsAndEndFrameAndIdWithSortMatchesRead()

	{
		writeMatchesRead(ResultsFileFormat.GDSC_BINARY, true, true, true, true);
	}

	// -=-=-=-=-

	@Test
	public void readWithScannerMatchesNonScanner()
	{
		readWithScannerMatchesNonScanner(false, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviations()

	{
		readWithScannerMatchesNonScanner(true, false, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithEndFrame()

	{
		readWithScannerMatchesNonScanner(false, true, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithId()

	{
		readWithScannerMatchesNonScanner(false, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrame()

	{
		readWithScannerMatchesNonScanner(true, true, false, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithId()

	{
		readWithScannerMatchesNonScanner(true, false, true, false);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrameWithId()

	{
		readWithScannerMatchesNonScanner(true, true, true, false);
	}

	// -=-=-=-=-

	@Test
	public void readWithScannerMatchesNonScannerWithSort()
	{
		readWithScannerMatchesNonScanner(false, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithSort()

	{
		readWithScannerMatchesNonScanner(true, false, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithEndFrameWithSort()

	{
		readWithScannerMatchesNonScanner(false, true, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithIdWithSort()

	{
		readWithScannerMatchesNonScanner(false, false, true, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrameWithSort()

	{
		readWithScannerMatchesNonScanner(true, true, false, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithIdWithSort()

	{
		readWithScannerMatchesNonScanner(true, false, true, true);
	}

	@Test
	public void readWithScannerMatchesNonScannerWithDeviationsWithEndFrameWithIdWithSort()

	{
		readWithScannerMatchesNonScanner(true, true, true, true);
	}

	// -=-=-=-=-

	@Test
	public void readWithScannerIsFasterThanNonScanner()
	{
		MemoryPeakResults out = createResults(20000, true, true, true);
		String filename = createFile();

		writeFile(ResultsFileFormat.GDSC_TEXT, true, true, true, false, out, filename);

		final int loops = 20;

		long time = getReadTime(filename, false, loops);
		long time2 = getReadTime(filename, true, loops);

		Assert.assertTrue(String.format("Scanner is slower: %d > %d", time2, time), time2 < time);
	}

	private void writeMatchesRead(ResultsFileFormat fileFormat, boolean showDeviations, boolean showEndFrame,
			boolean showId, boolean sort)
	{
		MemoryPeakResults out = createResults(200, showDeviations, showEndFrame, showId);
		String filename = createFile();

		writeFile(fileFormat, showDeviations, showEndFrame, showId, sort, out, filename);

		MemoryPeakResults in = readFile(filename, false);

		checkEqual(showDeviations, showEndFrame, showId, sort, out, in);
	}

	private void readWithScannerMatchesNonScanner(boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean sort)
	{
		MemoryPeakResults out = createResults(1000, showDeviations, showEndFrame, showId);
		String filename = createFile();

		writeFile(ResultsFileFormat.GDSC_TEXT, showDeviations, showEndFrame, showId, sort, out, filename);

		MemoryPeakResults in = readFile(filename, false);
		MemoryPeakResults in2 = readFile(filename, true);

		checkEqual(showDeviations, showEndFrame, showId, sort, in, in2);
	}

	private void checkEqual(boolean showDeviations, boolean showEndFrame, boolean showId, boolean sort,
			MemoryPeakResults expectedResults, MemoryPeakResults actualResults) throws ArrayComparisonFailure
	{
		Assert.assertNotNull("Input results are null", actualResults);
		Assert.assertEquals("Size differ", expectedResults.size(), actualResults.size());

		List<PeakResult> expected = expectedResults.getResults();
		List<PeakResult> actual = actualResults.getResults();
		if (sort)
		{
			// Results should be sorted by time
			Collections.sort(expected, new Comparator<PeakResult>()
			{
				public int compare(PeakResult o1, PeakResult o2)
				{
					return o1.peak - o2.peak;
				}
			});
		}
		for (int i = 0; i < actualResults.size(); i++)
		{
			PeakResult p1 = expected.get(i);
			PeakResult p2 = actual.get(i);

			Assert.assertEquals("Peak mismatch @ " + i, p1.peak, p2.peak);
			Assert.assertEquals("Orig X mismatch @ " + i, p1.origX, p2.origX);
			Assert.assertEquals("Orig Y mismatch @ " + i, p1.origY, p2.origY);
			Assert.assertEquals("Orig value mismatch @ " + i, p1.origValue, p2.origValue, 1e-6f);
			Assert.assertEquals("Error mismatch @ " + i, p1.error, p2.error, 1e-6);
			Assert.assertEquals("Noise mismatch @ " + i, p1.noise, p2.noise, 1e-6f);
			Assert.assertNotNull("Params is null @ " + i, p2.params);
			Assert.assertArrayEquals("Params mismatch @ " + i, p1.params, p2.params, 1e-6f);
			if (showDeviations)
			{
				Assert.assertNotNull(p2.paramsStdDev);
				Assert.assertArrayEquals("Params StdDev mismatch @ " + i, p1.paramsStdDev, p2.paramsStdDev, 1e-6f);
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
	}

	private MemoryPeakResults createResults(int i, boolean showDeviations, boolean showEndFrame, boolean showId)
	{
		MemoryPeakResults results = new MemoryPeakResults();
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
			if (showEndFrame || showId)
				results.add(new ExtendedPeakResult(startFrame, origX, origY, origValue, error, noise, params,
						paramsStdDev, startFrame + (int) (10 * rand.next()), i + 1));
			else
				results.add(startFrame, origX, origY, origValue, error, noise, params, paramsStdDev);
		}
		return results;
	}

	private float[] createData()
	{
		float[] data = new float[7];
		for (int i = 0; i < data.length; i++)
			data[i] = rand.next();
		return data;
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

	private void writeFile(ResultsFileFormat fileFormat, boolean showDeviations, boolean showEndFrame, boolean showId,
			boolean sort, MemoryPeakResults results, String filename)
	{
		PeakResults out;
		switch (fileFormat)
		{
			case GDSC_BINARY:
				out = new BinaryFilePeakResults(filename, showDeviations, showEndFrame, showId);
				break;
			case GDSC_TEXT:
				out = new FilePeakResults(filename, showDeviations, showEndFrame, showId);
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

		if (sort && out instanceof FilePeakResults)
		{
			((FilePeakResults) out).setSortAfterEnd(sort);
		}
		out.begin();
		
		// TODO - option to test adding using:
		// add(peak, origX, origY, origValue, chiSquared, noise, params, paramsStdDev);
		
		out.addAll(results.getResults());
		out.end();
	}

	private MemoryPeakResults readFile(String filename, boolean useScanner)
	{
		PeakResultsReader reader = new PeakResultsReader(filename);
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
