package gdsc.smlm.ij.plugins;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.core.utils.Random;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.results.IdPeakResult;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.TSFPeakResultsWriter;
import gdsc.smlm.tsf.TaggedSpotFile.FitMode;
import gdsc.smlm.tsf.TaggedSpotFile.FluorophoreType;
import gdsc.smlm.tsf.TaggedSpotFile.IntensityUnits;
import gdsc.smlm.tsf.TaggedSpotFile.LocationUnits;
import gdsc.smlm.tsf.TaggedSpotFile.Spot;
import gdsc.smlm.tsf.TaggedSpotFile.SpotList;
import ij.Macro;

/**
 * Test the ResultsManager functionality to load results from file when the file has options.
 */
public class ResultsManagerTest
{
	private gdsc.core.utils.Random rand = new Random();

	@Test
	public void writeTSFMatchesRead()
	{
		writeTSFMatchesRead(1, 1, 1, 1);
	}

	@Test
	public void writeTSFMatchesReadWithChannels()
	{
		writeTSFMatchesRead(3, 1, 1, 1);
	}

	@Test
	public void writeTSFMatchesReadWithSlices()
	{
		writeTSFMatchesRead(1, 3, 1, 1);
	}

	@Test
	public void writeTSFMatchesReadWithPositions()
	{
		writeTSFMatchesRead(1, 1, 3, 1);
	}

	@Test
	public void writeTSFMatchesReadWithTypes()
	{
		writeTSFMatchesRead(1, 1, 1, 3);
	}

	@Test
	public void writeTSFMatchesReadWithCombinations()
	{
		writeTSFMatchesRead(2, 2, 2, 2);
	}

	private void writeTSFMatchesRead(int channels, int slices, int positions, int types)
	{
		String filename = createFile();
		FileOutputStream out = null;
		try
		{
			out = new FileOutputStream(filename);
		}
		catch (Exception e)
		{
			closeOutput(out);
			e.printStackTrace();
			Assert.fail(e.getMessage());
		}

		// Write the offsets used in the TSF format
		try
		{
			DataOutputStream dos = new DataOutputStream(out);
			dos.writeInt(0);
			dos.writeLong(0);
		}
		catch (IOException e)
		{
			closeOutput(out);
			e.printStackTrace();
			Assert.fail(e.getMessage());
		}

		// Generate random spots
		int size = 1000;
		Spot[] spots = new Spot[size];
		for (int i = 1; i <= size; i++)
		{
			Spot.Builder builder = Spot.newBuilder();
			builder.setChannel(1 + rand.nextInt(channels));
			builder.setSlice(1 + rand.nextInt(slices));
			builder.setPos(1 + rand.nextInt(positions));
			builder.setFluorophoreType(rand.nextInt(1, types));

			builder.setMolecule(i); // This is a required field but is ignored when reading
			builder.setCluster(rand.nextInt(10));
			builder.setFrame(rand.nextInt(1, 100));
			builder.setXPosition(rand.nextInt(50));
			builder.setYPosition(rand.nextInt(50));
			builder.setBackground(rand.next());
			builder.setIntensity(rand.next());
			builder.setX(rand.next());
			builder.setY(rand.next());
			builder.setWidth(TSFPeakResultsWriter.SD_TO_FWHM_FACTOR * rand.next());

			Spot spot = builder.build();
			spots[i - 1] = spot;
			try
			{
				spot.writeDelimitedTo(out);
			}
			catch (IOException e)
			{
				closeOutput(out);
				e.printStackTrace();
				Assert.fail(e.getMessage());
			}
		}

		// Write the header
		// Get the offset to the SpotList message
		long offset = 0;
		try
		{
			// The offset is the amount to skip forward after reading the int 
			// magic number (4 bytes) and long offset (8 bytes)
			//out.flush();
			offset = out.getChannel().position() - 12;
		}
		catch (IOException e)
		{
			closeOutput(out);
			e.printStackTrace();
			Assert.fail(e.getMessage());
		}

		// Record the SpotList message
		SpotList.Builder builder = SpotList.newBuilder();

		builder.setApplicationId(1);
		builder.setNrSpots(size);
		builder.setLocationUnits(LocationUnits.PIXELS);
		builder.setIntensityUnits(IntensityUnits.COUNTS);
		builder.setFitMode(FitMode.ONEAXIS);

		builder.setNrChannels(channels);
		builder.setNrSlices(slices);
		builder.setNrPos(positions);
		for (int type = 1; type <= types; type++)
		{
			FluorophoreType.Builder typeBuilder = FluorophoreType.newBuilder();
			typeBuilder.setId(type);
			typeBuilder.setDescription("Type " + type);
			typeBuilder.setIsFiducial(rand.next() < 0.5f);
			builder.addFluorophoreTypes(typeBuilder.build());
		}

		SpotList spotList = builder.build();
		try
		{
			spotList.writeDelimitedTo(out);
		}
		catch (IOException e)
		{
			e.printStackTrace();
			Assert.fail(e.getMessage());
		}
		finally
		{
			closeOutput(out);
		}

		// Write the offset to the SpotList message into the offset position
		RandomAccessFile f = null;
		try
		{
			f = new RandomAccessFile(new File(filename), "rw");
			f.seek(4);
			f.writeLong(offset);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			Assert.fail(e.getMessage());
		}
		finally
		{
			if (f != null)
			{
				try
				{
					f.close();
				}
				catch (IOException e)
				{
				}
			}
		}

		// Read each combination 
		for (int channel = 1; channel <= channels; channel++)
			for (int slice = 1; slice <= slices; slice++)
				for (int position = 1; position <= positions; position++)
					for (int type = 1; type <= types; type++)
					{
						StringBuilder sb = new StringBuilder();
						sb.append(" channel=").append(channel);
						sb.append(" slice=").append(slice);
						sb.append(" position=").append(position);
						sb.append(" fluorophore_type=[").append(type).append(":Type ").append(type).append(']');
						// This is needed to trick the Macro class into returning the options for the thread
						Thread.currentThread().setName("Run$_");
						Macro.setOptions(sb.toString());

						ResultsManager.setInputFilename(filename);
						MemoryPeakResults in = ResultsManager.loadInputResults(ResultsManager.INPUT_FILE, false);
						checkEqual(spots, channel, slice, position, type, in);
					}
	}

	private void closeOutput(FileOutputStream out)
	{
		if (out == null)
			return;

		try
		{
			out.close();
		}
		catch (Exception e)
		{
			// Ignore exception
		}
		finally
		{
			out = null;
		}
	}

	private void checkEqual(Spot[] spots, int channel, int slice, int position, int type,
			MemoryPeakResults actualResults) throws ArrayComparisonFailure
	{
		Assert.assertNotNull("Input results are null", actualResults);

		MemoryPeakResults expectedResults = extract(spots, channel, slice, position, type);

		Assert.assertEquals("Size differ", expectedResults.size(), actualResults.size());

		final float delta = 0;

		List<PeakResult> expected = expectedResults.getResults();
		List<PeakResult> actual = actualResults.getResults();
		for (int i = 0; i < actualResults.size(); i++)
		{
			PeakResult p1 = expected.get(i);
			PeakResult p2 = actual.get(i);

			Assert.assertEquals("Peak mismatch @ " + i, p1.peak, p2.peak);

			Assert.assertEquals("Orig X mismatch @ " + i, p1.origX, p2.origX);
			Assert.assertEquals("Orig Y mismatch @ " + i, p1.origY, p2.origY);
			Assert.assertEquals("Orig value mismatch @ " + i, p1.origValue, p2.origValue, delta);
			Assert.assertEquals("Error mismatch @ " + i, p1.error, p2.error, 1e-6);
			Assert.assertEquals("Noise mismatch @ " + i, p1.noise, p2.noise, delta);
			Assert.assertNotNull("Params is null @ " + i, p2.params);

			Assert.assertEquals("Background mismatch @ " + i, p1.getBackground(), p2.getBackground(), delta);
			Assert.assertEquals("Signal mismatch @ " + i, p1.getSignal(), p2.getSignal(), delta);
			Assert.assertEquals("XPosition mismatch @ " + i, p1.getXPosition(), p2.getXPosition(), delta);
			Assert.assertEquals("YPosition mismatch @ " + i, p1.getYPosition(), p2.getYPosition(), delta);
			Assert.assertEquals("XSD mismatch @ " + i, p1.getXSD(), p2.getXSD(), 1e-6);
			Assert.assertEquals("YSD mismatch @ " + i, p1.getYSD(), p2.getYSD(), 1e-6);

			Assert.assertEquals("ID mismatch @ " + i, p1.getId(), p2.getId());
		}
	}

	private MemoryPeakResults extract(Spot[] spots, int channel, int slice, int position, int type)
	{
		MemoryPeakResults results = new MemoryPeakResults();
		for (Spot spot : spots)
		{
			if (spot.getChannel() == channel && spot.getSlice() == slice && spot.getPos() == position &&
					spot.getFluorophoreType() == type)
			{
				int id = spot.getCluster();
				int startFrame = spot.getFrame();
				int origX = spot.getXPosition();
				int origY = spot.getYPosition();
				float origValue = 0;
				double error = 0;
				float noise = 0;
				float[] params = new float[7];
				params[Gaussian2DFunction.BACKGROUND] = spot.getBackground();
				params[Gaussian2DFunction.SIGNAL] = spot.getIntensity();
				params[Gaussian2DFunction.X_POSITION] = spot.getX();
				params[Gaussian2DFunction.Y_POSITION] = spot.getY();
				params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = spot.getWidth() /
						TSFPeakResultsWriter.SD_TO_FWHM_FACTOR;
				float[] paramsStdDev = null;
				IdPeakResult peak = new IdPeakResult(startFrame, origX, origY, origValue, error, noise, params,
						paramsStdDev, id);
				results.add(peak);
			}
		}
		return results;
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
}
