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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.IdPeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.FitMode;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.FluorophoreType;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.IntensityUnits;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.LocationUnits;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.Spot;
import uk.ac.sussex.gdsc.smlm.tsf.TSFProtos.SpotList;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

import ij.Macro;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * Test the ResultsManager functionality to load results from file when the file has options.
 */
@SuppressWarnings({"javadoc"})
public class ResultsManagerTest {
  @SeededTest
  public void writeTSFMatchesRead(RandomSeed seed) {
    // This is redundant
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    writeTSFMatchesRead(seed, 1, 1, 1, 1);
  }

  @SeededTest
  public void writeTSFMatchesReadWithChannels(RandomSeed seed) {
    // Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    writeTSFMatchesRead(seed, 2, 1, 1, 1);
  }

  @SeededTest
  public void writeTSFMatchesReadWithSlices(RandomSeed seed) {
    // Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    writeTSFMatchesRead(seed, 1, 2, 1, 1);
  }

  @SeededTest
  public void writeTSFMatchesReadWithPositions(RandomSeed seed) {
    // Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    writeTSFMatchesRead(seed, 1, 1, 2, 1);
  }

  @SeededTest
  public void writeTSFMatchesReadWithTypes(RandomSeed seed) {
    // Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    writeTSFMatchesRead(seed, 1, 1, 1, 2);
  }

  @SeededTest
  public void writeTSFMatchesReadWithCombinations(RandomSeed seed) {
    // This takes longer
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    writeTSFMatchesRead(seed, 2, 2, 2, 2);
  }

  @SuppressWarnings("null")
  private static void writeTSFMatchesRead(RandomSeed seed, int channels, int slices, int positions,
      int types) {
    final String filename = createFile();
    FileOutputStream out = null;
    try {
      out = new FileOutputStream(filename);
    } catch (final Exception e) {
      closeOutput(out);
      e.printStackTrace();
      Assertions.fail(e.getMessage());
    }

    // Write the offsets used in the TSF format
    try {
      @SuppressWarnings("resource")
      final DataOutputStream dos = new DataOutputStream(out);
      dos.writeInt(0);
      dos.writeLong(0);
    } catch (final IOException e) {
      closeOutput(out);
      e.printStackTrace();
      Assertions.fail(e.getMessage());
    }

    // Generate random spots
    final UniformRandomProvider rand = RngUtils.create(seed.getSeedAsLong());
    final int size = 100;
    final Spot[] spots = new Spot[size];
    for (int i = 1; i <= size; i++) {
      final Spot.Builder builder = Spot.newBuilder();
      builder.setChannel(nextInt(rand, channels));
      builder.setSlice(nextInt(rand, slices));
      builder.setPos(nextInt(rand, positions));
      builder.setFluorophoreType(nextInt(rand, types));

      builder.setMolecule(i); // This is a required field but is ignored when reading
      builder.setCluster(rand.nextInt(10));
      builder.setFrame(nextInt(rand, 100));
      builder.setXPosition(rand.nextInt(50));
      builder.setYPosition(rand.nextInt(50));
      builder.setBackground(rand.nextFloat());
      builder.setIntensity(rand.nextFloat());
      builder.setX(rand.nextFloat());
      builder.setY(rand.nextFloat());
      builder.setZ(rand.nextFloat());
      builder.setWidth((float) (Gaussian2DFunction.SD_TO_FWHM_FACTOR * rand.nextDouble()));

      final Spot spot = builder.build();
      spots[i - 1] = spot;
      try {
        spot.writeDelimitedTo(out);
      } catch (final IOException e) {
        closeOutput(out);
        e.printStackTrace();
        Assertions.fail(e.getMessage());
      }
    }

    // Write the header
    // Get the offset to the SpotList message
    long offset = 0;
    try {
      // The offset is the amount to skip forward after reading the int
      // magic number (4 bytes) and long offset (8 bytes)
      // out.flush();
      offset = out.getChannel().position() - 12;
    } catch (final IOException e) {
      closeOutput(out);
      e.printStackTrace();
      Assertions.fail(e.getMessage());
    }

    // Record the SpotList message
    final SpotList.Builder builder = SpotList.newBuilder();

    builder.setApplicationId(1);
    builder.setNrSpots(size);
    builder.setLocationUnits(LocationUnits.PIXELS);
    builder.setIntensityUnits(IntensityUnits.COUNTS);
    builder.setFitMode(FitMode.ONEAXIS);

    builder.setNrChannels(channels);
    builder.setNrSlices(slices);
    builder.setNrPos(positions);
    for (int type = 1; type <= types; type++) {
      final FluorophoreType.Builder typeBuilder = FluorophoreType.newBuilder();
      typeBuilder.setId(type);
      typeBuilder.setDescription("Type " + type);
      typeBuilder.setIsFiducial(rand.nextDouble() < 0.5);
      builder.addFluorophoreTypes(typeBuilder.build());
    }

    final SpotList spotList = builder.build();
    try {
      spotList.writeDelimitedTo(out);
    } catch (final IOException e) {
      e.printStackTrace();
      Assertions.fail(e.getMessage());
    } finally {
      closeOutput(out);
    }

    // Write the offset to the SpotList message into the offset position
    try (RandomAccessFile f = new RandomAccessFile(new File(filename), "rw")) {
      f.seek(4);
      f.writeLong(offset);
    } catch (final Exception e) {
      e.printStackTrace();
      Assertions.fail(e.getMessage());
    }

    // Read each combination
    for (int channel = 1; channel <= channels; channel++) {
      for (int slice = 1; slice <= slices; slice++) {
        for (int position = 1; position <= positions; position++) {
          for (int type = 1; type <= types; type++) {
            final StringBuilder sb = new StringBuilder();
            sb.append(" channel=").append(channel);
            sb.append(" slice=").append(slice);
            sb.append(" position=").append(position);
            sb.append(" fluorophore_type=[").append(type).append(":Type ").append(type)
                .append(":fiducial=").append(builder.getFluorophoreTypes(type - 1).getIsFiducial())
                .append(']');
            // This is needed to trick the Macro class into returning the options
            // for the thread to the GenericDialog used in the ResultsManager
            Thread.currentThread().setName("Run$_");
            Macro.setOptions(sb.toString());

            ResultsManager.setInputFilename(filename);
            final MemoryPeakResults in =
                ResultsManager.loadInputResults(ResultsManager.INPUT_FILE, false, null, null);
            checkEqual(spots, channel, slice, position, type, in);
          }
        }
      }
    }
  }

  private static int nextInt(UniformRandomProvider rand, int n) {
    return (n == 1) ? 1 : 1 + rand.nextInt(n);
  }

  private static void closeOutput(FileOutputStream out) {
    if (out == null) {
      return;
    }

    try {
      out.close();
    } catch (final Exception e) {
      // Ignore exception
    } finally {
      out = null;
    }
  }

  private static void checkEqual(Spot[] spots, int channel, int slice, int position, int type,
      MemoryPeakResults actualResults) {
    Assertions.assertNotNull(actualResults, "Input results are null");

    final MemoryPeakResults expectedResults = extract(spots, channel, slice, position, type);

    Assertions.assertEquals(expectedResults.size(), actualResults.size(), "Size differ");

    final PeakResult[] expected = expectedResults.toArray();
    final PeakResult[] actual = actualResults.toArray();
    for (int i = 0; i < actualResults.size(); i++) {
      final int ii = i;
      final PeakResult p1 = expected[i];
      final PeakResult p2 = actual[i];

      Assertions.assertEquals(p1.getFrame(), p2.getFrame(), () -> "Peak mismatch @ " + ii);

      Assertions.assertEquals(p1.getOrigX(), p2.getOrigX(), () -> "Orig X mismatch @ " + ii);
      Assertions.assertEquals(p1.getOrigY(), p2.getOrigY(), () -> "Orig Y mismatch @ " + ii);
      Assertions.assertEquals(p1.getOrigValue(), p2.getOrigValue(),
          () -> "Orig value mismatch @ " + ii);
      Assertions.assertEquals(p1.getError(), p2.getError(), 1e-6, () -> "Error mismatch @ " + ii);
      Assertions.assertEquals(p1.getNoise(), p2.getNoise(), () -> "Noise mismatch @ " + ii);
      Assertions.assertNotNull(p2.getParameters(), () -> "Params is null @ " + ii);

      Assertions.assertEquals(p1.getBackground(), p2.getBackground(),
          () -> "Background mismatch @ " + ii);
      Assertions.assertEquals(p1.getIntensity(), p2.getIntensity(),
          () -> "Signal mismatch @ " + ii);
      Assertions.assertEquals(p1.getXPosition(), p2.getXPosition(),
          () -> "XPosition mismatch @ " + ii);
      Assertions.assertEquals(p1.getYPosition(), p2.getYPosition(),
          () -> "YPosition mismatch @ " + ii);
      Assertions.assertEquals(p1.getZPosition(), p2.getZPosition(),
          () -> "ZPosition mismatch @ " + ii);
      for (int j = PeakResult.STANDARD_PARAMETERS, size =
          p1.getNumberOfParameters(); j < size; j++) {
        Assertions.assertEquals(p1.getParameter(j), p2.getParameter(j), 1e-6,
            () -> "Parameter mismatch @ " + ii);
      }

      Assertions.assertEquals(p1.getId(), p2.getId(), () -> "ID mismatch @ " + ii);
    }
  }

  private static MemoryPeakResults extract(Spot[] spots, int channel, int slice, int position,
      int type) {
    final MemoryPeakResults results =
        new MemoryPeakResults(PSFHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D));
    for (final Spot spot : spots) {
      if (spot.getChannel() == channel && spot.getSlice() == slice && spot.getPos() == position
          && spot.getFluorophoreType() == type) {
        final int id = spot.getCluster();
        final int startFrame = spot.getFrame();
        final int origX = spot.getXPosition();
        final int origY = spot.getYPosition();
        final float origValue = 0;
        final double error = 0;
        final float noise = 0;
        final float[] params = Gaussian2DPeakResultHelper.createOneAxisParams(spot.getBackground(),
            spot.getIntensity(), spot.getX(), spot.getY(), spot.getZ(),
            (float) (spot.getWidth() / Gaussian2DFunction.SD_TO_FWHM_FACTOR));
        final float[] paramsStdDev = null;
        final IdPeakResult peak = new IdPeakResult(startFrame, origX, origY, origValue, error,
            noise, 0, params, paramsStdDev, id);
        results.add(peak);
      }
    }
    return results;
  }

  private static String createFile() {
    File file;
    try {
      file = File.createTempFile("test", null);
      file.deleteOnExit();
      final String filename = file.getPath();
      return filename;
    } catch (final IOException e) {
      Assertions.fail("Cannot create temp files for IO testing");
    }
    return null; // Allow compilation but the assert will stop the code
  }
}
