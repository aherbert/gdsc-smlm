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

import java.util.List;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.results.DynamicMultipleTargetTracing.DmttConfiguration;
import uk.ac.sussex.gdsc.smlm.results.DynamicMultipleTargetTracing.Trajectory;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class DynamicMultipleTargetTracingTest {
  @SeededTest
  void checkBuilderDefaults(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    final double diffusionCoefficientMaximum = 1 + rng.nextDouble();
    final DmttConfiguration.Builder b = DmttConfiguration.newBuilder(diffusionCoefficientMaximum);

    final DmttConfiguration c1 = b.build();

    // Check round-trip
    for (final DmttConfiguration config : new DmttConfiguration[] {c1, c1.toBuilder().build()}) {
      Assertions.assertEquals(diffusionCoefficientMaximum, config.getDiffusionCoefficientMaximum());
      Assertions.assertTrue(config.getTemporalWindow() > 1);
      Assertions.assertTrue(config.getLocalDiffusionWeight() >= 0);
      Assertions.assertTrue(config.getLocalDiffusionWeight() <= 1);
      Assertions.assertTrue(config.getOnIntensityWeight() >= 0);
      Assertions.assertTrue(config.getOnIntensityWeight() <= 1);
      Assertions.assertTrue(config.getDisappearanceDecayFactor() > 0);
      Assertions.assertTrue(config.getDisappearanceThreshold() > 0);
      Assertions.assertFalse(config.isDisableIntensityModel());
      Assertions.assertFalse(config.isDisableLocalDiffusionModel());
    }
  }

  @SeededTest
  void checkBuilder(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    final int temporalWindow = 2 + rng.nextInt(10);
    final double localDiffusionWeight = rng.nextDouble();
    final double diffusionCoefficientMaximum = 1 + rng.nextDouble();
    final double onIntensityWeight = rng.nextDouble();
    final double disappearanceDecayFactor = 1 + rng.nextDouble();
    final int disappearanceThreshold = 1 + rng.nextInt(10);
    final boolean disableIntensityModel = rng.nextBoolean();
    final boolean disableLocalDiffusionModel = rng.nextBoolean();

    final DmttConfiguration.Builder b = DmttConfiguration.newBuilder(45);

    final DmttConfiguration c1 =
        b.setTemporalWindow(temporalWindow).setLocalDiffusionWeight(localDiffusionWeight)
            .setDiffusionCoefficientMaximum(diffusionCoefficientMaximum)
            .setOnIntensityWeight(onIntensityWeight)
            .setDisappearanceDecayFactor(disappearanceDecayFactor)
            .setDisappearanceThreshold(disappearanceThreshold)
            .setDisableIntensityModel(disableIntensityModel)
            .setDisableLocalDiffusionModel(disableLocalDiffusionModel).build();

    // Check round-trip
    for (final DmttConfiguration config : new DmttConfiguration[] {c1, c1.toBuilder().build()}) {
      Assertions.assertEquals(diffusionCoefficientMaximum, config.getDiffusionCoefficientMaximum());
      Assertions.assertEquals(temporalWindow, config.getTemporalWindow());
      Assertions.assertEquals(localDiffusionWeight, config.getLocalDiffusionWeight());
      Assertions.assertEquals(onIntensityWeight, config.getOnIntensityWeight());
      Assertions.assertEquals(disappearanceDecayFactor, config.getDisappearanceDecayFactor());
      Assertions.assertEquals(disappearanceThreshold, config.getDisappearanceThreshold());
      Assertions.assertEquals(disableIntensityModel, config.isDisableIntensityModel());
      Assertions.assertEquals(disableLocalDiffusionModel, config.isDisableLocalDiffusionModel());
    }
  }

  @Test
  void checkBuilderThrows() {
    Assertions.assertThrows(IllegalArgumentException.class, () -> DmttConfiguration.newBuilder(0));

    final DmttConfiguration.Builder b = DmttConfiguration.newBuilder(1.0);
    // This is allowed: it will cause trajectories to be continuous (no gaps)
    b.setDisappearanceThreshold(0);

    Assertions.assertThrows(IllegalArgumentException.class, () -> b.setTemporalWindow(0));
    Assertions.assertThrows(IllegalArgumentException.class, () -> b.setLocalDiffusionWeight(-1));
    Assertions.assertThrows(IllegalArgumentException.class, () -> b.setLocalDiffusionWeight(1.01));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> b.setDiffusionCoefficientMaximum(0));
    Assertions.assertThrows(IllegalArgumentException.class, () -> b.setOnIntensityWeight(-1));
    Assertions.assertThrows(IllegalArgumentException.class, () -> b.setOnIntensityWeight(1.01));
    Assertions.assertThrows(IllegalArgumentException.class, () -> b.setDisappearanceDecayFactor(0));
    Assertions.assertThrows(IllegalArgumentException.class, () -> b.setDisappearanceThreshold(-1));
  }

  @Test
  void checkConstructorThrows() {
    Assertions.assertThrows(NullPointerException.class,
        () -> new DynamicMultipleTargetTracing(null), "null results");

    final MemoryPeakResults results = new MemoryPeakResults();
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> new DynamicMultipleTargetTracing(results), "Empty results");

    int frame = 0;
    final float x = 1;
    final float y = 2;
    final float intensity = 3;
    results.add(new PeakResult(frame++, x, y, intensity));
    Assertions.assertThrows(ConfigurationException.class,
        () -> new DynamicMultipleTargetTracing(results), "No calibration");

    final CalibrationWriter writer = results.getCalibrationWriterSafe();
    results.setCalibration(writer.getCalibration());
    Assertions.assertThrows(ConversionException.class,
        () -> new DynamicMultipleTargetTracing(results), "No distance calibration");

    writer.setDistanceUnit(DistanceUnit.PIXEL);
    writer.setNmPerPixel(100);
    results.setCalibration(writer.getCalibration());
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> new DynamicMultipleTargetTracing(results), "No time calibration");

    writer.setExposureTime(50);
    results.setCalibration(writer.getCalibration());

    // Should be OK
    new DynamicMultipleTargetTracing(results)
        .traceMolecules(DmttConfiguration.newBuilder(0.2).build());
  }

  @Test
  void testTrajectory() {
    final PeakResult p1 = new PeakResult(0, 1, 2);
    final PeakResult p2 = new PeakResult(1, 2, 3);
    final PeakResult p3 = new PeakResult(2, 3, 4);
    final PeakResult p4 = new PeakResult(3, 4, 5);
    final Trajectory t = new Trajectory(42, p1, true);
    Assertions.assertEquals(42, t.getId());
    Assertions.assertSame(p1, t.getLast(-1));
    Assertions.assertEquals(1, t.size());
    Assertions.assertEquals(1, t.onSize());
    t.add(p2, false);
    t.add(p3, false);
    t.add(p4, true);
    Assertions.assertSame(p4, t.getLast(-1));
    Assertions.assertSame(p3, t.getLast(-2));
    Assertions.assertSame(p2, t.getLast(-3));
    Assertions.assertSame(p1, t.getLast(-4));
    Assertions.assertEquals(4, t.size());
    Assertions.assertEquals(2, t.onSize());
    final LocalList<PeakResult> list = new LocalList<>();
    t.forLast(2, list::add);
    Assertions.assertEquals(2, list.size());
    Assertions.assertSame(p3, list.pop());
    Assertions.assertSame(p4, list.pop());
    t.forLastOn(2, list::add);
    Assertions.assertEquals(2, list.size());
    Assertions.assertSame(p1, list.pop());
    Assertions.assertSame(p4, list.pop());
    final Trace trace = t.toTrace();
    Assertions.assertEquals(4, trace.size());
    Assertions.assertEquals(42, trace.getId());

    // Test setting the gap
    final int frame = t.getLast(-1).getFrame();
    for (int i = 1; i <= 3; i++) {
      t.reset(frame + i);
      Assertions.assertEquals(i - 1, t.gap);
    }

    // Test setting the local statistics
    t.setLocalIntensity(10, 0);
    Assertions.assertEquals(10, t.meanI);
    Assertions.assertEquals(0, t.sdI);
    Assertions.assertFalse(t.isLocalIntensity);
    t.setLocalIntensity(9, 1);
    Assertions.assertEquals(9, t.meanI);
    Assertions.assertEquals(1, t.sdI);
    Assertions.assertTrue(t.isLocalIntensity);

    t.setLocalDiffusion(10, 9);
    Assertions.assertEquals(10, t.r2);
    Assertions.assertFalse(t.isLocalDiffusion);
    t.setLocalDiffusion(9, 10);
    Assertions.assertEquals(9, t.r2);
    Assertions.assertTrue(t.isLocalDiffusion);
  }

  @Test
  void testTrajectoryWithNoOnFrames() {
    final PeakResult p1 = new PeakResult(0, 1, 2);
    final PeakResult p2 = new PeakResult(1, 2, 3);
    final PeakResult p3 = new PeakResult(2, 3, 4);
    final PeakResult p4 = new PeakResult(3, 4, 5);
    final Trajectory t = new Trajectory(42, p1);
    Assertions.assertEquals(42, t.getId());
    Assertions.assertSame(p1, t.getLast(-1));
    Assertions.assertEquals(1, t.size());
    Assertions.assertThrows(NullPointerException.class, () -> t.onSize());
    t.add(p2);
    t.add(p3);
    t.add(p4);
    Assertions.assertSame(p4, t.getLast(-1));
    Assertions.assertSame(p3, t.getLast(-2));
    Assertions.assertSame(p2, t.getLast(-3));
    Assertions.assertSame(p1, t.getLast(-4));
    Assertions.assertEquals(4, t.size());
    final LocalList<PeakResult> list = new LocalList<>();
    t.forLast(2, list::add);
    Assertions.assertEquals(2, list.size());
    Assertions.assertSame(p3, list.pop());
    Assertions.assertSame(p4, list.pop());
    final Trace trace = t.toTrace();
    Assertions.assertEquals(4, trace.size());
    Assertions.assertEquals(42, trace.getId());

    // Test setting the gap
    final int frame = t.getLast(-1).getFrame();
    for (int i = 1; i <= 3; i++) {
      t.reset(frame + i);
      Assertions.assertEquals(i - 1, t.gap);
    }

    // Test setting the local statistics
    Assertions.assertFalse(t.isLocalIntensity);

    t.setLocalDiffusion(10, 9);
    Assertions.assertEquals(10, t.r2);
    Assertions.assertFalse(t.isLocalDiffusion);
    t.setLocalDiffusion(9, 10);
    Assertions.assertEquals(9, t.r2);
    Assertions.assertTrue(t.isLocalDiffusion);
  }

  @Test
  void testIsOn() {
    final double meanI = 1000;
    final double sdI = 100;
    Assertions.assertTrue(DynamicMultipleTargetTracing.isOn(1001, meanI, sdI));
    Assertions.assertTrue(DynamicMultipleTargetTracing.isOn(999, meanI, sdI));
    Assertions.assertFalse(DynamicMultipleTargetTracing.isOn(10, meanI, sdI));
  }

  /**
   * Test trace molecules using 2 molecules. One is fixed and the other moves across it. The tracing
   * should assign the fixed molecule correctly as it has a low local diffusion rate and different
   * intensity.
   */
  @Test
  void testTraceMolecules() {
    // The test is not very robust and fails 10% of the time. A fixed seed corrects this.

    final UniformRandomProvider rng = RngUtils.create(0x12345L);
    final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(rng);
    // localisation precision (in pixels)
    final double s = 0.1;
    final SharedStateContinuousSampler intensity1 =
        SamplerUtils.createGaussianSampler(rng, 1000, 100);
    final SharedStateContinuousSampler intensity2 =
        SamplerUtils.createGaussianSampler(rng, 500, 50);

    final MemoryPeakResults results = new MemoryPeakResults(100);
    final CalibrationWriter writer = results.getCalibrationWriterSafe();
    // 0.1 um pixels, 1 second exposure time
    writer.setDistanceUnit(DistanceUnit.PIXEL);
    writer.setNmPerPixel(100);
    writer.setExposureTime(1000);
    results.setCalibration(writer.getCalibration());

    // First molecule diffuses roughly across the field from top-left to bottom-right.
    // 5 frames is the default for local stats, 15 frames for trajectory removal.
    // Use 20 so we build local stats and can expire a trajectory.
    final int size = 20;
    for (int i = 0; i < size; i++) {
      results.add(new PeakResult(i, (float) (i + gauss.sample() * s),
          (float) (i + gauss.sample() * s), (float) (intensity1.sample())));
    }
    // Second molecule is fixed in the centre with a lower intensity (allow
    // correct matching when tracks overlap)
    final int x = size / 2;
    for (int i = 0; i < size; i++) {
      results.add(new PeakResult(i, (float) (x + gauss.sample() * s),
          (float) (x + gauss.sample() * s), (float) (intensity2.sample())));
    }
    // Add a single molecule that will not connect to anything in the second frame.
    // This should create a trajectory that will expire.
    results.add(new PeakResult(1, size, size, (float) (intensity1.sample())));

    // 1 diffuses top-left to bottom-right.
    // 2 is fixed in the centre.
    // 3 is in the bottom-right for 1 frame.
    //
    // 1
    // 1
    // 1
    // 12
    // 1
    // 1
    // 13
    //
    // Molecule 3 can sometimes connect to the long lifetime molecules once they have been alive
    // long enough to create a local probability model. The default lifetime is 5 frames.
    // Setting this to 10 frames allows a better local model to be created.

    // Move centre to centre each jump => sqrt(2 * 0.1^2) = 0.141 um or 0.02 um^2
    // MSD = 4D => D = 0.02 / 4 = 0.005
    final DmttConfiguration config =
        DmttConfiguration.newBuilder(0.005).setTemporalWindow(10).build();
    final List<Trace> traces = new DynamicMultipleTargetTracing(results).traceMolecules(config);

    // Should have 3 traces
    Assertions.assertEquals(3, traces.size());

    // Assert ids start from 1
    for (int i = 0; i < traces.size(); i++) {
      Assertions.assertEquals(i + 1, traces.get(i).getId());
    }

    // Traces should be 2 full length and 1 single peak
    Assertions.assertEquals(size, traces.get(0).size());
    Assertions.assertEquals(size, traces.get(1).size());
    Assertions.assertEquals(1, traces.get(2).size());

    // Do an analysis on the actual tracks.
    // One should be based in the centre and the other should have parts close to position (i,i)
    // for each frame i.
    final PeakResult[] peaks = results.toArray();
    // Assume traces are initially created using the input order of the results.
    final Trace t1 = traces.get(0);
    final Trace t2 = traces.get(1);
    for (int i = 0; i < size; i++) {
      Assertions.assertSame(peaks[i], t1.get(i));
      Assertions.assertSame(peaks[i + size], t2.get(i));
    }
  }

  /**
   * Test trace molecules using 2 molecules. One is fixed and the other moves past it. The tracing
   * should assign the fixed molecule correctly as it has a low local diffusion rate.
   */
  @Test
  void testTraceMoleculesDisableIntensityModel() {
    final UniformRandomProvider rng = RngUtils.create(125631236L);
    final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(rng);
    // localisation precision (in pixels)
    final double s = 0.1;
    final SharedStateContinuousSampler intensity =
        SamplerUtils.createGaussianSampler(rng, 1000, 100);

    final MemoryPeakResults results = new MemoryPeakResults(100);
    final CalibrationWriter writer = results.getCalibrationWriterSafe();
    // 0.1 um pixels, 1 second exposure time
    writer.setDistanceUnit(DistanceUnit.PIXEL);
    writer.setNmPerPixel(100);
    writer.setExposureTime(1000);
    results.setCalibration(writer.getCalibration());

    // First molecule diffuses roughly across the field from top-left to bottom-right.
    // 5 frames is the default for local stats, 15 frames for trajectory removal.
    // Use 20 so we build local stats and can expire a trajectory.
    final int size = 20;
    final float x1 = size / 2 + 0.5f;
    for (int i = 0; i < size; i++) {
      results.add(new PeakResult(i, (float) (x1 + gauss.sample() * s),
          (float) (i + gauss.sample() * s), (float) (intensity.sample())));
    }
    // Second molecule is fixed in the centre with a same intensity
    final int x = size / 2;
    for (int i = 0; i < size; i++) {
      results.add(new PeakResult(i, (float) (x + gauss.sample() * s),
          (float) (x + gauss.sample() * s), (float) (intensity.sample())));
    }
    // Add a single molecule that will not connect to anything in the second frame.
    // This should create a trajectory that will expire.
    results.add(new PeakResult(1, x1, size, (float) (intensity.sample())));

    // Move centre to centre each jump => 0.1 um or 0.01 um^2
    // MSD = 4D => D = 0.01 / 4 = 0.0025
    final DmttConfiguration config = DmttConfiguration.newBuilder(0.0025)
        .setDisableIntensityModel(true).setTemporalWindow(10).build();
    final List<Trace> traces = new DynamicMultipleTargetTracing(results).traceMolecules(config);

    // Should have 3 traces
    Assertions.assertEquals(3, traces.size());

    // Assert ids start from 1
    for (int i = 0; i < traces.size(); i++) {
      Assertions.assertEquals(i + 1, traces.get(i).getId());
    }

    // Traces should be 2 full length and 1 single peak
    Assertions.assertEquals(size, traces.get(0).size());
    Assertions.assertEquals(size, traces.get(1).size());
    Assertions.assertEquals(1, traces.get(2).size());

    // Do an analysis on the actual tracks.
    // One should be based in the centre and the other should have parts close to position (i,i)
    // for each frame i.
    final PeakResult[] peaks = results.toArray();
    // Assume traces are initially created using the input order of the results.
    final Trace t1 = traces.get(0);
    final Trace t2 = traces.get(1);
    for (int i = 0; i < size; i++) {
      Assertions.assertSame(peaks[i], t1.get(i));
      Assertions.assertSame(peaks[i + size], t2.get(i));
    }
  }

  /**
   * Test trace molecules using 2 molecules. One is fixed and the other moves past it. The tracing
   * should assign the fixed molecule correctly as it has a low local diffusion rate.
   */
  @Test
  void testTraceMoleculesDisableLocalDiffusionModel() {
    // The test is not very robust and fails 20% of the time. A fixed seed corrects this.

    final UniformRandomProvider rng = RngUtils.create(0x12345L);
    final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(rng);
    // localisation precision (in pixels)
    final double s = 0.1;
    final SharedStateContinuousSampler intensity1 =
        SamplerUtils.createGaussianSampler(rng, 1000, 100);
    final SharedStateContinuousSampler intensity2 =
        SamplerUtils.createGaussianSampler(rng, 500, 50);

    final MemoryPeakResults results = new MemoryPeakResults(100);
    final CalibrationWriter writer = results.getCalibrationWriterSafe();
    // 0.1 um pixels, 1 second exposure time
    writer.setDistanceUnit(DistanceUnit.PIXEL);
    writer.setNmPerPixel(100);
    writer.setExposureTime(1000);
    results.setCalibration(writer.getCalibration());

    // First molecule diffuses roughly across the field from top-left to bottom-right.
    // 5 frames is the default for local stats, 15 frames for trajectory removal.
    // Use 20 so we build local stats and can expire a trajectory.
    final int size = 20;
    for (int i = 0; i < size; i++) {
      results.add(new PeakResult(i, (float) (i + gauss.sample() * s),
          (float) (i + gauss.sample() * s), (float) (intensity1.sample())));
    }
    // Second molecule is fixed in the centre with a lower intensity (allow
    // correct matching when tracks overlap)
    final int x = size / 2;
    for (int i = 0; i < size; i++) {
      results.add(new PeakResult(i, (float) (x + gauss.sample() * s),
          (float) (x + gauss.sample() * s), (float) (intensity2.sample())));
    }
    // Add a single molecule that will not connect to anything in the second frame.
    // This should create a trajectory that will expire.
    results.add(new PeakResult(1, size, size, (float) (intensity1.sample())));

    // Move centre to centre each jump => sqrt(2 * 0.1^2) = 0.141 um or 0.02 um^2
    // MSD = 4D => D = 0.02 / 4 = 0.005
    final DmttConfiguration config = DmttConfiguration.newBuilder(0.005)
        .setDisableLocalDiffusionModel(true).setTemporalWindow(10).build();
    final List<Trace> traces = new DynamicMultipleTargetTracing(results).traceMolecules(config);

    // Should have 3 traces
    Assertions.assertEquals(3, traces.size());

    // Assert ids start from 1
    for (int i = 0; i < traces.size(); i++) {
      Assertions.assertEquals(i + 1, traces.get(i).getId());
    }

    // Traces should be 2 full length and 1 single peak
    Assertions.assertEquals(size, traces.get(0).size());
    Assertions.assertEquals(size, traces.get(1).size());
    Assertions.assertEquals(1, traces.get(2).size());

    // Do an analysis on the actual tracks.
    // One should be based in the centre and the other should have parts close to position (i,i)
    // for each frame i.
    final PeakResult[] peaks = results.toArray();
    // Assume traces are initially created using the input order of the results.
    final Trace t1 = traces.get(0);
    final Trace t2 = traces.get(1);
    for (int i = 0; i < size; i++) {
      Assertions.assertSame(peaks[i], t1.get(i));
      Assertions.assertSame(peaks[i + size], t2.get(i));
    }
  }
}
