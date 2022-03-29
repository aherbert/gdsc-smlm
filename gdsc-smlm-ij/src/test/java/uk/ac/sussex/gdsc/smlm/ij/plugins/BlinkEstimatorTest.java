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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.model.ActivationEnergyImageModel;
import uk.ac.sussex.gdsc.smlm.model.CompoundMoleculeModel;
import uk.ac.sussex.gdsc.smlm.model.DiffusionType;
import uk.ac.sussex.gdsc.smlm.model.FluorophoreSequenceModel;
import uk.ac.sussex.gdsc.smlm.model.ImageModel;
import uk.ac.sussex.gdsc.smlm.model.LocalisationModel;
import uk.ac.sussex.gdsc.smlm.model.MoleculeModel;
import uk.ac.sussex.gdsc.smlm.model.SpatialDistribution;
import uk.ac.sussex.gdsc.smlm.model.SpatialIllumination;
import uk.ac.sussex.gdsc.smlm.model.UniformDistribution;
import uk.ac.sussex.gdsc.smlm.model.UniformIllumination;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.FormatSupplier;

@SuppressWarnings({"javadoc"})
class BlinkEstimatorTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(BlinkEstimatorTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  // Set to sensible simulation parameters
  double diffusionRate = 0.25; // pixels^2/sec
  double pixelPitch = 107;
  int msPerFrame = 1000;
  double photons = 1000;
  float psfWidth = 1.2f;

  final double relativeError = 0.2;
  final double minPhotons = 20;
  final double probabilityDelete = 0;
  final double probabilityAdd = 0;

  static final int LOW = 0;
  static final int MEDIUM = 1;
  static final int HIGH = 2;

  static final int MIN_FITTED_POINTS = 5;
  static final int MAX_FITTED_POINTS = 15;
  double[] blinkingRate = {0.5, 1.5, 4};
  double[] ton = {1.5, 3, 8};
  double[] toff = {2.5, 5, 10};

  // If true then test against the real population statistics.
  // If false then test against the sampled statistics (i.e. using integer frames).
  // Note: When false the success rate is very low so the estimation method does actually
  // account for the integer frame sampling and get the population statistics.
  boolean usePopulationStatistics = true;

  @SeededTest
  void canEstimateBlinkingFromSimulationWithLowNBlinksAndMediumOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[LOW], ton[MEDIUM], toff[MEDIUM],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void
      canEstimateBlinkingFromSimulationWithMediumNBlinksAndMediumOnOffTimesWithFixedMolecules(
          RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[MEDIUM], ton[MEDIUM], toff[MEDIUM],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void canEstimateBlinkingFromSimulationWithHighNBlinksAndMediumOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[HIGH], ton[MEDIUM], toff[MEDIUM],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void canEstimateBlinkingFromSimulationWithLowNBlinksAndHighOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[LOW], ton[HIGH], toff[HIGH],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void canEstimateBlinkingFromSimulationWithMediumNBlinksAndHighOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[MEDIUM], ton[HIGH], toff[HIGH],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void canEstimateBlinkingFromSimulationWithHighNBlinksAndHighOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[HIGH], ton[HIGH], toff[HIGH],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void canEstimateBlinkingFromSimulationWithLowNBlinksAndLowOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[LOW], ton[LOW], toff[LOW],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void canEstimateBlinkingFromSimulationWithMediumNBlinksAndLowOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[MEDIUM], ton[LOW], toff[LOW],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void canEstimateBlinkingFromSimulationWithHighNBlinksAndLowOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngFactory.create(seed.get()), blinkingRate[HIGH], ton[LOW], toff[LOW],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  void findOptimalFittedPoints(RandomSeed seed) {
    // Skip this as it is slow
    Assumptions.assumeTrue(false);
    final UniformRandomProvider rg = RngFactory.create(seed.get());

    final int particles = 1000;
    final double fixedFraction = 1;

    for (final boolean timeAtLowerBound : new boolean[] {false}) {
      final int[] count = new int[MAX_FITTED_POINTS + 1];
      int tests = 0;
      for (int run = 0; run < 3; run++) {
        for (final double n : blinkingRate) {
          for (int i = 0; i < ton.length; i++) {
            tests++;
            final IntOpenHashSet ok = estimateBlinking(rg, n, ton[i], toff[i], particles,
                fixedFraction, timeAtLowerBound, false);
            ok.forEach((int value) -> {
              count[value]++;
            });
          }
        }
      }
      logger.info(FormatSupplier.getSupplier("Time@LowerBound = %b", timeAtLowerBound));
      for (int n = MIN_FITTED_POINTS; n <= MAX_FITTED_POINTS; n++) {
        if (logger.isLoggable(Level.INFO)) {
          final StringBuilder sb = new StringBuilder();
          TextUtils.formatTo(sb, "%2d = %2d/%2d |", n, count[n], tests);
          for (int i = 0; i < count[n]; i++) {
            sb.append('-');
          }
          logger.info(sb.toString());
        }
      }
    }
  }

  private IntOpenHashSet estimateBlinking(UniformRandomProvider rg, double blinkingRate, double ton,
      double toff, int particles, double fixedFraction, boolean timeAtLowerBound,
      boolean doAssert) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MAXIMUM));

    final SpatialIllumination activationIllumination = new UniformIllumination(100);
    int totalSteps = 100;
    final double eAct = totalSteps * 0.3 * activationIllumination.getAveragePhotons();

    final ImageModel imageModel = new ActivationEnergyImageModel(eAct, activationIllumination, ton,
        0, toff, 0, blinkingRate, rg);

    final double[] max = new double[] {256, 256, 32};
    final double[] min = new double[3];
    final SpatialDistribution distribution = new UniformDistribution(min, max, rg.nextInt());
    final List<CompoundMoleculeModel> compounds = new ArrayList<>(1);
    final CompoundMoleculeModel c =
        new CompoundMoleculeModel(1, 0, 0, 0, Arrays.asList(new MoleculeModel(0, 0, 0, 0)));
    c.setDiffusionRate(diffusionRate);
    c.setDiffusionType(DiffusionType.RANDOM_WALK);
    compounds.add(c);

    final List<CompoundMoleculeModel> molecules =
        imageModel.createMolecules(compounds, particles, distribution, false);

    // Activate fluorophores
    final List<? extends FluorophoreSequenceModel> fluorophores =
        imageModel.createFluorophores(molecules, totalSteps);

    totalSteps = checkTotalSteps(totalSteps, fluorophores);

    final List<LocalisationModel> localisations =
        imageModel.createImage(molecules, fixedFraction, totalSteps, photons, 0.5, false);

    // // Remove localisations to simulate missed counts.
    // List<LocalisationModel> newLocalisations = new
    // ArrayList<LocalisationModel>(localisations.size());
    // boolean[] id = new boolean[fluorophores.size() + 1];
    // Statistics photonStats = new Statistics();
    // for (LocalisationModel l : localisations)
    // {
    // photonStats.add(l.getIntensity());
    // // Remove by intensity threshold and optionally at random.
    // if (l.getIntensity() < minPhotons || rand.nextDouble() < pDelete)
    // continue;
    // newLocalisations.add(l);
    // id[l.getId()] = true;
    // }
    // localisations = newLocalisations;
    // logger.info("Photons = %f", photonStats.getMean());
    //
    // List<FluorophoreSequenceModel> newFluorophores = new
    // ArrayList<FluorophoreSequenceModel>(fluorophores.size());
    // for (FluorophoreSequenceModel f : fluorophores)
    // {
    // if (id[f.getId()])
    // newFluorophores.add(f);
    // }
    // fluorophores = newFluorophores;

    final MemoryPeakResults results = new MemoryPeakResults();
    final CalibrationWriter calibration = new CalibrationWriter();
    calibration.setNmPerPixel(pixelPitch);
    calibration.setExposureTime(msPerFrame);
    calibration.setCountPerPhoton(1);
    results.setCalibration(calibration.getCalibration());
    results.setPsf(PsfHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D));
    final float b = 0;
    float intensity;
    final float z = 0;
    for (final LocalisationModel l : localisations) {
      // Remove by intensity threshold and optionally at random.
      if (l.getIntensity() < minPhotons || rg.nextDouble() < probabilityDelete) {
        continue;
      }
      final int frame = l.getTime();
      intensity = (float) l.getIntensity();
      final float x = (float) l.getX();
      final float y = (float) l.getY();
      final float[] params =
          Gaussian2DPeakResultHelper.createParams(b, intensity, x, y, z, psfWidth);
      results.add(frame, 0, 0, 0, 0, 0, 0, params, null);
    }

    // Add random localisations
    // Intensity doesn't matter at the moment for tracing
    intensity = (float) photons;
    for (int i = (int) (localisations.size() * probabilityAdd); i-- > 0;) {
      final int frame = 1 + rg.nextInt(totalSteps);
      final float x = (float) (rg.nextDouble() * max[0]);
      final float y = (float) (rg.nextDouble() * max[1]);
      final float[] params =
          Gaussian2DPeakResultHelper.createParams(b, intensity, x, y, z, psfWidth);
      results.add(frame, 0, 0, 0, 0, 0, 0, params, null);
    }

    // Get actual simulated stats ...
    final Statistics statsNBlinks = new Statistics();
    final Statistics statsTOn = new Statistics();
    final Statistics statsTOff = new Statistics();
    final Statistics statsSampledNBlinks = new Statistics();
    final Statistics statsSampledTOn = new Statistics();
    final StoredDataStatistics statsSampledTOff = new StoredDataStatistics();
    for (final FluorophoreSequenceModel f : fluorophores) {
      statsNBlinks.add(f.getNumberOfBlinks());
      statsTOn.add(f.getOnTimes());
      statsTOff.add(f.getOffTimes());
      final int[] on = f.getSampledOnTimes();
      statsSampledNBlinks.add(on.length);
      statsSampledTOn.add(on);
      statsSampledTOff.add(f.getSampledOffTimes());
    }

    logger.info(
        FormatSupplier.getSupplier("N = %d (%d), N-blinks = %f, tOn = %f, tOff = %f, Fixed = %f",
            fluorophores.size(), localisations.size(), blinkingRate, ton, toff, fixedFraction));
    logger.info(FormatSupplier.getSupplier(
        "Actual N-blinks = %f (%f), tOn = %f (%f), tOff = %f (%f), 95%% = %f, max = %f",
        statsNBlinks.getMean(), statsSampledNBlinks.getMean(), statsTOn.getMean(),
        statsSampledTOn.getMean(), statsTOff.getMean(), statsSampledTOff.getMean(),
        statsSampledTOff.getStatistics().getPercentile(95),
        statsSampledTOff.getStatistics().getMax()));
    logger.info("-=-=--=-");

    final BlinkEstimator be = new BlinkEstimator();
    be.setMaxDarkTime((int) (toff * 10));
    be.setMsPerFrame(msPerFrame);
    be.setRelativeDistance(false);
    final double d = ImageModel.getRandomMoveDistance(diffusionRate);
    be.setSearchDistance((fixedFraction < 1) ? Math.sqrt(2 * d * d) * 3 : 0);
    be.setTimeAtLowerBound(timeAtLowerBound);

    // Assertions.assertTrue("Max dark time must exceed the dark time of the data (otherwise no
    // plateau)",
    // be.maxDarkTime > statsSampledTOff.getStatistics().getMax());

    final int nMolecules = fluorophores.size();
    if (usePopulationStatistics) {
      blinkingRate = statsNBlinks.getMean();
      toff = statsTOff.getMean();
    } else {
      blinkingRate = statsSampledNBlinks.getMean();
      toff = statsSampledTOff.getMean();
    }

    // See if any fitting regime gets a correct answer
    final IntOpenHashSet ok = new IntOpenHashSet();
    for (int numberOfFittedPoints = MIN_FITTED_POINTS; numberOfFittedPoints <= MAX_FITTED_POINTS;
        numberOfFittedPoints++) {
      be.setNumberOfFittedPoints(numberOfFittedPoints);
      be.computeBlinkingRate(results, true);

      final double moleculesError = DoubleEquality.relativeError(nMolecules, be.getNMolecules());
      final double blinksError = DoubleEquality.relativeError(blinkingRate, be.getNBlinks());
      final double offError = DoubleEquality.relativeError(toff * msPerFrame, be.getTOff());
      logger.info(FormatSupplier.getSupplier("Error %d: N = %f, blinks = %f, tOff = %f : %f",
          numberOfFittedPoints, moleculesError, blinksError, offError,
          (moleculesError + blinksError + offError) / 3));

      if (moleculesError < relativeError && blinksError < relativeError
          && offError < relativeError) {
        ok.add(numberOfFittedPoints);
        logger.info("-=-=--=-");
        logger.info(FormatSupplier.getSupplier("*** Correct at %d fitted points ***",
            numberOfFittedPoints));
        if (doAssert) {
          break;
        }
      }

      // if (!be.isIncreaseNFittedPoints())
      // break;
    }

    logger.info("-=-=--=-");

    if (doAssert) {
      Assertions.assertFalse(ok.isEmpty());
    }

    // Assertions.assertEquals("Invalid N-blinks", blinkingRate, be.getNBlinks(), blinkingRate *
    // relativeError);
    // Assertions.assertEquals("Invalid N-molecules", fluorophores.size(), be.getNMolecules(),
    // fluorophores.size() * relativeError);
    // Assertions.assertEquals("Invalid t-off", tOff * msPerFrame, be.getTOff(), tOff * msPerFrame *
    // relativeError);
    return ok;
  }

  private static int checkTotalSteps(int totalSteps,
      List<? extends FluorophoreSequenceModel> fluorophores) {
    for (final FluorophoreSequenceModel f : fluorophores) {
      if (totalSteps < f.getEndTime()) {
        totalSteps = (int) (f.getEndTime() + 1);
      }
    }
    return totalSteps;
  }
}
