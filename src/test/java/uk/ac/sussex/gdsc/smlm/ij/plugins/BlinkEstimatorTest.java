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

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
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
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

import gnu.trove.procedure.TIntProcedure;
import gnu.trove.set.hash.TIntHashSet;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class BlinkEstimatorTest {
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
  final double pDelete = 0;
  final double pAdd = 0;

  int LOW = 0;
  int MEDIUM = 1;
  int HIGH = 2;

  int MIN_FITTED_POINTS = 5;
  int MAX_FITTED_POINTS = 15;
  double[] nBlinks = {0.5, 1.5, 4};
  double[] tOn = {1.5, 3, 8};
  double[] tOff = {2.5, 5, 10};

  // If true then test against the real population statistics.
  // If false then test against the sampled statistics (i.e. using integer frames).
  // Note: When false the success rate is very low so the estimation method does actually
  // account for the integer frame sampling and get the population statistics.
  boolean usePopulationStatistics = true;

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithLowNBlinksAndMediumOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[LOW], tOn[MEDIUM], tOff[MEDIUM],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithMediumNBlinksAndMediumOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[MEDIUM], tOn[MEDIUM],
        tOff[MEDIUM], particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithHighNBlinksAndMediumOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[HIGH], tOn[MEDIUM],
        tOff[MEDIUM], particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithLowNBlinksAndHighOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[LOW], tOn[HIGH], tOff[HIGH],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithMediumNBlinksAndHighOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[MEDIUM], tOn[HIGH], tOff[HIGH],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithHighNBlinksAndHighOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[HIGH], tOn[HIGH], tOff[HIGH],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithLowNBlinksAndLowOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[LOW], tOn[LOW], tOff[LOW],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithMediumNBlinksAndLowOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[MEDIUM], tOn[LOW], tOff[LOW],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void canEstimateBlinkingFromSimulationWithHighNBlinksAndLowOnOffTimesWithFixedMolecules(
      RandomSeed seed) {
    final int particles = 1000;
    final double fixedFraction = 1;
    estimateBlinking(RngUtils.create(seed.getSeedAsLong()), nBlinks[HIGH], tOn[LOW], tOff[LOW],
        particles, fixedFraction, false, true);
  }

  @SeededTest
  public void findOptimalFittedPoints(RandomSeed seed) {
    // Skip this as it is slow
    Assumptions.assumeTrue(false);
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());

    final int particles = 1000;
    final double fixedFraction = 1;

    for (final boolean timeAtLowerBound : new boolean[] {false}) {
      final int[] count = new int[MAX_FITTED_POINTS + 1];
      int tests = 0;
      for (int run = 0; run < 3; run++) {
        for (final double n : nBlinks) {
          for (int i = 0; i < tOn.length; i++) {
            tests++;
            final TIntHashSet ok = estimateBlinking(rg, n, tOn[i], tOff[i], particles,
                fixedFraction, timeAtLowerBound, false);
            ok.forEach(new TIntProcedure() {
              @Override
              public boolean execute(int value) {
                count[value]++;
                return true;
              }
            });
          }
        }
      }
      logger.info(FunctionUtils.getSupplier("Time@LowerBound = %b", timeAtLowerBound));
      for (int nFittedPoints =
          MIN_FITTED_POINTS; nFittedPoints <= MAX_FITTED_POINTS; nFittedPoints++) {
        if (logger.isLoggable(Level.INFO)) {
          final StringBuilder sb = new StringBuilder();
          sb.append(String.format("%2d = %2d/%2d |", nFittedPoints, count[nFittedPoints], tests));
          for (int i = 0; i < count[nFittedPoints]; i++) {
            sb.append('-');
          }
          logger.info(sb.toString());
        }
      }
    }
  }

  private TIntHashSet estimateBlinking(UniformRandomProvider rg, double nBlinks, double tOn,
      double tOff, int particles, double fixedFraction, boolean timeAtLowerBound,
      boolean doAssert) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MAXIMUM));

    final SpatialIllumination activationIllumination = new UniformIllumination(100);
    int totalSteps = 100;
    final double eAct = totalSteps * 0.3 * activationIllumination.getAveragePhotons();

    final ImageModel imageModel =
        new ActivationEnergyImageModel(eAct, activationIllumination, tOn, 0, tOff, 0, nBlinks);
    imageModel.setRandomGenerator(new RandomGeneratorAdapter(rg));

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
    results.setPSF(PSFHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D));
    final float b = 0;
    float intensity;
    final float z = 0;
    for (final LocalisationModel l : localisations) {
      // Remove by intensity threshold and optionally at random.
      if (l.getIntensity() < minPhotons || rg.nextDouble() < pDelete) {
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
    for (int i = (int) (localisations.size() * pAdd); i-- > 0;) {
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
        FunctionUtils.getSupplier("N = %d (%d), N-blinks = %f, tOn = %f, tOff = %f, Fixed = %f",
            fluorophores.size(), localisations.size(), nBlinks, tOn, tOff, fixedFraction));
    logger.info(FunctionUtils.getSupplier(
        "Actual N-blinks = %f (%f), tOn = %f (%f), tOff = %f (%f), 95%% = %f, max = %f",
        statsNBlinks.getMean(), statsSampledNBlinks.getMean(), statsTOn.getMean(),
        statsSampledTOn.getMean(), statsTOff.getMean(), statsSampledTOff.getMean(),
        statsSampledTOff.getStatistics().getPercentile(95),
        statsSampledTOff.getStatistics().getMax()));
    logger.info("-=-=--=-");

    final BlinkEstimator be = new BlinkEstimator();
    be.maxDarkTime = (int) (tOff * 10);
    be.msPerFrame = msPerFrame;
    be.relativeDistance = false;
    final double d = ImageModel.getRandomMoveDistance(diffusionRate);
    be.searchDistance = (fixedFraction < 1) ? Math.sqrt(2 * d * d) * 3 : 0;
    be.timeAtLowerBound = timeAtLowerBound;
    be.showPlots = false;

    // Assertions.assertTrue("Max dark time must exceed the dark time of the data (otherwise no
    // plateau)",
    // be.maxDarkTime > statsSampledTOff.getStatistics().getMax());

    final int nMolecules = fluorophores.size();
    if (usePopulationStatistics) {
      nBlinks = statsNBlinks.getMean();
      tOff = statsTOff.getMean();
    } else {
      nBlinks = statsSampledNBlinks.getMean();
      tOff = statsSampledTOff.getMean();
    }

    // See if any fitting regime gets a correct answer
    final TIntHashSet ok = new TIntHashSet();
    for (int nFittedPoints =
        MIN_FITTED_POINTS; nFittedPoints <= MAX_FITTED_POINTS; nFittedPoints++) {
      be.nFittedPoints = nFittedPoints;
      be.computeBlinkingRate(results, true);

      final double moleculesError = DoubleEquality.relativeError(nMolecules, be.getNMolecules());
      final double blinksError = DoubleEquality.relativeError(nBlinks, be.getNBlinks());
      final double offError = DoubleEquality.relativeError(tOff * msPerFrame, be.getTOff());
      logger.info(FunctionUtils.getSupplier("Error %d: N = %f, blinks = %f, tOff = %f : %f",
          nFittedPoints, moleculesError, blinksError, offError,
          (moleculesError + blinksError + offError) / 3));

      if (moleculesError < relativeError && blinksError < relativeError
          && offError < relativeError) {
        ok.add(nFittedPoints);
        logger.info("-=-=--=-");
        logger
            .info(FunctionUtils.getSupplier("*** Correct at %d fitted points ***", nFittedPoints));
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

    // Assertions.assertEquals("Invalid N-blinks", nBlinks, be.getNBlinks(), nBlinks *
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
