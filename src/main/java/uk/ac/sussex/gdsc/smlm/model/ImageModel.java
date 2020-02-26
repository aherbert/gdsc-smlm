/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.PoissonSampler;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;

/**
 * Contains a model for an image of blinking fluorophores.
 *
 * <p>Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated
 * localization microscopy images for quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public abstract class ImageModel {
  /** Define the Z-axis. */
  private static final double[] Z_AXIS = new double[] {0, 0, 1};
  /** Average on-state time. */
  protected double onTime;
  /** Average off-state time for the first dark state. */
  protected double offTime1;
  /** Average off-state time for the second dark state. */
  protected double offTime2;
  /**
   * Average number of blinks int the first dark state (used for each burst between second dark
   * states).
   */
  protected double blinks1;
  /** Average number of blinks into the second dark state. */
  protected double blinks2;

  /**
   * Specifies the maximum number of frames that will be simulated in
   * {@link #createFluorophores(List, int)}. Any activation time above this limit returned by
   * {@link #createActivationTime(double[])} will be ignored and there will be no call to
   * {@link #createFluorophore(int, double[], double)}.
   */
  protected int frameLimit;

  private UniformRandomProvider random;
  private RealDistribution photonDistribution;
  private SpatialDistribution confinementDistribution;
  private int confinementAttempts = 5;
  private boolean useGeometricDistribution;
  private boolean photonBudgetPerFrame = true;
  private boolean diffusion2D;
  private boolean rotation2D;

  /**
   * Construct a new image model.
   *
   * @param onTime Average on-state time
   * @param offTime1 Average off-state time for the first dark state
   * @param offTime2 Average off-state time for the second dark state
   * @param blinks1 Average number of blinks in the first dark state (used for each burst between
   *        second dark states)
   * @param blinks2 Average number of blinks in the second dark state
   * @param rng the random generator for creating the image
   */
  public ImageModel(double onTime, double offTime1, double offTime2, double blinks1, double blinks2,
      UniformRandomProvider rng) {
    checkParameter("onTime", onTime);
    checkParameter("offTime1", offTime1);
    checkParameter("offTime2", offTime2);
    checkParameter("blinks1", blinks1);
    checkParameter("blinks2", blinks2);
    this.onTime = onTime;
    this.offTime1 = offTime1;
    this.offTime2 = offTime2;
    this.blinks1 = blinks1;
    this.blinks2 = blinks2;
    setUniformRandomProvider(rng);
  }

  /**
   * Check the parameter is not less than zero.
   *
   * @param param the parameter name
   * @param value the value
   */
  protected void checkParameter(String param, double value) {
    if (value < 0) {
      throw new IllegalArgumentException(param + " cannot be less than 0");
    }
  }

  /**
   * Get the average on-state time.
   *
   * @return the on time
   */
  public double getOnTime() {
    return onTime;
  }

  /**
   * Get the average off-state time for the first dark state.
   *
   * @return the off time 1
   */
  public double getOffTime1() {
    return offTime1;
  }

  /**
   * Get the average off-state time for the second dark state.
   *
   * @return the off time 2
   */
  public double getOffTime2() {
    return offTime2;
  }

  /**
   * Get the average number of blinks in the first dark state (used for each burst between second
   * dark states).
   *
   * @return the blinks
   */
  public double getBlinks1() {
    return blinks1;
  }

  /**
   * Get the average number of blinks in the second dark state.
   *
   * @return the blinks
   */
  public double getBlinks2() {
    return blinks2;
  }

  /**
   * Gets the random data generator.
   *
   * @return The random data generator.
   */
  protected UniformRandomProvider getRandom() {
    return random;
  }

  /**
   * Creates N particles by sampling from the distribution. If the distribution returns null (no
   * coordinates) then no more attempts to generate particles is made.
   *
   * <p>The compound's molecule coordinates will be updated by calling
   * {@link CompoundMoleculeModel#centre()}.
   *
   * @param compounds the compounds
   * @param particles the particles
   * @param distribution the distribution
   * @param rotate the rotate
   * @return the list
   */
  public List<CompoundMoleculeModel> createMolecules(List<CompoundMoleculeModel> compounds,
      int particles, SpatialDistribution distribution, boolean rotate) {
    if (compounds == null || compounds.isEmpty()) {
      throw new IllegalArgumentException("Compounds cannot be null or empty");
    }

    // Sum the fractions
    double total = 0;
    for (final CompoundMoleculeModel c : compounds) {
      total += c.getFraction();
      c.centre();
    }

    final ArrayList<CompoundMoleculeModel> molecules = new ArrayList<>(particles);
    for (int i = 1; i <= particles; i++) {
      // Get the next random compound
      double fraction = random.nextDouble() * total;
      for (int j = 0; j < compounds.size(); j++) {
        final CompoundMoleculeModel c = compounds.get(j);
        fraction -= c.getFraction();
        if (fraction <= 0) {
          final double[] xyz = distribution.next();
          if (xyz == null || xyz.length < 3) {
            return molecules;
          }

          final CompoundMoleculeModel m =
              new CompoundMoleculeModel(i, xyz, copyMolecules(c), false);
          m.setLabel(j);
          m.setDiffusionRate(c.getDiffusionRate());
          m.setDiffusionType(c.getDiffusionType());
          if (rotate) {
            rotate(m);
          }
          molecules.add(m);
          break;
        }
      }
    }
    return molecules;
  }

  private static List<MoleculeModel> copyMolecules(CompoundMoleculeModel compound) {
    int count = compound.getSize();
    final List<MoleculeModel> list = new ArrayList<>(count);
    while (count-- > 0) {
      final MoleculeModel molecule = compound.getMolecule(count);
      list.add(new MoleculeModel(count + 1, molecule.getMass(),
          Arrays.copyOf(molecule.getCoordinates(), 3)));
    }
    return list;
  }

  private void diffuse(CompoundMoleculeModel model, final double diffusionRate,
      final double[] axis) {
    final double z = model.xyz[2];
    switch (model.getDiffusionType()) {
      case GRID_WALK:
        model.walk(diffusionRate, random);
        break;

      case LINEAR_WALK:
        model.slide(diffusionRate, axis, random);
        break;

      case RANDOM_WALK:
        model.move(diffusionRate, random);
        break;

      default:
        throw new IllegalStateException("Unsupported diffusion type: " + model.getDiffusionType());
    }
    if (diffusion2D) {
      model.xyz[2] = z;
    }
  }

  private void rotate(CompoundMoleculeModel model) {
    if (rotation2D) {
      model.rotateRandomAngle(Z_AXIS, 180, random);
    } else {
      model.rotateRandom(180, random);
    }
  }

  /**
   * Generates fluorophores for the molecules spatial positions. Only simulate up to the given
   * maximum number of frames.
   *
   * <p>Replace the CompoundMoleculeModel objects in the list with new compounds containing the
   * generated FluorophoreSequenceModel objects. If no fluorophores can be generated for a compound
   * then it is removed. Since the fluorophores are part of a compound their coordinates are
   * relative to the compound centre of mass.
   *
   * <p>Note that the activation energy is sampled at the spatial position without movement. This is
   * therefore an approximation of the energy the molecule would receive if it were moving.
   *
   * @param molecules the molecules
   * @param frames the frames
   * @return the fluorophores
   */
  public List<FluorophoreSequenceModel> createFluorophores(List<CompoundMoleculeModel> molecules,
      int frames) {
    frameLimit = frames;
    final ArrayList<FluorophoreSequenceModel> list = new ArrayList<>(molecules.size());
    for (int i = 0; i < molecules.size();) {
      final CompoundMoleculeModel c = molecules.get(i);
      final List<FluorophoreSequenceModel> fluorophores = new ArrayList<>(c.getSize());
      final List<MoleculeModel> removed = new ArrayList<>(c.getSize());
      for (int n = c.getSize(); n-- > 0;) {
        // Molecule Id is ignored since we renumber the sorted collection at the end
        final double[] xyz = c.getRelativeCoordinates(n);
        final double tAct = createActivationTime(xyz);
        FluorophoreSequenceModel flouro = null;
        if (frames == 0 || tAct < frames) {
          flouro = createFluorophore(0, xyz, tAct);
        }
        if (flouro != null && flouro.getCoordinates() != null
            && flouro.getEndTime() > flouro.getStartTime()) {
          flouro.setLabel(c.getLabel());
          list.add(flouro);
          fluorophores.add(flouro);
        } else {
          // If we remove the molecule from the compound then the rotation will not be around
          // the true COM but the new COM of the reduced size compound. Keep these molecules
          // so the rotation is valid.
          // Note: Fluorophores have mass zero. If this changes then this will have to be updated.
          removed.add(new MoleculeModel(0, c.getRelativeCoordinates(n)));
        }
      }
      if (fluorophores.isEmpty()) {
        molecules.remove(i);
        // No increment as the list has been reduced in size
      } else {
        removed.addAll(fluorophores);
        final CompoundMoleculeModel newC =
            new CompoundMoleculeModel(c.getId(), c.getCoordinates(), removed, false);
        newC.setDiffusionRate(c.getDiffusionRate());
        newC.setDiffusionType(c.getDiffusionType());
        molecules.set(i, newC);

        // Debug the fluorophore position
        // c = molecules.get(i);
        // System.out.printf("fluorophores XYZ = %f,%f,%f (%d)\n", c.getX(), c.getY(), c.getZ(),
        // fluorophores.size());
        // for (int j = 0; j < c.getSize(); j++)
        // {
        // double[] xyz = c.getCoordinates(j);
        // System.out.printf("fluorophore %d XYZ = %f,%f,%f\n", j + 1, xyz[0], xyz[1], xyz[2]);
        // }

        i++;
      }
    }
    // Sort by time
    Collections.sort(list, FluorophoreSequenceModel::compare);
    // Renumber
    int id = 1;
    for (final FluorophoreSequenceModel f : list) {
      f.setId(id++);
    }
    return list;
  }

  /**
   * Create a fluorophore activation time for the given position. Derived classes can check the
   * {@link #frameLimit} variable to determine the limit of the simulation. If non-zero only
   * activation times below this will be used to generate fluorophores.
   *
   * @param xyz the xyz
   * @return the activation time
   */
  protected abstract double createActivationTime(double[] xyz);

  /**
   * Create a fluorophore with the given id, position and activation time.
   *
   * @param id the id
   * @param xyz the xyz
   * @param activationTime the activation time (generated using
   *        {@link #createActivationTime(double[])})
   * @return a fluorophore
   */
  protected abstract FluorophoreSequenceModel createFluorophore(int id, double[] xyz,
      double activationTime);

  /**
   * Add to the list but link up the continuous pulses with previous/next pointers.
   *
   * @param localisations the localisations
   * @param models the models
   */
  private static void linkLocalisations(List<LocalisationModel> localisations,
      LocalisationModel[] models) {
    LocalisationModel previous = null;
    for (int i = 0; i < models.length; i++) {
      if (models[i] != null) {
        models[i].setPrevious(previous);
        localisations.add(models[i]);
      }
      previous = models[i];
    }
  }

  private static void generateOnTimes(int maxFrames, final double frameInterval,
      List<double[]> bursts, int sequenceStartT, double[] onTime, int[] state) {
    for (int i = 0; i < bursts.size(); i++) {
      final double[] burst = bursts.get(i);
      final int on = (int) (burst[0] / frameInterval) - sequenceStartT;
      int off = (int) (burst[1] / frameInterval) - sequenceStartT;

      // If the off-time is truncated the modulus function used below is
      // not valid.
      boolean offTimeTrucated = false;

      if (on >= onTime.length) {
        break;
      }
      if (off >= onTime.length) {
        burst[1] = maxFrames * frameInterval;
        off = onTime.length - 1;
        offTimeTrucated = true;
      }

      if (on == off) {
        onTime[on] += (burst[1] - burst[0]);
        // This is a LocalisationModel.SINGLE (state value = 0)
      } else {
        for (int t = on; t <= off; t++) {
          if (t > on) {
            state[t] |= LocalisationModel.PREVIOUS;
          }
          if (t < off) {
            state[t] |= LocalisationModel.NEXT;
          }

          // Get signal on-time
          if (t == on) {
            // Signal turned on this frame
            onTime[t] += frameInterval - (burst[0] % frameInterval);
          } else if (t == off) {
            // Signal turned off this frame
            if (offTimeTrucated) {
              // The off-time was truncated so the burst is on for
              // the entire duration
              onTime[t] += 1;
            } else {
              onTime[t] += (burst[1] % frameInterval);
            }
          } else {
            // Continuous signal
            onTime[t] = frameInterval;
            state[t] |= LocalisationModel.CONTINUOUS;
          }
        }
      }
    }
  }

  private static void sortByTime(List<LocalisationModel> localisations) {
    // Sort by time-only (not signal as well which is the default).
    // This means the localisations are in the same order as the input
    // fluorophores.
    Collections.sort(localisations, new Comparator<LocalisationModel>() {
      @Override
      public int compare(LocalisationModel o1, LocalisationModel o2) {
        if (o1.getTime() < o2.getTime()) {
          return -1;
        }
        if (o1.getTime() > o2.getTime()) {
          return 1;
        }
        return 0;
      }
    });
  }

  /**
   * Simulate an image of fluorophores. The total amount of time a fluorophore is on (i.e. sum of
   * tOn) is used to create the signal strength using the specified correlation.
   *
   * <p>A second set of fluorophores, independent of the first, are generated using
   * {@link #createFluorophores(List, int)}. The correlated on-times will be created by combining
   * the times using the provided correlation (r):
   *
   * <pre>
   * // X1 : Fluorophore total on-times
   * // X2 : Independently generated fluorophore total on-times
   * // r  : correlation
   * a = sqrt(1 - r * r)
   * newX = (r * X1 + a * X2) / (r + a)
   * </pre>
   *
   * <p>The signal is proportional to newly generated on-times (newX) with an average of the
   * specified photon budget.
   *
   * <p>The photon budget can either be distributed evenly over the fluorophore lifetime or per
   * frame (see {@link #isPhotonBudgetPerFrame()}). Each frame signal output will be subject to
   * Poisson sampling.
   *
   * <p>If the input correlation is zero then the number of photons will be sampled from the
   * configured photon distribution or, if this is null, will be uniform.
   *
   * <p>A random fraction of the fluorophores will move using a random walk with the diffusion
   * coefficient defined in the compound.
   *
   * @param compoundFluorophores The compounds containing the fluorophores
   * @param fixedFraction The fraction of molecules that are fixed
   * @param maxFrames The maximum frame for the simulation
   * @param photonBudget the average number of photons per frame/lifetime; see
   *        {@link #isPhotonBudgetPerFrame()}
   * @param correlation The correlation between the number of photons and the total on-time of the
   *        fluorophore
   * @param rotate Rotate the molecule per frame
   * @return the localisations
   */
  public List<LocalisationModel> createImage(List<CompoundMoleculeModel> compoundFluorophores,
      double fixedFraction, int maxFrames, double photonBudget, double correlation,
      boolean rotate) {
    final List<LocalisationModel> localisations = new ArrayList<>();
    // Extract the fluorophores in all the compounds
    final ArrayList<FluorophoreSequenceModel> fluorophores =
        new ArrayList<>(compoundFluorophores.size());
    for (final CompoundMoleculeModel c : compoundFluorophores) {
      for (int i = c.getSize(); i-- > 0;) {
        if (c.getMolecule(i) instanceof FluorophoreSequenceModel) {
          fluorophores.add((FluorophoreSequenceModel) c.getMolecule(i));
        }
      }
    }
    final int nMolecules = fluorophores.size();

    // Check the correlation bounds.
    // Correlation for photons per frame verses total on time should be negative.
    // Correlation for total photons verses total on time should be positive.
    double boundedCorrelation;
    if (photonBudgetPerFrame) {
      boundedCorrelation = MathUtils.clip(-1, 0, correlation);
    } else {
      boundedCorrelation = MathUtils.clip(0, 1, correlation);
    }

    // Create a photon budget for each fluorophore
    final double[] photons = new double[nMolecules];

    // Generate a second set of on times using the desired correlation
    if (boundedCorrelation != 0) {
      // Observations on real data show:
      // - there is a weak positive correlation between total photons and time
      // - There is a weak negative correlation between photons per frame and total on-time

      // How to generate a random correlation:
      // http://www.uvm.edu/~dhowell/StatPages/More_Stuff/CorrGen.html
      // http://stats.stackexchange.com/questions/13382/how-to-define-a-distribution-such-that-draws-from-it-correlate-with-a-draw-from

      // Create a second set of fluorophores. This is used to generate the correlated photon data
      final List<FluorophoreSequenceModel> fluorophores2 = new ArrayList<>();
      while (fluorophores2.size() < fluorophores.size()) {
        final FluorophoreSequenceModel f = createFluorophore(0, new double[] {0, 0, 0}, maxFrames);
        if (f != null) {
          fluorophores2.add(f);
        }
      }

      final double a = Math.sqrt(1 - boundedCorrelation * boundedCorrelation);

      // Q - How to generate a negative correlation?
      // Generate a positive correlation then invert the data and shift to the same range
      final boolean invert = (boundedCorrelation < 0);
      boundedCorrelation = Math.abs(boundedCorrelation);

      StoredDataStatistics correlatedOnTime = new StoredDataStatistics();
      for (int i = 0; i < nMolecules; i++) {
        final double X = getTotalOnTime(fluorophores.get(i));
        final double Z = getTotalOnTime(fluorophores2.get(i));
        correlatedOnTime.add((boundedCorrelation * X + a * Z) / (boundedCorrelation + a));
      }

      if (invert) {
        final double min = correlatedOnTime.getStatistics().getMin();
        final double range = correlatedOnTime.getStatistics().getMax() - min;
        final StoredDataStatistics newCorrelatedOnTime = new StoredDataStatistics();
        for (final double d : correlatedOnTime.getValues()) {
          newCorrelatedOnTime.add(range - d + min);
        }
        correlatedOnTime = newCorrelatedOnTime;
      }

      // Get the average on time from the correlated sample.
      // Using the population value (tOn * (1+blinks1)) over-estimates the on time.
      final double averageTotalTOn = correlatedOnTime.getMean();

      // Now create the localisations
      final double[] correlatedOnTimes = correlatedOnTime.getValues();
      for (int i = 0; i < nMolecules; i++) {
        // Generate photons using the correlated on time data
        final double p = photonBudget * correlatedOnTimes[i] / averageTotalTOn;
        photons[i] = p;
      }

    } else if (photonDistribution != null) {
      // Sample from the provided distribution. Do not over-write the class level distribution
      // to allow running again with a different shape parameter / photon budget.

      // Ensure the custom distribution is scaled to the correct photon budget
      final double photonScale = photonBudget / photonDistribution.getNumericalMean();

      // Generate photons sampling from the photon budget
      for (int i = 0; i < nMolecules; i++) {
        photons[i] = photonDistribution.sample() * photonScale;
      }
    } else {
      // No distribution so use the same number for all
      Arrays.fill(photons, photonBudget);
    }

    int photonIndex = 0;
    if (fixedFraction > 0) {
      for (final CompoundMoleculeModel c : compoundFluorophores) {
        final boolean fixed = random.nextDouble() < fixedFraction;
        photonIndex += createLocalisation(c, localisations, fixed, maxFrames, photons, photonIndex,
            !fixed && rotate);
      }
    } else {
      // No fixed molecules
      for (final CompoundMoleculeModel c : compoundFluorophores) {
        photonIndex +=
            createLocalisation(c, localisations, false, maxFrames, photons, photonIndex, rotate);
      }
    }

    sortByTime(localisations);

    return localisations;
  }

  private static double getTotalOnTime(FluorophoreSequenceModel flouro) {
    return MathUtils.sum(flouro.getOnTimes());
  }

  /**
   * Create localisations for each time frame using the fluorophores in the compound. The compound
   * is moved using the specified diffusion rate.
   *
   * @param compound Compound containing one or more fluorophores
   * @param localisations Output list to put all localisations
   * @param fixed the fixed flag
   * @param maxFrames The maximum number of frames. Fluorophores still active beyond this frame will
   *        not be represented
   * @param photons Array of the photon budget for all the fluorophores
   * @param photonIndex The current index within the photon budget array
   * @param rotate Rotate at each step
   * @return The number of fluorophores that were drawn
   */
  private int createLocalisation(CompoundMoleculeModel compound,
      List<LocalisationModel> localisations, boolean fixed, int maxFrames, double[] photons,
      int photonIndex, boolean rotate) {
    // -=-=-
    // TODO:
    // Account for the variable nature of fluorophore intensities in real data?
    // Perhaps add quantum efficiency simulation using a Binomial distribution?
    // -=-=-
    final double frameInterval = 1;

    // Extract the fluorophores
    final ArrayList<FluorophoreSequenceModel> fluorophores = new ArrayList<>();
    for (int i = compound.getSize(); i-- > 0;) {
      if (compound.getMolecule(i) instanceof FluorophoreSequenceModel) {
        fluorophores.add((FluorophoreSequenceModel) compound.getMolecule(i));
      }
    }
    final int nFluorophores = fluorophores.size();

    // Get the duration of each fluorophore sequence
    final List<List<double[]>> bursts = new ArrayList<>(nFluorophores);
    final int[] sequenceStartT = new int[nFluorophores];
    final int[] sequenceEndT = new int[nFluorophores];

    for (int i = 0; i < nFluorophores; i++) {
      final FluorophoreSequenceModel fluorophore = fluorophores.get(i);
      bursts.add(fluorophore.getBurstSequence());
      sequenceStartT[i] = (int) (fluorophore.getStartTime() / frameInterval);
      sequenceEndT[i] = (int) (fluorophore.getEndTime() / frameInterval);
    }

    // Get the maximum and minimum start and end times
    final int sequenceStart = MathUtils.min(sequenceStartT);
    int sequenceEnd = MathUtils.max(sequenceEndT);

    // Time-range check. Note that the final frames are 1-based
    if (sequenceStart > maxFrames - 1) {
      return nFluorophores;
    }
    if (sequenceEnd > maxFrames - 1) {
      sequenceEnd = maxFrames - 1;
    }

    // For each fluorophore build a graph of on-time in frame intervals
    final int nSteps = sequenceEnd - sequenceStart + 1;
    final double[][] onTime = new double[nFluorophores][nSteps];
    final int[][] state = new int[nFluorophores][nSteps];
    final double[] totalOnTime = new double[nFluorophores];
    final double[] photonBudget = new double[nFluorophores];
    for (int i = 0; i < nFluorophores; i++) {
      generateOnTimes(maxFrames, frameInterval, bursts.get(i), sequenceStart, onTime[i], state[i]);
      totalOnTime[i] = MathUtils.sum(onTime[i]);
      photonBudget[i] = photons[i + photonIndex];
    }

    // Convert diffusion co-efficient into the standard deviation for the random walk
    final double diffusionRate;
    double[] axis = null;
    if (fixed) {
      diffusionRate = 0;
    } else if (compound.getDiffusionType() == DiffusionType.LINEAR_WALK) {
      diffusionRate = (diffusion2D) ? getRandomMoveDistance2D(compound.getDiffusionRate())
          : getRandomMoveDistance3D(compound.getDiffusionRate());
      // Create a random unit vector using 2/3 Gaussian variables
      axis = new double[3];
      double length = 0;
      final NormalizedGaussianSampler gauss =
          SamplerUtils.createNormalizedGaussianSampler(random);
      // Check the vector has a length
      while (length == 0) {
        axis[0] = gauss.sample();
        axis[1] = gauss.sample();
        if (!diffusion2D) {
          axis[2] = gauss.sample();
        }
        length = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
      }
      length = Math.sqrt(length);
      for (int i = 0; i < 3; i++) {
        axis[i] /= length;
      }
    } else {
      diffusionRate = getRandomMoveDistance(compound.getDiffusionRate());
    }

    // For the duration of the sequence, move the molecule and output
    // position and intensity
    final LocalisationModel[][] models = new LocalisationModel[nFluorophores][nSteps];

    if (confinementDistribution != null) {
      confinementDistribution.initialise(compound.getCoordinates());
    }

    for (int t = sequenceStart, step = 0; t <= sequenceEnd; t++, step++) {
      for (int i = 0; i < nFluorophores; i++) {
        if (onTime[i][step] > 0 && photonBudget[i] > 0) {
          final double intensity;
          if (photonBudgetPerFrame) {
            // Ignore the time that the fluorophore is on and simply use the random sample.
            // This allows the simulation to match observed experimental data.
            intensity = new PoissonSampler(getRandom(), photonBudget[i]).sample();
          } else {
            // -=-=-
            // When photon budget is for the entire lifetime then allocate proportionately
            // using the time fraction
            // -=-=-
            intensity =
                new PoissonSampler(getRandom(), photonBudget[i] * onTime[i][step] / totalOnTime[i])
                    .sample();
          }

          // TODO - Is this useful? Add orientation brightness
          // intensity *= randomGenerator.nextUniform(0.7, 1.3);

          // rawPhotons.add(intensity);
          final LocalisationModel model = new LocalisationModel(fluorophores.get(i).getId(), t + 1,
              compound.getCoordinates(i), intensity, state[i][step]);
          model.setLabel(compound.getLabel());
          models[i][step] = model;
        }
      }

      // Move the compound
      if (diffusionRate > 0) {
        if (confinementDistribution != null) {
          // TODO
          // Confined diffusion only asks if the molecule moved into a region that is allowed.
          // A better model would be to ask if the move went out of bounds. If yes then the
          // molecules
          // should be reflected back into the allowed region using the vector of the original move.

          // This only works because the coordinates are a reference
          final double[] xyz = compound.getCoordinates();
          final double[] originalXyz = Arrays.copyOf(xyz, 3);
          for (int n = confinementAttempts; n-- > 0;) {
            diffuse(compound, diffusionRate, axis);
            if (!confinementDistribution.isWithin(compound.getCoordinates())) {
              // Reset position
              for (int j = 0; j < 3; j++) {
                xyz[j] = originalXyz[j];
              }
            } else {
              // The move was allowed
              break;
            }
          }
        } else {
          diffuse(compound, diffusionRate, axis);
        }
      }

      if (rotate) {
        rotate(compound);
      }
    }

    for (int i = 0; i < nFluorophores; i++) {
      linkLocalisations(localisations, models[i]);
    }

    return nFluorophores;
  }

  /**
   * Set the random generator for creating the image.
   *
   * @param random the new random generator
   */
  public void setUniformRandomProvider(UniformRandomProvider random) {
    if (random == null) {
      throw new NullPointerException("Random generator must not be null");
    }
    this.random = random;
  }

  /**
   * Convert diffusion co-efficient (D) into the average step size required for a random diffusion.
   * The step size is per dimension.
   *
   * @param diffusionRateInPixelsPerStep the diffusion rate in pixels per step
   * @return The step size
   * @see MoleculeModel#move(double, UniformRandomProvider)
   * @see MoleculeModel#walk(double, UniformRandomProvider)
   */
  public static double getRandomMoveDistance(double diffusionRateInPixelsPerStep) {
    // Convert diffusion co-efficient into the standard deviation for the random move in each
    // dimension
    // For 1D diffusion: sigma^2 = 2D
    // sigma = sqrt(2D)
    return Math.sqrt(2 * diffusionRateInPixelsPerStep);
  }

  /**
   * Convert diffusion co-efficient (D) into the average step size required for a random diffusion.
   * The step size is for a movement along a random unit vector in the XY plane, i.e. 2 dimensions
   * together.
   *
   * @param diffusionRateInPixelsPerStep the diffusion rate in pixels per step
   * @return The step size
   * @see MoleculeModel#slide(double, double[], UniformRandomProvider)
   */
  public static double getRandomMoveDistance2D(double diffusionRateInPixelsPerStep) {
    // Convert diffusion co-efficient into the standard deviation for the random move
    // For 3D diffusion: sigma^2 = 4D
    // sigma = sqrt(4D)
    return Math.sqrt(4 * diffusionRateInPixelsPerStep);
  }

  /**
   * Convert diffusion co-efficient (D) into the average step size required for a random diffusion.
   * The step size is for a movement along a random unit vector, i.e. all 3 dimensions together.
   *
   * @param diffusionRateInPixelsPerStep the diffusion rate in pixels per step
   * @return The step size
   * @see MoleculeModel#slide(double, double[], UniformRandomProvider)
   */
  public static double getRandomMoveDistance3D(double diffusionRateInPixelsPerStep) {
    // Convert diffusion co-efficient into the standard deviation for the random move
    // For 3D diffusion: sigma^2 = 6D
    // sigma = sqrt(6D)
    return Math.sqrt(6 * diffusionRateInPixelsPerStep);
  }

  /**
   * Gets the photon distribution used for the fluorophore.
   *
   * @return the photon distribution used for the fluorophore.
   */
  public RealDistribution getPhotonDistribution() {
    return photonDistribution;
  }

  /**
   * Set the distribution used to generate the photon budget of a fluorophore.
   *
   * @param photonDistribution the photon distribution to set
   */
  public void setPhotonDistribution(RealDistribution photonDistribution) {
    this.photonDistribution = photonDistribution;
  }

  /**
   * Checks whether to use geometric distribution.
   *
   * @return true, if using a geometric distribution
   */
  public boolean isUseGeometricDistribution() {
    return useGeometricDistribution;
  }

  /**
   * Sets whether to use geometric distribution.
   *
   * @param useGeometricDistribution true, if using a geometric distribution
   */
  public void setUseGeometricDistribution(boolean useGeometricDistribution) {
    this.useGeometricDistribution = useGeometricDistribution;
  }

  /**
   * Checks if the image will be created using the photon budget per frame.
   *
   * @return if true the image will be created using the photon budget per frame.
   */
  public boolean isPhotonBudgetPerFrame() {
    return photonBudgetPerFrame;
  }

  /**
   * Set to true if the photon budget is per frame. The default is for the lifetime of the
   * fluorophore.
   *
   * @param photonBudgetPerFrame if true the image will be created using the photon budget per frame
   */
  public void setPhotonBudgetPerFrame(boolean photonBudgetPerFrame) {
    this.photonBudgetPerFrame = photonBudgetPerFrame;
  }

  /**
   * Set the distribution used to confine any diffusing molecules.
   *
   * @param confinementDistribution the new confinement distribution
   */
  public void setConfinementDistribution(SpatialDistribution confinementDistribution) {
    this.confinementDistribution = confinementDistribution;
  }

  /**
   * Gets the confinement attempts.
   *
   * @return the confinement attempts
   */
  public int getConfinementAttempts() {
    return confinementAttempts;
  }

  /**
   * Set the number of attempts to move a diffusing molecule. If the molecule cannot be successfully
   * moved within the confinement distribution then it remains fixed.
   *
   * @param confinementAttempts the confinementAttempts to set
   */
  public void setConfinementAttempts(int confinementAttempts) {
    this.confinementAttempts = confinementAttempts;
  }

  /**
   * Set to true if only diffusing in XY.
   *
   * @return True if only diffusing in XY
   */
  public boolean isDiffusion2D() {
    return diffusion2D;
  }

  /**
   * Set to true to only diffuse in XY.
   *
   * @param diffusion2d true to only diffuse in XY
   */
  public void setDiffusion2D(boolean diffusion2d) {
    diffusion2D = diffusion2d;
  }

  /**
   * Set to true if rotation is around the z-axis (i.e. rotation in the 2D plane)
   *
   * @return True if using 2D rotation
   */
  public boolean isRotation2D() {
    return rotation2D;
  }

  /**
   * Set to true to rotate around the z-axis (i.e. rotation in the 2D plane)
   *
   * @param rotation2D True if using 2D rotation
   */
  public void setRotation2D(boolean rotation2D) {
    this.rotation2D = rotation2D;
  }
}
