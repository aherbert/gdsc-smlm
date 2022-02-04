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

package uk.ac.sussex.gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.ga.Chromosome;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;
import uk.ac.sussex.gdsc.smlm.results.TraceManager.TraceMode;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Filter results using a precision threshold. Any results below the lower precision limit are
 * included. Any results above the upper precision limit are excluded. Any results between the
 * limits are included only if they can be traced through time, optionally via other candidates, to
 * a valid result.
 */
public abstract class HysteresisFilter extends Filter {
  /** The default absolute distance increment. */
  public static final double DEFAULT_ABSOLUTE_DISTANCE_INCREMENT = 5;
  /** The default relative distance increment. */
  public static final double DEFAULT_RELATIVE_DISTANCE_INCREMENT = 0.05;
  /** The default time increment (in seconds). */
  public static final double DEFAULT_SECONDS_TIME_INCREMENT = 0.05;
  /** The default time increment (in frames). */
  public static final double DEFAULT_FRAMES_TIME_INCREMENT = 1;
  /** The default range for absolute distance. */
  public static final double DEFAULT_ABSOLUTE_DISTANCE_RANGE = 200;
  /** The default range for relative distance. */
  public static final double DEFAULT_RELATIVE_DISTANCE_RANGE = 1;
  /** The default range for time (in seconds). */
  public static final double DEFAULT_SECONDS_TIME_RANGE = 5;
  /** The default range for time (in frames). */
  public static final double DEFAULT_FRAMES_TIME_RANGE = 10;

  /** The search distance. */
  @XStreamAsAttribute
  final double searchDistance;

  /** The search distance mode. */
  @XStreamAsAttribute
  final int searchDistanceMode;

  /** The time threshold. */
  @XStreamAsAttribute
  final double timeThreshold;

  /** The time threshold mode. */
  @XStreamAsAttribute
  final int timeThresholdMode;

  /** The ok. */
  @XStreamOmitField
  Set<PeakResult> ok;

  /**
   * The peak status.
   */
  protected enum PeakStatus {
    /** A valid peak. */
    OK,
    /** A candidate peak. */
    CANDIDATE,
    /** Marked for rejection. */
    REJECT
  }

  /**
   * Instantiates a new hysteresis filter.
   *
   * @param searchDistance the search distance
   * @param searchDistanceMode 0 = relative to the precision of the candidates; 1 = Absolute (in nm)
   * @param timeThreshold the time threshold
   * @param timeThresholdMode 0 = frames; 1 = seconds
   */
  public HysteresisFilter(double searchDistance, int searchDistanceMode, double timeThreshold,
      int timeThresholdMode) {
    this.searchDistance = Math.max(0, searchDistance);
    this.searchDistanceMode = searchDistanceMode;
    this.timeThreshold = Math.max(0, timeThreshold);
    this.timeThresholdMode = timeThresholdMode;
  }

  /**
   * Gets the search name.
   *
   * @return The name of the configured search distance mode
   */
  public String getSearchName() {
    switch (searchDistanceMode) {
      case 1:
        return "Absolute";

      case 0:
      default:
        return "Candidate precision";
    }
  }

  /**
   * Gets the default search range.
   *
   * @return the default search range
   */
  protected double getDefaultSearchRange() {
    switch (searchDistanceMode) {
      case 1:
        return DEFAULT_ABSOLUTE_DISTANCE_RANGE;

      case 0:
      default:
        return DEFAULT_RELATIVE_DISTANCE_RANGE;
    }
  }

  /**
   * Gets the time name.
   *
   * @return The name of the configured time threshold mode
   */
  public String getTimeName() {
    switch (timeThresholdMode) {
      case 1:
        return "Seconds";

      case 0:
      default:
        return "Frames";
    }
  }

  /**
   * Gets the default time range.
   *
   * @return the default time range
   */
  protected double getDefaultTimeRange() {
    switch (timeThresholdMode) {
      case 1:
        return DEFAULT_SECONDS_TIME_RANGE;

      case 0:
      default:
        return DEFAULT_FRAMES_TIME_RANGE;
    }
  }

  /**
   * Gets the trace parameters.
   *
   * @return the trace parameters
   */
  protected String getTraceParameters() {
    return String.format("@%.2f %s, %.2f %s", searchDistance, getSearchName(), timeThreshold,
        getTimeName());
  }

  @Override
  public int getNumberOfParameters() {
    return 4;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    switch (index) {
      case 0:
        return searchDistance;
      case 1:
        return searchDistanceMode;
      case 2:
        return timeThreshold;
      default:
        return timeThresholdMode;
    }
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    switch (index) {
      case 0:
        return (searchDistanceMode == 1) ? DEFAULT_RELATIVE_DISTANCE_INCREMENT
            : DEFAULT_ABSOLUTE_DISTANCE_INCREMENT;
      case 1:
        return 1;
      case 2:
        return (timeThresholdMode == 1) ? DEFAULT_SECONDS_TIME_INCREMENT
            : DEFAULT_FRAMES_TIME_INCREMENT;
      default:
        return 1;
    }
  }

  @Override
  public double getDisabledParameterValue(int index) {
    throw new NotImplementedException("Parameters in hysteresis filters cannot be disabled");
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    switch (index) {
      case 0:
        return ParameterType.DISTANCE_THRESHOLD;
      case 1:
        return ParameterType.DISTANCE_THRESHOLD_MODE;
      case 2:
        return ParameterType.TIME_THRESHOLD;
      default:
        return ParameterType.TIME_THRESHOLD_MODE;
    }
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMax(parameters, 0, searchDistance);
    setMax(parameters, 2, timeThreshold);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    ok = new HashSet<>();

    // Create a set of candidates and valid peaks
    final MemoryPeakResults traceResults = new MemoryPeakResults();

    // Initialise peaks to check
    final LinkedList<PeakResult> candidates = new LinkedList<>();
    peakResults.forEach((PeakResultProcedure) result -> {
      switch (getStatus(result)) {
        case OK:
          ok.add(result);
          traceResults.add(result);
          break;
        case CANDIDATE:
          candidates.add(result);
          traceResults.add(result);
          break;
        default:
          break;
      }
    });

    if (candidates.isEmpty()) {
      // No candidates for tracing so just return
      return;
    }

    double distanceThreshold;
    switch (searchDistanceMode) {
      case 1:
        distanceThreshold = searchDistance / peakResults.getNmPerPixel();
        break;

      case 0:
      default:
        distanceThreshold = getSearchDistanceUsingCandidates(peakResults, candidates);
    }

    if (distanceThreshold <= 0) {
      return;
    }

    // This must be in frames
    int myTimeThreshold;
    if (timeThresholdMode == 1) {
      // time threshold is in Seconds.
      // Default to 1 frame if not calibrated.
      myTimeThreshold = 1;
      if (peakResults.hasCalibration()) {
        // Convert time threshold in seconds to frames
        final CalibrationReader cr = peakResults.getCalibrationReader();
        final double et = cr.getExposureTime();
        if (et > 0) {
          myTimeThreshold = (int) Math.round((this.timeThreshold / et));
        }
      }
    } else {
      // frames
      myTimeThreshold = (int) this.timeThreshold;
    }

    if (myTimeThreshold <= 0) {
      return;
    }

    // Trace through candidates
    final TraceManager tm = new TraceManager(traceResults);
    tm.setTraceMode(TraceMode.LATEST_FORERUNNER);
    tm.traceMolecules(distanceThreshold, myTimeThreshold);
    final Trace[] traces = tm.getTraces();

    for (final Trace trace : traces) {
      if (trace.size() > 1) {
        // Check if the trace touches a valid point
        boolean isOk = false;
        for (int i = 0; i < trace.size(); i++) {
          if (ok.contains(trace.get(i))) {
            isOk = true;
            break;
          }
        }
        // Add the entire trace to the OK points
        if (isOk) {
          for (int i = 0; i < trace.size(); i++) {
            ok.add(trace.get(i));
          }
        }
      }
    }
  }

  /**
   * Find average precision of the candidates and use it for the search distance.
   *
   * @param peakResults the peak results
   * @param candidates the candidates
   * @return the search distance using candidates
   */
  private double getSearchDistanceUsingCandidates(MemoryPeakResults peakResults,
      LinkedList<PeakResult> candidates) {
    final Gaussian2DPeakResultCalculator calculator =
        Gaussian2DPeakResultHelper.create(peakResults.getPsf(), peakResults.getCalibration(),
            Gaussian2DPeakResultHelper.LSE_PRECISION);
    double sum = 0;
    for (final PeakResult peakResult : candidates) {
      sum += calculator.getLsePrecision(peakResult.getParameters(), peakResult.getNoise());
    }
    final double nmPerPixel = peakResults.getNmPerPixel();
    return (sum / candidates.size()) * searchDistance / nmPerPixel;
  }

  /**
   * Gets the status.
   *
   * @param result the result
   * @return the status
   */
  protected abstract PeakStatus getStatus(PeakResult result);

  /**
   * {@inheritDoc}
   *
   * @throws NullPointerException if not first initialised with a call to
   *         {@link #setup(MemoryPeakResults)}
   */
  @Override
  public boolean accept(PeakResult peak) {
    return ok.contains(peak);
  }

  @Override
  public void end() {
    ok.clear();
    ok = null;
  }

  @Override
  public String getDescription() {
    return "Any results between the limits (candidates) are included only if they can be traced "
        + "through time, potentially via other candidates, to a valid result. The distance used "
        + "for tracing is the search distance multiplied by the average precision of the "
        + "candidates.";
  }

  @Override
  public boolean subsetWithFailCount() {
    return false;
  }

  @Override
  public Chromosome<FilterScore> newChromosome(double[] sequence) {
    // Hysteresis filters remove their search and time mode parameters in their Chromosome sequence
    // so add it back
    final double[] parameters = new double[sequence.length];
    parameters[0] = sequence[0];
    parameters[1] = searchDistanceMode;
    parameters[2] = sequence[1];
    parameters[3] = timeThresholdMode;
    System.arraycopy(sequence, 2, parameters, 4, sequence.length - 2);
    return create(parameters);
  }

  @Override
  public int[] getChromosomeParameters() {
    // Hysteresis filters remove their search and time mode parameters in their Chromosome sequence
    // Skip the search mode [param 1]
    // Skip the time mode [param 3]
    final int[] indices = new int[length()];
    indices[1] = 2;
    for (int i = 2; i < indices.length; i++) {
      indices[i] = i + 2;
    }
    return indices;
  }

  @Override
  public int length() {
    // Hysteresis filters remove their search and time mode parameters in their Chromosome sequence
    return getNumberOfParameters() - 2;
  }

  @Override
  public double[] sequence() {
    // Implement a default version using the results of getChromosomeParameters() and
    // getParameters().
    final double[] sequence = new double[length()];
    final double[] params = getParameters();
    final int[] indices = getChromosomeParameters();
    for (int i = 0; i < indices.length; i++) {
      sequence[i] = params[indices[i]];
    }
    return sequence;
  }
}
