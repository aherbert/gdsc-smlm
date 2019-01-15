/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.smlm.data.NamedObject;

/**
 * Define the type of parameter.
 */
public enum ParameterType implements NamedObject {
  //@formatter:off
  /** Signal. */
    SIGNAL("Signal"),
    /** Signal Range. */
    SIGNAL_RANGE("Signal Range"),
    /** Signal-to-Noise Ratio */
    SNR("Signal-to-Noise Ratio","SNR"),
    /** Signal-to-Noise Ratio Range */
    SNR_RANGE("Signal-to-Noise Ratio Range","SNR Range"),
    /** Min Width. */
    MIN_WIDTH("Min Width"),
    /** Min Width Range. */
    MIN_WIDTH_RANGE("Min Width Range"),
    /** Max Width. */
    MAX_WIDTH("Max Width"),
    /** Max Width Range. */
    MAX_WIDTH_RANGE("Max Width Range"),
    /** Shift. */
    SHIFT("Shift"),
    /** Shift Range. */
    SHIFT_RANGE("Shift Range"),
    /** Euclidian Shift. */
    ESHIFT("Euclidian Shift","EShift"),
    /** Precision. */
    PRECISION("Precision"),
    /** Precision Range. */
    PRECISION_RANGE("Precision Range"),
    /** Precision using local background. */
    PRECISION2("Precision using local background","Precision2"),
    /** Precision using local background Range. */
    PRECISION2_RANGE("Precision using local background Range","Precision2 Range"),
    /** Precision using Cramér–Rao Lower Bounds (CRLB) */
    PRECISION_CRLB("Precision using CRLB","Precision CRLB"),
    /** Precision using Cramér–Rao Lower Bounds (CRLB) Range */
    PRECISION_CRLB_RANGE("Precision using CRLB Range","Precision CRLB Range"),
    /** Amplitude-to-Noise Ratio */
    ANR("Amplitude-to-Noise Ratio","ANR"),
    /** Signal-to-Background Ratio */
    SBR("Signal-to-Background Ratio","SBR"),
    /** Distance Threshold. */
    DISTANCE_THRESHOLD("Distance Threshold","D-threshold"),
    /** Distance Threshold Mode. */
    DISTANCE_THRESHOLD_MODE("Distance Threshold Mode","D-threshold mode"),
    /** Time Threshold. */
    TIME_THRESHOLD("Time Threshold","T-threshold"),
    /** Time Threshold Mode. */
    TIME_THRESHOLD_MODE("Time Threshold Mode","T-threshold mode"),
    /** Min X. */
    MIN_X("Min X"),
    /** Max X. */
    MAX_X("Max X"),
    /** Min Y. */
    MIN_Y("Min Y"),
    /** Max Y. */
    MAX_Y("Max Y"),
    /** Min Z. */
    MIN_Z("Min Z"),
    /** Max Z. */
    MAX_Z("Max Z"),
  ;

  //@formatter:on

  private final String name;
  private final String shortName;

  private ParameterType(String name) {
    this(name, name);
  }

  private ParameterType(String name, String sname) {
    this.name = name;
    this.shortName = sname;
  }

  @Override
  public String toString() {
    return shortName;
  }

  @Override
  public String getName() {
    return name;
  }

  @Override
  public String getShortName() {
    return shortName;
  }
}
