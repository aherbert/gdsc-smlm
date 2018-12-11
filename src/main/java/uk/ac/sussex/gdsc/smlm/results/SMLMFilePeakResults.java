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

package uk.ac.sussex.gdsc.smlm.results;

/**
 * Saves the fit results to file.
 */
public abstract class SMLMFilePeakResults extends FilePeakResults {
  /** The flag for the end frame. Used in the file format version string. */
  public static final int FLAG_END_FRAME = 0x0001;
  /** The flag for the id. Used in the file format version string. */
  public static final int FLAG_ID = 0x0002;
  /** The flag for the precision. Used in the file format version string. */
  public static final int FLAG_PRECISION = 0x0004;

  /**
   * The version.
   *
   * <ul> <li>V1 = Version 1 had signal and amplitude in the results. It did not have the version
   * string. <li>V2 = Version 2 has only signal in the results <li>V3 = Version 3 has an improved
   * calibration header with dynamic fields <li>V4 = Version 4 added mean intensity field to the
   * standard data </ul>
   */
  public static final int VERSION = 4;

  private final boolean showDeviations;
  private final boolean showEndFrame;
  private final boolean showId;
  private final boolean showPrecision;

  /** The peak id column name. */
  protected String peakIdColumnName = "Frame";

  /**
   * Instantiates a new SMLM file peak results.
   *
   * @param filename the filename
   */
  public SMLMFilePeakResults(String filename) {
    // showDeviations=true
    this(filename, true);
  }

  /**
   * Instantiates a new SMLM file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   */
  public SMLMFilePeakResults(String filename, boolean showDeviations) {
    this(filename, showDeviations, false);
  }

  /**
   * Instantiates a new SMLM file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   */
  public SMLMFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame) {
    this(filename, showDeviations, showEndFrame, false, false);
  }

  /**
   * Instantiates a new SMLM file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   */
  public SMLMFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId) {
    this(filename, showDeviations, showEndFrame, showId, false);
  }

  /**
   * Instantiates a new SMLM file peak results.
   *
   * @param filename the filename
   * @param showDeviations Set to true to show deviations
   * @param showEndFrame Set to true to show the end frame
   * @param showId Set to true to show the id
   * @param showPrecision Set to true to show the precision
   */
  public SMLMFilePeakResults(String filename, boolean showDeviations, boolean showEndFrame,
      boolean showId, boolean showPrecision) {
    super(filename);
    this.showDeviations = showDeviations;
    this.showEndFrame = showEndFrame;
    this.showId = showId;
    this.showPrecision = showPrecision;
  }

  /**
   * Check size of the parameter arrays. The array size must match the number of parameters.
   *
   * @param numberOfParams the number of params
   * @param params the params
   * @param paramsStdDevs the params std devs
   * @throws IllegalArgumentException if the parameter arrays are the wrong size
   */
  protected static void checkSize(int numberOfParams, float[] params, float[] paramsStdDevs)
      throws IllegalArgumentException {
    checkSize(numberOfParams, params);
    checkSize(numberOfParams, paramsStdDevs);
  }

  /**
   * Check size of the parameter arrays. The array size must match the number of parameters.
   *
   * @param numberOfParams the number of params
   * @param a the a
   * @throws IllegalArgumentException if the parameter array is the wrong size
   */
  protected static void checkSize(int numberOfParams, float[] a) throws IllegalArgumentException {
    if (a.length < numberOfParams) {
      throw new IllegalArgumentException("Incorrect number of parameters " + a.length);
    }
  }

  /**
   * Gets the file format version.
   *
   * @return A line containing the file format version
   */
  @Override
  protected String getVersion() {
    final StringBuilder sb = new StringBuilder();
    sb.append(isBinary() ? "Binary" : "Text");
    sb.append(".");
    sb.append(isShowDeviations() ? "D1" : "D0");
    sb.append(".E");
    int extended = 0;
    if (isShowEndFrame()) {
      extended += FLAG_END_FRAME;
    }
    if (isShowId()) {
      extended += FLAG_ID;
    }
    if (isShowPrecision()) {
      extended += FLAG_PRECISION;
    }
    sb.append(extended);
    sb.append(".V").append(VERSION);
    return sb.toString();
  }

  /**
   * Gets the name of the peak column.
   *
   * @return the name of the peak column
   */
  public String getPeakIdColumnName() {
    return peakIdColumnName;
  }

  /**
   * Sets the name of the peak column.
   *
   * @param peakIdColumnName the name of the peak column
   */
  public void setPeakIdColumnName(String peakIdColumnName) {
    this.peakIdColumnName = peakIdColumnName;
  }

  /**
   * Checks if the records contain the parameter deviations.
   *
   * @return True if the records contain the parameter deviations
   */
  public boolean isShowDeviations() {
    return showDeviations;
  }

  /**
   * Checks if the records contain the result end frame.
   *
   * @return True if the records contain the result end frame
   */
  public boolean isShowEndFrame() {
    return showEndFrame;
  }

  /**
   * Checks if the records contain a result Id.
   *
   * @return True if the records contain a result Id
   */
  public boolean isShowId() {
    return showId;
  }

  /**
   * Checks if the records contain the localisation precision.
   *
   * @return True if the records contain the localisation precision
   */
  public boolean isShowPrecision() {
    return showPrecision;
  }
}
