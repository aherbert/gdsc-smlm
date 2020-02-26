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

package uk.ac.sussex.gdsc.smlm.results.procedures;

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Contains functionality to obtain the standard calibrated data for results.
 */
//@formatter:off
public class StandardResultProcedure extends UnitResultProcedure implements
    BResultProcedure,
    BixyResultProcedure,
    BixyzResultProcedure,
    IResultProcedure,
    IxyResultProcedure,
    IxyrResultProcedure,
    IxyzResultProcedure,
    TResultProcedure,
    TxyResultProcedure,
    XyResultProcedure,
    XyrResultProcedure,
    XyzResultProcedure,
    ZResultProcedure {
  //@formatter:on

  /** The frame. */
  public int[] frame;

  /** The background. */
  public float[] background;

  /** The intensity. */
  public float[] intensity;

  /** The x. */
  public float[] x;

  /** The y. */
  public float[] y;

  /** The z. */
  public float[] z;

  /** The peak results. */
  public PeakResult[] peakResults;

  /**
   * Instantiates a new standard result procedure.
   *
   * @param results the results
   * @param distanceUnit the distance unit
   * @param intensityUnit the intensity unit
   */
  public StandardResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit,
      IntensityUnit intensityUnit) {
    super(results, distanceUnit, intensityUnit);
  }

  /**
   * Instantiates a new standard result procedure.
   *
   * @param results the results
   * @param distanceUnit the distance unit
   * @param intensityUnit the intensity unit
   */
  public StandardResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit,
      DistanceUnit distanceUnit) {
    super(results, distanceUnit, intensityUnit);
  }

  /**
   * Instantiates a new standard result procedure.
   *
   * @param results the results
   * @param distanceUnit the distance unit
   */
  public StandardResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit) {
    super(results, distanceUnit);
  }

  /**
   * Instantiates a new standard result procedure.
   *
   * @param results the results
   * @param intensityUnit the intensity unit
   */
  public StandardResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit) {
    super(results, intensityUnit);
  }

  /**
   * Instantiates a new standard result procedure.
   *
   * @param results the results
   */
  public StandardResultProcedure(MemoryPeakResults results) {
    super(results);
  }

  /**
   * Gets the B data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getB() {
    counter = 0;
    allocateB();
    results.forEach(getIntensityUnit(), (BResultProcedure) this);
  }

  @Override
  public void executeB(float background) {
    this.background[counter++] = background;
  }

  /**
   * Gets the BIXY data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getBixy() {
    counter = 0;
    allocateB();
    allocateI();
    allocateX();
    allocateY();
    results.forEach(getIntensityUnit(), getDistanceUnit(), (BixyResultProcedure) this);
  }

  @Override
  public void executeBixy(float background, float intensity, float x, float y) {
    this.background[counter] = background;
    this.intensity[counter] = intensity;
    this.x[counter] = x;
    this.y[counter] = y;
    counter++;
  }

  /**
   * Gets the BIXYZ data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getBixyz() {
    counter = 0;
    allocateB();
    allocateI();
    allocateX();
    allocateY();
    allocateZ();
    results.forEach(getIntensityUnit(), getDistanceUnit(), (BixyzResultProcedure) this);
  }

  @Override
  public void executeBixyz(float background, float intensity, float x, float y, float z) {
    this.background[counter] = background;
    this.intensity[counter] = intensity;
    this.x[counter] = x;
    this.y[counter] = y;
    this.z[counter] = z;
    counter++;
  }

  /**
   * Gets the I data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getI() {
    counter = 0;
    allocateI();
    results.forEach(getIntensityUnit(), (IResultProcedure) this);
  }

  @Override
  public void executeI(float intensity) {
    this.intensity[counter] = intensity;
    counter++;
  }

  /**
   * Gets the IXY data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getIxy() {
    counter = 0;
    allocateI();
    allocateX();
    allocateY();
    results.forEach(getIntensityUnit(), getDistanceUnit(), (IxyResultProcedure) this);
  }

  @Override
  public void executeIxy(float intensity, float x, float y) {
    this.intensity[counter] = intensity;
    this.x[counter] = x;
    this.y[counter] = y;
    counter++;
  }

  /**
   * Gets the IXYR data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getIxyr() {
    counter = 0;
    allocateI();
    allocateX();
    allocateY();
    allocateR();
    results.forEach(getIntensityUnit(), getDistanceUnit(), (IxyrResultProcedure) this);
  }

  @Override
  public void executeIxyr(float intensity, float x, float y, PeakResult result) {
    this.intensity[counter] = intensity;
    this.x[counter] = x;
    this.y[counter] = y;
    peakResults[counter] = result;
    counter++;
  }

  /**
   * Gets the IXYZ data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getIxyz() {
    counter = 0;
    allocateI();
    allocateX();
    allocateY();
    allocateZ();
    results.forEach(getIntensityUnit(), getDistanceUnit(), (IxyzResultProcedure) this);
  }

  @Override
  public void executeIxyz(float intensity, float x, float y, float z) {
    this.intensity[counter] = intensity;
    this.x[counter] = x;
    this.y[counter] = y;
    this.z[counter] = z;
    counter++;
  }

  /**
   * Gets the T data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getT() {
    counter = 0;
    allocateT();
    results.forEach(this);
  }

  @Override
  public void executeT(int frame) {
    this.frame[counter++] = frame;
  }

  /**
   * Gets the TXY data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getTxy() {
    counter = 0;
    allocateT();
    allocateX();
    allocateY();
    results.forEach(getDistanceUnit(), (TxyResultProcedure) this);
  }

  @Override
  public void executeTxy(int frame, float x, float y) {
    this.frame[counter] = frame;
    this.x[counter] = x;
    this.y[counter] = y;
    counter++;
  }

  /**
   * Gets the XY data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getXy() {
    counter = 0;
    allocateX();
    allocateY();
    results.forEach(getDistanceUnit(), (XyResultProcedure) this);
  }

  @Override
  public void executeXy(float x, float y) {
    this.x[counter] = x;
    this.y[counter] = y;
    counter++;
  }

  /**
   * Gets the XYR data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getXyr() {
    counter = 0;
    allocateX();
    allocateY();
    allocateR();
    results.forEach(getDistanceUnit(), (XyrResultProcedure) this);
  }

  @Override
  public void executeXyr(float x, float y, PeakResult result) {
    this.x[counter] = x;
    this.y[counter] = y;
    peakResults[counter] = result;
    counter++;
  }

  /**
   * Gets the XYZ data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getXyz() {
    counter = 0;
    allocateX();
    allocateY();
    allocateZ();
    results.forEach(getDistanceUnit(), (XyzResultProcedure) this);
  }

  @Override
  public void executeXyz(float x, float y, float z) {
    this.x[counter] = x;
    this.y[counter] = y;
    this.z[counter] = z;
    counter++;
  }

  /**
   * Gets the Z data in the configured units.
   *
   * @throws DataException if conversion to the required units is not possible
   */
  public void getZ() {
    counter = 0;
    allocateZ();
    results.forEach(getDistanceUnit(), (ZResultProcedure) this);
  }

  @Override
  public void executeZ(float z) {
    this.z[counter++] = z;
  }

  private void allocateT() {
    this.frame = allocate(this.frame);
  }

  private void allocateB() {
    this.background = allocate(this.background);
  }

  private void allocateI() {
    this.intensity = allocate(this.intensity);
  }

  private void allocateX() {
    this.x = allocate(this.x);
  }

  private void allocateY() {
    this.y = allocate(this.y);
  }

  private void allocateZ() {
    this.z = allocate(this.z);
  }

  private void allocateR() {
    this.peakResults = allocate(this.peakResults);
  }
}
