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
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using a z-coordinate range.
 */
public class ZCoordinateFilter extends DirectFilter {
  // Assuming units are in pixels with 100nm/px set the default range as +/- 1000nm

  /**
   * The default increment. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_INCREMENT = 0.1;
  /**
   * The default range. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double DEFAULT_RANGE = 10;
  /**
   * The default limit. Used for {@link uk.ac.sussex.gdsc.smlm.ga.Chromosome} interface.
   */
  public static final double UPPER_LIMIT = 50; // This may need to be changed

  @XStreamAsAttribute
  private final float minZ;
  @XStreamAsAttribute
  private final float maxZ;

  /**
   * Instantiates a new z coordinate filter.
   *
   * @param minZ the min Z
   * @param maxZ the max Z
   */
  public ZCoordinateFilter(float minZ, float maxZ) {
    if (maxZ < minZ) {
      final float f = maxZ;
      maxZ = minZ;
      minZ = f;
    }
    this.minZ = minZ;
    this.maxZ = maxZ;
  }

  @Override
  protected String generateName() {
    return "Z " + minZ + "-" + maxZ;
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    // Ignore
  }

  @Override
  public boolean accept(PeakResult peak) {
    return peak.getZPosition() >= minZ && peak.getZPosition() <= maxZ;
  }

  @Override
  public int getValidationFlags() {
    return FilterValidationFlag.Z;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    if (peak.getZ() < minZ || peak.getZ() > maxZ) {
      return getValidationFlags();
    }
    return 0;
  }

  @Override
  public String getDescription() {
    return "Filter results using a z-coordinate range.";
  }

  @Override
  public int getNumberOfParameters() {
    return 2;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    return (index == 0) ? minZ : maxZ;
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    return DEFAULT_INCREMENT;
  }

  @Override
  public double getDisabledParameterValue(int index) {
    checkIndex(index);
    return (index == 0) ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    return (index == 0) ? ParameterType.MIN_Z : ParameterType.MAX_Z;
  }

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    return (index == 0) ? new ZCoordinateFilter(updateParameter(minZ, delta, DEFAULT_RANGE), maxZ)
        : new ZCoordinateFilter(minZ, updateParameter(maxZ, delta, DEFAULT_RANGE));
  }

  @Override
  public Filter create(double... parameters) {
    return new ZCoordinateFilter((float) parameters[0], (float) parameters[1]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMin(parameters, 0, minZ);
    setMax(parameters, 1, maxZ);
  }

  @Override
  public int lowerBoundOrientation(int index) {
    return (index == 1) ? 1 : -1;
  }

  @Override
  public int length() {
    return 2;
  }

  @Override
  public double[] sequence() {
    // Ignore the mode parameters
    return new double[] {minZ, maxZ};
  }

  @Override
  public double[] mutationStepRange() {
    return new double[] {DEFAULT_RANGE, DEFAULT_RANGE};
  }
}
