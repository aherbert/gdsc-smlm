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

package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift, precision and
 * z-depth. Calculates the precision using the fitted X/Y parameter deviations.
 */
public class MultiFilterCRLB extends MultiFilter implements IMultiFilter {

  /**
   * Instantiates a new multi filter CRLB.
   *
   * @param signal the signal
   * @param snr the snr
   * @param minWidth the min width
   * @param maxWidth the max width
   * @param shift the shift
   * @param eshift the eshift
   * @param precision the precision
   * @param minZ the min Z
   * @param maxZ the max Z
   */
  public MultiFilterCRLB(double signal, float snr, double minWidth, double maxWidth, double shift,
      double eshift, double precision, float minZ, float maxZ) {
    super(signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
  }

  @Override
  protected String generateName() {
    return String.format(
        "MultiCRLB: Signal=%.1f, SNR=%.1f, Width=%.2f-%.2f, Shift=%.2f, EShift=%.2f, "
            + "Precision=%.1f, Width=%.2f-%.2f",
        signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
  }

  @Override
  protected void setupCalculator(MemoryPeakResults peakResults) {
    calculator =
        Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(), 0);
  }

  @Override
  protected MultiFilterComponent createPrecisionComponent() {
    return new MultiFilterVarianceCRLBComponent(precision);
  }

  @Override
  protected double getVariance(PeakResult peak) {
    if (peak.hasParameterDeviations()) {
      final float vx = peak.getParameterDeviation(PeakResult.X);
      final float vy = peak.getParameterDeviation(PeakResult.Y);
      return (vx * vx + vy * vy);
    }
    return variance; // Return the current limit
  }

  /** {@inheritDoc} */
  @Override
  public String getDescription() {
    return "Filter results using multiple thresholds: Signal, SNR, width, shift, "
        + "Euclidian shift, precision (uses fitted parameter variance) and Z-depth";
  }

  /** {@inheritDoc} */
  @Override
  public boolean requiresParameterDeviations() {
    return (precision != 0);
  }

  @Override
  protected ParameterType getPrecisionParamaterType() {
    return ParameterType.PRECISION_CRLB;
  }

  /** {@inheritDoc} */
  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    final double[] params =
        new double[] {signal, snr, minWidth, maxWidth, shift, eshift, precision};
    params[index] = updateParameter(params[index], delta, MultiFilter.defaultRange[index]);
    return new MultiFilterCRLB(params[0], (float) params[1], params[2], params[3], params[4],
        params[5], params[6], (float) params[7], (float) params[8]);
  }

  /** {@inheritDoc} */
  @Override
  public Filter create(double... parameters) {
    return new MultiFilterCRLB(parameters[0], (float) parameters[1], parameters[2], parameters[3],
        parameters[4], parameters[5], parameters[6], (float) parameters[7], (float) parameters[8]);
  }

  @Override
  public PrecisionType getPrecisionType() {
    return PrecisionType.CRLB;
  }
}
