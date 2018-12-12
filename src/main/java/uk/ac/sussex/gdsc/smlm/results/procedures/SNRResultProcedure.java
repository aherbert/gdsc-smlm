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

package uk.ac.sussex.gdsc.smlm.results.procedures;

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Contains functionality to obtain the Signal-to-Noise Ratio (SNR) for results.
 */
//@formatter:off
public class SNRResultProcedure extends AbstractResultProcedure implements
    PeakResultProcedure {
  //@formatter:on

  /** The Signal-to-Noise Ratio (SNR). */
  public float[] snr;

  /**
   * Instantiates a new precision result procedure.
   *
   * @param results the results
   * @throws DataException if the results have no noise
   */
  public SNRResultProcedure(MemoryPeakResults results) {
    super(results);
    if (!results.hasNoise()) {
      throw new DataException("Results do not have noise");
    }
    if (!results.hasMeanIntensity()) {
      throw new DataException("Results do not have mean intensity");
    }
  }

  /**
   * Gets the SNR for the results.
   *
   * <p>The SNR is computed using the mean signal divided by the noise.
   *
   * @return the snr
   */
  public float[] getSNR() {
    counter = 0;
    snr = allocate(snr);
    results.forEach(this);
    return snr;
  }

  @Override
  public void execute(PeakResult peakResult) {
    this.snr[counter++] = peakResult.getSNR();
  }
}
