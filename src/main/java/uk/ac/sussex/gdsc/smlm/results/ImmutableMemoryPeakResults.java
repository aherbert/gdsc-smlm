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

package uk.ac.sussex.gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.Collection;
import java.util.function.Predicate;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;

/**
 * Wraps peak results in memory and prevents modification of the results size.
 *
 * <p>Any method that modifies the size of the results set will throw a data exception.
 */
public class ImmutableMemoryPeakResults extends MemoryPeakResults {
  private static final String IMMUTABLE_MESSAGE = "This results set is immutable";

  private final boolean built;

  /**
   * Instantiates a new immutable memory peak results with the original results store.
   *
   * @param results the results
   */
  public ImmutableMemoryPeakResults(MemoryPeakResults results) {
    this(results, false);
  }

  /**
   * Instantiates a new immutable memory peak results with an optional copy of the results store.
   *
   * @param results the results
   * @param copy Set to true to copy the original results store
   */
  public ImmutableMemoryPeakResults(MemoryPeakResults results, boolean copy) {
    super((PeakResultStoreList) ((copy) ? results.results.copy() : results.results));
    copySettings(results);
    built = true;
  }

  @Override
  PeakResult getfX(int index) {
    return new ImmutablePeakResult(super.getfX(index));
  }

  @Override
  public void setSource(ImageSource source) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setSource(source);
  }

  @Override
  public void setBounds(Rectangle bounds) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setBounds(bounds);
  }

  @Override
  public Rectangle getBounds() {
    // Prevent modification
    return new Rectangle(super.getBounds());
  }

  @Override
  public void setCalibration(Calibration calibration) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setCalibration(calibration);
  }

  @Override
  public void setPsf(PSF psf) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setPsf(psf);
  }

  @Override
  public void setConfiguration(String configuration) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setConfiguration(configuration);
  }

  @Override
  public void setName(String name) {
    if (built) {
      throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
    }
    super.setName(name);
  }

  @Override
  public void add(PeakResult result) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void add(MemoryPeakResults results) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void addAll(Collection<PeakResult> results) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void addAll(PeakResult[] results) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void addAll(PeakResultStore results) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void removeNullResults() {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public boolean removeIf(Predicate<PeakResult> filter) {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void begin() {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public void end() {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public boolean isActive() {
    throw new UnsupportedOperationException(IMMUTABLE_MESSAGE);
  }

  @Override
  public PeakResultView getPeakResultView() {
    return new DynamicPeakResultView(new ImmutablePeakResultStore(results));
  }

  @Override
  public PeakResultView getSnapshotPeakResultView() {
    return new CachedPeakResultView(new ImmutablePeakResultStore(results));
  }
}
