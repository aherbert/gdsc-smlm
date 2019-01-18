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

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;

import java.awt.Rectangle;
import java.util.Collection;

/**
 * Wraps a peak results with synchronized methods.
 */
public class SynchronizedPeakResults implements ThreadSafePeakResults {
  private final PeakResults peakResults;
  private final Object lock = new Object();

  /**
   * Instantiates a new synchronized peak results.
   *
   * @param peakResults the peak results
   * @throws IllegalArgumentException if the results are null
   */
  public SynchronizedPeakResults(PeakResults peakResults) {
    if (peakResults == null) {
      throw new IllegalArgumentException("PeakResults must not be null");
    }
    this.peakResults = peakResults;
  }

  /**
   * Creates a PeakResults object that is synchronized if not already a thread-safe instance.
   *
   * <p>The input is unchanged if already a thread-safe instance.
   *
   * @param peakResults the peak results
   * @return the peak results
   * @throws IllegalArgumentException if the results are null
   */
  public static PeakResults create(PeakResults peakResults) {
    if (peakResults instanceof ThreadSafePeakResults) {
      return peakResults;
    }
    return new SynchronizedPeakResults(peakResults);
  }

  /**
   * Creates a PeakResults object that is synchronized if the thread count is above 1, otherwise the
   * input results are returned.
   *
   * <p>The input is unchanged if already a thread-safe instance.
   *
   * @param peakResults the peak results
   * @param threadCount the thread count
   * @return the peak results
   * @throws IllegalArgumentException if the results are null
   */
  public static PeakResults create(PeakResults peakResults, int threadCount) {
    if (threadCount <= 1 || peakResults instanceof ThreadSafePeakResults) {
      return peakResults;
    }
    return new SynchronizedPeakResults(peakResults);
  }

  @Override
  public void begin() {
    synchronized (lock) {
      peakResults.begin();
    }
  }

  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    synchronized (lock) {
      peakResults.add(peak, origX, origY, origValue, error, noise, meanIntensity, params,
          paramsStdDev);
    }
  }

  @Override
  public void add(PeakResult result) {
    synchronized (lock) {
      peakResults.add(result);
    }
  }

  @Override
  public void addAll(Collection<PeakResult> results) {
    synchronized (lock) {
      peakResults.addAll(results);
    }
  }

  @Override
  public void addAll(PeakResult[] results) {
    synchronized (lock) {
      peakResults.addAll(results);
    }
  }

  @Override
  public void addAll(PeakResultStore results) {
    synchronized (lock) {
      peakResults.addAll(results);
    }
  }

  @Override
  public int size() {
    synchronized (lock) {
      return peakResults.size();
    }
  }

  @Override
  public void end() {
    synchronized (lock) {
      peakResults.end();
    }
  }

  @Override
  public boolean isActive() {
    synchronized (lock) {
      return peakResults.isActive();
    }
  }

  @Override
  public void setSource(ImageSource source) {
    synchronized (lock) {
      peakResults.setSource(source);
    }
  }

  @Override
  public ImageSource getSource() {
    synchronized (lock) {
      return peakResults.getSource();
    }
  }

  @Override
  public void setBounds(Rectangle bounds) {
    synchronized (lock) {
      peakResults.setBounds(bounds);
    }
  }

  @Override
  public Rectangle getBounds() {
    synchronized (lock) {
      return peakResults.getBounds();
    }
  }

  @Override
  public void setCalibration(Calibration calibration) {
    synchronized (lock) {
      peakResults.setCalibration(calibration);
    }
  }

  @Override
  public Calibration getCalibration() {
    synchronized (lock) {
      return peakResults.getCalibration();
    }
  }

  @Override
  public void setPsf(PSF psf) {
    synchronized (lock) {
      peakResults.setPsf(psf);
    }
  }

  @Override
  public PSF getPsf() {
    synchronized (lock) {
      return peakResults.getPsf();
    }
  }

  @Override
  public void setConfiguration(String configuration) {
    synchronized (lock) {
      peakResults.setConfiguration(configuration);
    }
  }

  @Override
  public String getConfiguration() {
    synchronized (lock) {
      return peakResults.getConfiguration();
    }
  }

  @Override
  public String getName() {
    synchronized (lock) {
      return peakResults.getName();
    }
  }

  @Override
  public void setName(String name) {
    synchronized (lock) {
      peakResults.setName(name);
    }
  }

  @Override
  public void copySettings(PeakResults peakResults) {
    synchronized (lock) {
      peakResults.copySettings(peakResults);
    }
  }
}
