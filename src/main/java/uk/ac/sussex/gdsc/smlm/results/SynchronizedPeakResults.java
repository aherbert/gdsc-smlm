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

import java.awt.Rectangle;
import java.util.Collection;

import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;

/**
 * Wraps a peak results with synchronized methods.
 */
public class SynchronizedPeakResults implements ThreadSafePeakResults {
  private final PeakResults r;
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
    this.r = peakResults;
  }

  /**
   * Creates a PeakResults object that is synchronized if not already a thread-safe instance. <p>
   * The input is unchanged if already a thread-safe instance.
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
   * input results are returned. <p> The input is unchanged if already a thread-safe instance.
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

  //@formatter:off

  @Override
  public void begin()
  {
    synchronized (lock)  { r.begin(); }
  }

  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise, float meanIntensity, float[] params,
      float[] paramsStdDev)
  {
    synchronized (lock)  { r.add(peak, origX, origY, origValue, error, noise, meanIntensity, params, paramsStdDev); }
  }

  @Override
  public void add(PeakResult result)
  {
    synchronized (lock)  { r.add(result); }
  }

  @Override
  public void addAll(Collection<PeakResult> results)
  {
    synchronized (lock)  { r.addAll(results); }
  }

  @Override
  public void addAll(PeakResult[] results)
  {
    synchronized (lock)  { r.addAll(results); }
  }

  @Override
  public void addAll(PeakResultStore results)
  {
    synchronized (lock)  { r.addAll(results); }
  }

  @Override
  public int size()
  {
    synchronized (lock)  { return r.size(); }
  }

  @Override
  public void end()
  {
    synchronized (lock)  { r.end(); }
  }

  @Override
  public boolean isActive()
  {
    synchronized (lock)  { return r.isActive(); }
  }

  @Override
  public void setSource(ImageSource source)
  {
    synchronized (lock)  { r.setSource(source); }
  }

  @Override
  public ImageSource getSource()
  {
    synchronized (lock)  { return r.getSource(); }
  }

  @Override
  public void setBounds(Rectangle bounds)
  {
    synchronized (lock)  { r.setBounds(bounds); }
  }

  @Override
  public Rectangle getBounds()
  {
    synchronized (lock)  { return r.getBounds(); }
  }

  @Override
  public void setCalibration(Calibration calibration)
  {
    synchronized (lock)  { r.setCalibration(calibration); }
  }

  @Override
  public Calibration getCalibration()
  {
    synchronized (lock)  { return r.getCalibration(); }
  }

  @Override
  public void setPSF(PSF psf)
  {
    synchronized (lock)  { r.setPSF(psf); }
  }

  @Override
  public PSF getPSF()
  {
    synchronized (lock)  { return r.getPSF(); }
  }

  @Override
  public void setConfiguration(String configuration)
  {
    synchronized (lock)  { r.setConfiguration(configuration); }
  }

  @Override
  public String getConfiguration()
  {
    synchronized (lock)  { return r.getConfiguration(); }
  }

  @Override
  public String getName()
  {
    synchronized (lock)  { return r.getName(); }
  }

  @Override
  public void setName(String name)
  {
    synchronized (lock)  { r.setName(name); }
  }

  @Override
  public void copySettings(PeakResults peakResults)
  {
    synchronized (lock)  { r.copySettings(peakResults); }
  }

  //@formatter:on
}
