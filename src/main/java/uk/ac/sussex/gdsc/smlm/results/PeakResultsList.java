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

package uk.ac.sussex.gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;

/**
 * Wrapper class to output to multiple results destinations.
 */
public class PeakResultsList extends AbstractPeakResults {
  private final List<PeakResults> results = new ArrayList<>();

  /**
   * Add a result format to the output. If a PeakResultsList is passed then it will be separated
   * into the child PeakResults instances. This will break the size() function of any input
   * PeakResultsList since only the children will remain within this list.
   *
   * <p>Sets the settings (source and configuration) of the child to the same as this list
   *
   * @param peakResults the peak results
   */
  public void addOutput(PeakResults peakResults) {
    if (peakResults instanceof PeakResultsList) {
      for (final PeakResults r : ((PeakResultsList) peakResults).results) {
        addOutput(r);
      }
    } else {
      peakResults.copySettings(this);
      results.add(peakResults);
    }
  }

  /**
   * Get the number of outputs contained in the list.
   *
   * @return The number of outputs contained in the list.
   */
  public int numberOfOutputs() {
    return results.size();
  }

  /**
   * Gets the output.
   *
   * @param index the index
   * @return the output
   */
  public PeakResults getOutput(int index) {
    return results.get(index);
  }

  /**
   * Convert the list of outputs to an array.
   *
   * @return The outputs.
   */
  public PeakResults[] toArray() {
    return results.toArray(new PeakResults[0]);
  }

  @Override
  public void begin() {
    for (final PeakResults peakResults : results) {
      peakResults.begin();
    }
  }

  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsStdDev) {
    for (final PeakResults peakResults : results) {
      peakResults.add(peak, origX, origY, origValue, error, noise, meanIntensity, params,
          paramsStdDev);
    }
  }

  @Override
  public void add(PeakResult result) {
    for (final PeakResults peakResults : results) {
      peakResults.add(result);
    }
  }

  @Override
  public void addAll(PeakResult[] results) {
    for (final PeakResults peakResults : this.results) {
      peakResults.addAll(results);
    }
  }

  @Override
  public void addAll(Collection<PeakResult> results) {
    for (final PeakResults peakResults : this.results) {
      peakResults.addAll(results);
    }
  }

  @Override
  public int size() {
    return (results.isEmpty()) ? 0 : results.get(0).size();
  }

  @Override
  public void end() {
    for (final PeakResults peakResults : results) {
      peakResults.end();
    }
  }

  @Override
  public boolean isActive() {
    for (final PeakResults peakResults : this.results) {
      if (peakResults.isActive()) {
        return true;
      }
    }
    return false;
  }

  /**
   * Checks all the results in the list. If any are not thread safe then they are wrapped with a
   * SynchronizedPeakResults container.
   *
   * @return the thread safe list
   */
  public PeakResultsList getThreadSafeList() {
    final PeakResultsList newList = new PeakResultsList();
    newList.copySettings(this);
    for (final PeakResults peakResults : this.results) {
      // This assumes the settings are OK, i.e. the result was added
      // using addOutput(...).
      newList.results.add(SynchronizedPeakResults.create(peakResults));
    }
    return newList;
  }

  // Pass through all the modifications to the list objects as well as this
  // so that any new list objects can copy the settings.

  @Override
  public void setSource(ImageSource source) {
    super.setSource(source);
    for (final PeakResults peakResults : results) {
      peakResults.setSource(source);
    }
  }

  @Override
  public void setBounds(Rectangle bounds) {
    super.setBounds(bounds);
    for (final PeakResults peakResults : results) {
      peakResults.setBounds(bounds);
    }
  }

  @Override
  public void setCalibration(Calibration calibration) {
    super.setCalibration(calibration);
    for (final PeakResults peakResults : results) {
      peakResults.setCalibration(calibration);
    }
  }

  @Override
  public void setPsf(PSF psf) {
    super.setPsf(psf);
    for (final PeakResults peakResults : results) {
      peakResults.setPsf(psf);
    }
  }

  @Override
  public void setConfiguration(String configuration) {
    super.setConfiguration(configuration);
    for (final PeakResults peakResults : results) {
      peakResults.setConfiguration(configuration);
    }
  }

  @Override
  public void setName(String name) {
    super.setName(name);
    for (final PeakResults peakResults : results) {
      peakResults.setName(name);
    }
  }

  @Override
  public void copySettings(PeakResults results) {
    super.copySettings(results);
    for (final PeakResults peakResults : this.results) {
      peakResults.copySettings(results);
    }
  }
}
