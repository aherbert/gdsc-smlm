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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"deprecation", "javadoc"})
public class CalibrationTest {
  int minpoints = 3;
  int maxpoints = 20;

  @Test
  public void canGet() {
    final Calibration c = new Calibration();
    c.getNmPerPixel();
    c.getGain();
    c.getExposureTime();
    c.getReadNoise();
    c.getBias();
    c.isEmCcd();
    c.getAmplification();
  }

  @Test
  public void canGetWithNoException() {
    final Calibration c = new Calibration(true);
    c.setNmPerPixel(1);
    c.setGain(2);
    c.setExposureTime(3);
    c.setReadNoise(4);
    c.setBias(5);
    c.setEmCcd(true);
    c.setAmplification(6);

    c.getNmPerPixel();
    c.getGain();
    c.getExposureTime();
    c.getReadNoise();
    c.getBias();
    c.isEmCcd();
    c.getAmplification();
  }

  @Test
  public void getNmPerPixelThrowsException() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      new Calibration(true).getNmPerPixel();
    });
  }

  @Test
  public void getGainThrowsException() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      new Calibration(true).getGain();
    });
  }

  @Test
  public void getExposureTimeThrowsException() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      new Calibration(true).getExposureTime();
    });
  }

  @Test
  public void getReadNoiseThrowsException() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      new Calibration(true).getReadNoise();
    });
  }

  @Test
  public void getBiasThrowsException() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      new Calibration(true).getBias();
    });
  }

  @Test
  public void isEmCcdThrowsException() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      new Calibration(true).isEmCcd();
    });
  }

  @Test
  public void getAmplificationThrowsException() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      new Calibration(true).getAmplification();
    });
  }

  @Test
  public void getNmPerPixelThrowsExceptionWhenInvalid() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setNmPerPixel(0);
      c.getNmPerPixel();
    });
  }

  @Test
  public void getGainThrowsExceptionWhenInvalid() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setGain(0);
      c.getGain();
    });
  }

  @Test
  public void getExposureTimeThrowsExceptionWhenInvalid() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setExposureTime(0);
      c.getExposureTime();
    });
  }

  @Test
  public void getReadNoiseThrowsExceptionWhenInvalid() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setReadNoise(-1);
      c.getReadNoise();
    });
  }

  @Test
  public void getBiasThrowsExceptionWhenInvalid() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setBias(-1);
      c.getBias();
    });
  }

  @Test
  public void getAmplificationThrowsExceptionWhenInvalid() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setAmplification(0);
      c.getAmplification();
    });
  }

  @Test
  public void getNmPerPixelThrowsExceptionAfterClear() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setNmPerPixel(1);
      c.clearHasNmPerPixel();
      c.getNmPerPixel();
    });
  }

  @Test
  public void getGainThrowsExceptionAfterClear() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setGain(1);
      c.clearHasGain();
      c.getGain();
    });
  }

  @Test
  public void getExposureTimeThrowsExceptionAfterClear() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setExposureTime(1);
      c.clearHasExposureTime();
      c.getExposureTime();
    });
  }

  @Test
  public void getReadNoiseThrowsExceptionAfterClear() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setReadNoise(1);
      c.clearHasReadNoise();
      c.getReadNoise();
    });
  }

  @Test
  public void getBiasThrowsExceptionAfterClear() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setBias(1);
      c.clearHasBias();
      c.getBias();
    });
  }

  @Test
  public void getEmCcdThrowsExceptionAfterClear() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setEmCcd(true);
      c.clearHasEmCcd();
      c.isEmCcd();
    });
  }

  @Test
  public void getAmplificationThrowsExceptionAfterClear() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setAmplification(1);
      c.clearHasAmplification();
      c.getAmplification();
    });
  }

  @Test
  public void clearDoesNotResetFieldMissingFlag() {
    Assertions.assertThrows(IllegalStateException.class, () -> {
      final Calibration c = new Calibration(true);
      c.setAmplification(0);

      c.clear();
      c.getAmplification();
    });
  }
}
