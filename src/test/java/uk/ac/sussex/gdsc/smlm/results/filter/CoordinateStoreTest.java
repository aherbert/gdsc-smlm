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

package uk.ac.sussex.gdsc.smlm.results.filter;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.function.Supplier;

@SuppressWarnings({"javadoc"})
public class CoordinateStoreTest {
  // TODO - test for a crop store ...

  @Test
  public void canCreateStore() {
    CoordinateStore store;
    store = CoordinateStoreFactory.create(1, 2, 10, 11, -1, -1);
    Assertions.assertTrue(store instanceof NullCoordinateStore);
    store = CoordinateStoreFactory.create(1, 2, 10, 11, 0, -1);
    Assertions.assertTrue(store instanceof GridCoordinateStore);
    store = CoordinateStoreFactory.create(1, 2, 10, 11, 0.2, -1);
    Assertions.assertTrue(store instanceof GridCoordinateStore);
    store = CoordinateStoreFactory.create(1, 2, 10, 11, 0.5, -1);
    Assertions.assertTrue(store instanceof GridCoordinateStore1);
    store = CoordinateStoreFactory.create(1, 2, 10, 11, 1, -1);
    Assertions.assertTrue(store instanceof GridCoordinateStore1);
    store = CoordinateStoreFactory.create(1, 2, 10, 11, 1.5, -1);
    Assertions.assertTrue(store instanceof GridCoordinateStore);
    store = CoordinateStoreFactory.create(1, 2, 10, 11, 2, -1);
    Assertions.assertTrue(store instanceof GridCoordinateStore);
  }

  @Test
  public void canDetectXyDuplicates() {
    final double[] datax = {1.1, 4.1};
    final double[] datay = {3.1, 7.1};
    final double[] dataz = {0, 0.1};
    final double[] resolution = {0, 0.3, 0.5, 1.5};

    for (int i = 0; i < resolution.length; i++) {
      final int ii = i;
      final CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], -1);
      for (int j = 0; j < datax.length; j++) {
        s.add(datax[j], datay[j], dataz[j]);
      }

      for (int j = 0; j < datax.length; j++) {
        final int jj = j;
        final Supplier<String> msg = () -> resolution[ii] + " [" + jj + "]";
        //@formatter:off
        Assertions.assertTrue(s.contains(datax[j], datay[j], dataz[j]), msg);
        Assertions.assertTrue(s.contains(datax[j] + resolution[i] * 0.99, datay[j], dataz[j]), msg);
        Assertions.assertTrue(s.contains(datax[j], datay[j] + resolution[i] * 0.99, dataz[j]), msg);
        Assertions.assertFalse(s.contains(datax[j] + increase(resolution[i], 1.01), datay[j] + resolution[i], dataz[j]), msg);
        Assertions.assertFalse(s.contains(datax[j] + resolution[i], datay[j] + increase(resolution[i], 1.01), dataz[j]), msg);
        //@formatter:on
      }
    }
  }

  private static double increase(double value, double delta) {
    return (value == 0) ? delta : value * delta;
  }

  @Test
  public void canDetectZDuplicates() {
    final double x = 3.1;
    final double y = 4.3;
    final double z = 1.1;
    final double[] resolution = {0.3, 0.5, 1.5};

    CoordinateStore store;
    Supplier<String> msg;

    for (int i = 0; i < resolution.length; i++) {
      final int ii = i;

      // No 3D
      final double zResolution_1 = -1;
      store = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], zResolution_1);
      store.add(x, y, z);

      msg = () -> resolution[ii] + "," + zResolution_1;
      Assertions.assertTrue(store.contains(x, y, z), msg);

      // 3D exact match
      final double zResolution0 = 0;
      store = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], zResolution0);
      store.add(x, y, z);

      msg = () -> resolution[ii] + "," + zResolution0;
      Assertions.assertTrue(store.contains(x, y, z), msg);
      Assertions.assertFalse(store.contains(x, y, z + 0.01), msg);
      Assertions.assertFalse(store.contains(x, y, z - 0.01), msg);

      // 3D match within z-resolution
      final double zResolution1 = 1;
      store = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], zResolution1);
      store.add(x, y, z);

      msg = () -> resolution[ii] + "," + zResolution1;
      Assertions.assertTrue(store.contains(x, y, z), msg);
      Assertions.assertTrue(store.contains(x, y, z + zResolution1), msg);
      Assertions.assertTrue(store.contains(x, y, z - zResolution1), msg);
      Assertions.assertFalse(store.contains(x, y, z + zResolution1 * 1.01), msg);
      Assertions.assertFalse(store.contains(x, y, z - zResolution1 * 1.01), msg);
    }
  }

  @Test
  public void canQueueToGrid() {
    final double[] datax = {1.1, 4.1};
    final double[] datay = {3.1, 7.1};
    final double[] dataz = {0, 0.1};
    final double[] resolution = {0.3, 0.5, 1.5};

    for (int i = 0; i < resolution.length; i++) {
      final CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], -1);
      for (int j = 0; j < datax.length; j++) {
        s.addToQueue(datax[j], datay[j], dataz[j]);
        Assertions.assertFalse(s.contains(datax[j], datay[j], dataz[j]));
        for (int k = 0; k < j; k++) {
          Assertions.assertTrue(s.contains(datax[k], datay[k], dataz[j]));
        }
        s.flush();
        for (int k = 0; k <= j; k++) {
          Assertions.assertTrue(s.contains(datax[k], datay[k], dataz[j]));
        }
      }
    }
  }

  @Test
  public void canClearGrid() {
    final double[] datax = {1.1, 4.1};
    final double[] datay = {3.1, 7.1};
    final double[] dataz = {0, 0.1};
    final double[] resolution = {0.3, 0.5, 1.5};

    for (int i = 0; i < resolution.length; i++) {
      final CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, resolution[i], -1);

      // Add then clear
      for (int j = 0; j < datax.length; j++) {
        s.add(datax[j], datay[j], dataz[j]);
        Assertions.assertTrue(s.contains(datax[j], datay[j], dataz[j]));
      }

      s.clear();

      for (int j = 0; j < datax.length; j++) {
        Assertions.assertFalse(s.contains(datax[j], datay[j], dataz[j]));
      }

      // Queue then flush then clear
      for (int j = 0; j < datax.length; j++) {
        s.addToQueue(datax[j], datay[j], dataz[j]);
        Assertions.assertFalse(s.contains(datax[j], datay[j], dataz[j]));
      }
      s.flush();
      for (int j = 0; j < datax.length; j++) {
        Assertions.assertTrue(s.contains(datax[j], datay[j], dataz[j]));
      }

      s.clear();

      for (int j = 0; j < datax.length; j++) {
        Assertions.assertFalse(s.contains(datax[j], datay[j], dataz[j]));
      }
    }
  }

  @Test
  public void cannotAddOutsideGrid1XLow() {
    Assertions.assertThrows(IndexOutOfBoundsException.class, () -> {
      final CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
      s.add(s.getMinX() - 1, s.getMinY(), 0);
    });
  }

  @Test
  public void cannotAddOutsideGridXHigh() {
    Assertions.assertThrows(IndexOutOfBoundsException.class, () -> {
      final CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
      s.add(s.getMinX() + s.getWidth() + 1, s.getMinY(), 0);
    });
  }

  @Test
  public void cannotAddOutsideGrid1YLow() {
    Assertions.assertThrows(IndexOutOfBoundsException.class, () -> {
      final CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
      s.add(s.getMinX(), s.getMinY() - 1, 0);
    });
  }

  @Test
  public void cannotAddOutsideGridYHigh() {
    Assertions.assertThrows(IndexOutOfBoundsException.class, () -> {
      final CoordinateStore s = CoordinateStoreFactory.create(1, 2, 10, 11, 1, 0);
      s.add(s.getMinX(), s.getMinY() + s.getHeight() + 1, 0);
    });
  }

  @Test
  public void canSafeAddOutsideGrid() {
    final GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 1, 0.0);
    s.safeAdd(s.getMinX() - 1, s.getMinY(), 0);
    s.safeAdd(s.getMinX() + s.getWidth() + 1, s.getMinY(), 0);
    s.safeAdd(s.getMinX(), s.getMinY() - 1, 0);
    s.safeAdd(s.getMinX(), s.getMinY() + s.getHeight() + 1, 0);
  }

  @Test
  public void containsOutsideGridIsFalse() {
    final GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 1, 0.0);
    Assertions.assertFalse(addAndFind(s, s.getMinX() - 1, s.getMinY(), 0));
    Assertions.assertFalse(addAndFind(s, s.getMinX() + s.getWidth() + 1, s.getMinY(), 0));
    Assertions.assertFalse(addAndFind(s, s.getMinX(), s.getMinY() - 1, 0));
    Assertions.assertFalse(addAndFind(s, s.getMinX(), s.getMinY() + s.getHeight() + 1, 0));
  }

  private static boolean addAndFind(GridCoordinateStore store, double x, double y, double z) {
    store.safeAdd(x, y, z);
    return store.contains(x, y, z);
  }

  @Test
  public void findOutsideGridIsNull() {
    final GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 1, 0.0);
    s.safeAdd(-1, 0, 0);
    s.safeAdd(11, 0, 0);
    s.safeAdd(-1, 0, 0);
    s.safeAdd(0, 11, 0);
    Assertions.assertNull(s.find(-1, 0, 0));
    Assertions.assertNull(s.find(11, 0, 0));
    Assertions.assertNull(s.find(-1, 0, 0));
    Assertions.assertNull(s.find(0, 11, 0));
  }

  @Test
  public void canChangeXyResolution() {
    final double[] datax = {1.1, 4.1};
    final double[] datay = {3.1, 7.1};
    final double[] dataz = {0, 0.1};
    final double[] resolution = {0, 0.3, 0.5, 1.5};

    final GridCoordinateStore s = new GridCoordinateStore(1, 2, 10, 11, 0, 0.0);
    for (int i = 0; i < resolution.length; i++) {
      final int ii = i;
      s.changeXyResolution(resolution[i]);

      for (int j = 0; j < datax.length; j++) {
        s.add(datax[j], datay[j], dataz[j]);
      }

      for (int j = 0; j < datax.length; j++) {
        final int jj = j;
        final Supplier<String> msg = () -> resolution[ii] + " [" + jj + "]";
        //@formatter:off
        Assertions.assertTrue(s.contains(datax[j], datay[j], dataz[j]), msg);
        Assertions.assertTrue(s.contains(datax[j] + resolution[i] * 0.99, datay[j], dataz[j]), msg);
        Assertions.assertTrue(s.contains(datax[j], datay[j] + resolution[i] * 0.99, dataz[j]), msg);
        Assertions.assertFalse(s.contains(datax[j] + increase(resolution[i], 1.01), datay[j] + resolution[i], dataz[j]), msg);
        Assertions.assertFalse(s.contains(datax[j] + resolution[i], datay[j] + increase(resolution[i], 1.01), dataz[j]), msg);
        //@formatter:on
      }
    }
  }

  @Test
  public void canChangeXyResolutionOnFixedStore() {
    final double[] datax = {1.1, 4.1};
    final double[] datay = {3.1, 7.1};
    final double[] dataz = {0, 0.1};
    final double[] resolution = {0, 0.3, 0.5};

    final GridCoordinateStore1 s = new GridCoordinateStore1(1, 2, 10, 11, 0, 0.0);
    for (int i = 0; i < resolution.length; i++) {
      final int ii = i;
      s.changeXyResolution(resolution[i]);

      for (int j = 0; j < datax.length; j++) {
        s.add(datax[j], datay[j], dataz[j]);
      }

      for (int j = 0; j < datax.length; j++) {
        final int jj = j;
        final Supplier<String> msg = () -> resolution[ii] + " [" + jj + "]";
        //@formatter:off
        Assertions.assertTrue(s.contains(datax[j], datay[j], dataz[j]), msg);
        Assertions.assertTrue(s.contains(datax[j] + resolution[i] * 0.99, datay[j], dataz[j]), msg);
        Assertions.assertTrue(s.contains(datax[j], datay[j] + resolution[i] * 0.99, dataz[j]), msg);
        Assertions.assertFalse(  s.contains(datax[j] + increase(resolution[i], 1.01), datay[j] + resolution[i], dataz[j]), msg);
        Assertions.assertFalse(  s.contains(datax[j] + resolution[i], datay[j] + increase(resolution[i], 1.01), dataz[j]), msg);
        //@formatter:on
      }
    }
  }

  @Test
  public void cannotChangeToBadXyResolutionOnFixedStore() {
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      final GridCoordinateStore1 s = new GridCoordinateStore1(1, 2, 10, 11, 0, 0.0);
      s.changeXyResolution(1.1);
    });
  }
}
