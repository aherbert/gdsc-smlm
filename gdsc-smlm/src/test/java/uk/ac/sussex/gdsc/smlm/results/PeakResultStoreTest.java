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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Objects;
import java.util.function.BiPredicate;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.smlm.results.sort.FrameIdPeakResultComparator;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class PeakResultStoreTest {
  int capacity = 1;

  /**
   * Subset of specific tests for {@link ArrayPeakResultStore} to ensure coverage. Note: The
   * behaviour may be different from other stores that use a backing collection, for example the
   * methods may accept null for a PeakResult[] argument where a collection may not due to use of
   * Arrays.asList to convert the raw array to a collection.
   */
  static class ArrayPeakResultStoreTest {
    @Test
    void canAddEmptyArray() {
      Assertions.assertThrows(NullPointerException.class,
          () -> new ArrayPeakResultStore((PeakResult[]) null));
      final ArrayPeakResultStore store = new ArrayPeakResultStore(10);
      Assertions.assertFalse(store.addArray(null));
      Assertions.assertFalse(store.addArray(new PeakResult[0]));
      Assertions.assertEquals(0, store.size());
    }

    @SeededTest
    void canRemoveEmptyArray(RandomSeed seed) {
      final ArrayPeakResultStore store = new ArrayPeakResultStore(10);
      store.add(create(RngFactory.create(seed.get())));
      Assertions.assertFalse(store.removeArray(null));
      Assertions.assertFalse(store.removeArray(new PeakResult[0]));
      Assertions.assertEquals(1, store.size());
    }

    @SeededTest
    void canRetainEmptyArray(RandomSeed seed) {
      final ArrayPeakResultStore store = new ArrayPeakResultStore(10);
      final UniformRandomProvider r = RngFactory.create(seed.get());

      store.add(create(r));
      Assertions.assertTrue(store.retainArray(null));
      Assertions.assertEquals(0, store.size());
      Assertions.assertFalse(store.retainArray(null));
      Assertions.assertEquals(0, store.size());

      store.add(create(r));
      Assertions.assertTrue(store.retainArray(new PeakResult[0]));
      Assertions.assertEquals(0, store.size());
      Assertions.assertFalse(store.retainArray(new PeakResult[0]));
      Assertions.assertEquals(0, store.size());
    }

    @SeededTest
    void canRemoveIfWhenPredicateThrows(RandomSeed seed) {
      final ArrayPeakResultStore store = new ArrayPeakResultStore(10);
      final UniformRandomProvider r = RngFactory.create(seed.get());

      final int n = 20;
      for (int i = 0; i < n; i++) {
        store.add(create(r));
      }
      final PeakResultStore store2 = store.copy();
      final int[] count = {n / 2};
      Assertions.assertThrows(IllegalStateException.class, () -> store.removeIf(x -> {
        if (count[0]-- < 0) {
          throw new IllegalStateException();
        }
        // Frames will be random integers. Pick some to remove using the sign.
        if (x.getFrame() < 0) {
          store2.remove(x);
          return true;
        }
        return false;
      }));
      // The store should have removed the correct number and not corrupted the internal state
      // when an exception occurred.
      Assertions.assertArrayEquals(store2.toArray(), store.toArray());
    }
  }

  @SeededTest
  void canStoreResultsUsingArrayList(RandomSeed seed) {
    canStoreResults(seed, new ArrayListPeakResultStore(capacity));
  }

  @SeededTest
  void canStoreResultsUsingArray(RandomSeed seed) {
    canStoreResults(seed, new ArrayPeakResultStore(capacity));
  }

  @SeededTest
  void canStoreResultsUsingHashSet(RandomSeed seed) {
    canStoreResults(seed, new SetPeakResultStore(capacity));
  }

  @SuppressWarnings("null")
  private static void canStoreResults(RandomSeed seed, PeakResultStore store) {
    final boolean isList = store instanceof PeakResultStoreList;
    final PeakResultStoreList storeList = (isList) ? (PeakResultStoreList) store : null;
    PeakResult result;
    final UniformRandomProvider r = RngFactory.create(seed.get());

    PeakResult[] list = new PeakResult[20];
    int size = 0;

    Assertions.assertEquals(size, store.size(), "Not empty on construction");
    Assertions.assertEquals(size, store.toArray().length, "Not empty list");

    // Can store data in order
    for (int i = 0; i < 10; i++) {
      result = create(r);
      list[size++] = result;
      store.add(result);
    }

    assertEquals(list, size, store);

    // Can sort
    if (isList) {
      Arrays.sort(list, 0, size, FrameIdPeakResultComparator.INSTANCE);
      storeList.sort();

      for (int i = 0; i < size; i++) {
        Assertions.assertTrue(list[i] == storeList.get(i), "List entry not same reference");
      }

      final Comparator<PeakResult> c = new Comparator<PeakResult>() {
        @Override
        public int compare(PeakResult o1, PeakResult o2) {
          return o1.getOrigX() - o2.getOrigX();
        }
      };
      Arrays.sort(list, 0, size, c);
      storeList.sort(c);

      for (int i = 0; i < size; i++) {
        Assertions.assertTrue(list[i] == storeList.get(i),
            "List entry not same reference after sort");
      }
      Assertions.assertThrows(IndexOutOfBoundsException.class, () -> storeList.get(-1));
      // This may return stale data instead of throwing as a code optimisation
      // Assertions.assertThrows(IndexOutOfBoundsException.class, () ->
      // storeList.get(store.size()));
    }

    // Can trim to size
    store.trimToSize();
    assertEquals(list, size, store);
    store.trimToSize();
    assertEquals(list, size, store);

    // Can remove null results
    store.add(null);
    result = create(r);
    list[size++] = result;
    store.add(result);
    store.add(null);
    for (int i = 0; i < 3; i++) {
      result = create(r);
      list[size++] = result;
      store.add(result);
    }
    Assertions.assertTrue(size != store.size(), "Same size after adding null results");
    store.removeIf(Objects::isNull);
    assertEquals(list, size, store);
    store.removeIf(Objects::isNull);
    assertEquals(list, size, store);

    // Can clear
    size = 0;
    store.clear();
    assertEquals(list, size, store);

    // Can shuffle
    if (isList) {
      // Test with a random walk of +1/-1 steps
      final PeakResult r0 = new PeakResult(-1, 0, 0, 0);
      final PeakResult r1 = new PeakResult(1, 0, 0, 0);
      final int n = 1024;
      for (int i = n / 2; i-- != 0;) {
        store.add(r0);
      }
      for (int i = n / 2; i-- != 0;) {
        store.add(r1);
      }
      // Fixed seed for stability
      storeList.shuffle(RngFactory.create(123678621384682L));
      // The test is equivalent to the NIST Monobit test with a fixed p-value of 0.01.
      // The reference distribution is Normal with a standard deviation of sqrt(n).
      // Check the absolute position is not too far from the mean of 0 with a fixed
      // p-value of 0.01 taken from a 2-tailed Normal distribution. Computation of
      // the p-value requires the complimentary error function.
      int sum = 0;
      for (int i = n; i-- != 0;) {
        sum += storeList.get(i).getFrame();
      }
      final double absSum = Math.abs(sum);
      final double max = Math.sqrt(n) * 2.5758293035489004;
      Assertions.assertTrue(absSum <= max, () -> "Walked too far astray: " + absSum + " > " + max
          + " (test will fail randomly about 1 in 100 times)");
      store.clear();
    }

    // Can add collection
    for (int i = 0; i < 10; i++) {
      result = create(r);
      list[size++] = result;
    }
    store.addCollection(Arrays.asList(Arrays.copyOf(list, size)));
    assertEquals(list, size, store);

    // Can add array
    size = 0;
    store.clear();
    for (int i = 0; i < 10; i++) {
      result = create(r);
      list[size++] = result;
    }
    Assertions.assertTrue(store.addArray(Arrays.copyOf(list, size)));
    Assertions.assertFalse(store.addArray(new PeakResult[0]));
    assertEquals(list, size, store);
    store.clear();
    store.trimToSize();
    Assertions.assertTrue(store.addArray(Arrays.copyOf(list, size)));
    assertEquals(list, size, store);

    // Can add store
    size = 0;
    store.clear();
    for (int i = 0; i < 10; i++) {
      result = create(r);
      list[size++] = result;
    }
    // Make sure it is not the same type
    final PeakResultStore store2 =
        store instanceof PeakResultStoreCollection ? new ArrayPeakResultStore(10)
            : new ArrayListPeakResultStore(10);
    store2.addArray(Arrays.copyOf(list, size));
    // Ensure the store is not at full capacity
    store2.remove(list[--size]);
    store.addStore(store2);
    assertEquals(list, size, store);

    // Can copy
    final PeakResultStore copy = store.copy(false);
    // Assertions.assertNotEquals(copy, store);
    Assertions.assertTrue(copy != store, "Copy is the same reference");
    assertEquals(list, size, copy);

    PeakResultStore copy2 = store.copy(true);
    Assertions.assertTrue(copy2 != store, "Copy is the same reference");
    assertEquals(list, size, copy2, PeakResultStoreTest::equals);

    // Can add an instance of itself
    store.clear();
    store.addStore(copy);
    assertEquals(list, size, store);

    // Test contains
    store.forEach(x -> {
      Assertions.assertTrue(copy.contains(x));
    });
    Assertions.assertFalse(copy.contains(null));

    // Test index of
    if (isList) {
      size = store.size();
      for (int i = 0; i < size; i++) {
        Assertions.assertEquals(i, storeList.indexOf(storeList.get(i)));
        Assertions.assertEquals(i, storeList.lastIndexOf(storeList.get(i)));
      }
      final PeakResult peakResult = create(r);
      Assertions.assertEquals(-1, storeList.indexOf(peakResult));
      Assertions.assertEquals(-1, storeList.lastIndexOf(peakResult));
      Assertions.assertEquals(-1, storeList.indexOf(null));
      Assertions.assertEquals(-1, storeList.lastIndexOf(null));
      store.add(null);
      Assertions.assertEquals(size, storeList.indexOf(null));
      Assertions.assertEquals(size, storeList.lastIndexOf(null));
      store.removeIf(Objects::isNull);
      Assertions.assertEquals(size, store.size());
    }

    // Can remove single
    for (int j = 0; j < 5; j++) {
      size = 0;
      store.clear();
      PeakResult toRemove = null;
      for (int i = 0; i < 5; i++) {
        result = create(r);
        if (j != i) {
          list[size++] = result;
        } else {
          toRemove = result;
        }
        store.add(result);
      }
      Assertions.assertTrue(size != store.size(), "Same size after adding extra single result");
      Assertions.assertTrue(store.remove(toRemove));
      assertEquals(list, size, store);
      Assertions.assertFalse(store.remove(toRemove));
    }

    if (isList) {
      // Can remove index
      for (int j = 0; j < 5; j++) {
        size = 0;
        store.clear();
        PeakResult toRemove = null;
        for (int i = 0; i < 5; i++) {
          result = create(r);
          if (j != i) {
            list[size++] = result;
          } else {
            toRemove = result;
          }
          store.add(result);
        }
        Assertions.assertTrue(size != store.size(), "Same size after adding extra single result");
        Assertions.assertSame(toRemove, storeList.remove(j));
        assertEquals(list, size, store);
      }
    }

    // Can remove/retain many results
    for (int j = 0; j < 5; j++) {
      size = 0;
      store.clear();
      for (int i = 0; i < 20; i++) {
        result = create(r);
        list[size++] = result;
        store.add(result);
      }

      final int[] indices = RandomUtils.sample(3, size, r);
      Arrays.sort(indices);
      final PeakResult[] list1 = new PeakResult[size - indices.length];
      final PeakResult[] toRemove = new PeakResult[indices.length];
      for (int i = 0; i < indices.length; i++) {
        toRemove[i] = list[indices[i]];
      }
      for (int i = 0, k = 0; i < size; i++) {
        if (Arrays.binarySearch(indices, i) < 0) {
          list1[k++] = list[i];
        }
      }

      final ArrayListPeakResultStore s = new ArrayListPeakResultStore(toRemove.length);
      s.addArray(toRemove);

      // Remove
      copy2 = store.copy();
      copy2.removeCollection(Arrays.asList(toRemove));
      assertEquals(list1, list1.length, copy2);
      copy2 = store.copy();
      copy2.removeArray(toRemove);
      assertEquals(list1, list1.length, copy2);
      copy2 = store.copy();
      copy2.removeStore(new ArrayPeakResultStore(toRemove));
      assertEquals(list1, list1.length, copy2);
      copy2 = store.copy();
      copy2.removeStore(s);
      assertEquals(list1, list1.length, copy2);

      // Retain
      copy2 = store.copy();
      copy2.retainCollection(Arrays.asList(toRemove));
      assertEquals(toRemove, toRemove.length, copy2);
      copy2 = store.copy();
      copy2.retainArray(toRemove);
      assertEquals(toRemove, toRemove.length, copy2);
      copy2 = store.copy();
      copy2.retainStore(new ArrayPeakResultStore(toRemove));
      assertEquals(toRemove, toRemove.length, copy2);
      copy2 = store.copy();
      copy2.retainStore(s);
      assertEquals(toRemove, toRemove.length, copy2);

      // Can subset results
      final PeakResult[] subset = store.subset(x -> {
        for (final PeakResult peakResult : toRemove) {
          if (x == peakResult) {
            return true;
          }
        }
        return false;
      });
      // Use the same store to allow assertion
      store.clear();
      store.addArray(subset);
      assertEquals(toRemove, toRemove.length, store);
    }

    // Can subset results
    for (int j = 0; j < 5; j++) {
      size = 0;
      store.clear();
      for (int i = 0; i < 20; i++) {
        result = create(r);
        list[size++] = result;
        store.add(result);
      }

      final int[] indices = RandomUtils.sample(3, size, r);
      Arrays.sort(indices);
      final PeakResult[] list1 = new PeakResult[size - indices.length];
      final PeakResult[] toRemove = new PeakResult[indices.length];
      for (int i = 0; i < indices.length; i++) {
        toRemove[i] = list[indices[i]];
      }
      for (int i = 0, k = 0; i < size; i++) {
        if (Arrays.binarySearch(indices, i) < 0) {
          list1[k++] = list[i];
        }
      }

      final ArrayListPeakResultStore s = new ArrayListPeakResultStore(toRemove.length);
      s.addArray(toRemove);

      // Remove
      copy2 = store.copy();
      copy2.removeCollection(Arrays.asList(toRemove));
      assertEquals(list1, list1.length, copy2);
      copy2 = store.copy();
      copy2.removeArray(toRemove);
      assertEquals(list1, list1.length, copy2);
      copy2 = store.copy();
      copy2.removeStore(new ArrayPeakResultStore(toRemove));
      assertEquals(list1, list1.length, copy2);
      copy2 = store.copy();
      copy2.removeStore(s);
      assertEquals(list1, list1.length, copy2);

      // Retain
      copy2 = store.copy();
      copy2.retainCollection(Arrays.asList(toRemove));
      assertEquals(toRemove, toRemove.length, copy2);
      copy2 = store.copy();
      copy2.retainArray(toRemove);
      assertEquals(toRemove, toRemove.length, copy2);
      copy2 = store.copy();
      copy2.retainStore(new ArrayPeakResultStore(toRemove));
      assertEquals(toRemove, toRemove.length, copy2);
      copy2 = store.copy();
      copy2.retainStore(s);
      assertEquals(toRemove, toRemove.length, copy2);
    }

    if (store instanceof PeakResultStoreCollection) {
      final PeakResultStoreCollection storeCollection = (PeakResultStoreCollection) store;
      final Collection<PeakResult> col1 = storeCollection.getCollectionReference();
      final Collection<PeakResult> col2 = storeCollection.getCollection();
      Assertions.assertNotSame(col1, col2, "Reference and copy are the same object");
      Assertions.assertEquals(col1, col2, "Reference and copy are not the same collection");
    }

    // Can remove range
    if (!isList) {
      return;
    }

    Assertions.assertThrows(IllegalArgumentException.class, () -> storeList.remove(2, 1));
    Assertions.assertThrows(IndexOutOfBoundsException.class, () -> storeList.remove(-1, -1));
    Assertions.assertThrows(IndexOutOfBoundsException.class,
        () -> storeList.remove(store.size(), store.size()));
    Assertions.assertThrows(IndexOutOfBoundsException.class,
        () -> storeList.remove(0, store.size()));

    size = 0;
    store.clear();
    for (int i = 0; i < 20; i++) {
      result = create(r);
      list[size++] = result;
      store.add(result);
    }

    // From the end
    storeList.remove(size - 1, size - 1);
    assertEquals(list, --size, store);

    storeList.remove(size - 2, size - 1);
    size -= 2;
    assertEquals(list, size, store);

    // From the start
    storeList.remove(0, 0);
    list = Arrays.copyOfRange(list, 1, size);
    size = list.length;
    assertEquals(list, size, store);

    storeList.remove(0, 1);
    list = Arrays.copyOfRange(list, 2, size);
    size = list.length;
    assertEquals(list, size, store);

    // From the middle
    storeList.remove(3, 3);
    System.arraycopy(list, 4, list, 3, size - 4);
    size--;
    assertEquals(list, size, store);

    storeList.remove(3, 4);
    System.arraycopy(list, 5, list, 3, size - 5);
    size -= 2;
    assertEquals(list, size, store);
  }

  private static PeakResult create(UniformRandomProvider rng) {
    return new PeakResult(rng.nextInt(), rng.nextInt(), rng.nextInt(), rng.nextFloat(),
        rng.nextDouble(), rng.nextFloat(), rng.nextFloat(), PeakResult.createParams(rng.nextFloat(),
            rng.nextFloat(), rng.nextFloat(), rng.nextFloat(), rng.nextFloat()),
        null);
  }

  private static void assertEquals(PeakResult[] list, int size, PeakResultStore store) {
    assertEquals(list, size, store, (a, b) -> a == b);
  }

  private static void assertEquals(PeakResult[] list, int size, PeakResultStore store,
      BiPredicate<PeakResult, PeakResult> equality) {
    Assertions.assertEquals(size, store.size(), "Not the same size");

    final boolean isList = store instanceof PeakResultStoreList;
    if (isList) {
      final PeakResultStoreList storeList = (PeakResultStoreList) store;
      for (int i = 0; i < size; i++) {
        Assertions.assertTrue(equality.test(list[i], storeList.get(i)),
            "Not the same list index reference");
        Assertions.assertEquals(i, storeList.indexOf(storeList.get(i)), "indexOf finds wrong item");
      }
    }
    // Set equals
    final PeakResult[] list2 = store.toArray();
    Assertions.assertEquals(size, list2.length, "toArray() creates wrong size");
    for (int i = 0; i < size; i++) {
      Assertions.assertTrue(contains(list, size, list2[i], equality),
          "Cannot find item in the array");
    }
  }

  private static boolean contains(PeakResult[] list, int size, PeakResult peakResult,
      BiPredicate<PeakResult, PeakResult> equality) {
    for (int i = 0; i < size; i++) {
      if (equality.test(list[i], peakResult)) {
        return true;
      }
    }
    return false;
  }

  private static boolean equals(PeakResult r1, PeakResult r2) {
    return r1.getFrame() == r2.getFrame() && r1.getOrigX() == r2.getOrigX()
        && r1.getOrigY() == r2.getOrigY() && r1.getOrigValue() == r2.getOrigValue()
        && r1.getError() == r2.getError() && r1.getNoise() == r2.getNoise()
        && r1.getMeanIntensity() == r2.getMeanIntensity()
        && Arrays.equals(r1.getParameters(), r2.getParameters());
  }
}
