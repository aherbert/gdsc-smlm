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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.BitSet;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TrackPopulationImporter.CustomLong2IntOpenHashMap;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.IdCategoryPeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class TrackPopulationImporterTest {
  private static final float[] PARAMS = new float[5];

  @SeededTest
  void canEncodeAndDecodeKey(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    for (int i = 0; i < 10; i++) {
      final int frame = rng.nextInt();
      final int id = rng.nextInt();
      final long key = TrackPopulationImporter.getKey(frame, id);
      Assertions.assertEquals(frame, TrackPopulationImporter.getFrame(key));
      Assertions.assertEquals(id, TrackPopulationImporter.getId(key));
    }
  }

  @Test
  void createResultMapIgnoresZeroId() {
    final PeakResult[] data = {create(0, 0), create(0, 0)};
    final Long2IntMap map = TrackPopulationImporter.createResultMap(data);
    Assertions.assertEquals(0, map.size());
  }

  @Test
  void createResultMapDuplicateFrameIdThrows() {
    final PeakResult[] data = {create(0, 1), create(0, 1)};
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> TrackPopulationImporter.createResultMap(data));
  }

  @Test
  void createResultMapCanMapFrameIdToResult() {
    final PeakResult[] data = {create(13, 1), create(13, 2), create(14, 2)};
    final Long2IntMap map = TrackPopulationImporter.createResultMap(data);
    Assertions.assertEquals(3, map.size());
    for (int i = 0; i < data.length; i++) {
      Assertions.assertEquals(i, map.get(getKey(data[i])));
    }
  }

  @Test
  void createCatgoryMapIgnoresZeroId() throws IOException {
    final int[][] data = {{0, 0, 1}, {1, 0, 1}};
    final Long2IntMap map;
    try (BufferedReader br = toReader(data)) {
      map = TrackPopulationImporter.createCategoryMap(br);
    }
    Assertions.assertEquals(0, map.size());
  }

  @Test
  void createCatgoryMapIgnoresZeroCategory() throws IOException {
    final int[][] data = {{0, 1, 0}, {1, 1, 0}};
    final Long2IntMap map;
    try (BufferedReader br = toReader(data)) {
      map = TrackPopulationImporter.createCategoryMap(br);
    }
    Assertions.assertEquals(0, map.size());
  }

  @Test
  void createCatgoryMapWith2FieldsThrows() throws IOException {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    try (BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out))) {
      bw.write(String.format("%d,%d", 42, 13));
      bw.newLine();
    }
    try (BufferedReader br =
        new BufferedReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())))) {
      Assertions.assertThrows(IllegalArgumentException.class,
          () -> TrackPopulationImporter.createCategoryMap(br));
    }
  }

  @Test
  void createCatgoryMapWithNonIntegerFieldsThrows() throws IOException {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    try (BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out))) {
      bw.write(String.format("%d,a,%d", 42, 13));
      bw.newLine();
    }
    try (BufferedReader br =
        new BufferedReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())))) {
      Assertions.assertThrows(IllegalArgumentException.class,
          () -> TrackPopulationImporter.createCategoryMap(br));
    }
  }

  @Test
  void createCatgoryMapCanMapFrameIdToCategory() throws IOException {
    final int[][] data = {{0, 1, 13}, {1, 1, 13}, {2, 1, 7}, {3, 2, 7}};
    final Long2IntMap map;
    try (BufferedReader br = toReader(data)) {
      map = TrackPopulationImporter.createCategoryMap(br);
    }
    Assertions.assertEquals(4, map.size());
    for (int i = 0; i < data.length; i++) {
      Assertions.assertEquals(data[i][2],
          map.get(TrackPopulationImporter.getKey(data[i][0], data[i][1])));
    }
  }

  @Test
  void assignCategoryWithMissingResultKeyThrows() {
    final PeakResult[] data = {create(1, 1), create(2, 1), create(2, 2), create(3, 2)};
    final Long2IntMap resultMap = TrackPopulationImporter.createResultMap(data);
    final CustomLong2IntOpenHashMap categoryMap = new CustomLong2IntOpenHashMap(4);
    categoryMap.put(TrackPopulationImporter.getKey(3, 15), 42);
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> TrackPopulationImporter.assignCategory(data, resultMap, categoryMap));
  }

  @Test
  void assignCategoryCanMarkAssigned() {
    final PeakResult[] data = {create(1, 1), create(2, 1), create(2, 2), create(3, 2)};
    final Long2IntMap resultMap = TrackPopulationImporter.createResultMap(data);
    final CustomLong2IntOpenHashMap categoryMap = new CustomLong2IntOpenHashMap(4);
    categoryMap.put(getKey(data[0]), 1);
    categoryMap.put(getKey(data[2]), 4);
    categoryMap.put(getKey(data[3]), 3);
    final BitSet assigned = TrackPopulationImporter.assignCategory(data, resultMap, categoryMap);
    Assertions.assertEquals(3, assigned.cardinality());
    Assertions.assertTrue(assigned.get(0));
    Assertions.assertFalse(assigned.get(1));
    Assertions.assertTrue(assigned.get(2));
    Assertions.assertTrue(assigned.get(3));
  }

  @Test
  void canChangeCategory() {
    final IdCategoryPeakResult r1 =
        new IdCategoryPeakResult(0, 0, 0, 0, 0, 0, 0, PARAMS, null, 0, 1);
    Assertions.assertSame(r1, TrackPopulationImporter.changeCategory(r1, 1));
    final PeakResult r1b = TrackPopulationImporter.changeCategory(r1, 2);
    Assertions.assertEquals(2, r1b.getCategory());

    final AttributePeakResult r2 = new AttributePeakResult(0, 0, 0);
    r2.setCategory(3);
    Assertions.assertSame(r2, TrackPopulationImporter.changeCategory(r2, 3));
    final PeakResult r2b = TrackPopulationImporter.changeCategory(r2, 42);
    Assertions.assertSame(r2, r2b);
    Assertions.assertEquals(42, r2b.getCategory());
  }

  private static PeakResult create(int frame, int id) {
    return new IdCategoryPeakResult(frame, 0, 0, 0, 0, 0, 0, PARAMS, null, id, 0);
  }

  private static long getKey(PeakResult r) {
    return TrackPopulationImporter.getKey(r.getFrame(), r.getId());
  }

  private static BufferedReader toReader(int[][] data) throws IOException {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    try (BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out))) {
      for (final int[] d : data) {
        bw.write(String.format("%d,%d,%d", d[0], d[1], d[2]));
        bw.newLine();
      }
    }
    return new BufferedReader(new InputStreamReader(new ByteArrayInputStream(out.toByteArray())));
  }
}
