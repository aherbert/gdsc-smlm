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

package uk.ac.sussex.gdsc.smlm.data.config;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
class TestProtosTest {
  @Test
  void canChangeFieldNameAndDeserialise() throws IOException {
    final int value = 35;
    // These two messages have the same type for field 1 but a different name
    final TestProtos.Message1.Builder b1 = TestProtos.Message1.newBuilder();
    final TestProtos.Message2.Builder b2 = TestProtos.Message2.newBuilder();
    b1.setField1(value);
    b2.setEntry1(value);

    final TestProtos.Message1 m1 = b1.build();
    final TestProtos.Message2 m2 = b2.build();

    // Write as message 1
    try (ByteArrayOutputStream output = new ByteArrayOutputStream()) {
      m1.writeDelimitedTo(output);
      output.close();
      final byte[] buf = output.toByteArray();

      // Read as message 2
      try (ByteArrayInputStream input = new ByteArrayInputStream(buf)) {
        final TestProtos.Message2 o2 = TestProtos.Message2.parseDelimitedFrom(input);

        // They should match as the field value is the same
        Assertions.assertEquals(m2, o2);
      }
    }
  }

  @Test
  void cannotChangeFieldTypeAndDeserialise() throws IOException {
    final int value = 35;
    // These two messages have a different type (and byte size) for field 1
    final TestProtos.Message1.Builder b1 = TestProtos.Message1.newBuilder();
    final TestProtos.Message3.Builder b3 = TestProtos.Message3.newBuilder();
    b1.setField1(value);
    // b1.setField2(true);
    b3.setField1(value);
    // b3.setField2(1);

    final TestProtos.Message1 m1 = b1.build();
    final TestProtos.Message3 m3 = b3.build();

    // Write as message 1
    try (ByteArrayOutputStream output = new ByteArrayOutputStream()) {
      m1.writeDelimitedTo(output);
      output.close();
      final byte[] buf = output.toByteArray();

      // Read as message 3
      try (ByteArrayInputStream input = new ByteArrayInputStream(buf)) {
        final TestProtos.Message3 o3 = TestProtos.Message3.parseDelimitedFrom(input);

        // They should not match as the field size has changed
        Assertions.assertNotEquals(m3, o3);
      }
    }
  }
}
