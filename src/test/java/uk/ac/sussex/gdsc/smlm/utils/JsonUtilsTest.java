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

package uk.ac.sussex.gdsc.smlm.utils;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
public class JsonUtilsTest {
  @Test
  public void canSimplify() {
    assertSimplify(null, "");
    assertSimplify("", "");
    assertSimplify("hello", "hello");
    assertSimplify("\"hello\"", "hello");
    assertSimplify("\"hello", "\"hello");
    assertSimplify("\"hello world\"", "\"hello world\"");
    assertSimplify("\"hello.world\"", "\"hello.world\"");
    assertSimplify("\"hello_world\"", "hello_world");
    assertSimplify("\"hello-world\"", "hello-world");
    assertSimplify("\"Say \\\"hello\\\".", "\"Say \\\"hello\\\".");
  }

  private static void assertSimplify(String json, String expected) {
    Assertions.assertEquals(expected, JsonUtils.simplify(json));
  }

  @Test
  public void canPrettyPrint() {
    assertPrettyPrint(null, "");
    assertPrettyPrint("", "");
    assertPrettyPrint("hello", "hello");

    // Real example data:
    //@formatter:off
    // names:[
    //   {
    //     first:"john",
    //     last:"doe"
    //   },
    //   {
    //     first:"bill, ben",
    //     last:"blogs"
    //   }
    // ]
    //@formatter:on

    assertPrettyPrint(
        "names:[{first: \"john\", last: \"doe\"}, {first:\"bill, ben\", last:\"blogs\"}]",
        "names:[\n  {\n    first:\"john\",\n    last:\"doe\"\n  },\n  {\n    first:\"bill, ben\","
            + "\n    last:\"blogs\"\n  }\n]");
  }

  private static void assertPrettyPrint(String json, String expected) {
    Assertions.assertEquals(expected, JsonUtils.prettyPrint(json));
  }
}
