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

import java.io.IOException;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;
import java.util.Set;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;
import org.opentest4j.TestAbortedException;

@SuppressWarnings({"javadoc"})
class HelpUrlsTest {
  @Test
  void canLoadMissingUrls() {
    Assertions.assertTrue(HelpUrls.loadUrls(null).isEmpty());
    Assertions.assertTrue(HelpUrls.loadUrls("").isEmpty());
    Assertions.assertTrue(HelpUrls.loadUrls("  ").isEmpty());
    Assertions.assertTrue(HelpUrls.loadUrls("path/does/not/exist").isEmpty());
  }

  @Test
  void canLoadTestUrls() {
    final Properties props = HelpUrls.loadUrls("/uk/ac/sussex/gdsc/smlm/ij/test.help.config");
    Assertions.assertEquals(1, props.size());
    Assertions.assertEquals("path/to/docs/page.html", props.get("testplugin"));

    // url with suffix
    String url = HelpUrls.getUrl("testplugin", props);
    Assertions.assertTrue(url.endsWith("path/to/docs/page.html"));
    Assertions.assertTrue(url.startsWith("http"));

    // Default base url
    url = HelpUrls.getUrl("testplugin", props);
    Assertions.assertTrue(url.startsWith("http"));
  }

  @Test
  void canGetUrls() {
    final String url = HelpUrls.getUrl("non-existing key");
    Assertions.assertTrue(url.startsWith("http"));
    Assertions.assertEquals(url, HelpUrls.getUrl(), "Non-existent key should match base URL");
  }

  /**
   * Connects to the online documentation, downloads the pages and tests the links and anchors are
   * reachable.
   *
   * @throws TestAbortedException the test aborted exception
   */
  @Test
  void canReachUrls() throws TestAbortedException {
    // Test we are online
    Assumptions.assumeTrue(testInternetAvailable());

    // Iterate all the URLs and test we can connect.
    // URLs may have anchor tags so build the pages and all their anchors.
    final HashMap<String, Set<String>> pages = new HashMap<>();
    HelpUrls.forEach((k, v) -> {
      // Strip anchor tags
      final int index = v.indexOf('#');
      if (index < 0) {
        pages.computeIfAbsent(v, key -> new HashSet<>());
      } else {
        final String page = v.substring(0, index);
        final String anchor = v.substring(index + 1);
        pages.computeIfAbsent(page, key -> new HashSet<>()).add(anchor);
      }
    });

    // Obtain each page. Then download the page html and check the named anchor tag exists.
    pages.forEach((page, anchors) -> {
      // Use Jsoup library to get the HTML as a document...
      final Document doc = connectToUrl(page);
      anchors.forEach(anchor -> {
        final Elements e = doc.select("a[href=#" + anchor + "]");
        Assertions.assertTrue(e.size() > 0, () -> "Missing " + page + "#" + anchor);
      });
    });
  }

  private static Document connectToUrl(String address) {
    try {
      return Jsoup.connect(address).get();
    } catch (final IOException ex) {
      Assertions.fail("Failed to connect to " + address, ex);
      return null;
    }
  }

  private static boolean testInternetAvailable() {
    return isHostAvailable("google.com") || isHostAvailable("amazon.com")
        || isHostAvailable("facebook.com") || isHostAvailable("apple.com");
  }

  private static boolean isHostAvailable(String hostName) {
    try (Socket socket = new Socket()) {
      final int port = 80;
      final InetSocketAddress socketAddress = new InetSocketAddress(hostName, port);
      socket.connect(socketAddress, 3000);
      return true;
    } catch (final IOException exHost) {
      return false;
    }
  }
}
