
/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

import ij.IJ;
import ij.ImageJ;
import ij.plugin.PlugIn;
import java.io.File;
import java.net.URISyntaxException;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmTools;

/**
 * Default ImageJ plugin (no Java package) to run the {@link SmlmTools} plugin.
 *
 * <p><strong>This class is not included in the packaged uk.ac.sussex.gdsc.smlm jar.</strong>
 *
 * <p>This class can run ImageJ and load all the GDSC plugins using the {@link #main(String[])}
 * function to launch a Java application.
 *
 * <p>Alternatively this class allows the project to be run in from an IDE (e.g. Eclipse) using
 * ImageJ's detection of classes in the plugins folder. The Maven output directory will be
 * target/classes. Create a symbolic link to that directory from the project root and name it
 * plugins. Optionally create a link to the macros directory to allow the toolset to be loaded:
 *
 * <pre>
 * ${root}/plugins -&gt; ${root}/target/classes
 * ${root}/macros -&gt; ${root}/target/classes/macros
 * </pre>
 *
 * <p>Set the project to run ij.ImageJ as the main class and use the root directory as the ImageJ
 * path:
 *
 * <pre>
 * ij.ImageJ -ijpath ${root}
 * </pre>
 *
 * <p>ImageJ will load this class from the plugins directory. This class can call all other plugins.
 */
public class Smlm_PlugIn implements PlugIn {
  @SuppressWarnings("unused")
  @Override
  public void run(String arg) {
    // Create the SMLM Tools plugin.
    new SmlmTools();
  }


  /**
   * Main method for debugging.
   *
   * <p>For debugging, it is convenient to have a method that starts ImageJ and calls the plugin,
   * e.g. after setting breakpoints.
   *
   * @param args unused
   * @throws URISyntaxException if the URL cannot be converted to a URI
   */
  public static void main(String[] args) throws URISyntaxException {
    // Set the base directory for plugins
    // see: https://stackoverflow.com/a/7060464/1207769
    Class<Smlm_PlugIn> clazz = Smlm_PlugIn.class;
    java.net.URL url = clazz.getProtectionDomain().getCodeSource().getLocation();
    File file = new File(url.toURI());
    // Note: This returns the base path. ImageJ will find plugins in here that have an
    // underscore in the name. But it will not search recursively through the
    // package structure to find plugins. Adding this at least puts it on ImageJ's
    // classpath so plugins not satisfying these requirements can be loaded.
    System.setProperty("plugins.dir", file.getAbsolutePath());

    // Start ImageJ and exit when closed
    ImageJ imagej = new ImageJ();
    imagej.exitWhenQuitting(true);

    // Run the plugin
    IJ.runPlugIn(clazz.getName(), "");
  }
}
