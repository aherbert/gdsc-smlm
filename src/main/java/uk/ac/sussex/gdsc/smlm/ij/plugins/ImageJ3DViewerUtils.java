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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij3d.ImageJ_3D_Viewer;

import java.security.AccessController;
import java.security.PrivilegedAction;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Utility class for the {@link ImageJ_3D_Viewer}.
 */
final class ImageJ3DViewerUtils {

  // To debug this from Eclipse relies on being able to find the native
  // runtime libraries for Open GL. See the README in the eclipse project folder.

  /** The Java 3D version if available, otherwise null. */
  static final String JAVA_3D_VERSION;

  static {
    // Try setting -Dj3d.sortShape3DBounds for faster centroid computation.
    // See org.scijava.java3d.MasterControl.sortShape3DBounds.
    // This only works if the VirtualUniverse has not been created, i.e. no java3d has been loaded.
    AccessController.doPrivileged((PrivilegedAction<
        String>) () -> System.setProperty("j3d.sortShape3DBounds", Boolean.toString(true)));

    // Support gracefully handling missing dependencies for the 3D viewer
    String version = null;
    try {
      version = ImageJ_3D_Viewer.getJava3DVersion();
    } catch (final Throwable thrown) {
      Logger.getLogger(ImageJ3DViewerUtils.class.getName()).log(Level.SEVERE,
          "Java 3D is not available", thrown);
    } finally {
      JAVA_3D_VERSION = version;
    }
  }

  /**
   * No public constructor.
   */
  private ImageJ3DViewerUtils() {}
}
