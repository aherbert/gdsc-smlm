/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import customnode.CustomMesh;
import customnode.CustomTriangleMesh;
import customnode.MeshMaker;

import ij.gui.GUI;
import ij.plugin.PlugIn;

import ij3d.DefaultUniverse;
import ij3d.Image3DUniverse;
import ij3d.ImageWindow3D;

import org.scijava.java3d.PolygonAttributes;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

import java.util.List;

/**
 * This is a demo to demonstrate the issues with transparency in the ImageJ 3D Viewer.
 */
public class Transparency_Demo implements PlugIn {
  @Override
  public void run(String arg) {
    createUniverse("Transparency_Demo Backface Normal Flip (Default)", false, false);
    createUniverse("Transparency_Demo No Backface Normal Flip", true, false);
    createUniverse("Transparency_Demo No Backface Normal Flip + Backface Cull", true, true);
  }

  private static void createUniverse(String title, boolean disableBackfaceNormalFlip,
      boolean backfaceCull) {
    final Image3DUniverse univ = new Image3DUniverse();
    univ.showAttribute(DefaultUniverse.ATTRIBUTE_SCALEBAR, false);
    univ.show();
    final ImageWindow3D w = univ.getWindow();
    GUI.center(w);
    w.setTitle(title);

    addPoint(univ, disableBackfaceNormalFlip, backfaceCull, -1.5f, 0, 0, new Color3f(1, 0, 0));
    addPoint(univ, disableBackfaceNormalFlip, backfaceCull, 1.5f, 0, 0, new Color3f(0, 1, 0));
    addPoint(univ, disableBackfaceNormalFlip, backfaceCull, -1.8f, -0.3f, -1.5f,
        new Color3f(0, 0, 1));
    univ.sync(true);
  }

  private static void addPoint(Image3DUniverse univ, boolean disableBackfaceNormalFlip,
      boolean backfaceCull, float x, float y, float z, Color3f c) {
    final List<Point3f> points = MeshMaker.createIcosahedron(0, 1f);
    for (final Point3f p : points) {
      p.x += x;
      p.y += y;
      p.z += z;
    }
    final CustomMesh mesh = new CustomTriangleMesh(points, c, 0.5f);
    mesh.getAppearance().getPolygonAttributes().setBackFaceNormalFlip(!disableBackfaceNormalFlip);
    if (backfaceCull) {
      mesh.getAppearance().getPolygonAttributes().setCullFace(PolygonAttributes.CULL_BACK);
    }
    univ.addCustomMesh(mesh, "Point " + x);
  }
}
