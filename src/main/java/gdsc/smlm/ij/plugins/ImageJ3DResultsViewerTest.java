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
package gdsc.smlm.ij.plugins;

import java.util.List;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.Shape3D;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.View;
import org.scijava.java3d.utils.geometry.Primitive;
import org.scijava.java3d.utils.geometry.Sphere;
import org.scijava.vecmath.Point3f;

import customnode.CustomMesh;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.ij.ij3d.CustomContent;
import gdsc.smlm.ij.ij3d.CustomContentHelper;
import gdsc.smlm.ij.ij3d.ItemGeometryGroup;
import gdsc.smlm.ij.ij3d.ItemGroupNode;
import gdsc.smlm.ij.ij3d.ItemIndexedTriangleMesh;
import gdsc.smlm.ij.ij3d.ItemMesh;
import gdsc.smlm.ij.ij3d.ItemTriangleMesh;
import gdsc.smlm.ij.ij3d.ReferenceItemMesh;
import gdsc.smlm.ij.ij3d.Shape3DHelper;
import gdsc.smlm.ij.ij3d.Shape3DHelper.Rendering;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.utils.Pair;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GUI;
import ij.plugin.PlugIn;
import ij3d.ContentInstant;
import ij3d.DefaultUniverse;
import ij3d.Image3DUniverse;
import ij3d.ImageJ_3D_Viewer;
import ij3d.ImageWindow3D;

/**
 * Tests for drawing a localisation results set using an ImageJ 3D image
 *
 * @see <A href="https://imagej.net/3D_Viewer">https://imagej.net/3D_Viewer</a>
 */
public class ImageJ3DResultsViewerTest implements PlugIn
{
	private final static String TITLE = "ImageJ 3D Results Viewer";

	// To debug this from Eclipse relies on being able to find the native
	// runtime libraries for Open GL. See the README in the eclipse project folder.

	private static String version = "";
	static
	{
		// Try setting -Dj3d.sortShape3DBounds for faster centroid computation
		// See MasterControl.sortShape3DBounds. This only works if The VirtualUniverse
		// has not been created.
		java.security.AccessController.doPrivileged(new java.security.PrivilegedAction<String>()
		{
			@Override
			public String run()
			{
				return System.setProperty("j3d.sortShape3DBounds", Boolean.toString(true));
			}
		});

		// Support gracefully handling missing dependencies for the 3D viewer
		try
		{
			version = ImageJ_3D_Viewer.getJava3DVersion();
		}
		catch (final Throwable t)
		{
			t.printStackTrace();
			version = null;
		}
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		if (version == null)
		{
			IJ.error(TITLE, "Java 3D is not available");
			return;
		}

		final boolean renderingTest = true;
		if (renderingTest)
		{
			// Put all the shapes from rendering on a view and see how they look.
			final Image3DUniverse univ = new Image3DUniverse();
			univ.showAttribute(DefaultUniverse.ATTRIBUTE_SCALEBAR, false);
			univ.show();
			final ImageWindow3D w = univ.getWindow();
			GUI.center(w);

			// Test how many vertices a primitive sphere has
			float x = 0;
			final float y = 0;
			final float space = 2.5f;

			//for (Rendering rendering : new Rendering[] { Rendering.SQUARE, Rendering.OCTAGON, Rendering.ICOSAHEDRON })
			for (final Rendering rendering : Rendering.values())
			{
				final Shape3D shape = Shape3DHelper.createShape(rendering, 3);

				final Appearance a = shape.getAppearance();

				final ItemMesh mesh = new ReferenceItemMesh(new Point3f[] { new Point3f(x, y, 0) },
						(GeometryArray) shape.getGeometry(), a, null, null, 0f);

				System.out.printf("R=%s %s  Vc=%d  V=%d  T=%d\n", rendering,
						shape.getGeometry().getClass().getSimpleName(), mesh.getVerticesCountPerItem(),
						mesh.getVerticesPerItem(), Shape3DHelper.getNumberOfTriangles(rendering));

				if (rendering == Rendering.POINT)
					a.getPointAttributes().setPointSize(10);
				univ.addCustomMesh(mesh, x + "," + y);
				x += space;
			}

			return;
		}

		final boolean sphereTest = true;
		if (sphereTest)
		{
			// Sphere test
			// Put lots of spheres on a view and see how they look.
			// Is it worth supporting an ItemTriangleStripMesh so we can control
			// the sphere better than the icosahedron?
			final Image3DUniverse univ = new Image3DUniverse();
			univ.showAttribute(DefaultUniverse.ATTRIBUTE_SCALEBAR, false);
			univ.show();
			final ImageWindow3D w = univ.getWindow();
			GUI.center(w);

			// Test how many vertices a primitive sphere has
			float x = 0, y = 0;
			final float space = 2.5f;

			Appearance a = null;
			for (int d = 0; d < 4; d++)
			{
				final List<Point3f> points = customnode.MeshMaker.createIcosahedron(d, 1f);
				final Pair<Point3f[], int[]> p = CustomContentHelper.createIndexedObject(points);
				final int v = points.size();
				final int t = v / 3;
				System.out.printf("Icosahedron divisions = %d, V=%d, T=%d, Vi=%d (%.2f), i=%d\n", d, v, t, p.a.length,
						v / (double) p.a.length, p.b.length);

				CustomMesh mesh = new ItemTriangleMesh(points.toArray(new Point3f[0]),
						new Point3f[] { new Point3f(x, y, 0) }, null, null, 0);

				a = mesh.getAppearance();
				univ.addCustomMesh(mesh, x + "," + y + "," + t);

				final float y2 = y + space;
				mesh = new ItemIndexedTriangleMesh(p.a, p.b, new Point3f[] { new Point3f(x, y2, 0) }, null, null, 0);
				univ.addCustomMesh(mesh, x + "," + y2 + "," + t);

				x += space;
			}

			// Avoid null pointer warnings
			if (a == null)
				throw new NullPointerException();

			// The T=800 sphere looks about the same as the Icosahedron(div=3) T=1280
			// This may be a better super-high resolution option.

			x = 0;
			y += 2 * space;
			a = (Appearance) a.cloneNodeComponent(true);
			//a.getColoringAttributes().setColor(0, 1, 0);
			a.getMaterial().setDiffuseColor(0, 1, 0);
			for (int d = 4; d < 50; d += 4)
			{
				// This is a triangle strip array so is more space efficient
				final Sphere s = new Sphere(1, Primitive.GENERATE_NORMALS, d);
				final int t = s.getNumTriangles();
				System.out.printf("Sphere divisions = %d, V=%d, T=%d\n", d, s.getNumVertices(), t);

				final ItemGeometryGroup g = new ItemGeometryGroup(new Point3f[] { new Point3f(x, y, 0) },
						(GeometryArray) s.getShape().getGeometry(), a, null, null, null);
				final String name = x + "," + y + "," + t;
				final CustomContent content = new CustomContent(name, true);
				content.getCurrent().display(new ItemGroupNode(g));
				univ.addContent(content);

				x += space;
			}

			return;
		}

		TurboList<Point3f> pointList;
		float scale;

		if (MemoryPeakResults.isMemoryEmpty())
		{
			pointList = new TurboList<>();
			int range;
			range = 1; // 9 points
			//range = 12; // 25^3 = 15625 points
			//range = 17; // 35^3 = 42875 points
			//range = 22; // 45^3 = 91125 points
			//range = 49; // 99^3 = 970299 points
			for (int x = -range; x <= range; x++)
				for (int y = -range; y <= range; y++)
					for (int z = -range; z <= range; z++)
						pointList.add(new Point3f(x, y, z));
			scale = 0.25f;
		}
		else
		{
			final ImageJ3DResultsViewerSettings.Builder settings = SettingsManager.readImageJ3DResultsViewerSettings(0)
					.toBuilder();
			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select a dataset to display");
			ResultsManager.addInput(gd, settings.getInputOption(), InputSource.MEMORY);
			gd.showDialog();
			if (gd.wasCanceled())
				return;
			final String name = ResultsManager.getInputSource(gd);
			settings.setInputOption(name);

			SettingsManager.writeSettings(settings);
			final MemoryPeakResults results = ResultsManager.loadInputResults(name, false, null, null);
			if (results == null || results.size() == 0)
			{
				IJ.error(TITLE, "No results could be loaded");
				IJ.showStatus("");
				return;
			}
			pointList = ImageJ3DResultsViewer.getPoints(results, settings);
			if (pointList == null)
				return;
			scale = 10f;
		}

		final Image3DUniverse univ = new Image3DUniverse();
		univ.showAttribute(DefaultUniverse.ATTRIBUTE_SCALEBAR, false);
		univ.show();
		final ImageWindow3D w = univ.getWindow();
		GUI.center(w);

		final View view = univ.getViewer().getView();
		view.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);
		// I am not sure if this is required if objects are sorted.
		//view.setDepthBufferFreezeTransparent(false);

		IJ.showStatus("Creating points ...");
		final Point3f[] points = pointList.toArray(new Point3f[pointList.size()]);
		final Point3f[] sizes = new Point3f[] { new Point3f(scale, scale, scale) };
		ItemGeometryGroup pointGroup;
		final Appearance appearance = new Appearance();
		final TransparencyAttributes ta = new TransparencyAttributes();
		ta.setTransparency(0.5f);
		ta.setTransparencyMode(TransparencyAttributes.FASTEST);
		appearance.setTransparencyAttributes(ta);
		pointGroup = new ItemGeometryGroup(points, null, appearance, sizes, null, null);
		//pointGroup = new OrderedItemGeometryGroup(points, null, appearance, sizes, null, null);

		//		// This supports transparency sorting
		//		BranchGroup bg = new BranchGroup();
		//		bg.addChild(pointGroup);
		//		bg.compile();
		//		univ.getScene().addChild(bg);

		//		// This does not since ContentInstant uses an OrderedPath to show:
		//		// the content; the object bounding box, coordinate system and point list
		//		final Content c = new Content("Test");
		//		final ContentInstant content = c.getCurrent();
		//		content.display(new PointGroupNode(pointGroup));
		//		univ.addContent(c);

		// This does since ItemGeometryNode uses a Group to show the points as individual shapes
		// and the CustomContentInstant uses a group not an ordered group so show all the adornments.
		IJ.showStatus("Displaying points ...");
		final CustomContent c = new CustomContent("Test", false);
		final ContentInstant content = c.getCurrent();
		content.display(new ItemGroupNode(pointGroup));
		univ.addContent(c);

		IJ.showStatus("Done");
	}
}
