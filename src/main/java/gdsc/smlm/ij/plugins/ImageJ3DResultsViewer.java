package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuItem;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.BranchGroup;
import org.scijava.java3d.Canvas3D;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.IndexedGeometryArray;
import org.scijava.java3d.PickInfo;
import org.scijava.java3d.PickInfo.IntersectionInfo;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.SceneGraphPath;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.View;
import org.scijava.java3d.utils.pickfast.PickCanvas;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

import customnode.CustomMesh;
import customnode.CustomMeshNode;
import customnode.CustomPointMesh;
import customnode.CustomTriangleMesh;
import gdsc.core.data.DataException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.core.utils.TurboList;
import gdsc.smlm.data.NamedObject;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.GUIProtos.Image3DDrawingMode;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings.Builder;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettingsOrBuilder;
import gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.ij.ij3d.CustomMeshHelper;
import gdsc.smlm.ij.ij3d.RepeatedIndexedTriangleMesh;
import gdsc.smlm.ij.ij3d.RepeatedTriangleMesh;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.XYZResultProcedure;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TObjectIntHashMap;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionListener;
import ij.gui.GUI;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;
import ij3d.Content;
import ij3d.ContentInstant;
import ij3d.Image3DMenubar;
import ij3d.Image3DUniverse;
import ij3d.ImageCanvas3D;
import ij3d.ImageJ_3D_Viewer;
import ij3d.ImageWindow3D;
import ij3d.UniverseListener;

/**
 * Draws a localisation results set using an ImageJ 3D image
 * 
 * @see <A href="https://imagej.net/3D_Viewer">https://imagej.net/3D_Viewer</a>
 */
public class ImageJ3DResultsViewer implements PlugIn, ActionListener, UniverseListener
{
	private final static String TITLE = "ImageJ 3D Results Viewer";

	//@formatter:off
	private enum Rendering implements NamedObject
	{
		POINT { public String getName() { return "Point"; } 
				public boolean is2D() { return true; }},
		TRIANGLE { public String getName() { return "Triangle"; }
				public boolean is2D() { return true; }},
		SQUARE { public String getName() { return "Square"; }
				public boolean is2D() { return true; }},
        TETRAHEDRON	{ public String getName() { return "Tetrahedron"; }},
        OCTAHEDRON{ public String getName() { return "Octahedron"; }},
        ICOSAHEDRON	{ public String getName() { return "Icosahedron"; }},
        LOW_RES_SPHERE { public String getName() { return "Low Resolution Sphere"; }
        		public boolean isHighResolution() { return true; }},
        HIGH_RES_SPHERE	{ public String getName() { return "High Resolution Sphere"; }
        		public boolean isHighResolution() { return true; }},
        ;

		public String getShortName()
		{
			return getName();
		}		
		
		public boolean is2D() { return false; }
		
		public boolean isHighResolution() { return false; }

		public static Rendering forNumber(int rendering)
		{
			Rendering[] values = Rendering.values();
			if (rendering < 0 || rendering >= values.length)
				throw new IllegalArgumentException();
			return values[rendering];
		}
	};
	//@formatter:on

	private final static String[] RENDERING = SettingsManager.getNames((Object[]) Rendering.values());

	// To debug this from Eclipse relies on being able to find the native 
	// runtime libraries for Open GL. See the README in the eclipse project folder.

	private static String version = "";
	static
	{
		// Support gracefully handling missing dependencies for the 3D viewer
		try
		{
			version = ImageJ_3D_Viewer.getJava3DVersion();
		}
		catch (Throwable t)
		{
			t.printStackTrace();
			version = null;
		}
	}

	private static class ResultsMetaData
	{
		ImageJ3DResultsViewerSettings settings;

		/** The results size when the object mesh was constructed. */
		final int resultsSize;

		public ResultsMetaData(ImageJ3DResultsViewerSettings settings, int resultsSize)
		{
			this.settings = settings;
			this.resultsSize = resultsSize;
		}
	}

	// No ned to store this in settings as when the plugin is first run there are no windows 
	private static String lastWindow = "";

	private static HashMap<String, IJTablePeakResults> resultsTables = new HashMap<String, IJTablePeakResults>();
	private static ResultsTableSettings.Builder resultsTableSettings;

	private Image3DUniverse univ;
	private JMenuItem resetRotation;
	private JMenuItem resetTranslation;
	private JMenuItem resetZoom;
	private JMenuItem resetAll;
	private JMenuItem changeColour;
	private JMenuItem toggleShading;
	private JMenuItem resetSelectedView;
	private JMenuItem colourSurface;
	private JMenuItem updateTableSettings;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (version == null)
		{
			IJ.error(TITLE, "Java 3D is not available");
			return;
		}

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		final ImageJ3DResultsViewerSettings.Builder settings = SettingsManager.readImageJ3DResultsViewerSettings(0)
				.toBuilder();

		// Get a list of the window titles available. Allow the user to select 
		// an existing window or a new one.
		String title = TITLE;
		List<Image3DUniverse> univList = new TurboList<Image3DUniverse>();
		List<String> titleList = new TurboList<String>();
		titleList.add("New window");
		buildWindowList(title, univList, titleList);
		String[] titles = titleList.toArray(new String[titleList.size()]);

		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Select a dataset to display");
		ResultsManager.addInput(gd, settings.getInputOption(), InputSource.MEMORY);
		// The behaviour is to allow the settings to store if the user prefers a new window
		// or to reuse an existing window. If 'new window' then a new window should always
		// be selected. Otherwise open in the same window as last time. If there was no last
		// window then the settings will carried over from the last ImageJ session.
		String window = (settings.getNewWindow()) ? "" : lastWindow;
		gd.addChoice("Window", titles, window);
		gd.addSlider("Transparancy", 0, 0.9, settings.getTransparency());
		gd.addChoice("Colour", LUTHelper.luts, settings.getLut());
		gd.addChoice("Rendering", RENDERING, settings.getRendering(), new OptionListener<Integer>()
		{
			public boolean collectOptions(Integer value)
			{
				settings.setRendering(value);
				return collectOptions(false);
			}

			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog("Drawing mode options", null);
				int rendering = settings.getRendering();
				if (rendering == 0)
				{
					egd.addNumericField("Pixel_size", settings.getPixelSize(), 2, 6, "px");
				}
				else
				{
					// Other modes do not require options
					return false;
				}
				egd.setSilent(silent);
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				settings.setPixelSize(egd.getNextNumber());
				return true;
			}
		});
		gd.addMessage("3D options");
		gd.addCheckbox("Shaded", settings.getShaded());
		gd.addChoice("Drawing_mode", SettingsManager.getImage3DDrawingModeNames(), settings.getDrawingModeValue(),
				new OptionListener<Integer>()
				{
					public boolean collectOptions(Integer value)
					{
						settings.setDrawingModeValue(value);
						return collectOptions(false);
					}

					public boolean collectOptions()
					{
						return collectOptions(true);
					}

					private boolean collectOptions(boolean silent)
					{
						ExtendedGenericDialog egd = new ExtendedGenericDialog("Drawing mode options", null);
						int mode = settings.getDrawingModeValue();
						if (mode == Image3DDrawingMode.DRAW_3D_FIXED_SIZE_VALUE)
						{
							egd.addNumericField("Size", settings.getSize(), 2, 6, "nm");
						}
						else
						{
							// Other modes do not require options
							return false;
						}
						egd.setSilent(silent);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return false;
						settings.setSize(egd.getNextNumber());
						return true;
					}
				});

		gd.showDialog();
		if (gd.wasCanceled())
			return;
		final String name = ResultsManager.getInputSource(gd);
		int windowChoice = gd.getNextChoiceIndex();
		lastWindow = titles[windowChoice];
		settings.setInputOption(name);
		settings.setTransparency(gd.getNextNumber());
		settings.setLut(gd.getNextChoiceIndex());
		settings.setRendering(gd.getNextChoiceIndex());
		settings.setShaded(gd.getNextBoolean());
		settings.setDrawingModeValue(gd.getNextChoiceIndex());
		gd.collectOptions();

		if (windowChoice == 0)
		{
			// Store if the user chose a new window when they had a choice of an existing window
			if (titleList.size() > 1)
				settings.setNewWindow(true);
			// Otherwise they had no choice so leave the preferences as they are.
		}
		else
		{
			// This was not a new window
			settings.setNewWindow(false);
		}

		SettingsManager.writeSettings(settings);
		MemoryPeakResults results = ResultsManager.loadInputResults(name, false, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		// Test object rendering with a small results set
		//		MemoryPeakResults results2 = results;
		//		results = new MemoryPeakResults();
		//		results.copySettings(results2);
		//		results.add(results2.getFirst());
		//		results.add(results2.getFirst().clone());
		//		results.get(1).setZPosition(results.get(0).getZPosition() + 2);
		//		results.add(results2.getFirst().clone());
		//		results.get(2).setZPosition(results.get(0).getZPosition() - 2);

		// Determine if the drawing mode is supported and compute the point size
		final Point3f[] sphereSize = createSphereSize(results, settings);
		if (sphereSize == null)
			return;

		// Cache the table settings
		resultsTableSettings = settings.getResultsTableSettingsBuilder();

		float transparency = getTransparency(settings);

		// Create a 3D viewer.
		if (windowChoice == 0)
			univ = createImage3DUniverse(title, titleList);
		else
			univ = univList.get(windowChoice - 1); // Ignore the new window

		lastWindow = univ.getWindow().getTitle();

		//		{
		//			// Testing what mesh to use
		//			//@formatter:off
        //			CustomPointMesh m1 = new CustomPointMesh(Arrays.asList(new Point3f(-2, -2, -2)));
        //			m1.setPointSize(10);
        //			CustomIndexedTriangleMesh m2 = new CustomIndexedTriangleMesh(
        //					new Point3f[] { 
        //							new Point3f(0, 0, 0), 
        //							new Point3f(1, 0, 0), 
        //							new Point3f(0, 1, 0) },
        //					new int[] { 0, 1, 2 });
        //			CustomQuadMesh m3 = new CustomQuadMesh(Arrays.asList(
        //					new Point3f(2, 2, 2), 
        //					new Point3f(3, 2, 2), 
        //					new Point3f(3, 3, 2),
        //					new Point3f(2, 3, 2)
        //					));
        //			RepeatedIndexedTriangleMesh m4 = new RepeatedIndexedTriangleMesh(
        //					new Point3f[] {
        //							new Point3f(-1, -1, 0), 
        //							new Point3f(1, -1, 0), 
        //							new Point3f(1, 1, 0),
        //							new Point3f(-1, 1, 0)
        //					}, new int[] {0,1,2,0,2,3}, 
        //					new Point3f[]{new Point3f(5, 5, 0)}, null, null, 0);
        //			//@formatter:on
		//			univ.addContent(univ.createContent(m1, "m1"));
		//			univ.addContent(univ.createContent(m2, "m2"));
		//			univ.addContent(univ.createContent(m3, "m3"));
		//			univ.addContent(univ.createContent(m4, "m4"));
		//			if (true)
		//				return;
		//		}

		CustomMesh mesh = createMesh(settings, results, sphereSize, transparency);
		if (mesh == null)
			return;

		updateAppearance(mesh, settings);

		ResultsMetaData data = new ResultsMetaData(settings.build(), results.size());

		mesh.setShaded(settings.getShaded());

		IJ.showStatus("Creating 3D mesh colour ...");
		changeColour(mesh, results, data.settings);

		IJ.showStatus("Creating 3D content ...");
		Content content = univ.createContent(mesh, name);

		content.setUserData(data);

		IJ.showStatus("Drawing 3D content ...");

		// Preserve orientation on the content
		boolean auto = univ.getAutoAdjustView();
		if (univ.getContent(name) == null)
		{
			univ.setAutoAdjustView(true);
		}
		else
		{
			univ.removeContent(name);
			univ.setAutoAdjustView(false);
		}
		univ.addContent(content);
		univ.setAutoAdjustView(auto);

		//		// Debug vertex normals
		//		GeometryArray ga = (GeometryArray) mesh.getGeometry();
		//		TurboList<Point3f> lines = new TurboList<Point3f>();
		//		float[] coordinates = new float[ga.getValidVertexCount() * 3];
		//		float[] normals = new float[coordinates.length];
		//		ga.getCoordinates(0, coordinates);
		//		ga.getNormals(0, normals);
		//		for (int i = 0; i < coordinates.length; i += 3)
		//		{
		//			Point3f p = new Point3f(coordinates[i], coordinates[i + 1], coordinates[i + 2]);
		//			lines.add(p);
		//			p = (Point3f) p.clone();
		//			p.x += normals[i];
		//			p.y += normals[i + 1];
		//			p.z += normals[i + 2];
		//			lines.add(p);
		//		}
		//		univ.removeContent("Normals");
		//		univ.addLineMesh(lines, null, "Normals", false);

		IJ.showStatus("");
	}

	private static Point3f[] createSphereSize(MemoryPeakResults results, Builder settings)
	{
		// Special case for point rendering
		if (settings.getRendering() == 0)
		{
			final float size = getFixedSize(settings.getPixelSize());
			return new Point3f[] { new Point3f(size, size, size) };
		}

		// Store XYZ size for each localisation
		switch (settings.getDrawingMode())
		{
			case DRAW_3D_FIXED_SIZE:
				final float size = getFixedSize(settings.getSize());
				Point3f[] sizes = new Point3f[results.size()];
				Arrays.fill(sizes, new Point3f(size, size, size));
				return sizes;
			case DRAW_3D_XYZ_DEVIATIONS:
				return createSphereSizeFromDeviations(results);
			case DRAW_3D_XY_PRECISION:
				return createSphereSizeFromPrecision(results);
			case UNRECOGNIZED:
			default:
				break;
		}
		return null;
	}

	private static float getFixedSize(double size)
	{
		return (size > 0) ? (float) size : 1f;
	}

	private static Point3f[] createSphereSizeFromDeviations(MemoryPeakResults results)
	{
		if (!results.hasDeviations())
		{
			IJ.error(TITLE, "The results have no deviations");
			return null;
		}
		// Currently the rendering is in nm
		final TypeConverter<DistanceUnit> dc = results.getDistanceConverter(DistanceUnit.NM);

		final Point3f[] size = new Point3f[results.size()];
		boolean failed = results.forEach(new PeakResultProcedureX()
		{
			int i = 0;

			public boolean execute(PeakResult peakResult)
			{
				float x = peakResult.getParameterDeviation(PeakResult.X);
				float y = peakResult.getParameterDeviation(PeakResult.Y);
				float z = peakResult.getParameterDeviation(PeakResult.Z);
				// Check x & y are not zero. 
				// This should be OK as 2D fitting should provide these.
				if (x == 0 || y == 0)
					return true;
				if (z == 0)
				{
					z = (float) Math.sqrt((x * x + y * y) / 2); // Mean variance
					//z = (x + y) / 2; // Mean Std Dev
				}
				size[i++] = new Point3f(dc.convert(x), dc.convert(y), dc.convert(z));
				return false;
			}
		});
		return (failed) ? null : size;
	}

	private static Point3f[] createSphereSizeFromPrecision(MemoryPeakResults results)
	{
		PrecisionResultProcedure p = new PrecisionResultProcedure(results);
		try
		{
			PrecisionMethod m = p.getPrecision();
			IJ.log("Using precision method " + FitProtosHelper.getName(m));
			final Point3f[] size = new Point3f[results.size()];
			for (int i = 0, j = 0; i < p.precision.length; i++)
			{
				// Precision is in NM which matches the rendering
				final float v = (float) p.precision[i];
				size[j++] = new Point3f(v, v, v);
			}
			return size;
		}
		catch (DataException e)
		{
			IJ.error(TITLE, "The results have no precision: " + e.getMessage());
			return null;
		}
	}

	private static float getTransparency(ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		float transparency = Maths.clip(0, 1, (float) settings.getTransparency());
		// Do not allow fully transparent objects
		if (transparency == 1)
			transparency = 0;
		return transparency;
	}

	private static CustomMesh createMesh(final ImageJ3DResultsViewerSettingsOrBuilder settings,
			MemoryPeakResults results, final Point3f[] sphereSize, float transparency)
	{
		// Support drawing as square pixels ...
		if (settings.getRendering() == 0)
		{
			final TurboList<Point3f> points = getPoints(results);
			CustomPointMesh mesh = new CustomPointMesh(points, null, transparency);
			mesh.setPointSize(sphereSize[0].x);
			return mesh;
		}

		// Repeated mesh creation is much faster as the normals are cached.
		// There does not appear to be a difference in the speed the image responds
		// to user interaction between indexed or standard triangles.

		// Currently the RepeatedIndexedTriangleMesh computes the normals a different way to 
		// the super class to preserve the orientation of the normals. So if the coordinates 
		// are modified through the mesh then the appearance will change. For now just use
		// the RepeatedTriangleMesh.

		// Also the IndexedTriangleMesh has one normal per vertex and this causes a colour fall-off
		// on the triangle plane towards the edges. The TriangleMesh colours the entire surface
		// of each triangle the same which looks 'normal'.

		int mode = 2;
		if (mode == 1)
		{
			Pair<Point3f[], int[]> pair = createIndexedObject(settings.getRendering());
			Point3f[] objectVertices = pair.s;
			int[] objectFaces = pair.r;
			long size = (long) results.size() * objectVertices.length;
			long size2 = (long) results.size() * objectFaces.length;
			if (size > 10000000L)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
				egd.addMessage("The results will generate a large mesh of " + size + " vertices and " + size2 +
						" faces.\nThis may take a long time to render and may run out of memory.");
				egd.setOKLabel("Continue");
				egd.showDialog();
				if (egd.wasCanceled())
					return null;
			}
			final TurboList<Point3f> points = getPoints(results);

			IJ.showStatus("Creating 3D mesh ...");
			return new RepeatedIndexedTriangleMesh(objectVertices, objectFaces,
					points.toArray(new Point3f[points.size()]), sphereSize, null, transparency);
		}

		final List<Point3f> point = createLocalisationObject(settings.getRendering());

		final int singlePointSize = point.size();
		long size = (long) results.size() * singlePointSize;
		if (size > 10000000L)
		{
			ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
			egd.addMessage("The results will generate a large mesh of " + size +
					" vertices.\nThis may take a long time to render and may run out of memory.");
			egd.setOKLabel("Continue");
			egd.showDialog();
			if (egd.wasCanceled())
				return null;
		}
		final TurboList<Point3f> points = getPoints(results);

		if (mode == 2)
		{
			IJ.showStatus("Creating 3D mesh ...");
			Rendering r = Rendering.forNumber(settings.getRendering());
			double creaseAngle = (r.isHighResolution()) ? 44 : 0;
			IJTrackProgress progress = null; // Used for debugging construction time
			return new RepeatedTriangleMesh(point.toArray(new Point3f[singlePointSize]),
					points.toArray(new Point3f[points.size()]), sphereSize, null, transparency, creaseAngle, progress);
		}

		// Old method:		
		// Adapted from Image3DUniverse.addIcospheres:
		// We create a grainy unit sphere an the origin. 
		// This is then scaled and translated for each localisation.
		final List<Point3f> allPoints = new TurboList<Point3f>();

		IJ.showStatus("Creating 3D objects ...");
		for (int i = 0; i < sphereSize.length; i++)
		{
			final Point3f p = points.getf(i);
			allPoints.addAll(
					copyScaledTranslated(point, sphereSize[i].x, sphereSize[i].y, sphereSize[i].z, p.x, p.y, p.z));
		}

		IJ.showStatus("Creating 3D mesh ...");

		//mesh = new CustomTransparentTriangleMesh(allPoints, color, transparency);
		// This avoids computing the volume
		CustomTriangleMesh mesh = new CustomTriangleMesh(null, null, transparency);
		mesh.setMesh(allPoints);

		return mesh;
	}

	private static TurboList<Point3f> getPoints(MemoryPeakResults results)
	{
		final TurboList<Point3f> points = new TurboList<Point3f>(results.size());
		results.forEach(DistanceUnit.NM, new XYZResultProcedure()
		{
			public void executeXYZ(float x, float y, float z)
			{
				points.addf(new Point3f(x, y, z));
			}
		});
		return points;
	}

	private static void updateAppearance(CustomMesh mesh, final ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		Appearance appearance = mesh.getAppearance();
		PolygonAttributes pa = appearance.getPolygonAttributes();

		// For all 3D polygons we want to support a true face orientation so transparency works
		Rendering r = Rendering.forNumber(settings.getRendering());
		pa.setBackFaceNormalFlip(r.is2D());
	}

	private static void changeColour(CustomMesh mesh, MemoryPeakResults results,
			ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		// Colour by z
		//if (results.is3D())
		//{
		GeometryArray ga = (GeometryArray) mesh.getGeometry();

		final int vertices = ga.getValidVertexCount();

		Color3f[] allColors = new Color3f[vertices];
		StandardResultProcedure p = new StandardResultProcedure(results);
		p.getZ();
		float[] limits = Maths.limits(p.z);
		final float minimum = limits[0], maximum = limits[1];
		final float range = maximum - minimum;
		LUT lut = LUTHelper.createLUT(LutColour.forNumber(settings.getLut()), false);

		if (range == 0)
		{
			mesh.setColor(new Color3f(new Color(lut.getRGB(255))));
		}
		else
		{
			// Create 256 Colors
			final float scale = 255f / range;
			Color3f[] colors = new Color3f[256];
			for (int i = 0; i < 256; i++)
			{
				Color c = new Color(lut.getRGB(i));
				colors[i] = new Color3f(c);
			}

			final int verticesPerLocalisation = vertices / results.size();
			for (int i = 0, j = 0, size = results.size(); i < size; i++)
			{
				float value = p.z[i];
				value = value - minimum;
				if (value < 0f)
					value = 0f;
				int ivalue = (int) ((value * scale) + 0.5f);
				if (ivalue > 255)
					ivalue = 255;
				for (int k = verticesPerLocalisation; k-- > 0;)
					allColors[j++] = colors[ivalue];
			}
			mesh.setColor(Arrays.asList(allColors));
		}
		mesh.setTransparency(getTransparency(settings));
		//}
	}

	/**
	 * Builds the window list of all visible windows starting with the title prefix.
	 *
	 * @param titlePrefix
	 *            the title prefix
	 * @param univList
	 *            the univ list
	 * @param titleList
	 *            the title list
	 */
	private void buildWindowList(String titlePrefix, List<Image3DUniverse> univList, List<String> titleList)
	{
		for (Image3DUniverse univ : Image3DUniverse.universes)
		{
			ImageWindow3D w = univ.getWindow();
			if (w != null && w.isVisible() && w.getTitle().startsWith(titlePrefix))
			{
				univList.add(univ);
				titleList.add(w.getTitle());
			}
		}
	}

	/**
	 * Creates the image 3D universe with a unique name
	 *
	 * @param title
	 *            the title
	 * @param titleList
	 *            the title list (of titles to ignore)
	 * @return the image 3D universe
	 */
	private Image3DUniverse createImage3DUniverse(String title, List<String> titleList)
	{
		// Get a unique name by appending numbers to the end
		String title2 = title;
		int counter = 2;
		while (titleList.contains(title2))
		{
			title2 = title + " " + (counter++);
		}

		final Image3DUniverse univ = new Image3DUniverse();
		univ.show();
		ImageWindow3D w = univ.getWindow();
		GUI.center(w);
		w.setTitle(title2);

		// See what these are initialised as ...
		//System.out.println(univ.getViewer().getView().getBackClipDistance());
		//System.out.println(univ.getViewer().getView().getFrontClipDistance());

		// Add a new menu for SMLM functionality
		Image3DMenubar menubar = (Image3DMenubar) univ.getMenuBar();
		menubar.add(createSMLMMenuBar());
		univ.setMenubar(menubar);

		univ.addUniverseListener(this);

		// Capture a canvas mouse click/region and identify the coordinates.
		final ImageCanvas3D canvas = (ImageCanvas3D) univ.getCanvas();
		final BranchGroup scene = univ.getScene();

		// We want to be the first mouse listener for the canvas to prevent
		// the content selected event.
		MouseListener[] l = canvas.getMouseListeners();
		for (int i = 0; i < l.length; i++)
			canvas.removeMouseListener(l[i]);

		canvas.addMouseListener(new MouseAdapter()
		{
			@Override
			public void mouseClicked(final MouseEvent e)
			{
				//System.out.println("plugin mouseClicked");
				if (!consumeEvent(e))
					return;

				// This finds the vertex indices of the rendered object.
				Pair<Content, IntersectionInfo> pair = getPickedContent(canvas, scene, e.getX(), e.getY());
				if (pair == null)
					return;

				// Only process content added from localisations
				Content c = pair.s;
				if (!(c.getUserData() instanceof ResultsMetaData))
					return;

				ResultsMetaData data = (ResultsMetaData) c.getUserData();

				MemoryPeakResults results = ResultsManager.loadInputResults(data.settings.getInputOption(), false, null,
						null);
				// Results must be the same size
				if (results == null || results.size() != data.resultsSize)
				{
					IJ.error(TITLE, "Results are not the same size as when the mesh was constructed: " +
							data.settings.getInputOption());
					return;
				}

				// Look up the localisation from the clicked vertex
				final ContentInstant content = c.getInstant(0);
				CustomMeshNode node = (CustomMeshNode) content.getContent();
				CustomMesh mesh = node.getMesh();
				int nVertices;
				GeometryArray ga = (GeometryArray) mesh.getGeometry();
				if (ga instanceof IndexedGeometryArray)
					// An indexed mesh has the correct number of vertex indices
					nVertices = ((IndexedGeometryArray) ga).getValidIndexCount();
				else
					// Default to the number of vertices
					nVertices = ga.getValidVertexCount();

				int nPerLocalisation = nVertices / data.resultsSize;

				// Determine the localisation
				int vertexIndex = pair.r.getVertexIndices()[0];
				int index = vertexIndex / nPerLocalisation;
				//System.out.printf("n=%d [%d]  %s  %s\n", nPerLocalisation, index,
				//		Arrays.toString(pair.r.getVertexIndices()), pair.r.getIntersectionPoint());

				PeakResult p = results.get(index);

				if (e.getClickCount() > 1)
				{
					// Centre on the localisation
					Point3d coordinate = new Point3d();
					ga.getCoordinate(vertexIndex, coordinate);
					univ.centerAt(coordinate);
				}
				else
				{
					// Output the result to a table.
					// Just create a table and add to it.
					IJTablePeakResults table = createTable(results);
					table.add(p);
					table.flush();
				}

				// Consume the event so the item is not selected
				e.consume();
			}

			private IJTablePeakResults createTable(MemoryPeakResults results)
			{
				String name = results.getName();
				IJTablePeakResults table = resultsTables.get(name);
				if (table == null || !table.getResultsWindow().isVisible())
				{
					// Have table settings in the settings.
					// Allow it to be set in the GDSC SMLM menu.
					table = new IJTablePeakResults(results.hasDeviations());
					table.setTableTitle(TITLE + " " + name);
					table.copySettings(results);
					table.setDistanceUnit(resultsTableSettings.getDistanceUnit());
					table.setIntensityUnit(resultsTableSettings.getIntensityUnit());
					table.setAngleUnit(resultsTableSettings.getAngleUnit());
					table.setShowPrecision(resultsTableSettings.getShowPrecision());
					if (resultsTableSettings.getShowPrecision())
						table.setComputePrecision(true);
					table.setShowEndFrame(results.hasEndFrame());
					table.setRoundingPrecision(resultsTableSettings.getRoundingPrecision());
					table.setShowZ(results.is3D());
					table.setShowFittingData(resultsTableSettings.getShowFittingData());
					table.setShowNoiseData(resultsTableSettings.getShowNoiseData());
					table.setShowId(results.hasId());
					table.setAddCounter(true);
					table.setHideSourceText(false);
					table.begin();
					resultsTables.put(name, table);
				}
				return table;
			}

			@Override
			public void mousePressed(final MouseEvent e)
			{
				//System.out.println("plugin mousePressed");
				if (consumeEvent(e))
				{
					// Image3DUniverse adds a canvas mouse listener that calls showPopup.
					// This selects the content. It ignores the e.isConsumed() flag so
					// we translate the point offscreen.
					e.consume();
					e.translatePoint(1000000, 1000000);
				}
			}

			@Override
			public void mouseReleased(final MouseEvent e)
			{
				//System.out.println("plugin mouseReleased");
				if (consumeEvent(e))
				{
					// Image3DUniverse adds a canvas mouse listener that calls showPopup.
					// This selects the content. It ignores the e.isConsumed() flag so
					// we translate the point offscreen.
					e.consume();
					e.translatePoint(1000000, 1000000);
				}
			}

			private boolean consumeEvent(final MouseEvent e)
			{
				// Consume left-mouse clicks with the Ctrl or Alt key down.
				// Single clicks only if showing the results table.
				// Double clicks for centring the universe.

				if (e.isConsumed() || e.getButton() != MouseEvent.BUTTON1 || !(e.isControlDown() || e.isAltDown()))
					return false;
				if (resultsTableSettings.getShowTable() && e.getClickCount() == 1)
					return true;
				if (e.getClickCount() == 2)
					return true;
				return false;
			}
		});

		for (int i = 0; i < l.length; i++)
			canvas.addMouseListener(l[i]);

		return univ;
	}

	/**
	 * Get the Content and closest intersection point at the specified canvas position
	 * <p>
	 * Adapted from Picker.getPickedContent(...).
	 * 
	 * @param x
	 * @param y
	 * @return the Content and closest intersection point
	 */
	private static Pair<Content, IntersectionInfo> getPickedContent(Canvas3D canvas, BranchGroup scene, final int x,
			final int y)
	{
		final PickCanvas pickCanvas = new PickCanvas(canvas, scene);
		pickCanvas.setMode(PickInfo.PICK_GEOMETRY);
		pickCanvas.setFlags(PickInfo.SCENEGRAPHPATH |
				//PickInfo.CLOSEST_INTERSECTION_POINT | 
				PickInfo.CLOSEST_GEOM_INFO);
		pickCanvas.setTolerance(3);
		pickCanvas.setShapeLocation(x, y);
		try
		{
			final PickInfo[] result = pickCanvas.pickAllSorted();
			if (result == null)
				return null;
			for (int i = 0; i < result.length; i++)
			{
				final SceneGraphPath path = result[i].getSceneGraphPath();
				Content c = null;
				for (int j = path.nodeCount(); j-- > 0;)
					if (path.getNode(j) instanceof Content)
						c = (Content) path.getNode(j);
				if (c == null)
					continue;
				return new Pair<Content, IntersectionInfo>(c, result[i].getIntersectionInfos()[0]);
			}
			return null;
		}
		catch (final Exception ex)
		{
			return null;
		}
	}

	/**
	 * Creates the SMLM menu bar.
	 *
	 * @return the menu bar
	 */
	private JMenu createSMLMMenuBar()
	{
		final JMenu add = new JMenu("GDSC SMLM");

		resetRotation = new JMenuItem("Reset Global Rotation");
		resetRotation.addActionListener(this);
		add.add(resetRotation);

		resetTranslation = new JMenuItem("Reset Global Translation");
		resetTranslation.addActionListener(this);
		add.add(resetTranslation);

		resetZoom = new JMenuItem("Reset Global Zoom");
		resetZoom.addActionListener(this);
		add.add(resetZoom);

		add.addSeparator();

		resetAll = new JMenuItem("Reset All");
		resetAll.addActionListener(this);
		add.add(resetAll);

		add.addSeparator();

		changeColour = new JMenuItem("Change Colour");
		changeColour.addActionListener(this);
		add.add(changeColour);

		toggleShading = new JMenuItem("Toggle Shading");
		toggleShading.addActionListener(this);
		add.add(toggleShading);

		resetSelectedView = new JMenuItem("Reset Selected");
		resetSelectedView.addActionListener(this);
		add.add(resetSelectedView);

		add.addSeparator();

		colourSurface = new JMenuItem("Colour surface from 2D image");
		colourSurface.addActionListener(this);
		add.add(colourSurface);

		updateTableSettings = new JMenuItem("Update results table settings");
		updateTableSettings.addActionListener(this);
		add.add(updateTableSettings);

		return add;
	}

	private interface ContentAction
	{
		/**
		 * Run the action
		 * 
		 * @param c
		 *            The content
		 * @return 0 negative for error. No further content can be processed.
		 */
		public int run(Content c);
	}

	private static class ToggleShadedContentAction implements ContentAction
	{
		public int run(Content c)
		{
			c.setShaded(!c.isShaded());
			return 0;
		}
	}

	private static class ChangeColourContentAction implements ContentAction
	{
		ImageJ3DResultsViewerSettings.Builder settings = null;

		public int run(Content c)
		{
			if (!(c.getUserData() instanceof ResultsMetaData))
				return 0;

			ResultsMetaData data = (ResultsMetaData) c.getUserData();

			MemoryPeakResults results = ResultsManager.loadInputResults(data.settings.getInputOption(), false, null,
					null);
			// Results must be the same size
			if (results == null || results.size() != data.resultsSize)
			{
				IJ.error(TITLE, "Results are not the same size as when the mesh was constructed: " +
						data.settings.getInputOption());
				return 1;
			}

			// Change the colour
			if (settings == null)
			{
				// Use the latest settings
				settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
				ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
				gd.addSlider("Transparancy", 0, 0.9, settings.getTransparency());
				gd.addChoice("Colour", LUTHelper.luts, settings.getLut());
				gd.showDialog();
				if (gd.wasCanceled())
					return -1;
				settings.setTransparency(gd.getNextNumber());
				settings.setLut(gd.getNextChoiceIndex());
				SettingsManager.writeSettings(settings);
			}

			final ContentInstant content = c.getInstant(0);
			CustomMeshNode node = (CustomMeshNode) content.getContent();
			CustomMesh mesh = node.getMesh();
			changeColour(mesh, results, settings);
			return 0;
		}
	}

	private static class ResetViewContentAction implements ContentAction
	{
		public int run(Content c)
		{
			final Transform3D t = new Transform3D();
			c.setTransform(t);
			return 0;
		}
	}

	private static class ColourSurfaceContentAction implements ContentAction
	{
		static String title = "";
		static boolean resetTransparency = true;
		ImagePlus imp = null;

		public int run(Content c)
		{
			if (!(c.getUserData() instanceof ResultsMetaData))
				return 0;

			if (imp == null)
			{
				ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
				String[] list = Utils.getImageList(Utils.SINGLE);
				if (list.length == 0)
					return -1;
				gd.addChoice("Image", list, title);
				gd.addCheckbox("Reset_transparency", resetTransparency);
				gd.showDialog();
				if (gd.wasCanceled())
					return -1;
				title = gd.getNextChoice();
				resetTransparency = gd.getNextBoolean();
				imp = WindowManager.getImage(title);
				if (imp == null)
					return -1;
			}

			final ContentInstant content = c.getInstant(0);
			CustomMeshNode node = (CustomMeshNode) content.getContent();
			CustomMesh mesh = node.getMesh();
			CustomMeshHelper helper = new CustomMeshHelper(mesh);
			helper.loadSurfaceColorsFromImage2D(imp);
			// Remove the transparency
			if (resetTransparency)
				mesh.setTransparency(0);
			return 0;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	@SuppressWarnings("unchecked")
	public void actionPerformed(ActionEvent e)
	{
		final Object src = e.getSource();

		// Universe actions
		// Adapted from univ.resetView();
		if (src == resetRotation)
		{
			univ.fireTransformationStarted();
			// rotate so that y shows downwards
			final Transform3D t = new Transform3D();
			final AxisAngle4d aa = new AxisAngle4d(1, 0, 0, Math.PI);
			t.set(aa);
			univ.getRotationTG().setTransform(t);
			univ.fireTransformationUpdated();
			univ.fireTransformationFinished();
			return;
		}
		if (src == resetTranslation)
		{
			univ.fireTransformationStarted();
			final Transform3D t = new Transform3D();
			univ.getTranslateTG().setTransform(t);
			univ.recalculateGlobalMinMax();
			univ.getViewPlatformTransformer().centerAt(univ.getGlobalCenterPoint());
			univ.fireTransformationUpdated();
			univ.fireTransformationFinished();
			return;
		}
		if (src == resetZoom)
		{
			univ.fireTransformationStarted();
			final Transform3D t = new Transform3D();
			univ.getZoomTG().setTransform(t);
			Point3d max = new Point3d();
			Point3d min = new Point3d();
			univ.getGlobalMaxPoint(max);
			univ.getGlobalMinPoint(min);
			float range = (float) (max.x - min.x);
			final double d = (range) / Math.tan(Math.PI / 8);
			univ.getViewPlatformTransformer().zoomTo(d);
			univ.fireTransformationUpdated();
			univ.fireTransformationFinished();
			return;
		}
		if (src == updateTableSettings)
		{
			ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			ResultsSettings.Builder s = ResultsSettings.newBuilder();
			s.setResultsTableSettings(resultsTableSettings); // This is from the cache
			gd.addMessage("Click on the image to view localisation data.\nCtrl/Alt key must be pressed.");
			ResultsManager.addTableResultsOptions(gd, s, ResultsManager.FLAG_NO_SECTION_HEADER);
			gd.showDialog();
			if (gd.wasCanceled())
				return;
			resultsTableSettings = s.getResultsTableSettingsBuilder();
			resultsTableSettings.setShowTable(gd.getNextBoolean());

			// Save updated settings
			final ImageJ3DResultsViewerSettings.Builder settings = SettingsManager.readImageJ3DResultsViewerSettings(0)
					.toBuilder();
			settings.setResultsTableSettings(resultsTableSettings);
			SettingsManager.writeSettings(settings);
			return;
		}

		// Actions to perform on content
		ContentAction action = null;
		if (src == toggleShading)
		{
			action = new ToggleShadedContentAction();
		}
		else if (src == changeColour)
		{
			action = new ChangeColourContentAction();
		}
		else if (src == resetAll)
		{
			univ.resetView();
			univ.select(null);
			action = new ResetViewContentAction();
		}
		else if (src == resetSelectedView)
		{
			action = new ResetViewContentAction();
		}
		else if (src == colourSurface)
		{
			action = new ColourSurfaceContentAction();
		}
		if (action == null)
			return;

		if (univ.getSelected() != null)
		{
			action.run(univ.getSelected());
			return;
		}

		for (Iterator<Content> it = univ.contents(); it.hasNext();)
		{
			if (action.run(it.next()) < 0)
				return;
		}
	}

	/**
	 * Creates the object used to draw a single localisation.
	 * 
	 * @param rendering
	 *
	 * @return the list of triangle vertices for the object
	 */
	private static List<Point3f> createLocalisationObject(int rendering)
	{
		switch (rendering)
		{
			case 1:
				return createTriangle();
			case 2:
				return createSquare();
			case 3:
				return createTetrahedron();
			case 4:
				return createOctahedron();
		}
		// All spheres based on icosahedron for speed
		final int subdivisions = rendering - 5;
		return customnode.MeshMaker.createIcosahedron(subdivisions, 1);
	}

	// Note: The triangles are rendered on both sides so the handedness 
	// does not matter for the vertices to define the faces. However 
	// the handedness does matter for the vertices if using transparency.
	// The orders below have been worked out by trial and error.

	private static float sqrt(double d)
	{
		return (float) Math.sqrt(d);
	}

	static final private float[][] triVertices = { { sqrt(8d / 9), 0, 0 }, { -sqrt(2d / 9), sqrt(2d / 3), 0 },
			{ -sqrt(2d / 9), -sqrt(2d / 3), 0 } };
	static final private int[][] triFaces = { { 0, 1, 2 } };

	static final private float[][] squareVertices = { { 1, 1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, { -1, -1, 0 } };
	static final private int[][] squareFaces = { { 0, 1, 2 }, { 2, 1, 3 } };

	// https://en.m.wikipedia.org/wiki/Tetrahedron
	// based on alternated cube
	static final private float[][] tetraVertices = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };
	static final private int[][] tetraFaces = { { 0, 1, 2 }, { 0, 1, 3 }, { 1, 2, 3 }, { 0, 2, 3 } };

	// https://en.m.wikipedia.org/wiki/Octahedron
	static final private float[][] octaVertices = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 },
			{ 0, 0, -1 }, };
	static final private int[][] octaFaces = { { 0, 3, 4 }, { 3, 1, 4 }, { 1, 2, 4 }, { 2, 0, 4 }, { 3, 0, 5 },
			{ 1, 3, 5 }, { 2, 1, 5 }, { 0, 2, 5 }, };

	/**
	 * Creates the triangle with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createTriangle()
	{
		return createSolid(triVertices, triFaces, true);
	}

	/**
	 * Creates the square with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createSquare()
	{
		return createSolid(squareVertices, squareFaces, false);
	}

	/**
	 * Creates the tetrahedron with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createTetrahedron()
	{
		return createSolid(tetraVertices, tetraFaces, true);
	}

	/**
	 * Creates the octahedron with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createOctahedron()
	{
		// This is already normalised
		return createSolid(octaVertices, octaFaces, false);
	}

	/**
	 * Creates the solid with the defined faces and vertices on a unit sphere.
	 *
	 * @param vertices
	 *            the vertices
	 * @param faces
	 *            the faces
	 * @param normalise
	 *            the normalise
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createSolid(float[][] vertices, int[][] faces, boolean normalise)
	{
		List<Point3f> ps = new ArrayList<Point3f>();
		for (int i = 0; i < faces.length; i++)
		{
			for (int k = 0; k < 3; k++)
			{
				ps.add(new Point3f(vertices[faces[i][k]]));
			}
		}
		// Project all vertices to the surface of a sphere of radius 1
		if (normalise)
		{
			final Vector3f v = new Vector3f();
			for (final Point3f p : ps)
			{
				v.set(p);
				v.normalize();
				p.set(v);
			}
		}
		return ps;
	}

	/**
	 * Creates the object used to draw a single localisation.
	 * 
	 * @param rendering
	 *
	 * @return the vertices and faces of the the object
	 */
	private static Pair<Point3f[], int[]> createIndexedObject(int rendering)
	{
		List<Point3f> list = createLocalisationObject(rendering);

		// Compact the vertices to a set of vertices and faces
		final TObjectIntHashMap<Point3f> m = new TObjectIntHashMap<Point3f>(list.size(), 0.5f, -1);
		TurboList<Point3f> vertices = new TurboList<Point3f>(list.size());
		TIntArrayList faces = new TIntArrayList(list.size());
		int index = 0;
		// Process triangles
		for (int i = 0; i < list.size(); i += 3)
		{
			index = addFace(m, vertices, faces, list.get(i), index);
			index = addFace(m, vertices, faces, list.get(i + 1), index);
			index = addFace(m, vertices, faces, list.get(i + 2), index);
		}

		return new Pair<Point3f[], int[]>(vertices.toArray(new Point3f[vertices.size()]), faces.toArray());
	}

	private static int addFace(TObjectIntHashMap<Point3f> m, TurboList<Point3f> vertices, TIntArrayList faces,
			Point3f p, int index)
	{
		// Add the point if it is not in the set of vertices.
		// Get the index associated with the vertex.
		int value = m.putIfAbsent(p, index);
		if (value == -1)
		{
			// Store the points in order
			vertices.add(p);
			value = index++;
		}
		faces.add(value);
		return index;
	}

	/**
	 * Copy scaled and translated.
	 * <p>
	 * Adapted from customnode.MeshMaker.copyTranslated
	 *
	 * @param ps
	 *            the ps
	 * @param sx
	 *            the sx
	 * @param sy
	 *            the sy
	 * @param sz
	 *            the sz
	 * @param dx
	 *            the dx
	 * @param dy
	 *            the dy
	 * @param dz
	 *            the dz
	 * @return the list
	 */
	private static List<Point3f> copyScaledTranslated(final List<Point3f> ps, final float sx, final float sy,
			final float sz, final float dx, final float dy, final float dz)
	{
		final TurboList<Point3f> verts = new TurboList<Point3f>(ps.size());

		//		// Reuse points to save memory and scaling time
		//		final HashMap<Point3f, Point3f> m = new HashMap<Point3f, Point3f>();
		//		for (final Point3f p : ps)
		//		{
		//			Point3f p2 = m.get(p);
		//			if (null == p2)
		//			{
		//				p2 = new Point3f(p.x * sx + dx, p.y * sy + dy, p.z * sz + dz);
		//				m.put(p, p2);
		//			}
		//			verts.addf(p2);
		//		}

		// Duplicate all points
		for (final Point3f p : ps)
		{
			verts.addf(new Point3f(p.x * sx + dx, p.y * sy + dy, p.z * sz + dz));
		}

		return verts;
	}

	public void transformationStarted(View view)
	{
	}

	public void transformationUpdated(View view)
	{
		// This is called when the zoom is adjusted. We can update clipping
		//System.out.println(univ.getViewer().getView().getBackClipDistance());
		//System.out.println(univ.getViewer().getView().getFrontClipDistance());
	}

	public void transformationFinished(View view)
	{
	}

	public void contentAdded(Content c)
	{
	}

	public void contentRemoved(Content c)
	{
	}

	public void contentChanged(Content c)
	{
	}

	public void contentSelected(Content c)
	{
	}

	public void canvasResized()
	{
	}

	public void universeClosed()
	{
	}
}
