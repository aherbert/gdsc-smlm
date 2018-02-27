package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuItem;

import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.View;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

import customnode.CustomMeshNode;
import customnode.CustomTransparentTriangleMesh;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.logging.Ticker;
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
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettingsOrBuilder;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.XYZRResultProcedure;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GUI;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;
import ij3d.Content;
import ij3d.ContentInstant;
import ij3d.Image3DUniverse;
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
	private final static String[] RENDERING = { 
		"Triangle (F=1)",
		"Square (F=2)",
		"Tetrahedron (F=4)",
		"Octahedron (F=8)",
		"Icosahedron (F=20)",
		"Low Res Sphere (F=80)",
		"High Res Sphere (F=320)",
	};

	/**
	 * Gets the number of vertices for each rendering shape.
	 *
	 * @param rendering
	 *            the rendering
	 * @return the number of vertices
	 */
	private static int getNumberOfVertices(int rendering)
	{
		// This is just the F number (faces) multiplied by 3
		switch (rendering)
		{
			case 0:	return 4;
			case 1:	return 8;
			case 2:	return 12;
			case 3:	return 24;
			case 4:	return 60;
			case 5:	return 240;
			case 6:	return 960;
			default:
				throw new IllegalStateException("Unknown rendering");
		}
	}
	//@formatter:on

	// Prevent processing massive mesh objects. 
	// Arrays are used to store vertices that may be 2x larger than the number 
	// of mesh coords input so we could use 1<<30. Make it lower just in case.
	private final static int MAX_SIZE = 1 << 28;

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
		final ImageJ3DResultsViewerSettings settings;
		final int resultsSize;

		public ResultsMetaData(ImageJ3DResultsViewerSettings settings, int resultsSize)
		{
			this.settings = settings;
			this.resultsSize = resultsSize;
		}
	}

	private Image3DUniverse univ;
	private JMenuItem changeColour;
	private JMenuItem changeShading;

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

		ImageJ3DResultsViewerSettings.Builder settings = SettingsManager.readImageJ3DResultsViewerSettings(0)
				.toBuilder();

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Select a dataset to display");
		ResultsManager.addInput(gd, settings.getInputOption(), InputSource.MEMORY);
		gd.addNumericField("Size", settings.getSize(), 2, 6, "nm");
		gd.addSlider("Transparancy", 0, 0.9, settings.getTransparency());
		gd.addChoice("Colour", LUTHelper.luts, settings.getLut());
		gd.addChoice("Rendering", RENDERING, settings.getRendering());
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		final String name = ResultsManager.getInputSource(gd);
		settings.setInputOption(name);
		settings.setSize(gd.getNextNumber());
		settings.setTransparency(gd.getNextNumber());
		settings.setLut(gd.getNextChoiceIndex());
		settings.setRendering(gd.getNextChoiceIndex());
		SettingsManager.writeSettings(settings);
		MemoryPeakResults results = ResultsManager.loadInputResults(name, false, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		// Create a 3D viewer ...

		// TODO:
		// Configure the units used for display.
		// Make the sphere size depend on:
		// 1. the localisation precision.
		// 2. use configured input.
		// 3. any other (e.g. intensity, local density, etc)

		univ = getImage3DUniverse(TITLE);

		// Adapted from Image3DUniverse.addIcospheres.

		// We create a grainy unit sphere. This is then scaled and translated for each localisation.
		final List<Point3f> point = createLocalisationObject(settings.getRendering());
		final List<Point3f> allPoints = new TurboList<Point3f>();

		final int singlePointSize = point.size();
		long size = (long) results.size() * singlePointSize;
		if (size > 10000000L)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("The results will generate a large mesh: " + size);
			gd.setOKLabel("Continue");
			gd.showDialog();
			if (gd.wasCanceled())
				return;
		}

		IJ.showStatus("Creating 3D objects ...");
		final Ticker ticker = Ticker.createStarted(new IJTrackProgress(), results.size(), false);
		final float pointSize = (settings.getSize() > 0) ? (float) settings.getSize() : 1f;
		results.forEach(DistanceUnit.NM, new XYZRResultProcedure()
		{
			//int MAX_SIZE = 1; // For debugging
			public void executeXYZR(float x, float y, float z, PeakResult result)
			{
				// Note sure what the limits are for the graphics library.
				// Assume it can support a max array.
				if (allPoints.size() > MAX_SIZE)
					return;
				allPoints.addAll(copyScaledTranslated(point, pointSize, pointSize, pointSize, x, y, z));
				ticker.tick();
			}
		});
		ticker.stop();

		IJ.showStatus("Creating 3D mesh ...");

		// Default color
		final Color3f color = new Color3f(1, 1, 1);

		float transparency = getTransparency(settings);
		CustomTransparentTriangleMesh mesh;
		mesh = new CustomTransparentTriangleMesh(allPoints, color, transparency);
		//mesh = new CustomTransparentTriangleMesh(null, color, transparency);
		//mesh.setMesh(allPoints);

		// Better colour options.
		// call CustomTransparentTriangleMesh.setTransparentColor(). 
		// Each vertex of the triangles can have a colour. 
		// So we can create a colour for each localisation and duplicate it across the 
		// size of the single point.

		ResultsMetaData data = new ResultsMetaData(settings.build(), results.size());

		changeColour(mesh, results, data.settings);

		mesh.setShaded(true);

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

		IJ.showStatus("");
	}

	private float getTransparency(ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		float transparency = Maths.clip(0, 1, (float) settings.getTransparency());
		// Do not allow fully transparent objects
		if (transparency == 1)
			transparency = 0;
		return transparency;
	}

	private void changeColour(CustomTransparentTriangleMesh mesh, MemoryPeakResults results,
			ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		// Colour by z
		if (results.is3D())
		{
			final GeometryArray ga = (GeometryArray) mesh.getGeometry();
			final int N = ga.getValidVertexCount();

			Color4f[] allColors = new Color4f[N];
			StandardResultProcedure p = new StandardResultProcedure(results);
			p.getZ();
			float[] limits = Maths.limits(p.z);
			final float minimum = limits[0], maximum = limits[1];
			final float scale = 255f / (maximum - minimum);
			LUT lut = LUTHelper.createLUT(LutColour.forNumber(settings.getLut()), false);
			// Create 256 Colors
			Color4f[] colors = new Color4f[256];
			final float w = 1 - getTransparency(settings);
			for (int i = 0; i < 256; i++)
			{
				Color c = new Color(lut.getRGB(i));
				colors[i] = new Color4f(c);
				colors[i].setW(w);
			}

			int vertices = getNumberOfVertices(settings.getRendering());
			for (int i = 0, j = 0, localisations = N / vertices; i < localisations; i++)
			{
				float value = p.z[i];
				value = value - minimum;
				if (value < 0f)
					value = 0f;
				int ivalue = (int) ((value * scale) + 0.5f);
				if (ivalue > 255)
					ivalue = 255;
				for (int k = vertices; k-- > 0;)
					allColors[j++] = colors[ivalue];
			}
			mesh.setTransparentColor(Arrays.asList(allColors));
		}
	}

	/**
	 * Gets the image 3D universe, reusing if possible the same window.
	 *
	 * @param title
	 *            the title
	 * @return the image 3D universe
	 */
	private Image3DUniverse getImage3DUniverse(String title)
	{
		for (Image3DUniverse univ : Image3DUniverse.universes)
		{
			ImageWindow3D w = univ.getWindow();
			if (w != null && w.isVisible() && title.equals(w.getTitle()))
			{
				//univ.removeAllContents();
				//univ.resetView();
				return univ;
			}
		}
		Image3DUniverse univ = new Image3DUniverse();
		univ.show();
		ImageWindow3D w = univ.getWindow();
		GUI.center(w);
		w.setTitle(title);

		// See what these are initialised as ...
		//System.out.println(univ.getViewer().getView().getBackClipDistance());
		//System.out.println(univ.getViewer().getView().getFrontClipDistance());

		// Add a new menu for SMLM functionality
		univ.getMenuBar().add(createSMLMMenuBar());

		univ.addUniverseListener(this);

		return univ;
	}

	/**
	 * Creates the SMLM menu bar.
	 *
	 * @return the menu bar
	 */
	private JMenu createSMLMMenuBar()
	{
		final JMenu add = new JMenu("GDSC SMLM");

		changeColour = new JMenuItem("Change Colour");
		changeColour.addActionListener(this);
		add.add(changeColour);

		changeShading = new JMenuItem("Toggle Shading");
		changeShading.addActionListener(this);
		add.add(changeShading);

		//add.addSeparator();

		return add;
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

		if (src == changeShading)
		{
			for (Iterator<Content> it = univ.contents(); it.hasNext();)
			{
				Content c = it.next();
				Object o = c.getUserData();
				if (o instanceof ResultsMetaData)
				{
					c.setShaded(!c.isShaded());
				}
			}
		}
		else if (src == changeColour)
		{
			ImageJ3DResultsViewerSettings.Builder settings = null;

			for (Iterator<Content> it = univ.contents(); it.hasNext();)
			{
				Content c = it.next();
				Object o = c.getUserData();
				if (o instanceof ResultsMetaData)
				{
					ResultsMetaData data = (ResultsMetaData) o;

					MemoryPeakResults results = ResultsManager.loadInputResults(data.settings.getInputOption(), false,
							null, null);
					// Results must be the same size
					if (results == null || results.size() != data.resultsSize)
					{
						IJ.error(TITLE, "Results are not the same size as when the mesh was constructed: " +
								data.settings.getInputOption());
						continue;
					}

					final ContentInstant content = c.getInstant(0);
					CustomMeshNode node = (CustomMeshNode) content.getContent();
					CustomTransparentTriangleMesh mesh = (CustomTransparentTriangleMesh) node.getMesh();

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
							return;
						settings.setTransparency(gd.getNextNumber());
						settings.setLut(gd.getNextChoiceIndex());
						SettingsManager.writeSettings(settings);
					}

					// Restore original object rendering
					settings.setRendering(data.settings.getRendering());
					changeColour(mesh, results, settings);
				}
			}
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
			case 0:
				return createTriangle();
			case 1:
				return createSquare();
			case 2:
				return createTetrahedron();
			case 3:
				return createOctahedron();
		}
		// All spheres based on icosahedron for speed
		final int subdivisions = rendering - 4;
		return customnode.MeshMaker.createIcosahedron(subdivisions, 1);
	}

	// Note: The triangles are rendered on both sides so the handedness 
	// does not matter for the vertices to define the faces.

	private static float sqrt(double d)
	{
		return (float) Math.sqrt(d);
	}

	static final private float[][] triVertices = { { sqrt(8d / 9), 0, 0 }, { -sqrt(2d / 9), sqrt(2d / 3), 0 },
			{ -sqrt(2d / 9), -sqrt(2d / 3), 0 } };
	static final private int[][] triFaces = { { 0, 1, 2 } };

	static final private float[][] squareVertices = { { 1, 1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, { -1, -1, 0 } };
	static final private int[][] squareFaces = { { 0, 1, 2 }, { 1, 2, 3 } };

	// https://en.m.wikipedia.org/wiki/Tetrahedron
	// based on alternated cube
	static final private float[][] tetraVertices = { { 1, 1, 1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 } };
	// The triangles are rendered on both sides so the handedness does not matter
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
		return createSolid(triVertices, triFaces);
	}

	/**
	 * Creates the square with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createSquare()
	{
		return createSolid(squareVertices, squareFaces);
	}

	/**
	 * Creates the tetrahedron with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createTetrahedron()
	{
		return createSolid(tetraVertices, tetraFaces);
	}

	/**
	 * Creates the octahedron with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createOctahedron()
	{
		return createSolid(octaVertices, octaFaces);
	}

	/**
	 * Creates the solid with the defined faces and vertices on a unit sphere.
	 *
	 * @param vertices
	 *            the vertices
	 * @param faces
	 *            the faces
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createSolid(float[][] vertices, int[][] faces)
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
		final Vector3f v = new Vector3f();
		for (final Point3f p : ps)
		{
			v.set(p);
			v.normalize();
			p.set(v);
		}
		return ps;
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
