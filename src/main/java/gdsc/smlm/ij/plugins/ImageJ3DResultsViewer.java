package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuItem;

import org.scijava.java3d.View;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3f;

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
	// Prevent processing massive mesh objects. 
	// Arrays are used to store vertices that may be 2x larger than the number 
	// of mesh coords input so we could use 1<<30. Make it lower just in case.
	private final static int MAX_SIZE = 1 << 28;

	// To debug this from Eclipse relies on being able to find the native 
	// runtime libraries for Open GL.

	// Linux
	// sudo aptitude install libgluegen2-rt-java libgluegen2-build-java libjogl2-java libjogl2-jni libjogl2-toolkits

	// I have tried adding the /usr/lib/jni folder to the Eclipse project properties
	// build path. This did not work.	

	// This did work:
	// Created [GDSC SMLM]natives/linux-amd64 and put the following links in:
	// ln -s /usr/lib/jni/libgluegen2-rt.so libgluegen-rt.so
	// ln -s /usr/lib/jni/libnativewindow_awt.so

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

	private Image3DUniverse univ;
	private JMenuItem reset;

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
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		final String name = ResultsManager.getInputSource(gd);
		settings.setInputOption(name);
		settings.setSize(gd.getNextNumber());
		settings.setTransparency(gd.getNextNumber());
		settings.setLut(gd.getNextChoiceIndex());
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
		// Q. Can spheres be translucent so higher density makes it more opaque?

		univ = getImage3DUniverse(TITLE);

		// Adapted from Image3DUniverse.addIcospheres.

		// We create a grainy unit sphere. This is then scaled and translated for each localisation.
		final List<Point3f> point = createSphere();
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

		univ.removeContent(name);

		float transparency = Maths.clip(0, 1, (float) settings.getTransparency());
		// Do not allow fully transparent objects
		if (transparency == 1)
			transparency = 0;
		CustomTransparentTriangleMesh mesh;
		mesh = new CustomTransparentTriangleMesh(allPoints, color, transparency);
		//mesh = new CustomTransparentTriangleMesh(null, color, transparency);
		//mesh.setMesh(allPoints);

		// Better colour options.
		// call CustomTransparentTriangleMesh.setTransparentColor(). 
		// Each vertex of the triangles can have a colour. 
		// So we can create a colour for each localisation and duplicate it across the 
		// size of the single point.

		// Color by z
		if (results.is3D())
		{
			Color4f[] allColors = new Color4f[allPoints.size()];
			StandardResultProcedure p = new StandardResultProcedure(results);
			p.getZ();
			float[] limits = Maths.limits(p.z);
			final float minimum = limits[0], maximum = limits[1];
			final float scale = 255f / (maximum - minimum);
			LUT lut = LUTHelper.createLUT(LutColour.forNumber(settings.getLut()), false);
			// Create 256 Colors
			Color4f[] colors = new Color4f[256];
			for (int i = 0; i < 256; i++)
			{
				Color c = new Color(lut.getRGB(i));
				colors[i] = new Color4f(c);
				colors[i].setW(transparency);
			}

			for (int i = 0, j = 0, total = p.z.length; i < total; i++)
			{
				float value = p.z[i];
				value = value - minimum;
				if (value < 0f)
					value = 0f;
				int ivalue = (int) ((value * scale) + 0.5f);
				if (ivalue > 255)
					ivalue = 255;
				for (int k = singlePointSize; k-- > 0; )
					allColors[j++] = colors[ivalue];
			}
			mesh.setTransparentColor(Arrays.asList(allColors));
		}

		IJ.showStatus("Creating 3D content ...");
		Content content = univ.createContent(mesh, name);

		IJ.showStatus("Drawing 3D content ...");
		boolean auto = univ.getAutoAdjustView();
		univ.setAutoAdjustView(true);
		univ.addContent(content);
		univ.setAutoAdjustView(auto);

		IJ.showStatus("");
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
				univ.removeAllContents();
				univ.resetView();
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

		reset = new JMenuItem("Reset");
		reset.addActionListener(this);
		add.add(reset);

		//add.addSeparator();
		
		return add;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e)
	{
		final Object src = e.getSource();

		if (src == reset)
			univ.resetView();
	}

	/**
	 * Creates the sphere used to draw a single localisation.
	 *
	 * @return the list of triangle vertices for the sphere
	 */
	private static List<Point3f> createSphere()
	{
		final int subdivisions = 0;
		final float radius = 1;
		return customnode.MeshMaker.createIcosahedron(subdivisions, radius);
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
