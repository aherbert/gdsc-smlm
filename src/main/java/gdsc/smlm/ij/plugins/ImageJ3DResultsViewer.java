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
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.View;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

import customnode.CustomMeshNode;
import customnode.CustomTransparentTriangleMesh;
import gdsc.core.data.DataException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.logging.Ticker;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;

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
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.GUIProtos.Image3DDrawingMode;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings.Builder;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettingsOrBuilder;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.XYZRResultProcedure;
import ij.IJ;
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
		"Low Resolution Sphere (F=80)",
		"High Resolution Sphere (F=320)",
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
			case 0:	return 3;
			case 1:	return 6;
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

	// No ned to store this in settings as when the plugin is first run there are no windows 
	private static String lastWindow = "";

	private Image3DUniverse univ;
	private JMenuItem resetRotation;
	private JMenuItem resetTranslation;
	private JMenuItem resetZoom;
	private JMenuItem resetAll;
	private JMenuItem changeColour;
	private JMenuItem toggleShading;
	private JMenuItem resetSelectedView;

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
		gd.addChoice("Window", titles, lastWindow);
		gd.addSlider("Transparancy", 0, 0.9, settings.getTransparency());
		gd.addChoice("Colour", LUTHelper.luts, settings.getLut());
		gd.addChoice("Rendering", RENDERING, settings.getRendering());
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
						if (mode == Image3DDrawingMode.DRAW_3D_FIXED_SIZE_VALUE)
						{
							settings.setSize(egd.getNextNumber());
						}
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

		SettingsManager.writeSettings(settings);
		MemoryPeakResults results = ResultsManager.loadInputResults(name, false, null, null);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		// Determine if the drawing mode is supported and compute the point size
		final float[] sphereSize = createSphereSize(results, settings);
		if (sphereSize == null)
			return;

		// Create a 3D viewer.
		if (windowChoice == 0)
			univ = createImage3DUniverse(title, titleList);
		else
			univ = univList.get(windowChoice - 1); // Ignore the new window

		// Adapted from Image3DUniverse.addIcospheres:
		// We create a grainy unit sphere an the origin. 
		// This is then scaled and translated for each localisation.
		final List<Point3f> point = createLocalisationObject(settings.getRendering());
		final List<Point3f> allPoints = new TurboList<Point3f>();

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
				return;
		}

		IJ.showStatus("Creating 3D objects ...");
		final Ticker ticker = Ticker.createStarted(new IJTrackProgress(), results.size(), false);
		results.forEach(DistanceUnit.NM, new XYZRResultProcedure()
		{
			int i = 0; // For the sphere size array

			//int MAX_SIZE = 1; // For debugging
			public void executeXYZR(float x, float y, float z, PeakResult result)
			{
				// Note sure what the limits are for the graphics library.
				// Assume it can support a max array.
				if (allPoints.size() > MAX_SIZE)
					return;
				allPoints.addAll(
						copyScaledTranslated(point, sphereSize[i], sphereSize[i + 1], sphereSize[i + 2], x, y, z));
				i += 3;
				ticker.tick();
			}
		});
		ticker.stop();

		IJ.showStatus("Creating 3D mesh ...");

		// Default color
		final Color3f color = new Color3f(1, 1, 1);

		float transparency = getTransparency(settings);
		CustomTransparentTriangleMesh mesh;
		//mesh = new CustomTransparentTriangleMesh(allPoints, color, transparency);
		// This avoids computing the volume
		mesh = new CustomTransparentTriangleMesh(null, color, transparency);
		mesh.setMesh(allPoints);

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

		IJ.showStatus("");
	}

	private float[] createSphereSize(MemoryPeakResults results, Builder settings)
	{
		// Store XYZ size for each localisation
		switch (settings.getDrawingMode())
		{
			case DRAW_3D_FIXED_SIZE:
				final float size = (settings.getSize() > 0) ? (float) settings.getSize() : 1f;
				return SimpleArrayUtils.newFloatArray(results.size() * 3, size);
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

	private float[] createSphereSizeFromDeviations(MemoryPeakResults results)
	{
		if (!results.hasDeviations())
		{
			IJ.error(TITLE, "The results have no deviations");
			return null;
		}
		// Currently the rendering is in nm
		final TypeConverter<DistanceUnit> dc = results.getDistanceConverter(DistanceUnit.NM);

		final float[] size = new float[results.size() * 3];
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
				size[i++] = dc.convert(x);
				size[i++] = dc.convert(y);
				size[i++] = dc.convert(z);
				return false;
			}
		});
		return (failed) ? null : size;
	}

	private float[] createSphereSizeFromPrecision(MemoryPeakResults results)
	{
		PrecisionResultProcedure p = new PrecisionResultProcedure(results);
		try
		{
			PrecisionMethod m = p.getPrecision();
			IJ.log("Using precision method " + FitProtosHelper.getName(m));
			final float[] size = new float[results.size() * 3];
			for (int i = 0, j = 0; i < p.precision.length; i++)
			{
				// Precision is in NM which matches the rendering
				final float v = (float) p.precision[i];
				size[j++] = v;
				size[j++] = v;
				size[j++] = v;
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

	private static void changeColour(CustomTransparentTriangleMesh mesh, MemoryPeakResults results,
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

		Image3DUniverse univ = new Image3DUniverse();
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
					return -1;
				settings.setTransparency(gd.getNextNumber());
				settings.setLut(gd.getNextChoiceIndex());
				SettingsManager.writeSettings(settings);
			}

			// Restore original object rendering
			settings.setRendering(data.settings.getRendering());
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
