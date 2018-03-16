package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.scijava.java3d.Appearance;
import org.scijava.java3d.BranchGroup;
import org.scijava.java3d.Canvas3D;
import org.scijava.java3d.ColoringAttributes;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.IndexedGeometryArray;
import org.scijava.java3d.LineAttributes;
import org.scijava.java3d.PickInfo;
import org.scijava.java3d.PickInfo.IntersectionInfo;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.SceneGraphPath;
import org.scijava.java3d.Shape3D;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.TriangleArray;
import org.scijava.java3d.View;
import org.scijava.java3d.utils.pickfast.PickCanvas;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point2d;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3d;
import org.scijava.vecmath.Vector3f;

import customnode.CustomLineMesh;
import customnode.CustomMesh;
import customnode.CustomMeshNode;
import customnode.CustomPointMesh;
import gdsc.core.data.DataException;
import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.ij.roi.RoiTest;
import gdsc.core.ij.roi.RoiTestFactory;
import gdsc.core.logging.Ticker;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Sort;
import gdsc.core.utils.TextUtils;

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
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings.Builder;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettingsOrBuilder;
import gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.ij.ij3d.CustomContent;
import gdsc.smlm.ij.ij3d.CustomMeshHelper;
import gdsc.smlm.ij.ij3d.ItemGeometryGroup;
import gdsc.smlm.ij.ij3d.ItemGeometryNode;
import gdsc.smlm.ij.ij3d.ItemPointMesh;
import gdsc.smlm.ij.ij3d.ItemShape;
import gdsc.smlm.ij.ij3d.ItemTriangleMesh;
import gdsc.smlm.ij.ij3d.OrderedItemGeometryGroup;
import gdsc.smlm.ij.ij3d.TransparentItemPointMesh;
import gdsc.smlm.ij.ij3d.TransparentItemShape;
import gdsc.smlm.ij.ij3d.TransparentItemTriangleMesh;
import gdsc.smlm.ij.ij3d.UpdateableItemShape;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.results.procedures.PrecisionResultProcedure;
import gdsc.smlm.results.procedures.StandardResultProcedure;
import gdsc.smlm.results.procedures.XYResultProcedure;
import gdsc.smlm.results.procedures.XYZResultProcedure;
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
import ij3d.DefaultUniverse;
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
	private enum SizeMode implements NamedObject
	{
		FIXED_SIZE { public String getName() { return "Fixed"; }},
		XY_PRECISION { public String getName() { return "XY Precision"; }},
		XYZ_DEVIATIONS { public String getName() { return "XYZ Deviations"; }},
        ;

		public String getShortName()
		{
			return getName();
		}		
		
		public static SizeMode forNumber(int number)
		{
			SizeMode[] values = SizeMode.values();
			if (number < 0 || number >= values.length)
				throw new IllegalArgumentException();
			return values[number];
		}
	};
	
	private final static String[] SIZE_MODE = SettingsManager.getNames((Object[]) SizeMode.values());
	
	private enum Rendering implements NamedObject
	{
		POINT { public String getName() { return "Point"; } 
				public boolean is2D() { return true; }},
		TRIANGLE { public String getName() { return "Triangle"; }
				public boolean is2D() { return true; }},
		SQUARE { public String getName() { return "Square"; }
				public boolean is2D() { return true; }},
		OCTAGON{ public String getName() { return "Octagon"; }
			public boolean is2D() { return true; }},
		LOW_RES_CIRCLE { public String getName() { return "Low resolution circle"; }
			public boolean is2D() { return true; }},
		HIGH_RES_CIRCLE { public String getName() { return "High resolution circle"; }
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

		public static Rendering forNumber(int number)
		{
			Rendering[] values = Rendering.values();
			if (number < 0 || number >= values.length)
				throw new IllegalArgumentException();
			return values[number];
		}
	};

	private final static String[] RENDERING = SettingsManager.getNames((Object[]) Rendering.values());

	private enum DepthMode implements NamedObject
	{
		NONE { public String getName() { return "None"; }},
		INTENSITY { public String getName() { return "Intensity"; }},
		DITHER { public String getName() { return "Dither"; }},
        ;

		public String getShortName()
		{
			return getName();
		}		
		
		public static DepthMode forNumber(int number)
		{
			DepthMode[] values = DepthMode.values();
			if (number < 0 || number >= values.length)
				throw new IllegalArgumentException();
			return values[number];
		}
	};
	
	private final static String[] DEPTH_MODE = SettingsManager.getNames((Object[]) DepthMode.values());
	
	private enum TransparencyMode implements NamedObject
	{
		NONE { public String getName() { return "None"; }},
		SIZE { public String getName() { return "Size"; }},
		XY_PRECISION { public String getName() { return "XY Precision"; }},
		XYZ_DEVIATIONS { public String getName() { return "XYZ Deviations"; }},
        ;

		public String getShortName()
		{
			return getName();
		}		
		
		public static TransparencyMode forNumber(int number)
		{
			TransparencyMode[] values = TransparencyMode.values();
			if (number < 0 || number >= values.length)
				throw new IllegalArgumentException();
			return values[number];
		}
	};
	
	private final static String[] TRANSPARENCY_MODE = SettingsManager.getNames((Object[]) TransparencyMode.values());	

	private enum SortMode implements NamedObject
	{
		NONE { public String getName() { return "None"; }
		public String getDescription() { return ""; }},
		XYZ { public String getName() { return "XYZ"; }
		public String getDescription() { return "Sort using XYZ. The order is defined by the direction with the major component used first, e.g. 1,2,3 for zyx ascending, -3,-2,-1 for xyz descending."; }},
		OTHOGRAPHIC { public String getName() { return "Othographic"; }
		public String getDescription() { return "Project all points to the plane defined by the direction. Rank by distance to the plane."; }},
		PERSPECTIVE { public String getName() { return "Perspective"; }
		public String getDescription() { return "Rank by distance to the eye position for true depth perspective rendering."; }},
        ;

		public String getShortName()
		{
			return getName();
		}		
		
		public static SortMode forNumber(int number)
		{
			SortMode[] values = SortMode.values();
			if (number < 0 || number >= values.length)
				throw new IllegalArgumentException();
			return values[number];
		}

		public abstract String getDescription();
		
		public String getDetails()
		{
			return getName() + ": " + getDescription();
		}
	};
	
	private final static String[] SORT_MODE = SettingsManager.getNames((Object[]) SortMode.values());	
	//@formatter:on	

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
		catch (Throwable t)
		{
			t.printStackTrace();
			version = null;
		}
	}

	private static class ResultsMetaData
	{
		final ImageJ3DResultsViewerSettings settings;

		/** The results when the object mesh was constructed. */
		final MemoryPeakResults results;

		/**
		 * The actual positions of each item in the object mesh. The coordinates may differ from the results as 2D
		 * datasets can be dithered.
		 */
		final TurboList<Point3f> points;

		/** The sizes of each rendered point. */
		final Point3f[] sizes;

		/** The rendering mode. */
		final Rendering rendering;

		/** The point outline for the rendering. */
		final Point3f[] pointOutline;

		public ResultsMetaData(ImageJ3DResultsViewerSettings settings, MemoryPeakResults results,
				TurboList<Point3f> points, Point3f[] sizes)
		{
			this.settings = settings;
			this.results = results;
			this.points = points;
			this.sizes = sizes;

			// Create the point outline for the rendering
			rendering = Rendering.forNumber(settings.getRendering());

			if (settings.getRendering() == 0)
			{
				pointOutline = new Point3f[1];
				pointOutline[0] = new Point3f();
				return;
			}

			List<Point3f> outline;
			if (rendering.is2D())
			{
				// Handle all the 2D objects to create an outline. 
				outline = createLocalisationObjectOutline(rendering);
			}
			else
			{
				// 3D objects can use the same rendering but then post-process to a set of lines
				// and create a line mesh.
				outline = createLocalisationObject(rendering);

				// For outlining
				//Pair<Point3f[], int[]> pair = CustomMeshHelper.createIndexedObject(outline);
				//Point3f[] vertices = pair.s;
				//int[] faces = pair.r;
				//outline.clear();
				//// Only add lines not yet observed. Use an array since 
				//int max = Maths.max(faces) + 1;
				//boolean[] observed = new boolean[max * max];
				//for (int i = 0; i < faces.length; i += 3)
				//{
				//	int t1 = faces[i];
				//	int t2 = faces[i + 1];
				//	int t3 = faces[i + 2];
				//	add(observed, max, t1, t2, vertices, outline);
				//	add(observed, max, t2, t3, vertices, outline);
				//	add(observed, max, t3, t1, vertices, outline);
				//}
			}
			this.pointOutline = outline.toArray(new Point3f[outline.size()]);
		}

		private static void add(boolean[] observed, int max, int t1, int t2, Point3f[] vertices, List<Point3f> point)
		{
			int index = (t1 < t2) ? t1 * max + t2 : t2 * max + t1;
			if (observed[index])
				return;
			observed[index] = true;
			point.add(vertices[t1]);
			point.add(vertices[t2]);
		}

		public CustomMesh createOutline(int index)
		{
			Point3f centre = points.get(index);
			TurboList<Point3f> outline = new TurboList<Point3f>(pointOutline.length);

			if (rendering == Rendering.POINT)
			{
				// The scaling does not matter for the point rendering as the point was (0,0,0)
				Point3f p = pointOutline[0];
				float x = (float) (p.x + centre.x);
				float y = (float) (p.y + centre.y);
				float z = (float) (p.z + centre.z);
				outline.addf(new Point3f(x, y, z));
				ItemPointMesh mesh = new ItemPointMesh(outline, highlightColor, 0);
				//mesh.setPointSize(sizes[0].x); // Handled later
				mesh.getAppearance().getPolygonAttributes().setPolygonMode(PolygonAttributes.POLYGON_LINE);
				return mesh;
			}

			Point3f scale = (sizes.length == 0) ? sizes[0] : sizes[index];
			if (rendering.is2D())
			{
				// Translate and scale the outline
				for (int i = 0; i < pointOutline.length; i++)
				{
					Point3f p = pointOutline[i];
					float x = (float) (p.x * scale.x + centre.x);
					float y = (float) (p.y * scale.y + centre.y);
					float z = (float) (p.z * scale.z + centre.z);
					outline.addf(new Point3f(x, y, z));
				}

				CustomLineMesh mesh = new CustomLineMesh(outline, CustomLineMesh.CONTINUOUS, highlightColor, 0);
				mesh.setAntiAliasing(true);
				mesh.setPattern(LineAttributes.PATTERN_SOLID);
				return mesh;
			}

			// 3D polygon so use a non-filled polygon
			float enlarge = 1.1f;
			scale = new Point3f(scale);
			scale.x *= enlarge;
			scale.y *= enlarge;
			scale.z *= enlarge;
			ItemTriangleMesh mesh = new ItemTriangleMesh(pointOutline, new Point3f[] { centre },
					new Point3f[] { scale }, highlightColor, 0);

			//updateAppearance(mesh, settings);

			// Outline
			mesh.setShaded(false);

			Appearance appearance = mesh.getAppearance();
			PolygonAttributes pa = appearance.getPolygonAttributes();
			pa.setCullFace(PolygonAttributes.CULL_BACK);
			pa.setBackFaceNormalFlip(false);
			final ColoringAttributes ca = appearance.getColoringAttributes();
			ca.setShadeModel(ColoringAttributes.SHADE_FLAT);

			return mesh;
		}
	}

	final static Transform3D IDENTITY = new Transform3D();

	// No ned to store this in settings as when the plugin is first run there are no windows 
	private static String lastWindow = "";

	private static HashMap<String, IJTablePeakResults> resultsTables = new HashMap<String, IJTablePeakResults>();
	private static ResultsTableSettings.Builder resultsTableSettings;
	private static Color3f highlightColor = null;

	private Image3DUniverse univ;
	private JMenuItem resetRotation;
	private JMenuItem resetTranslation;
	private JMenuItem resetZoom;
	private JMenuItem resetAll;
	private JMenuItem changeColour;
	private JMenuItem changePointSize;
	private JMenuItem resetSelectedView;
	private JMenuItem findEyePoint;
	private JMenuItem sortBackToFront;
	private JMenuItem sortFrontToBack;
	private JMenuItem colourSurface;
	private JMenuItem toggleTransparent;
	private JMenuItem toggleShaded;
	private JMenuItem updateSettings;
	private JMenuItem cropResults;
	private JMenuItem toggleDynamicTransparency;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		// For testing
		if (Utils.isExtraOptions())
		{
			new ImageJ3DResultsViewerTest().run(arg);
			return;
		}

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
		gd.addSlider("Transparancy", 0, 0.9, settings.getTransparency(), new OptionListener<Double>()
		{
			public boolean collectOptions(Double value)
			{
				return collectOptions(false);
			}

			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog("Transparancy options", null);
				egd.addCheckbox("Support_dynamic_transparency", settings.getSupportDynamicTransparency());
				egd.addCheckbox("Enable_dynamic_transparency", settings.getEnableDynamicTransparency());
				egd.setSilent(silent);
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				settings.setSupportDynamicTransparency(egd.getNextBoolean());
				settings.setEnableDynamicTransparency(egd.getNextBoolean());
				return true;
			}
		});
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
				if (rendering != 0)
					return false;
				egd.addNumericField("Pixel_size", settings.getPixelSize(), 2, 6, "px");
				egd.setSilent(silent);
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				settings.setPixelSize(egd.getNextNumber());
				return true;
			}
		});
		gd.addCheckbox("Shaded", settings.getShaded());
		gd.addChoice("Size_mode", SIZE_MODE, settings.getSizeMode(), new OptionListener<Integer>()
		{
			public boolean collectOptions(Integer value)
			{
				settings.setSizeMode(value);
				return collectOptions(false);
			}

			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog("Size mode options", null);
				SizeMode mode = SizeMode.forNumber(settings.getSizeMode());
				if (mode == SizeMode.FIXED_SIZE)
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
		gd.addChoice("Sort_mode", SORT_MODE, settings.getSortMode(), new OptionListener<Integer>()
		{
			public boolean collectOptions(Integer value)
			{
				settings.setSortMode(value);
				return collectOptions(false);
			}

			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog("Sort mode options", null);
				SortMode mode = SortMode.forNumber(settings.getSortMode());
				if (mode == SortMode.NONE)
					return false;
				egd.addMessage(TextUtils.wrap(
						"Note: The sort mode is used to correctly render transparent objects. For non-transparent objects " +
								"faster rendering is achieved with a reverse sort to put close objects at the front.",
						80));
				egd.addMessage(TextUtils.wrap(mode.getDetails(), 80));
				egd.addMessage("Define the direction of the view");
				egd.addNumericField("Direction_x", settings.getSortDirectionX(), 3, 10, "");
				egd.addNumericField("Direction_y", settings.getSortDirectionY(), 3, 10, "");
				egd.addNumericField("Direction_z", settings.getSortDirectionZ(), 3, 10, "");
				if (mode == SortMode.PERSPECTIVE)
				{
					egd.addMessage("Define the view eye position");
					egd.addNumericField("Eye_x", settings.getSortEyeX(), 3, 10, "nm");
					egd.addNumericField("Eye_y", settings.getSortEyeY(), 3, 10, "nm");
					egd.addNumericField("Eye_z", settings.getSortEyeZ(), 3, 10, "nm");
				}
				egd.setSilent(silent);
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				settings.setSortDirectionX(egd.getNextNumber());
				settings.setSortDirectionY(egd.getNextNumber());
				settings.setSortDirectionZ(egd.getNextNumber());
				if (mode == SortMode.PERSPECTIVE)
				{
					settings.setSortEyeX(egd.getNextNumber());
					settings.setSortEyeY(egd.getNextNumber());
					settings.setSortEyeZ(egd.getNextNumber());
				}
				return true;
			}
		});
		gd.addChoice("Transparency_mode", TRANSPARENCY_MODE, settings.getTransparencyMode(),
				new OptionListener<Integer>()
				{
					public boolean collectOptions(Integer value)
					{
						settings.setTransparencyMode(value);
						return collectOptions(false);
					}

					public boolean collectOptions()
					{
						return collectOptions(true);
					}

					private boolean collectOptions(boolean silent)
					{
						ExtendedGenericDialog egd = new ExtendedGenericDialog("Transparency mode options", null);
						TransparencyMode mode = TransparencyMode.forNumber(settings.getTransparencyMode());
						if (mode == TransparencyMode.NONE)
							return false;
						egd.addSlider("Min_transparancy", 0, 0.95, settings.getMinTransparency());
						egd.addSlider("Max_transparancy", 0, 0.95, settings.getMaxTransparency());
						egd.setSilent(silent);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return false;
						settings.setMinTransparency(egd.getNextNumber());
						settings.setMaxTransparency(egd.getNextNumber());
						return true;
					}
				});
		gd.addMessage("2D options");
		gd.addChoice("Depth_mode", DEPTH_MODE, settings.getDepthMode(), new OptionListener<Integer>()
		{
			public boolean collectOptions(Integer value)
			{
				settings.setDepthMode(value);
				return collectOptions(false);
			}

			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog("Depth mode options", null);
				DepthMode mode = DepthMode.forNumber(settings.getDepthMode());
				if (mode == DepthMode.NONE)
					return false;
				egd.addNumericField("Depth_range", settings.getDepthRange(), 2, 6, "nm");
				if (mode == DepthMode.DITHER)
					egd.addNumericField("Dither_seed", settings.getDitherSeed(), 0);
				egd.setSilent(silent);
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				settings.setDepthRange(egd.getNextNumber());
				if (mode == DepthMode.DITHER)
					settings.setDitherSeed((int) egd.getNextNumber());
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
		settings.setSizeMode(gd.getNextChoiceIndex());
		settings.setSortMode(gd.getNextChoiceIndex());
		settings.setTransparencyMode(gd.getNextChoiceIndex());
		settings.setDepthMode(gd.getNextChoiceIndex());
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

		// Create a 3D viewer.
		if (windowChoice == 0)
			univ = createImage3DUniverse(title, titleList);
		else
			univ = univList.get(windowChoice - 1); // Ignore the new window

		lastWindow = univ.getWindow().getTitle();

		results = results.copy();

		final TurboList<Point3f> points = getPoints(results, settings);

		ResultsMetaData data = new ResultsMetaData(settings.build(), results, points, sphereSize);

		sort(data);

		float[] alpha = createAlpha(results, settings, sphereSize);

		float transparency = getTransparency(settings);

		Color3f[] colors = createColour(results, settings);

		Content content;

		// Build to support true transparency (depends on settings).
		// Currently this is not supported for PointArrays as they require colouring
		// in the coordinate data.
		if (settings.getSupportDynamicTransparency())
		{
			ItemGeometryGroup pointGroup = createItemGroup(settings, sphereSize, points, alpha, transparency, colors);
			if (pointGroup == null)
				return;

			long total = points.size() + getTotalTransparentObjects(univ);

			activateDynamicTransparency(univ, total, settings.getEnableDynamicTransparency());

			IJ.showStatus("Creating 3D content ...");
			content = new CustomContent(name);
			content.getCurrent().display(new ItemGeometryNode(pointGroup));
		}
		else
		{
			CustomMesh mesh = createMesh(settings, points, sphereSize, transparency, alpha);
			if (mesh == null)
				return;

			updateAppearance(mesh, settings);

			IJ.showStatus("Creating 3D mesh colour ...");
			setColour((ItemShape) mesh, colors);

			IJ.showStatus("Creating 3D content ...");
			content = univ.createContent(mesh, name);
		}

		createHighlightColour(settings.getHighlightColour());

		content.setUserData(data);
		// Prevent relative rotation
		content.setLocked(true);

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

	private static Point3f[] createSphereSize(MemoryPeakResults results, Builder settings)
	{
		// Special case for point rendering
		if (settings.getRendering() == 0)
		{
			final float size = getFixedSize(settings.getPixelSize());
			return new Point3f[] { new Point3f(size, size, size) };
		}

		// Store XYZ size for each localisation
		SizeMode mode = SizeMode.forNumber(settings.getSizeMode());
		switch (mode)
		{
			case FIXED_SIZE:
				final float size = getFixedSize(settings.getSize());
				Point3f[] sizes = new Point3f[results.size()];
				Arrays.fill(sizes, new Point3f(size, size, size));
				return sizes;
			case XYZ_DEVIATIONS:
				return createSphereSizeFromDeviations(results);
			case XY_PRECISION:
				return createSphereSizeFromPrecision(results);
			default:
				throw new IllegalStateException("Unknown drawing mode: " + mode);
		}
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

	private static float[] createAlpha(MemoryPeakResults results, Builder settings, Point3f[] sphereSize)
	{
		if (settings.getTransparencyMode() == 0)
			return null;

		double min = Maths.clip(0, 1, settings.getMinTransparency());
		double max = Maths.clip(0, 1, settings.getMaxTransparency());
		if (min == max)
		{
			// No per item transparency
			Utils.log("No per-item transparency as min == max");
			return null;
		}
		// Convert to alpha
		double minA = 1 - max;
		double maxA = 1 - min;

		SizeMode sizeMode = SizeMode.forNumber(settings.getSizeMode());
		TransparencyMode mode = TransparencyMode.forNumber(settings.getTransparencyMode());
		switch (mode)
		{
			case XYZ_DEVIATIONS:
				if (sizeMode != SizeMode.XYZ_DEVIATIONS)
					sphereSize = createSphereSizeFromDeviations(results);
				return createAlpha(minA, maxA, sphereSize);
			case XY_PRECISION:
				if (sizeMode != SizeMode.XY_PRECISION)
					sphereSize = createSphereSizeFromPrecision(results);
				return createAlpha(minA, maxA, sphereSize);
			case SIZE:
				return createAlphaFromSize(results, settings, minA, maxA, sphereSize);
			default:
				throw new IllegalStateException("Unknown transparency mode: " + mode);
		}
	}

	private static float[] createAlphaFromSize(MemoryPeakResults results, Builder settings, double minA, double maxA,
			Point3f[] sphereSize)
	{
		SizeMode sizeMode = SizeMode.forNumber(settings.getSizeMode());
		if (sizeMode == SizeMode.FIXED_SIZE)
		{
			Utils.log("No per-item transparency as size is fixed");
			return null;
		}

		if (settings.getRendering() == 0)
		{
			// No size was created for fixed point rendering so create it now
			switch (sizeMode)
			{
				case XYZ_DEVIATIONS:
					sphereSize = createSphereSizeFromDeviations(results);
					break;
				case XY_PRECISION:
					sphereSize = createSphereSizeFromPrecision(results);
					break;
				default:
					throw new IllegalStateException("Unknown drawing mode: " + sizeMode);
			}
		}

		return createAlpha(minA, maxA, sphereSize);
	}

	private static float[] createAlpha(double minA, double maxA, Point3f[] sphereSize)
	{
		if (sphereSize == null)
			return null;

		double[] d = new double[sphereSize.length];
		for (int i = 0; i < d.length; i++)
		{
			Point3f p = sphereSize[i];
			if (p.x == p.y && p.y == p.z)
			{
				d[i] = 3.0 * p.x * p.x;
			}
			else
			{
				// Use the squared distance. This is the equivalent of the area of the shape projected to 2D.
				d[i] = (double) p.x * p.x + p.y * p.y + p.z * p.z;
			}

			// Use the average radius. This is the equivalent of the mean radius of an enclosing ellipsoid.
			//d[i] = Math.sqrt(d[i] / 3);
		}
		double[] limits = Maths.limits(d);
		double min = limits[0];
		double max = limits[1];
		if (min == max)
		{
			Utils.log("No per-item transparency as size is fixed");
			return null;
		}

		double range = (maxA - minA) / (max - min);
		float[] alpha = new float[d.length];
		for (int i = 0; i < alpha.length; i++)
		{
			// Largest distance has lowest alpha (more transparent) 
			alpha[i] = (float) (minA + range * (max - d[i]));
		}
		return alpha;
	}

	private static float getTransparency(ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		float transparency = Maths.clip(0, 1, (float) settings.getTransparency());
		// Do not allow fully transparent objects
		if (transparency == 1)
			transparency = 0;
		return transparency;
	}

	static TurboList<Point3f> getPoints(MemoryPeakResults results, ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		final TurboList<Point3f> points = new TurboList<Point3f>(results.size());
		if (results.is3D())
		{
			results.forEach(DistanceUnit.NM, new XYZResultProcedure()
			{
				public void executeXYZ(float x, float y, float z)
				{
					points.addf(new Point3f(x, y, z));
				}
			});
		}
		else
		{
			results.forEach(DistanceUnit.NM, new XYResultProcedure()
			{

				public void executeXY(float x, float y)
				{
					points.addf(new Point3f(x, y, 0));
				}
			});

			final double range = settings.getDepthRange();
			if (range > 0 && results.size() > 1)
			{
				DepthMode mode = DepthMode.forNumber(settings.getDepthMode());
				final double min = -settings.getDepthRange() / 2;
				switch (mode)
				{
					case DITHER:
						final RandomGenerator r = new Well19937c(settings.getDitherSeed());
						for (int i = points.size(); i-- > 0;)
						{
							points.getf(i).z += (min + r.nextDouble() * range);
						}
						break;
					case INTENSITY:
						// Rank by intensity, highest first
						StandardResultProcedure p = new StandardResultProcedure(results);
						p.getI();
						int[] indices = SimpleArrayUtils.newArray(results.size(), 0, 1);
						Sort.sort(indices, p.intensity);
						double inc = range / indices.length;
						for (int i = 0; i < indices.length; i++)
						{
							// The standard rendering has +z going away so put the highest rank at min
							points.getf(indices[i]).z += (min + i * inc);
						}
						break;
					case NONE:
						break;
					default:
						throw new IllegalStateException("Unknown depth mode: " + mode);
				}
			}
		}
		return points;
	}

	private void sort(ResultsMetaData data)
	{
		SortMode mode = SortMode.forNumber(data.settings.getSortMode());
		switch (mode)
		{
			case NONE:
				return;
			case PERSPECTIVE:
				sortPerspective(data);
				break;
			case OTHOGRAPHIC:
				sortOrthographic(data);
				break;
			case XYZ:
				sortXYZ(data);
				break;
			default:
				throw new IllegalStateException("Unknown sort mode " + mode);
		}
	}

	private void sortPerspective(ResultsMetaData data)
	{
		ImageJ3DResultsViewerSettingsOrBuilder settings = data.settings;
		Vector3d direction = getViewDirection(data.settings);
		if (direction == null)
			throw new IllegalStateException("The view direction is not valid");
		Point3d eye = new Point3d(settings.getSortEyeX(), settings.getSortEyeY(), settings.getSortEyeZ());

		double[] d = getDistance(data.points, direction, eye);

		int[] indices = SimpleArrayUtils.newArray(d.length, 0, 1);
		Sort.sort(indices, d);

		reorder(indices, data);
	}

	private double[] getDistance(TurboList<Point3f> points, Vector3d direction, Point3d eye)
	{
		//System.out.printf("Dir %s : Eye %s\n", v, eye);
		double[] d = new double[points.size()];
		for (int i = 0; i < d.length; i++)
		{
			Point3f p = points.getf(i);

			Vector3d v2 = new Vector3d(p.x - eye.x, p.y - eye.y, p.z - eye.z);

			// Compute distance of all points from the eye.
			d[i] = v2.length();

			// We need to know if point is in-front of the eye or behind.
			// Compute dot product (if positive then this is an acute angle).
			// We use a descending sort so acute angles (in front of the eye)
			// should be higher (ranked first) and obtuse angles (behind the 
			// eye) ranked later.
			if (v2.dot(direction) < 0)
				d[i] = -d[i];

			//System.out.printf("[%d] %s %s %g = %g\n", i, p, v2, v2.dot(v), d[i]);
		}
		return d;
	}

	private Vector3d getViewDirection(ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		Vector3d dir = new Vector3d(settings.getSortDirectionX(), settings.getSortDirectionY(),
				settings.getSortDirectionZ());
		double l1 = dir.lengthSquared();
		if (!Maths.isFinite(l1))
			return null;
		return dir;
	}

	private void reorder(int[] indices, ResultsMetaData data)
	{
		MemoryPeakResults results = data.results;
		TurboList<Point3f> points = data.points;

		PeakResult[] originalPeakResults = results.toArray();
		Point3f[] originalPoints = points.toArray(new Point3f[points.size()]);

		// We need another array to store the output 
		PeakResult[] peakResults = new PeakResult[originalPeakResults.length];

		// Rewrite order
		for (int i = 0; i < indices.length; i++)
		{
			int index = indices[i];
			points.setf(i, originalPoints[index]);
			peakResults[i] = originalPeakResults[index];
		}

		// Bulk update the results
		results.setSortAfterEnd(false);
		results.begin();
		results.addAll(peakResults);
		results.end();

		Point3f[] sizes = data.sizes;
		if (sizes.length == indices.length)
		{
			Point3f[] originalSizes = sizes.clone();
			// Rewrite order
			for (int i = 0; i < indices.length; i++)
			{
				int index = indices[i];
				sizes[i] = originalSizes[index];
			}
		}
	}

	private void sortOrthographic(ResultsMetaData data)
	{
		Vector3d v = getViewDirection(data.settings);
		if (v == null)
		{
			// Default
			v = new Vector3d(0, 0, -1);
		}
		v.normalize();
		final double a = v.x;
		final double b = v.y;
		final double c = v.z;

		TurboList<Point3f> points = data.points;

		double[] d = new double[points.size()];
		for (int i = 0; i < d.length; i++)
		{
			Point3f p = points.getf(i);

			// Compute signed distance of all points from the plane
			// defined by normal v and point (0,0,0)
			d[i] = a * p.x + b * p.y + c * p.z;
		}

		int[] indices = SimpleArrayUtils.newArray(d.length, 0, 1);
		Sort.sort(indices, d);

		reorder(indices, data);
	}

	private class CustomSortObject implements Comparable<CustomSortObject>
	{
		final float f1, f2, f3;
		final int index;

		CustomSortObject(int i, float f1, float f2, float f3)
		{
			this.f1 = f1;
			this.f2 = f2;
			this.f3 = f3;
			this.index = i;
		}

		public int compareTo(CustomSortObject o)
		{
			if (f1 < o.f1)
				return -1;
			if (f1 > o.f1)
				return 1;
			if (f2 < o.f2)
				return -1;
			if (f2 > o.f2)
				return 1;
			if (f3 < o.f3)
				return -1;
			if (f3 > o.f3)
				return 1;
			return 0;
		}
	}

	private void sortXYZ(ResultsMetaData data)
	{
		Vector3d v = getViewDirection(data.settings);
		if (v == null)
		{
			// Default to z axis (away), then y then x in ascending order
			v = new Vector3d(-1, -2, 3);
		}
		v.normalize();

		// Use the vector lengths in each dimension to set the order
		int[] indices = SimpleArrayUtils.newArray(3, 0, 1);
		double[] values = new double[] { v.x, v.y, v.z };
		double[] absValues = new double[] { Math.abs(v.x), Math.abs(v.y), Math.abs(v.z) };
		Sort.sort(indices, absValues);

		final int ix = search(indices, 0);
		final int iy = search(indices, 1);
		final int iz = search(indices, 2);

		// Use the vector signs to set the direction
		final int sx = (values[0] <= 0) ? 1 : -1;
		final int sy = (values[1] <= 0) ? 1 : -1;
		final int sz = (values[2] <= 0) ? 1 : -1;

		// Sort using the points since these have dithered positions for 2D results.
		TurboList<Point3f> points = data.points;
		CustomSortObject[] toSort = new CustomSortObject[points.size()];
		float[] f = new float[3];
		for (int i = 0; i < toSort.length; i++)
		{
			Point3f p = points.getf(i);
			f[ix] = sx * p.x;
			f[iy] = sy * p.y;
			f[iz] = sz * p.z;
			toSort[i] = new CustomSortObject(i, f[0], f[1], f[2]);
		}

		Arrays.sort(toSort);

		indices = new int[toSort.length];
		for (int i = 0; i < toSort.length; i++)
			indices[i] = toSort[i].index;

		reorder(indices, data);
	}

	private static int search(int[] values, int key)
	{
		for (int i = 0; i < values.length; i++)
			if (values[i] == key)
				return i;
		return -1;
	}

	private static void updateAppearance(CustomMesh mesh, final ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		mesh.setShaded(settings.getShaded());

		Appearance appearance = mesh.getAppearance();
		PolygonAttributes pa = appearance.getPolygonAttributes();

		// For all 3D polygons we want to support a true face orientation so transparency works
		Rendering r = Rendering.forNumber(settings.getRendering());
		if (r.is2D())
		{
			pa.setCullFace(PolygonAttributes.CULL_NONE);
			pa.setBackFaceNormalFlip(true);
		}
		else
		{
			pa.setCullFace(PolygonAttributes.CULL_BACK);
			pa.setBackFaceNormalFlip(false);
		}

		//TransparencyAttributes ta = appearance.getTransparencyAttributes();
		//ta.setSrcBlendFunction(TransparencyAttributes.BLEND_SRC_ALPHA);
		//ta.setDstBlendFunction(TransparencyAttributes.BLEND_ONE);
		//ta.setDstBlendFunction(TransparencyAttributes.BLEND_ONE_MINUS_SRC_ALPHA); // Default

		ItemTriangleMesh.setTransparencyMode(TransparencyAttributes.FASTEST);
		//ItemTriangleMesh.setTransparencyMode(TransparencyAttributes.SCREEN_DOOR);
		//ItemTriangleMesh.setTransparencyMode(TransparencyAttributes.BLENDED);

		final ColoringAttributes ca = appearance.getColoringAttributes();
		if (r.isHighResolution() || r.is2D())
			// Smooth across vertices. Required to show 2D surfaces smoothly
			ca.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
		else
			// Faster polygon rendering with flat shading
			ca.setShadeModel(ColoringAttributes.SHADE_FLAT);
	}

	private static HashMap<String, Color3f> colours;
	static
	{
		colours = new HashMap<String, Color3f>();
		Field[] fields = Color.class.getFields();
		for (Field field : fields)
		{
			if (Modifier.isStatic(field.getModifiers()) && field.getType() == Color.class)
			{
				try
				{
					Color c = (Color) field.get(null);
					colours.put(field.getName().toLowerCase(), new Color3f(c));
				}
				catch (IllegalArgumentException e1)
				{
				}
				catch (IllegalAccessException e1)
				{
				}
			}
		}
	}

	private static void createHighlightColour(String highlightColour)
	{
		highlightColor = null;
		for (int i = 0; i < highlightColour.length(); i++)
		{
			if (Character.isDigit(highlightColour.charAt(i)))
			{
				// Try and extract RGB
				String[] split = highlightColour.split("[\\s,:]+");
				if (split.length > 2)
				{
					try
					{
						// RGB
						int red = Integer.parseInt(split[0]);
						int green = Integer.parseInt(split[1]);
						int blue = Integer.parseInt(split[2]);
						highlightColor = new Color3f(new Color(red, green, blue));
						return;
					}
					catch (NumberFormatException e)
					{

					}
				}
			}
		}
		highlightColor = colours.get(highlightColour.toLowerCase());
	}

	private static Color3f[] createColour(MemoryPeakResults results, ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		// Colour by z
		LUT lut = LUTHelper.createLUT(LutColour.forNumber(settings.getLut()), false);

		StandardResultProcedure p = null;
		float range = 0;
		float[] limits = null;
		if (results.is3D())
		{
			p = new StandardResultProcedure(results);
			p.getZ();
			limits = Maths.limits(p.z);
			range = limits[1] - limits[0];
		}

		if (range == 0)
		{
			return new Color3f[] { new Color3f(new Color(lut.getRGB(255))) };
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

			final float minimum = limits[0];
			Color3f[] allColors = new Color3f[results.size()];
			for (int i = 0, size = results.size(); i < size; i++)
			{
				float value = p.z[i];
				value = value - minimum;
				if (value < 0f)
					value = 0f;
				int ivalue = (int) ((value * scale) + 0.5f);
				if (ivalue > 255)
					ivalue = 255;
				allColors[i] = colors[ivalue];
			}
			return allColors;
		}
	}

	private static void changeColour(ItemShape itemShape, MemoryPeakResults results,
			ImageJ3DResultsViewerSettingsOrBuilder settings)
	{
		setColour(itemShape, createColour(results, settings));
	}

	private static void setColour(ItemShape itemShape, Color3f[] colours)
	{
		if (colours.length == 1)
			itemShape.setItemColor(colours[0]);
		else
			itemShape.setItemColor(colours);
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

		univ.addUniverseListener(this);

		univ.setShowBoundingBoxUponSelection(false);
		univ.showAttribute(DefaultUniverse.ATTRIBUTE_SCALEBAR, false);

		// Capture a canvas mouse click/region and identify the coordinates.
		final ImageCanvas3D canvas = (ImageCanvas3D) univ.getCanvas();
		final BranchGroup scene = univ.getScene();

		MouseListener mouseListener = new MouseAdapter()
		{
			@Override
			public void mouseClicked(final MouseEvent e)
			{
				//System.out.println("plugin mouseClicked");
				if (!consumeEvent(e))
					return;

				// Consume the event
				e.consume();

				String name = "Clicked localisation";

				// This finds the vertex indices of the rendered object.
				Pair<Content, IntersectionInfo> pair = getPickedContent(canvas, scene, e.getX(), e.getY());
				if (pair == null)
				{
					univ.select(null); // Do the same as the mouseClicked in Image3DUniverse
					univ.removeContent(name);
					return;
				}

				// Only process content added from localisations
				Content c = pair.s;
				if (!(c.getUserData() instanceof ResultsMetaData))
				{
					univ.select(c); // Do the same as the mouseClicked in Image3DUniverse
					return;
				}

				ResultsMetaData data = (ResultsMetaData) c.getUserData();

				MemoryPeakResults results = data.results;

				// Look up the localisation from the clicked vertex
				final ContentInstant content = c.getCurrent();
				int index = -1;
				if (content.getContent() instanceof CustomMeshNode)
				{
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

					int nPerLocalisation = nVertices / results.size();

					// Determine the localisation
					int vertexIndex = pair.r.getVertexIndices()[0];
					index = vertexIndex / nPerLocalisation;
					//System.out.printf("n=%d [%d]  %s  %s\n", nPerLocalisation, index,
					//		Arrays.toString(pair.r.getVertexIndices()), pair.r.getIntersectionPoint());
				}
				else if (content.getContent() instanceof ItemGeometryNode)
				{
					//ItemGeometryNode node = (ItemGeometryNode) content.getContent();
					//ItemGeometryGroup g = node.getItemGeometry();
					// All shapes have the index as the user data
					Object o = pair.r.getGeometry().getUserData();
					if (o instanceof Integer)
						index = (Integer) pair.r.getGeometry().getUserData();
				}
				if (index == -1)
					return;

				PeakResult p = results.get(index);

				if (e.getClickCount() > 1)
				{
					// Centre on the localisation
					Point3d coordinate = new Point3d();
					//ga.getCoordinate(vertexIndex, coordinate);
					coordinate.set(data.points.get(index));

					// Handle the local transform of the content ...
					final Transform3D vWorldToLocal = getVworldToLocal(content);
					//vWorldToLocal.invert();
					vWorldToLocal.transform(coordinate);

					univ.centerAt(coordinate);
				}
				else
				{
					// Output the result to a table.
					// Just create a table and add to it.
					IJTablePeakResults table = createTable(results);
					table.add(p);
					table.flush();

					// Highlight the localisation using an outline.
					CustomMesh pointOutline = data.createOutline(index);
					if (pointOutline instanceof CustomPointMesh)
					{
						// Make the point size correct as points are resizeable
						float pointSize = 1;
						if (content.getContent() instanceof CustomMeshNode)
						{
							CustomMeshNode node = (CustomMeshNode) content.getContent();
							CustomMesh mesh = node.getMesh();
							pointSize = ((CustomPointMesh) mesh).getPointSize();
						}
						else if (content.getContent() instanceof ItemGeometryNode)
						{
							ItemGeometryNode node = (ItemGeometryNode) content.getContent();
							ItemGeometryGroup g = node.getItemGeometry();
							pointSize = g.getPointSize();
						}
						((CustomPointMesh) pointOutline).setPointSize(pointSize);
					}

					univ.setAutoAdjustView(false);

					pointOutline.setPickable(false);

					Content outline = univ.createContent(pointOutline, name);
					outline.getCurrent().setTransform(getVworldToLocal(content));
					outline.setLocked(true);

					univ.removeContent(name);
					univ.addContent(outline);
					univ.setAutoAdjustView(true);
				}
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
				// Not needed now the content is locked to prevent annoying local rotation
				//if (consumeEvent(e))
				//{
				//	// Image3DUniverse adds a canvas mouse listener that calls showPopup.
				//	// This selects the content. It ignores the e.isConsumed() flag so
				//	// we translate the point offscreen.
				//	e.consume();
				//	e.translatePoint(1000000, 1000000);
				//}
			}

			@Override
			public void mouseReleased(final MouseEvent e)
			{
				//System.out.println("plugin mouseReleased");
				// Not needed now the content is locked to prevent annoying local rotation
				//if (consumeEvent(e))
				//{
				//	// Image3DUniverse adds a canvas mouse listener that calls showPopup.
				//	// This selects the content. It ignores the e.isConsumed() flag so
				//	// we translate the point offscreen.
				//	e.consume();
				//	e.translatePoint(1000000, 1000000);
				//}
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
		};

		// 0 = ImageCanvas3D
		// 1 = DefaultUniverse
		// 2 = Image3DUniverse
		MouseListener[] l = canvas.getMouseListeners();
		for (int i = 0; i < l.length; i++)
		{
			if (l[i].getClass().getName().contains("Image3DUniverse"))
			{
				// We want to be before the Image3DUniverse to allow consuming the click event.
				// Only allow the click event. 
				// This disables the right-click pop-up menu. 
				// It doesn't have anything of use for localisations anyway. 
				canvas.removeMouseListener(l[i]);
				canvas.addMouseListener(mouseListener);
				canvas.addMouseListener(new MouseListenerWrapper(l[i], MouseListenerWrapper.MOUSE_CLICKED));
			}
		}

		// 0 = ImageCanvas3D
		// 1 = DefaultUniverse
		// 2 = Image3DUniverse
		// 3 = EventCatcher	(from scijava)	
		MouseMotionListener[] ml = canvas.getMouseMotionListeners();
		for (int i = 0; i < ml.length; i++)
		{
			if (ml[i].getClass().getName().contains("Image3DUniverse"))
			{
				// Ignore this as it just shows the name in the IJ status bar 
				canvas.removeMouseMotionListener(ml[i]);
			}
		}

		// Finally display the window
		univ.show();
		ImageWindow3D w = univ.getWindow();
		GUI.center(w);
		w.setTitle(title2);

		// Add a new menu for SMLM functionality
		Image3DMenubar menubar = (Image3DMenubar) univ.getMenuBar();
		JMenu menu = createSMLMMenuBar();
		menubar.add(menu);
		// Add back so it is redrawn
		univ.setMenubar(menubar);

		return univ;
	}

	private static class MouseListenerWrapper implements MouseListener
	{
		static final int MOUSE_CLICKED = 0x01;
		static final int MOUSE_PRESSED = 0x02;
		static final int MOUSE_RELEASED = 0x04;
		static final int MOUSE_ENTERED = 0x08;
		static final int MOUSE_EXITED = 0x10;

		final MouseListener l;
		final int flags;

		MouseListenerWrapper(MouseListener l, int flags)
		{
			this.l = l;
			this.flags = flags;
		}

		private boolean run(MouseEvent e, int flag)
		{
			return (flags & flag) != 0 && !e.isConsumed();
		}

		public void mouseClicked(MouseEvent e)
		{
			if (run(e, MOUSE_CLICKED))
				l.mouseClicked(e);
		}

		public void mousePressed(MouseEvent e)
		{
			if (run(e, MOUSE_PRESSED))
				l.mousePressed(e);
		}

		public void mouseReleased(MouseEvent e)
		{
			if (run(e, MOUSE_RELEASED))
				l.mouseReleased(e);
		}

		public void mouseEntered(MouseEvent e)
		{
			if (run(e, MOUSE_ENTERED))
				l.mouseEntered(e);
		}

		public void mouseExited(MouseEvent e)
		{
			if (run(e, MOUSE_EXITED))
				l.mouseExited(e);
		}
	}

	@SuppressWarnings("unused")
	private static class MouseMotionListenerWrapper implements MouseMotionListener
	{
		static final int MOUSE_DRAGGED = 0x01;
		static final int MOUSE_MOVED = 0x02;

		final MouseMotionListener l;
		final int flags;

		MouseMotionListenerWrapper(MouseMotionListener l, int flags)
		{
			this.l = l;
			this.flags = flags;
		}

		private boolean run(MouseEvent e, int flag)
		{
			return (flags & flag) != 0 && !e.isConsumed();
		}

		public void mouseDragged(MouseEvent e)
		{
			if (run(e, MOUSE_DRAGGED))
				l.mouseDragged(e);
		}

		public void mouseMoved(MouseEvent e)
		{
			if (run(e, MOUSE_MOVED))
				l.mouseMoved(e);
		}
	}

	private static long getTotalTransparentObjects(Image3DUniverse univ)
	{
		long size = 0;
		for (@SuppressWarnings("unchecked")
		Iterator<Content> it = univ.contents(); it.hasNext();)
		{
			Content c = it.next();
			if (!(c.getUserData() instanceof ResultsMetaData))
				return 0;
			final ContentInstant content = c.getCurrent();
			if (content.getContent() instanceof ItemGeometryNode)
			{
				ItemGeometryNode node = (ItemGeometryNode) content.getContent();
				ItemGeometryGroup g = node.getItemGeometry();
				if (g instanceof OrderedItemGeometryGroup)
					continue;
				size += g.size();
			}
		}
		return size;
	}

	private static void activateDynamicTransparency(Image3DUniverse univ, long size, boolean enable)
	{
		View view = univ.getViewer().getView();
		boolean isEnabled = view.getTransparencySortingPolicy() == View.TRANSPARENCY_SORT_GEOMETRY;
		if (enable == isEnabled)
			return;

		if (enable)
		{
			if (size > 20000L)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
				egd.addMessage("The results contain " + size +
						" transparent objects.\nDynamic transparency may take a long time to render.");
				egd.setOKLabel("Continue");
				egd.showDialog();
				if (egd.wasCanceled())
					return;
			}

			view.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);
			// I am not sure if this is required if objects are sorted.
			//view.setDepthBufferFreezeTransparent(false);
		}
		else
		{
			view.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_NONE);
		}
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
		add.setMnemonic(KeyEvent.VK_G);

		resetRotation = new JMenuItem("Reset global rotation", KeyEvent.VK_R);
		resetRotation.addActionListener(this);
		add.add(resetRotation);

		resetTranslation = new JMenuItem("Reset global translation", KeyEvent.VK_T);
		resetTranslation.addActionListener(this);
		add.add(resetTranslation);

		resetZoom = new JMenuItem("Reset global zoom", KeyEvent.VK_Z);
		resetZoom.addActionListener(this);
		add.add(resetZoom);

		add.addSeparator();

		resetAll = new JMenuItem("Reset all transformations", KeyEvent.VK_A);
		resetAll.addActionListener(this);
		add.add(resetAll);

		resetSelectedView = new JMenuItem("Reset selected transformation", KeyEvent.VK_S);
		resetSelectedView.addActionListener(this);
		add.add(resetSelectedView);

		findEyePoint = new JMenuItem("Find eye point", KeyEvent.VK_E);
		findEyePoint.addActionListener(this);
		add.add(findEyePoint);

		sortBackToFront = new JMenuItem("Sort Back-to-Front", KeyEvent.VK_B);
		sortBackToFront.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed B"));
		sortBackToFront.addActionListener(this);
		add.add(sortBackToFront);

		sortFrontToBack = new JMenuItem("Sort Front-to-Back", KeyEvent.VK_F);
		sortFrontToBack.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed R"));
		sortFrontToBack.addActionListener(this);
		add.add(sortFrontToBack);

		add.addSeparator();

		changeColour = new JMenuItem("Change colour", KeyEvent.VK_O);
		changeColour.addActionListener(this);
		add.add(changeColour);

		changePointSize = new JMenuItem("Change point size", KeyEvent.VK_H);
		changePointSize.addActionListener(this);
		add.add(changePointSize);

		toggleTransparent = new JMenuItem("Toggle transparent", KeyEvent.VK_P);
		toggleTransparent.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed E"));
		toggleTransparent.addActionListener(this);
		add.add(toggleTransparent);

		toggleShaded = new JMenuItem("Toggle shaded", KeyEvent.VK_S);
		toggleShaded.addActionListener(this);
		add.add(toggleShaded);

		toggleDynamicTransparency = new JMenuItem("Toggle dynamic transparency", KeyEvent.VK_D);
		toggleDynamicTransparency.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed D"));
		toggleDynamicTransparency.addActionListener(this);
		add.add(toggleDynamicTransparency);

		colourSurface = new JMenuItem("Colour surface from 2D image", KeyEvent.VK_I);
		colourSurface.addActionListener(this);
		add.add(colourSurface);

		add.addSeparator();

		cropResults = new JMenuItem("Crop results", KeyEvent.VK_C);
		cropResults.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed X"));
		cropResults.addActionListener(this);
		add.add(cropResults);

		add.addSeparator();

		updateSettings = new JMenuItem("Update settings", KeyEvent.VK_U);
		updateSettings.addActionListener(this);
		add.add(updateSettings);

		return add;
	}

	private interface ContentAction
	{
		/**
		 * Run the action
		 * 
		 * @param c
		 *            The content
		 * @return negative for error. No further content can be processed.
		 */
		public int run(Content c);

		public void finish();
	}

	private static abstract class BaseContentAction implements ContentAction
	{
		public void finish()
		{
		}
	}

	private static class ChangeColourContentAction extends BaseContentAction
	{
		ImageJ3DResultsViewerSettings.Builder settings = null;

		public int run(Content c)
		{
			if (!(c.getUserData() instanceof ResultsMetaData))
				return 0;

			ResultsMetaData data = (ResultsMetaData) c.getUserData();

			MemoryPeakResults results = data.results;

			// Change the colour
			if (settings == null)
			{
				// Use the latest settings
				settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
				ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
				// Transparency can be set interactively using: Edit > Change Transparency
				//gd.addSlider("Transparancy", 0, 0.9, settings.getTransparency());
				gd.addChoice("Colour", LUTHelper.luts, settings.getLut());
				gd.showDialog();
				if (gd.wasCanceled())
					return -1;
				//settings.setTransparency(gd.getNextNumber());
				settings.setLut(gd.getNextChoiceIndex());
				SettingsManager.writeSettings(settings);
			}

			final ContentInstant content = c.getCurrent();
			ItemShape itemShape = null;
			if (content.getContent() instanceof CustomMeshNode)
			{
				CustomMeshNode node = (CustomMeshNode) content.getContent();
				CustomMesh mesh = node.getMesh();
				if (mesh instanceof ItemShape)
					itemShape = (ItemShape) mesh;
			}
			else if (content.getContent() instanceof ItemGeometryNode)
			{
				ItemGeometryNode node = (ItemGeometryNode) content.getContent();
				itemShape = node.getItemGeometry();
			}
			if (itemShape != null)
				changeColour(itemShape, results, settings);
			return 0;
		}
	}

	private static class ChangePointSizeContentAction extends BaseContentAction
	{
		ImageJ3DResultsViewerSettings.Builder settings = null;

		public int run(Content c)
		{
			final ContentInstant content = c.getCurrent();
			if (content.getContent() instanceof CustomMeshNode)
			{
				CustomMeshNode node = (CustomMeshNode) content.getContent();
				CustomMesh mesh = node.getMesh();
				if (mesh instanceof CustomPointMesh)
				{
					// Change the point size
					if (!getSettings())
						return -1;
					((CustomPointMesh) mesh).setPointSize((float) settings.getPixelSize());
				}
			}
			else if (content.getContent() instanceof ItemGeometryNode)
			{
				if (!getSettings())
					return -1;
				ItemGeometryNode node = (ItemGeometryNode) content.getContent();
				ItemGeometryGroup g = node.getItemGeometry();
				g.setPointSize((float) settings.getPixelSize());
			}

			return 0;
		}

		private boolean getSettings()
		{
			if (settings == null)
			{
				// Use the latest settings
				settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
				ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
				gd.addNumericField("Pixel_size", settings.getPixelSize(), 2, 6, "px");
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				settings.setPixelSize(gd.getNextNumber());
				SettingsManager.writeSettings(settings);
			}
			return true;
		}
	}

	private static class ResetViewContentAction extends BaseContentAction
	{
		final boolean error;

		public ResetViewContentAction(boolean error)
		{
			this.error = error;
		}

		public int run(Content c)
		{
			if (c.isLocked())
			{
				if (error)
				{
					IJ.error(TITLE, c.getName() + " is locked");
					return -1;
				}
				return 0;
			}
			final Transform3D t = new Transform3D();
			c.setTransform(t);
			return 0;
		}
	}

	private static class ToggleTransparentAction extends BaseContentAction
	{
		private static class TransparencyData
		{
			float transparency;
			float[] alpha;

			public void save(CustomMesh mesh)
			{
				transparency = mesh.getTransparency();
				mesh.setTransparency(0);
				if (mesh instanceof TransparentItemShape)
				{
					TransparentItemShape t = (TransparentItemShape) mesh;
					int size = t.size();
					if (alpha == null || alpha.length != size)
						alpha = new float[size];
					t.getItemAlpha(alpha);
					t.setItemAlpha(1);
				}
			}

			public void restore(CustomMesh mesh)
			{
				mesh.setTransparency(transparency);
				if (mesh instanceof TransparentItemShape)
				{
					TransparentItemShape t = (TransparentItemShape) mesh;
					int size = t.size();
					if (alpha != null && alpha.length == size)
						t.setItemAlpha(alpha);
				}
			}

			public void save(ItemGeometryGroup group)
			{
				int size = group.size();
				if (alpha == null || alpha.length != size)
					alpha = new float[size];
				group.getItemAlpha(alpha);
				transparency = group.getTransparency();
				group.setItemAlpha(1f, 0f);
			}

			public void restore(ItemGeometryGroup group)
			{
				int size = group.size();
				if (alpha != null && alpha.length == size)
					group.setItemAlpha(alpha, transparency);
			}
		}

		public int run(Content c)
		{
			if (!(c.getUserData() instanceof ResultsMetaData))
				return 0;
			final ContentInstant content = c.getCurrent();
			if (content.getContent() instanceof CustomMeshNode)
			{
				CustomMeshNode node = (CustomMeshNode) content.getContent();
				CustomMesh mesh = node.getMesh();

				// Polygons can just switch the transparency mode
				TransparencyAttributes ta = mesh.getAppearance().getTransparencyAttributes();
				boolean off;
				if (ta.getTransparencyMode() == TransparencyAttributes.NONE)
				{
					ta.setTransparencyMode(ItemTriangleMesh.getTransparencyMode());
					off = false;
				}
				else
				{
					ta.setTransparencyMode(TransparencyAttributes.NONE);
					off = true;
				}

				// The point mesh does not support the transparency mode switching off.
				// So switch the actual transparency.
				if (mesh instanceof CustomPointMesh)
				{
					TransparencyData d;
					if (mesh.getUserData() instanceof TransparencyData)
					{
						d = (TransparencyData) mesh.getUserData();
					}
					else
					{
						d = new TransparencyData();
						mesh.setUserData(d);
					}

					if (off)
					{
						d.save(mesh);
					}
					else
					{
						d.restore(mesh);
					}
				}
			}
			else if (content.getContent() instanceof ItemGeometryNode)
			{
				ItemGeometryNode node = (ItemGeometryNode) content.getContent();
				ItemGeometryGroup g = node.getItemGeometry();

				TransparencyData d;
				if (g.getUserData() instanceof TransparencyData)
				{
					d = (TransparencyData) g.getUserData();
				}
				else
				{
					d = new TransparencyData();
					g.setUserData(d);
				}

				// All shapes have their own transparency to switch
				if (g.isTransparent())
				{
					d.save(g);
				}
				else
				{
					d.restore(g);
				}
			}

			return 0;
		}
	}

	private static class ToggleShadedAction extends BaseContentAction
	{
		public int run(Content c)
		{
			if (!(c.getUserData() instanceof ResultsMetaData))
				return 0;
			final ContentInstant content = c.getCurrent();
			if (content.getContent() instanceof CustomMeshNode)
			{
				CustomMeshNode node = (CustomMeshNode) content.getContent();
				CustomMesh mesh = node.getMesh();
				mesh.setShaded(!mesh.isShaded());
			}
			else if (content.getContent() instanceof ItemGeometryNode)
			{
				ItemGeometryNode node = (ItemGeometryNode) content.getContent();
				ItemGeometryGroup g = node.getItemGeometry();
				g.setShaded(!g.isShaded());
			}

			return 0;
		}
	}

	private class FindEyePointContentAction extends BaseContentAction
	{
		ImageJ3DResultsViewerSettings.Builder settings;
		final Point3d eyePtInVWorld = new Point3d();
		final Point3d dir0InVWorld = new Point3d();
		final Point3d dir1InVWorld = new Point3d(0, 0, -1);

		Point3d eye;
		Vector3d direction;

		FindEyePointContentAction()
		{
			this.settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
			if (!settings.getSaveEyePoint())
				settings = null;

			final Transform3D ipToVWorld = new Transform3D();
			univ.getCanvas().getImagePlateToVworld(ipToVWorld);
			//ipToVWorldInverse.invert(ipToVWorld);
			//System.out.printf("ipToVWorld\n%s", ipToVWorld);

			univ.getCanvas().getCenterEyeInImagePlate(eyePtInVWorld);
			ipToVWorld.transform(eyePtInVWorld);

			// Work out where the camera is looking in the virtual world
			final Transform3D cameraToVWorld = new Transform3D();
			univ.getVworldToCameraInverse(cameraToVWorld);
			cameraToVWorld.transform(dir0InVWorld);
			cameraToVWorld.transform(dir1InVWorld);
		}

		public int run(Content c)
		{
			final Transform3D vWorldToLocal = getVworldToLocal(c.getCurrent());

			boolean identity = vWorldToLocal.equals(IDENTITY);

			eye = new Point3d(eyePtInVWorld);

			final Point3d dir0InLocalWorld = new Point3d(dir0InVWorld);
			final Point3d dir1InLocalWorld = new Point3d(dir1InVWorld);

			if (!identity)
			{
				// Since we require the eye-point relative to the local world 
				// the matrix is inverted
				vWorldToLocal.invert();
				vWorldToLocal.transform(eye);
				vWorldToLocal.transform(dir0InLocalWorld);
				vWorldToLocal.transform(dir1InLocalWorld);
			}

			direction = new Vector3d();
			direction.sub(dir1InLocalWorld, dir0InLocalWorld);

			// Print the eye coords and direction in the virtual world. 
			// This can be used for a custom sort.
			//Utils.log("%s : Eye point = %s : Direction = %s", c.getName(), eye, direction);
			Rounder rounder = RounderFactory.create(4);
			Utils.log("%s : Eye point = (%s,%s,%s) : Direction = (%s,%s,%s)", c.getName(), rounder.round(eye.x),
					rounder.round(eye.y), rounder.round(eye.z), rounder.round(direction.x), rounder.round(direction.y),
					rounder.round(direction.z));

			if (settings != null)
			{
				settings.setSortEyeX(eye.x);
				settings.setSortEyeY(eye.y);
				settings.setSortEyeZ(eye.z);
				settings.setSortDirectionX(direction.x);
				settings.setSortDirectionY(direction.y);
				settings.setSortDirectionZ(direction.z);
			}

			return 0;
		}

		@Override
		public void finish()
		{
			if (settings != null)
				SettingsManager.writeSettings(settings);
		}
	}

	private class SortContentAction extends FindEyePointContentAction
	{
		final boolean reverse;

		SortContentAction(boolean reverse)
		{
			this.reverse = reverse;
		}

		public int run(Content c)
		{
			int result = super.run(c);
			if (result != 0)
				return result;

			final ContentInstant content = c.getCurrent();
			if (c.getUserData() == null)
				return 0;
			UpdateableItemShape updateable = null;
			boolean reorderData = false;
			if (content.getContent() instanceof CustomMeshNode)
			{
				CustomMeshNode node = (CustomMeshNode) content.getContent();
				CustomMesh mesh = node.getMesh();
				if (!(mesh instanceof UpdateableItemShape))
					return 0;
				updateable = (UpdateableItemShape) mesh;
				reorderData = true;
			}
			else if (content.getContent() instanceof ItemGeometryNode)
			{
				ItemGeometryNode node = (ItemGeometryNode) content.getContent();
				ItemGeometryGroup g = node.getItemGeometry();
				if (!(g instanceof UpdateableItemShape))
					return 0;
				updateable = (UpdateableItemShape) g;
			}

			ResultsMetaData data = (ResultsMetaData) c.getUserData();

			if (reverse)
				direction.negate();

			double[] d = getDistance(data.points, direction, eye);

			int[] indices = SimpleArrayUtils.newArray(d.length, 0, 1);
			Sort.sort(indices, d);

			// Switch to fast mode when not debugging
			//updateable.reorder(indices);
			updateable.reorderFast(indices);

			// We reorder the data that is used to create colours and clicked point size.
			// This is not needed for the ItemGeometryNode as it uses the indices directly.
			// Do this second as points is updated in-line so may break reordering the mesh
			// if it has a reference to the points list (e.g. ItemPointMesh initially uses 
			// the points list but will create a new internal list when it is re-ordered).
			if (reorderData)
				reorder(indices, data);

			return 0;
		}
	}

	private static class ColourSurfaceContentAction extends BaseContentAction
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

			final ContentInstant content = c.getCurrent();
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

	private class CropResultsAction extends BaseContentAction
	{
		RoiTest shape = null;
		private final Point2d p2d = new Point2d();
		ImageJ3DResultsViewerSettings.Builder settings = null;

		public int run(Content c)
		{
			if (!(c.getUserData() instanceof ResultsMetaData))
				return 0;

			Canvas3D canvas = univ.getCanvas();
			if (shape == null)
			{
				shape = RoiTestFactory.create(((ImageCanvas3D) canvas).getRoi());
				if (shape == null)
					return -1;
				settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
			}

			final ContentInstant content = c.getCurrent();

			// Adapted from CustomTriangleMesh.retain
			final Transform3D vWorldToIP = new Transform3D();
			canvas.getImagePlateToVworld(vWorldToIP);
			vWorldToIP.invert();

			final Transform3D vWorldToLocal = getVworldToLocal(content);
			if (!vWorldToIP.equals(IDENTITY))
			{
				vWorldToLocal.invert(); // To make localToVworld
				vWorldToIP.mul(vWorldToLocal); // To make localToIP
			}

			ResultsMetaData data = (ResultsMetaData) c.getUserData();

			// Transform all points to the image plate and test if they are in the ROI 
			MemoryPeakResults results = data.results;
			TurboList<Point3f> points = data.points;
			MemoryPeakResults newResults = new MemoryPeakResults();
			newResults.copySettings(results);
			// Get the output name
			String outputName;
			if (settings.getNameOption() == CropResults.NAME_OPTION_NAME)
			{
				ExtendedGenericDialog egd = new ExtendedGenericDialog("Crop results");
				String name = (TextUtils.isNullOrEmpty(settings.getOutputName())) ? (results.getName() + " Cropped")
						: settings.getOutputName();
				egd.addStringField("Output_name", name, Maths.clip(60, 120, name.length()));
				egd.showDialog();
				if (egd.wasCanceled())
					return -1;
				settings.setOutputName(egd.getNextString());
				outputName = settings.getOutputName();
				if (TextUtils.isNullOrEmpty(outputName))
				{
					IJ.error(TITLE, "No output name");
					return -1;
				}
			}
			else if (settings.getNameOption() == CropResults.NAME_OPTION_SUFFIX)
			{
				String suffix = settings.getNameSuffix();
				if (TextUtils.isNullOrEmpty(suffix))
				{
					IJ.error(TITLE, "No output suffix");
					return -1;
				}
				outputName = results.getName() + suffix;
			}
			else if (settings.getNameOption() == CropResults.NAME_OPTION_SEQUENCE)
			{
				outputName = results.getName();
				String suffix = settings.getNameSuffix();
				if (!TextUtils.isNullOrEmpty(suffix))
				{
					outputName += suffix;
				}
				int counter = settings.getNameCounter();
				outputName += counter;
				settings.setNameCounter(counter + 1); // Increment for next time
			}
			else
			{
				IJ.error(TITLE, "No output name");
				return -1;
			}
			newResults.setName(outputName);
			newResults.begin();
			Utils.showStatus("Cropping " + results.getName());
			Ticker ticker = Ticker.createStarted(new IJTrackProgress(), results.size(), false);
			PeakResult[] allResults = results.toArray();
			for (int i = 0, size = results.size(); i < size; i++)
			{
				Point3d locInImagePlate = new Point3d(points.getf(i));
				vWorldToIP.transform(locInImagePlate);
				canvas.getPixelLocationFromImagePlate(locInImagePlate, p2d);
				if (shape.contains(p2d.x, p2d.y))
					newResults.add(allResults[i]);
				ticker.tick();
			}

			int size = newResults.size();
			IJ.showStatus("Cropped " + TextUtils.pleural(size, "result"));
			ticker.stop();
			newResults.end();
			if (size != 0)
				MemoryPeakResults.addResults(newResults);

			return 0;
		}

		@Override
		public void finish()
		{
			SettingsManager.writeSettings(settings, 0);
		}
	}

	private static Transform3D getVworldToLocal(ContentInstant content)
	{
		final Transform3D vWorldToLocal = new Transform3D();
		Transform3D translate = new Transform3D();
		Transform3D rotate = new Transform3D();
		content.getLocalTranslate(translate);
		content.getLocalRotate(rotate);
		vWorldToLocal.mul(translate, rotate);
		return vWorldToLocal;
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
		if (src == updateSettings)
		{
			final ImageJ3DResultsViewerSettings.Builder settings = SettingsManager.readImageJ3DResultsViewerSettings(0)
					.toBuilder();

			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			ResultsSettings.Builder s = ResultsSettings.newBuilder();
			s.setResultsTableSettings(resultsTableSettings); // This is from the cache
			gd.addMessage("Click on the image to view localisation data.\nCtrl/Alt key must be pressed.");
			final TextField[] tf = new TextField[1];
			gd.addStringField("Highlight_colour", settings.getHighlightColour(), new OptionListener<String>()
			{
				public boolean collectOptions(String value)
				{
					createHighlightColour(value);
					int r, g, b;
					if (highlightColor == null)
					{
						r = b = 0;
						g = 255;
					}
					else
					{
						r = (int) (highlightColor.x * 255);
						g = (int) (highlightColor.y * 255);
						b = (int) (highlightColor.z * 255);
					}
					ExtendedGenericDialog egd = new ExtendedGenericDialog("Highlight colour", null);
					egd.addSlider("Red", 0, 255, r);
					egd.addSlider("Green", 0, 255, g);
					egd.addSlider("Blue", 0, 255, b);
					egd.showDialog(true, gd);
					if (egd.wasCanceled())
						return false;
					r = (int) egd.getNextNumber();
					g = (int) egd.getNextNumber();
					b = (int) egd.getNextNumber();
					Color c = new Color(r, g, b);
					//highlightColor = new Color3f();
					String cvalue = c.getRed() + "," + c.getGreen() + "," + c.getBlue();
					//settings.setHighlightColour(cvalue);
					tf[0].setText(cvalue);
					return true;
				}

				public boolean collectOptions()
				{
					return false;
				}
			});
			tf[0] = gd.getLastTextField();
			ResultsManager.addTableResultsOptions(gd, s, ResultsManager.FLAG_NO_SECTION_HEADER);
			gd.addMessage("Allow the 'Find Eye Point' command to save to settings");
			gd.addCheckbox("Save_eye_point", settings.getSaveEyePoint());
			// Same as CropResults
			gd.addChoice("Crop_name_option", CropResults.NAME_OPTIONS, settings.getNameOption(),
					new OptionListener<Integer>()
					{
						public boolean collectOptions(Integer value)
						{
							settings.setNameOption(value);
							ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
							if (settings.getNameOption() == CropResults.NAME_OPTION_NAME)
							{
								return false;
							}
							else if (settings.getNameOption() == CropResults.NAME_OPTION_SUFFIX)
							{
								String name = (TextUtils.isNullOrEmpty(settings.getNameSuffix())) ? " Cropped"
										: settings.getNameSuffix();
								egd.addStringField("Name_suffix", name, Maths.clip(20, 60, name.length()));
							}
							else if (settings.getNameOption() == CropResults.NAME_OPTION_SEQUENCE)
							{
								String name = settings.getNameSuffix();
								egd.addStringField("Name_suffix", name, Maths.clip(20, 60, name.length()));
								int c = settings.getNameCounter();
								if (c < 1)
									c = 1;
								egd.addNumericField("Name_counter", c, 0);
							}
							else
							{
								throw new IllegalStateException("Unknown name option: " + settings.getNameOption());
							}
							egd.showDialog(true, gd);
							if (egd.wasCanceled())
								return false;
							if (settings.getNameOption() == CropResults.NAME_OPTION_SUFFIX)
							{
								settings.setNameSuffix(egd.getNextString());
							}
							else if (settings.getNameOption() == CropResults.NAME_OPTION_SEQUENCE)
							{
								settings.setNameSuffix(egd.getNextString());
								settings.setNameCounter(Math.max(1, (int) egd.getNextNumber()));
							}

							return true;
						}

						public boolean collectOptions()
						{
							return false;
						}
					});
			gd.showDialog();
			if (gd.wasCanceled())
				return;
			settings.setHighlightColour(gd.getNextString());
			resultsTableSettings = s.getResultsTableSettingsBuilder();
			resultsTableSettings.setShowTable(gd.getNextBoolean());
			settings.setSaveEyePoint(gd.getNextBoolean());
			settings.setNameOption(gd.getNextChoiceIndex());

			createHighlightColour(settings.getHighlightColour());

			// Save updated settings
			settings.setResultsTableSettings(resultsTableSettings);
			SettingsManager.writeSettings(settings);
			return;
		}
		if (src == toggleDynamicTransparency)
		{
			View view = univ.getViewer().getView();
			boolean activate = view.getTransparencySortingPolicy() == View.TRANSPARENCY_SORT_NONE;
			long total = (activate) ? getTotalTransparentObjects(univ) : 0;
			activateDynamicTransparency(univ, total, activate);
			return;
		}

		// Actions to perform on content
		ContentAction action = null;
		if (src == changeColour)
		{
			action = new ChangeColourContentAction();
		}
		else if (src == resetAll)
		{
			univ.resetView();
			univ.select(null);
			action = new ResetViewContentAction(false);
		}
		else if (src == resetSelectedView)
		{
			action = new ResetViewContentAction(true);
		}
		else if (src == findEyePoint)
		{
			action = new FindEyePointContentAction();
		}
		else if (src == sortBackToFront)
		{
			action = new SortContentAction(false);
		}
		else if (src == sortFrontToBack)
		{
			action = new SortContentAction(true);
		}
		else if (src == colourSurface)
		{
			action = new ColourSurfaceContentAction();
		}
		else if (src == toggleTransparent)
		{
			action = new ToggleTransparentAction();
		}
		else if (src == toggleShaded)
		{
			action = new ToggleShadedAction();
		}
		else if (src == changePointSize)
		{
			action = new ChangePointSizeContentAction();
		}
		else if (src == cropResults)
		{
			action = new CropResultsAction();
		}
		if (action == null)
			return;

		if (univ.getSelected() != null)
		{
			action.run(univ.getSelected());
		}
		else
		{
			for (Iterator<Content> it = univ.contents(); it.hasNext();)
			{
				if (action.run(it.next()) < 0)
					break;
			}
		}

		action.finish();
	}

	private ItemGeometryGroup createItemGroup(final ImageJ3DResultsViewerSettings.Builder settings,
			final Point3f[] sphereSize, final TurboList<Point3f> points, float[] alpha, float transparency,
			Color3f[] colors)
	{
		// Create the single localisation shape
		Shape3D shape = createShape(settings);

		GeometryArray ga = (GeometryArray) shape.getGeometry();

		long size = (long) points.size() * ga.getValidVertexCount();
		if (size > 10000000L)
		{
			ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
			egd.addMessage("The results will generate a large dataset of " + size +
					" vertices.\nThis may take a long time to render and may run out of memory.");
			egd.setOKLabel("Continue");
			egd.showDialog();
			if (egd.wasCanceled())
				return null;
		}

		Appearance appearance = shape.getAppearance();
		TransparencyAttributes ta = new TransparencyAttributes();
		ta.setTransparency(transparency);
		ta.setTransparencyMode((transparency == 0) ? TransparencyAttributes.NONE : TransparencyAttributes.FASTEST);
		appearance.setTransparencyAttributes(ta);
		return new ItemGeometryGroup(points.toArray(new Point3f[points.size()]), ga, appearance, sphereSize, colors,
				alpha);
	}

	private static Shape3D createShape(Builder settings)
	{
		TurboList<Point3f> points = new TurboList<Point3f>(1);
		points.addf(new Point3f());

		// We try and match the geometry and appearance of the standard mesh.
		// Do this by creating a mesh with a single point and get the Geometry and Appearance.
		GeometryArray ga;
		CustomMesh mesh;

		float transparency = getTransparency(settings);

		// Support drawing as points ...
		if (settings.getRendering() == 0)
		{
			mesh = new TransparentItemPointMesh(points, null, transparency);
			((ItemPointMesh) mesh).setPointSize((float) settings.getPixelSize());

			updateAppearance(mesh, settings);

			// Simple point array
			//ga = new PointArray(1, GeometryArray.COORDINATES);
			//ga.setCoordinates(0, new float[3]);
			//ga.setValidVertexCount(1);
			// Assume the TransparentItemPointMesh sets COLOR_4
			ga = (GeometryArray) mesh.getGeometry();
		}
		else
		{
			Rendering r = Rendering.forNumber(settings.getRendering());

			final List<Point3f> point = createLocalisationObject(r);
			Point3f[] vertices = point.toArray(new Point3f[1]);

			// Correct the direction
			ItemTriangleMesh.checkFacets(vertices);

			double creaseAngle = (r.isHighResolution()) ? 44 : 0;
			mesh = new ItemTriangleMesh(vertices, points.toArray(new Point3f[1]), null, null, transparency, creaseAngle,
					null);

			updateAppearance(mesh, settings);

			int nVertices = vertices.length;
			ga = new TriangleArray(nVertices, GeometryArray.COORDINATES | GeometryArray.NORMALS);

			// Copy the coords and normals. We don't require the vertex colours.
			float[] coords = new float[nVertices * 3];
			float[] normals = new float[nVertices * 3];
			GeometryArray gaToCopy = (GeometryArray) mesh.getGeometry();
			gaToCopy.getCoordinates(0, coords);
			gaToCopy.getNormals(0, normals);
			ga.setCoordinates(0, coords);
			ga.setNormals(0, normals);
			ga.setValidVertexCount(nVertices);
		}

		return new Shape3D(ga, mesh.getAppearance());
	}

	private static CustomMesh createMesh(final ImageJ3DResultsViewerSettingsOrBuilder settings,
			TurboList<Point3f> points, final Point3f[] sphereSize, float transparency, float[] alpha)
	{
		// Support drawing as points ...
		if (settings.getRendering() == 0)
		{
			CustomPointMesh mesh;
			if (alpha != null)
			{
				TransparentItemPointMesh mesh2 = new TransparentItemPointMesh(points, null, transparency);
				mesh = mesh2;
				mesh2.setItemAlpha(alpha);
			}
			else
			{
				mesh = new ItemPointMesh(points, null, transparency);
			}
			mesh.setPointSize(sphereSize[0].x);
			return mesh;
		}

		Rendering r = Rendering.forNumber(settings.getRendering());

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

		final List<Point3f> point = createLocalisationObject(r);

		final int singlePointSize = point.size();
		long size = (long) points.size() * singlePointSize;
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

		IJ.showStatus("Creating 3D mesh ...");
		double creaseAngle = (r.isHighResolution()) ? 44 : 0;
		IJTrackProgress progress = null; // Used for debugging construction time
		if (alpha != null)
		{
			TransparentItemTriangleMesh mesh = new TransparentItemTriangleMesh(
					point.toArray(new Point3f[singlePointSize]), points.toArray(new Point3f[points.size()]), sphereSize,
					null, transparency, creaseAngle, progress);
			mesh.setItemAlpha(alpha);
			return mesh;
		}

		return new ItemTriangleMesh(point.toArray(new Point3f[singlePointSize]),
				points.toArray(new Point3f[points.size()]), sphereSize, null, transparency, creaseAngle, progress);
	}

	/**
	 * Creates the object used to draw a single localisation.
	 * 
	 * @param rendering
	 *
	 * @return the list of triangle vertices for the object
	 */
	private static List<Point3f> createLocalisationObject(Rendering rendering)
	{
		int subdivisions = 0;
		switch (rendering)
		{
			case OCTAHEDRON:
				return createOctahedron();
			case SQUARE:
				return createSquare();
			case TETRAHEDRON:
				return createTetrahedron();
			case TRIANGLE:
				return createTriangle();
			case OCTAGON:
				return createDisc(0, 0, 0, 0, 0, 1, 1, 8);
			case LOW_RES_CIRCLE:
				return createDisc(0, 0, 0, 0, 0, 1, 1, 12);
			case HIGH_RES_CIRCLE:
				return createDisc(0, 0, 0, 0, 0, 1, 1, 20);

			// All handle the same way
			case HIGH_RES_SPHERE:
				subdivisions++;
			case LOW_RES_SPHERE:
				subdivisions++;
			case ICOSAHEDRON:
				break;

			case POINT:
			default:
				throw new IllegalStateException("Unknown rendering " + rendering);
		}

		// All spheres based on icosahedron for speed
		return customnode.MeshMaker.createIcosahedron(subdivisions, 1f);
	}

	/**
	 * Creates the disc. This is copied from MeshMaker but the duplication of the vertices for both sides on the
	 * disc is
	 * removed.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param nx
	 *            the nx
	 * @param ny
	 *            the ny
	 * @param nz
	 *            the nz
	 * @param radius
	 *            the radius
	 * @param edgePoints
	 *            the edge points
	 * @return the list
	 */
	static public List<Point3f> createDisc(final double x, final double y, final double z, final double nx,
			final double ny, final double nz, final double radius, final int edgePoints)
	{
		double ax, ay, az;

		if (Math.abs(nx) >= Math.abs(ny))
		{
			final double scale = 1 / Math.sqrt(nx * nx + nz * nz);
			ax = -nz * scale;
			ay = 0;
			az = nx * scale;
		}
		else
		{
			final double scale = 1 / Math.sqrt(ny * ny + nz * nz);
			ax = 0;
			ay = nz * scale;
			az = -ny * scale;
		}

		/*
		 * Now to find the other vector in that plane, do the
		 * cross product of (ax,ay,az) with (nx,ny,nz)
		 */

		double bx = (ay * nz - az * ny);
		double by = (az * nx - ax * nz);
		double bz = (ax * ny - ay * nx);
		final double bScale = 1 / Math.sqrt(bx * bx + by * by + bz * bz);
		bx *= bScale;
		by *= bScale;
		bz *= bScale;

		final double[] circleX = new double[edgePoints + 1];
		final double[] circleY = new double[edgePoints + 1];
		final double[] circleZ = new double[edgePoints + 1];

		for (int i = 0; i < edgePoints + 1; ++i)
		{
			final double angle = (i * 2 * Math.PI) / edgePoints;
			final double c = Math.cos(angle);
			final double s = Math.sin(angle);
			circleX[i] = x + radius * c * ax + radius * s * bx;
			circleY[i] = y + radius * c * ay + radius * s * by;
			circleZ[i] = z + radius * c * az + radius * s * bz;
		}
		final ArrayList<Point3f> list = new ArrayList<Point3f>();
		final Point3f centre = new Point3f((float) x, (float) y, (float) z);
		for (int i = 0; i < edgePoints; ++i)
		{
			final Point3f t2 = new Point3f((float) circleX[i], (float) circleY[i], (float) circleZ[i]);
			final Point3f t3 = new Point3f((float) circleX[i + 1], (float) circleY[i + 1], (float) circleZ[i + 1]);
			list.add(centre);
			list.add(t2);
			list.add(t3);

			// We do not duplicate the triangle for both sides as we render the object as 2D
			// with setBackFaceNormalFlip(true)
			//list.add(centre);
			//list.add(t3);
			//list.add(t2);
		}
		return list;
	}

	// Note: The triangles are rendered using a right-hand coordinate system. 
	// However for 2D shapes the handedness does matter as we set back-face cull off.
	// For polygons we use ItemTriangleMesh which checks the handedness is 
	// facing away from the centre so back-face cull can be on.

	private static float sqrt(double d)
	{
		return (float) Math.sqrt(d);
	}

	static final private float[][] triVertices = { { sqrt(8d / 9), 0, 0 }, { -sqrt(2d / 9), sqrt(2d / 3), 0 },
			{ -sqrt(2d / 9), -sqrt(2d / 3), 0 } };
	static final private int[][] triFaces = { { 0, 1, 2 } };

	static final private float[][] squareVertices = { { 1, 1, 0 }, { -1, 1, 0 }, { -1, -1, 0 }, { 1, -1, 0 } };
	static final private int[][] squareFaces = { { 0, 1, 3 }, { 3, 1, 2 } };

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
	 * Creates the object used to outline a single localisation.
	 * 
	 * @param rendering
	 *
	 * @return the list of triangle vertices for the object
	 */
	private static List<Point3f> createLocalisationObjectOutline(Rendering rendering)
	{
		switch (rendering)
		{
			case SQUARE:
				return createSquareOutline();
			case TRIANGLE:
				return createTriangleOutline();
			case OCTAGON:
				return createDiscOutline(0, 0, 0, 0, 0, 1, 1, 8);
			case LOW_RES_CIRCLE:
				return createDiscOutline(0, 0, 0, 0, 0, 1, 1, 12);
			case HIGH_RES_CIRCLE:
				return createDiscOutline(0, 0, 0, 0, 0, 1, 1, 20);

			default:
				return createLocalisationObject(rendering);
		}
	}

	/**
	 * Creates the triangle with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createTriangleOutline()
	{
		return createSolidOutline(triVertices, true);
	}

	/**
	 * Creates the square with vertices on a unit sphere.
	 *
	 * @return the list of vertices for the triangles
	 */
	private static List<Point3f> createSquareOutline()
	{
		return createSolidOutline(squareVertices, false);
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
	private static List<Point3f> createSolidOutline(float[][] vertices, boolean normalise)
	{
		List<Point3f> ps = new ArrayList<Point3f>();
		for (int i = 0; i < vertices.length; i++)
		{
			ps.add(new Point3f(vertices[i]));
		}
		// Make continuous 
		ps.add(new Point3f(vertices[0]));
		return ps;
	}

	/**
	 * Creates the disc. This is copied from MeshMaker but the duplication of the vertices for both sides on the
	 * disc is
	 * removed.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param nx
	 *            the nx
	 * @param ny
	 *            the ny
	 * @param nz
	 *            the nz
	 * @param radius
	 *            the radius
	 * @param edgePoints
	 *            the edge points
	 * @return the list
	 */
	static public List<Point3f> createDiscOutline(final double x, final double y, final double z, final double nx,
			final double ny, final double nz, final double radius, final int edgePoints)
	{
		double ax, ay, az;

		if (Math.abs(nx) >= Math.abs(ny))
		{
			final double scale = 1 / Math.sqrt(nx * nx + nz * nz);
			ax = -nz * scale;
			ay = 0;
			az = nx * scale;
		}
		else
		{
			final double scale = 1 / Math.sqrt(ny * ny + nz * nz);
			ax = 0;
			ay = nz * scale;
			az = -ny * scale;
		}

		/*
		 * Now to find the other vector in that plane, do the
		 * cross product of (ax,ay,az) with (nx,ny,nz)
		 */

		double bx = (ay * nz - az * ny);
		double by = (az * nx - ax * nz);
		double bz = (ax * ny - ay * nx);
		final double bScale = 1 / Math.sqrt(bx * bx + by * by + bz * bz);
		bx *= bScale;
		by *= bScale;
		bz *= bScale;

		final ArrayList<Point3f> list = new ArrayList<Point3f>();
		for (int i = 0; i < edgePoints + 1; ++i)
		{
			final double angle = (i * 2 * Math.PI) / edgePoints;
			final double c = Math.cos(angle);
			final double s = Math.sin(angle);
			float px = (float) (x + radius * c * ax + radius * s * bx);
			float py = (float) (y + radius * c * ay + radius * s * by);
			float pz = (float) (z + radius * c * az + radius * s * bz);
			list.add(new Point3f(px, py, pz));
		}
		return list;
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
