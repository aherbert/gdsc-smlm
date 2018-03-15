package gdsc.smlm.ij.plugins;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.View;
import org.scijava.vecmath.Point3f;

import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.ij.ij3d.CustomContent;
import gdsc.smlm.ij.ij3d.ItemGeometryGroup;
import gdsc.smlm.ij.ij3d.ItemGeometryNode;
import gdsc.smlm.ij.ij3d.OrderedItemGeometryGroup;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (version == null)
		{
			IJ.error(TITLE, "Java 3D is not available");
			return;
		}

		TurboList<Point3f> pointList;
		float scale;

		if (MemoryPeakResults.isMemoryEmpty())
		{
			pointList = new TurboList<Point3f>();
			int range;
			range = 1; // 9 points
			//range = 12; // 25^3 = 15625 points
			//range = 17; // 35^3 = 42875 points
			//range = 22; // 45^3 = 91125 points
			//range = 49; // 99^3 = 970299 points
			for (int x = -range; x <= range; x++)
			{
				for (int y = -range; y <= range; y++)
				{
					for (int z = -range; z <= range; z++)
					{
						pointList.add(new Point3f(x, y, z));
					}
				}
			}
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
			MemoryPeakResults results = ResultsManager.loadInputResults(name, false, null, null);
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
		ImageWindow3D w = univ.getWindow();
		GUI.center(w);

		View view = univ.getViewer().getView();
		view.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);
		// I am not sure if this is required if objects are sorted.
		//view.setDepthBufferFreezeTransparent(false);

		IJ.showStatus("Creating points ...");
		Point3f[] points = pointList.toArray(new Point3f[pointList.size()]);
		Point3f[] sizes = new Point3f[] { new Point3f(scale, scale, scale) };
		ItemGeometryGroup pointGroup;
		Appearance appearance = new Appearance();
		TransparencyAttributes ta = new TransparencyAttributes();
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
		final CustomContent c = new CustomContent("Test");
		final ContentInstant content = c.getCurrent();
		content.display(new ItemGeometryNode(pointGroup));
		univ.addContent(c);

		IJ.showStatus("Done");
	}
}
