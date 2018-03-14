package gdsc.smlm.ij.plugins;

import org.scijava.vecmath.Point3f;

import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import gdsc.smlm.ij.ij3d.CustomContent;
import gdsc.smlm.ij.ij3d.PointGroup;
import gdsc.smlm.ij.ij3d.PointGroupNode;
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
		float radius;

		if (MemoryPeakResults.isMemoryEmpty())
		{
			pointList = new TurboList<Point3f>();
			int range;
			//range = 1; // 9 points
			range = 12; // 25^3 = 15625 points
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
			radius = 0.25f;
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
			radius = 10f;
		}

		final Image3DUniverse univ = new Image3DUniverse();
		univ.showAttribute(DefaultUniverse.ATTRIBUTE_SCALEBAR, false);
		univ.show();
		ImageWindow3D w = univ.getWindow();
		GUI.center(w);

		//View view = univ.getViewer().getView();
		//view.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);
		//view.setDepthBufferFreezeTransparent(false);

		Point3f[] points = pointList.toArray(new Point3f[pointList.size()]);
		PointGroup pointGroup = new PointGroup(points, radius, null);
		//pointGroup.setRadius(0.25f);

		//		// This supports transparency
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

		// This does since CustomContentInstant uses a Group to show the nodes 
		final CustomContent c = new CustomContent("Test");
		final ContentInstant content = c.getCurrent();
		content.display(new PointGroupNode(pointGroup));
		univ.addContent(c);
	}
}
