package gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.io.File;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import gdsc.core.ij.Utils;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationProtos.CameraModelResource;
import gdsc.smlm.data.config.CalibrationProtos.CameraModelSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.camera.PerPixelCameraModel;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;

/**
 * This plugin handle the save and load of per-pixel camera models
 */
public class CameraModelManager implements PlugIn
{
	private static final String TITLE = "Camera Model Manager";
	private static final String INFO_TAG = "Per-pixel camera model data";

	private static CameraModelSettings.Builder settings = null;
	// Cache camera models for speed
	private static LinkedHashMap<String, PerPixelCameraModel> map = new LinkedHashMap<String, PerPixelCameraModel>();

	private static CameraModelSettings.Builder getSettings()
	{
		return getSettings(0);
	}

	private static CameraModelSettings.Builder getSettings(int flags)
	{
		if (settings == null)
			settings = SettingsManager.readCameraModelSettings(flags).toBuilder();
		return settings;
	}

	/**
	 * Save the camera model. The model will be named in the resources using the filename without the extension or
	 * leading path entries.
	 *
	 * @param cameraModel
	 *            the camera model
	 * @param filename
	 *            the filename
	 * @return true, if successful
	 */
	public static boolean save(PerPixelCameraModel cameraModel, String filename)
	{
		if (cameraModel == null || TextUtils.isNullOrEmpty(filename))
			return false;

		// Try to save to file
		//filename = Utils.replaceExtension(filename, ".tif");
		File file = new File(filename);
		String name = Utils.removeExtension(file.getName());

		ImageStack stack = new ImageStack(cameraModel.getWidth(), cameraModel.getHeight());
		stack.addSlice("Bias", cameraModel.getBias());
		stack.addSlice("Gain", cameraModel.getGain());
		stack.addSlice("Normalised Variance", cameraModel.getNormalisedVariance());
		ImagePlus imp = new ImagePlus(name, stack);
		imp.setIgnoreGlobalCalibration(true);
		Calibration cal = imp.getCalibration();
		cal.xOrigin = cameraModel.getXOrigin();
		cal.yOrigin = cameraModel.getYOrigin();
		imp.setProperty("Info", INFO_TAG);
		// Do this to allow the filename to be something other than .tif
		boolean ok = new FileSaver(imp).saveAsTiffStack(filename);

		if (ok)
		{
			CameraModelResource.Builder resource = CameraModelResource.newBuilder();
			resource.setX(cameraModel.getXOrigin());
			resource.setY(cameraModel.getYOrigin());
			resource.setWidth(cameraModel.getWidth());
			resource.setHeight(cameraModel.getHeight());
			resource.setFilename(filename);

			CameraModelSettings.Builder settings = getSettings();
			settings.getCameraModelResourcesMap().put(name, resource.build());
			SettingsManager.writeSettings(settings.build());

			// Cache this
			map.put(name, cameraModel);
		}

		return ok;
	}

	/**
	 * Load the camera model. Returns null if the named model does not exist. Writes to the ImageJ log if a problems
	 * occurred loading the model.
	 *
	 * @param name
	 *            the name
	 * @return the per pixel camera model (or null)
	 */
	public static PerPixelCameraModel load(String name)
	{
		PerPixelCameraModel model = map.get(name);
		if (model == null)
		{
			CameraModelSettings.Builder settings = getSettings();
			// Try and get the named resource
			CameraModelResource resource = settings.getCameraModelResourcesMap().get(name);
			if (resource == null)
				return null;
			// Try and load the resource
			ImagePlus imp = IJ.openImage(resource.getFilename());
			if (imp == null)
			{
				IJ.log("Failed to load camera data for model: " + name);
				return null;
			}
			// Check stack size
			ImageStack stack = imp.getImageStack();
			if (stack.getSize() != 3)
			{
				IJ.log("Camera model requires 3 image stack for model: " + name);
				return null;
			}
			// Get the origin
			imp.setIgnoreGlobalCalibration(true);
			Calibration cal = imp.getCalibration();
			Rectangle bounds = new Rectangle((int) cal.xOrigin, (int) cal.yOrigin, stack.getWidth(), stack.getHeight());
			try
			{
				float[] bias = (float[]) stack.getPixels(1);
				float[] gain = (float[]) stack.getPixels(2);
				float[] normalisedVariance = (float[]) stack.getPixels(3);
				model = PerPixelCameraModel.create(bounds, bias, gain, normalisedVariance);
			}
			catch (Exception e)
			{
				IJ.log("Failed to load camera model: " + name + ". " + e.getMessage());
			}

			// Cache this
			map.put(name, model);
		}
		return model;
	}

	/**
	 * List the camera models.
	 *
	 * @param includeNone
	 *            Set to true to include an empty string
	 * @return the list
	 */
	public static String[] listCameraModels(boolean includeNone)
	{
		CameraModelSettings.Builder settings = getSettings();
		List<String> list = createList(includeNone);
		list.addAll(settings.getCameraModelResourcesMap().keySet());
		return list.toArray(new String[list.size()]);
	}

	private static List<String> createList(boolean includeNone)
	{
		List<String> list = new TurboList<String>();
		if (includeNone)
			list.add("");
		return list;
	}

	/**
	 * List the camera models with the correct dimensions.
	 *
	 * @param includeNone
	 *            Set to true to include an empty string
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @return the list
	 */
	public static String[] listCameraModels(boolean includeNone, int width, int height)
	{
		CameraModelSettings.Builder settings = getSettings();
		List<String> list = createList(includeNone);
		for (Map.Entry<String, CameraModelResource> entry : settings.getCameraModelResourcesMap().entrySet())
		{
			CameraModelResource resource = entry.getValue();
			if (resource.getWidth() == width && resource.getHeight() == height)
				list.add(entry.getKey());
		}
		return list.toArray(new String[list.size()]);
	}

	private static String[] OPTIONS = { "Print model details", "View a camera model" };
	private static int option = 0;
	private static String selected = "";

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		CameraModelSettings.Builder settings = getSettings(SettingsManager.FLAG_SILENT);
		if (settings.getCameraModelResourcesCount() == 0)
		{
			IJ.error(TITLE, "No camera models found");
			return;
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addChoice("Option", OPTIONS, OPTIONS[option]);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		option = gd.getNextChoiceIndex();

		switch (option)
		{
			case 1:
				viewCameraModel();
				break;
			default:
				printCameraModels();
		}
	}

	private void viewCameraModel()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		String[] MODELS = listCameraModels(false);
		gd.addChoice("Model", MODELS, selected);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		String name = selected = gd.getNextChoice();
		
		// Try and get the named resource
		CameraModelResource resource = settings.getCameraModelResourcesMap().get(name);
		if (resource == null)
		{
			IJ.log("Failed to find camera data for model: " + name);
			return;
		}
		// Try and load the resource
		ImagePlus imp = IJ.openImage(resource.getFilename());
		if (imp == null)
		{
			IJ.log("Failed to load camera data for model: " + name);
			return;
		}
		Utils.log("Camera model: %s\n%s", name, resource);
		imp.show();
	}

	private void printCameraModels()
	{
		IJ.log(settings.toString());		
	}
}
