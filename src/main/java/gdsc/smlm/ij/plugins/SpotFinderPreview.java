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

import java.awt.AWTEvent;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Label;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Vector;

import gdsc.core.ij.Utils;
import gdsc.core.match.AUCCalculator;
import gdsc.core.match.BasePoint;
import gdsc.core.match.ClassificationResult;
import gdsc.core.match.Coordinate;
import gdsc.core.match.FractionalAssignment;
import gdsc.core.match.ImmutableFractionalAssignment;
import gdsc.core.match.RankedScoreCalculator;
import gdsc.core.utils.Maths;
import gdsc.core.utils.RampedScore;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.StoredData;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import gdsc.smlm.data.config.FitProtos.DataFilterType;
import gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.filters.SpotFilterHelper;
import gdsc.smlm.ij.plugins.PeakFit.FitConfigurationProvider;
import gdsc.smlm.ij.plugins.PeakFit.RelativeParameterProvider;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.model.camera.FakePerPixelCameraModel;
import gdsc.smlm.results.MemoryPeakResults;
import gnu.trove.map.hash.TIntObjectHashMap;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionCollectedEvent;
import ij.gui.ExtendedGenericDialog.OptionCollectedListener;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.WindowOrganiser;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;

/**
 * Runs the candidate maxima identification on the image and provides a preview using an overlay
 */
public class SpotFinderPreview implements ExtendedPlugInFilter, DialogListener, ImageListener, ItemListener,
		FitConfigurationProvider, OptionCollectedListener
{
	private final static String TITLE = "Spot Finder Preview";

	private static DataFilterMethod defaultDataFilterMethod;
	private static double defaultSmooth;
	static
	{
		FitEngineConfiguration c = new FitEngineConfiguration();
		defaultDataFilterMethod = c.getDataFilterMethod(0);
		defaultSmooth = c.getDataFilterParameterValue(0);
	}

	private int flags = DOES_16 | DOES_8G | DOES_32;
	private FitEngineConfiguration config = null;
	private FitConfiguration fitConfig = null;
	private Overlay o = null;
	private ImagePlus imp = null;
	private boolean preview = false;
	private Label label = null;
	private TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates = null;
	private static double distance = 1.5;
	private static double lowerDistance = 50;
	private static boolean multipleMatches = false;
	private static boolean showTP = true;
	private static boolean showFP = true;
	private static int topN = 100;
	private static int select = 1;
	private static int neighbourRadius = 4;

	private int currentSlice = 0;
	private MaximaSpotFilter filter = null;

	// All the fields that will be updated when reloading the configuration file
	private Choice textCameraModelName;
	private Choice textPSF;
	private Choice textDataFilterType;
	private Choice textDataFilterMethod;
	private TextField textSmooth;
	private Choice textDataFilterMethod2;
	private TextField textSmooth2;
	private TextField textSearch;
	private TextField textBorder;

	// For adjusting the selction sliders
	private Scrollbar topNScrollBar, selectScrollBar;

	private boolean refreshing = false;
	private NonBlockingExtendedGenericDialog gd;

	private SpotFilterHelper spotFilterHelper = new SpotFilterHelper();

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}

		Roi roi = imp.getRoi();
		if (roi != null && roi.getType() != Roi.RECTANGLE)
		{
			IJ.error("Rectangular ROI required");
			return DONE;
		}

		return flags;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#showDialog(ij.ImagePlus, java.lang.String,
	 * ij.plugin.filter.PlugInFilterRunner)
	 */
	@Override
	@SuppressWarnings("unchecked")
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		this.o = imp.getOverlay();
		this.imp = imp;

		config = SettingsManager.readFitEngineConfiguration(0);
		fitConfig = config.getFitConfiguration();

		gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Preview candidate maxima");

		String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);

		String[] models = CameraModelManager.listCameraModels(true);
		gd.addChoice("Camera_model_name", models, fitConfig.getCameraModelName());

		PeakFit.addPSFOptions(gd, this);
		PeakFit.SimpleFitEngineConfigurationProvider provider = new PeakFit.SimpleFitEngineConfigurationProvider(
				config);
		PeakFit.addDataFilterOptions(gd, provider);
		gd.addChoice("Spot_filter_2", SettingsManager.getDataFilterMethodNames(),
				config.getDataFilterMethod(1, defaultDataFilterMethod).ordinal());
		//gd.addSlider("Smoothing_2", 2.5, 4.5, config.getDataFilterParameterValue(1, defaultSmooth));
		PeakFit.addRelativeParameterOptions(gd, new RelativeParameterProvider(2.5, 4.5, "Smoothing_2", provider)
		{
			@Override
			void setAbsolute(boolean absolute)
			{
				FitEngineConfiguration c = fitEngineConfigurationProvider.getFitEngineConfiguration();
				DataFilterMethod m = c.getDataFilterMethod(1, defaultDataFilterMethod);
				double smooth = c.getDataFilterParameterValue(1, defaultSmooth);
				c.setDataFilter(m, smooth, absolute, 1);
			}

			@Override
			boolean isAbsolute()
			{
				return fitEngineConfigurationProvider.getFitEngineConfiguration().getDataFilterParameterAbsolute(1,
						false);
			}

			@Override
			double getValue()
			{
				return fitEngineConfigurationProvider.getFitEngineConfiguration().getDataFilterParameterValue(1,
						defaultSmooth);
			}
		});

		PeakFit.addSearchOptions(gd, provider);
		PeakFit.addBorderOptions(gd, provider);
		//gd.addNumericField("Top_N", topN, 0);
		gd.addSlider("Top_N", 0, 100, topN);
		topNScrollBar = gd.getLastScrollbar();
		gd.addSlider("Select", 0, 100, select);
		selectScrollBar = gd.getLastScrollbar();
		gd.addSlider("Neigbour_radius", 0, 10, neighbourRadius);

		// Find if this image was created with ground truth data
		if (imp.getID() == CreateData.getImageId())
		{
			MemoryPeakResults results = CreateData.getResults();
			if (results != null)
			{
				gd.addSlider("Match_distance", 0, 2.5, distance);
				gd.addSlider("Lower_match_distance (%)", 0, 100, lowerDistance);
				gd.addCheckbox("Multiple_matches", multipleMatches);
				gd.addCheckbox("Show_TP", showTP);
				gd.addCheckbox("Show_FP", showFP);
				gd.addMessage("");
				label = (Label) gd.getMessage();
				boolean integerCoords = false;
				actualCoordinates = ResultsMatchCalculator.getCoordinates(results, integerCoords);
			}
		}

		if (Utils.isShowGenericDialog())
		{
			// Listen for changes in the dialog options
			gd.addOptionCollectedListener(this);
			// Listen for changes to an image
			ImagePlus.addImageListener(this);

			// Support template settings
			Vector<TextField> numerics = gd.getNumericFields();
			Vector<Choice> choices = gd.getChoices();

			int n = 0;
			int ch = 0;

			Choice textTemplate = choices.get(ch++);
			textTemplate.removeItemListener(gd);
			textTemplate.removeKeyListener(gd);
			textTemplate.addItemListener(this);

			textCameraModelName = choices.get(ch++);
			textPSF = choices.get(ch++);
			textDataFilterType = choices.get(ch++);
			textDataFilterMethod = choices.get(ch++);
			textSmooth = numerics.get(n++);
			textDataFilterMethod2 = choices.get(ch++);
			textSmooth2 = numerics.get(n++);
			textSearch = numerics.get(n++);
			textBorder = numerics.get(n++);
		}

		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.setOKLabel("Save");
		gd.setCancelLabel("Close");
		gd.showDialog();

		if (!(IJ.isMacro() || java.awt.GraphicsEnvironment.isHeadless()))
			ImagePlus.removeImageListener(this);

		if (!gd.wasCanceled())
		{
			if (!SettingsManager.writeSettings(config, SettingsManager.FLAG_SILENT))
				IJ.error(TITLE, "Failed to save settings");
		}

		// Reset
		imp.setOverlay(o);

		return DONE;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		if (refreshing)
			return false;

		gd.getNextChoice(); // Ignore template

		// Set a camera model
		fitConfig.setCameraModelName(gd.getNextChoice());
		CameraModel model = CameraModelManager.load(fitConfig.getCameraModelName());
		if (model == null)
			model = new FakePerPixelCameraModel(0, 1, 1);
		fitConfig.setCameraModel(model);

		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);

		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 1);
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		topN = (int) gd.getNextNumber();
		select = (int) gd.getNextNumber();
		neighbourRadius = (int) gd.getNextNumber();

		if (label != null)
		{
			distance = gd.getNextNumber();
			lowerDistance = gd.getNextNumber();
			multipleMatches = gd.getNextBoolean();
			showTP = gd.getNextBoolean();
			showFP = gd.getNextBoolean();
		}
		preview = gd.getNextBoolean();

		((ExtendedGenericDialog) gd).collectOptions();

		boolean result = !gd.invalidNumber();
		if (!preview)
		{
			setLabel("");
			this.imp.setOverlay(o);
		}
		// For astigmatism PSF.
		// TODO - See if this is slowing the preview down. If so only do if the PSF type changes.
		if (!PeakFit.configurePSFModel(config, PeakFit.FLAG_NO_SAVE))
			return false;
		return result;
	}

	private void setLabel(String message)
	{
		if (label == null)
			return;
		label.setText(message);
	}

	private Calibration lastCalibration = null;
	private FitEngineSettings lastFitEngineSettings = null;
	private PSF lastPSF = null;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip)
	{
		if (refreshing)
			return;

		Rectangle bounds = ip.getRoi();

		// Only do this if the settings changed
		Calibration calibration = fitConfig.getCalibration();
		FitEngineSettings fitEngineSettings = config.getFitEngineSettings();
		PSF psf = fitConfig.getPSF();

		boolean newCameraModel = false;
		if (!calibration.equals(lastCalibration))
		{
			newCameraModel = true;
			// Set a camera model.
			// We have to set the camera type too to avoid configuration errors.
			CameraModel cameraModel = CameraModelManager.load(fitConfig.getCameraModelName());
			if (cameraModel == null)
			{
				cameraModel = new FakePerPixelCameraModel(0, 1, 1);
				fitConfig.setCameraType(CameraType.EMCCD);
			}
			else
			{
				fitConfig.setCameraType(CameraType.SCMOS);
			}
			fitConfig.setCameraModel(cameraModel);
		}

		if (newCameraModel || !fitEngineSettings.equals(lastFitEngineSettings) || !psf.equals(lastPSF))
		{
			// Configure a jury filter
			if (config.getDataFilterType() == DataFilterType.JURY)
			{
				if (!PeakFit.configureDataFilter(config, PeakFit.FLAG_NO_SAVE))
					return;
			}

			try
			{
				filter = config.createSpotFilter();
			}
			catch (Exception e)
			{
				filter = null;
				this.imp.setOverlay(o);
				throw new RuntimeException(e); // Required for ImageJ to disable the preview
				//Utils.log("ERROR: " + e.getMessage());			
				//return;
			}
			Utils.log(filter.getDescription());
		}

		lastCalibration = calibration;
		lastFitEngineSettings = fitEngineSettings;
		lastPSF = psf;

		// Note: 
		// Do not support origin selection when there is a width/height mismatch 
		// as this will break the live preview with multiple pop-ups. Origin 
		// could be added to the input dialog.
		//cameraModel = PeakFit.cropCameraModel(cameraModel, ip.getWidth(), ip.getHeight(), 0, 0, true);		

		// Instead just warn if the roi cannot be extracted from the selected model 
		// or there is a mismatch
		Rectangle modelBounds = fitConfig.getCameraModel().getBounds();
		if (modelBounds != null)
		{
			if (!modelBounds.contains(bounds))
			{
			//@formatter:off
			Utils.log("WARNING: Camera model bounds [x=%d,y=%d,width=%d,height=%d] does not contain image target bounds [x=%d,y=%d,width=%d,height=%d]",
					modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, 
					bounds.x, bounds.y, bounds.width, bounds.height 
					);
			//@formatter:on
			}
			else
			// Warn if the model bounds are mismatched than the image as this may be an incorrect
			// selection for the camera model
			if (modelBounds.x != 0 || modelBounds.y != 0 || modelBounds.width > ip.getWidth() ||
					modelBounds.height > ip.getHeight())
			{
			//@formatter:off
			Utils.log("WARNING: Probably an incorrect camera model!\nModel bounds [x=%d,y=%d,width=%d,height=%d]\ndo not match the image target bounds [width=%d,height=%d].",
					modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, 
					ip.getWidth(),  ip.getHeight()
					);
			//@formatter:on
			}
		}

		run(ip, filter);
	}

	private void run(ImageProcessor ip, MaximaSpotFilter filter)
	{
		if (refreshing)
			return;

		currentSlice = imp.getCurrentSlice();

		Rectangle bounds = ip.getRoi();

		// Crop to the ROI
		FloatProcessor fp = ip.crop().toFloat(0, null);

		float[] data = (float[]) fp.getPixels();

		int width = fp.getWidth();
		int height = fp.getHeight();

		// Set weights
		CameraModel cameraModel = fitConfig.getCameraModel();
		if (!(cameraModel instanceof FakePerPixelCameraModel))
		{
			// This should be done on the normalised data
			float[] w = cameraModel.getNormalisedWeights(bounds);
			filter.setWeights(w, width, height);
			data = data.clone();
			cameraModel.removeBiasAndGain(data);
		}

		Spot[] spots = filter.rank(data, width, height);
		data = filter.getPreprocessedData();

		int size = spots.length;
		topNScrollBar.setMaximum(size);
		selectScrollBar.setMaximum(size);

		fp = new FloatProcessor(width, height, data);
		FloatProcessor out = new FloatProcessor(ip.getWidth(), ip.getHeight());
		out.copyBits(ip, 0, 0, Blitter.COPY);
		out.insert(fp, bounds.x, bounds.y);
		//ip.resetMinAndMax();
		double min = fp.getMin();
		double max = fp.getMax();
		out.setMinAndMax(min, max);

		Overlay o = new Overlay();
		o.add(new ImageRoi(0, 0, out));

		if (label != null)
		{
			// Get results for frame
			Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, imp.getCurrentSlice());

			Coordinate[] predicted = new Coordinate[size];
			for (int i = 0; i < size; i++)
			{
				predicted[i] = new BasePoint(spots[i].x + bounds.x, spots[i].y + bounds.y);
			}

			// Compute assignments
			TurboList<FractionalAssignment> fractionalAssignments = new TurboList<FractionalAssignment>(
					3 * predicted.length);
			double matchDistance = distance * fitConfig.getInitialPeakStdDev();
			final RampedScore score = new RampedScore(matchDistance * lowerDistance / 100, matchDistance);
			final double dmin = matchDistance * matchDistance;
			final int nActual = actual.length;
			final int nPredicted = predicted.length;
			for (int j = 0; j < nPredicted; j++)
			{
				// Centre in the middle of the pixel
				final float x = predicted[j].getX() + 0.5f;
				final float y = predicted[j].getY() + 0.5f;
				// Any spots that match 
				for (int i = 0; i < nActual; i++)
				{
					final double dx = (x - actual[i].getX());
					final double dy = (y - actual[i].getY());
					final double d2 = dx * dx + dy * dy;
					if (d2 <= dmin)
					{
						final double d = Math.sqrt(d2);
						double s = score.score(d);

						if (s == 0)
							continue;

						double distance = 1 - s;
						if (distance == 0)
						{
							// In the case of a match below the distance thresholds
							// the distance will be 0. To distinguish between candidates all below 
							// the thresholds just take the closest.
							// We know d2 is below dmin so we subtract the delta.
							distance -= (dmin - d2);
						}

						// Store the match
						fractionalAssignments.add(new ImmutableFractionalAssignment(i, j, distance, s));
					}
				}
			}

			FractionalAssignment[] assignments = fractionalAssignments
					.toArray(new FractionalAssignment[fractionalAssignments.size()]);

			// Compute matches
			RankedScoreCalculator calc = new RankedScoreCalculator(assignments, nActual - 1, nPredicted - 1);
			boolean save = showTP || showFP;
			double[] calcScore = calc.score(nPredicted, multipleMatches, save);
			ClassificationResult result = RankedScoreCalculator.toClassificationResult(calcScore, nActual);

			// Compute AUC and max jaccard (and plot)
			double[][] curve = RankedScoreCalculator.getPrecisionRecallCurve(assignments, nActual, nPredicted);
			double[] precision = curve[0];
			double[] recall = curve[1];
			double[] jaccard = curve[2];
			double auc = AUCCalculator.auc(precision, recall);

			// Show scores
			String label = String.format("Slice=%d, AUC=%s, R=%s, Max J=%s", imp.getCurrentSlice(), Utils.rounded(auc),
					Utils.rounded(result.getRecall()), Utils.rounded(Maths.maxDefault(0, jaccard)));
			setLabel(label);

			// Plot
			String title = TITLE + " Performance";
			Plot2 plot = new Plot2(title, "Spot Rank", "");
			double[] rank = SimpleArrayUtils.newArray(precision.length, 0, 1.0);
			plot.setLimits(0, nPredicted, 0, 1.05);
			plot.setColor(Color.blue);
			plot.addPoints(rank, precision, Plot.LINE);
			plot.setColor(Color.red);
			plot.addPoints(rank, recall, Plot.LINE);
			plot.setColor(Color.black);
			plot.addPoints(rank, jaccard, Plot.LINE);
			plot.setColor(Color.black);
			plot.addLabel(0, 0, label);

			WindowOrganiser windowOrganiser = new WindowOrganiser();
			Utils.display(title, plot, 0, windowOrganiser);

			title = TITLE + " Precision-Recall";
			plot = new Plot2(title, "Recall", "Precision");
			plot.setLimits(0, 1, 0, 1.05);
			plot.setColor(Color.red);
			plot.addPoints(recall, precision, Plot.LINE);
			plot.drawLine(recall[recall.length - 1], precision[recall.length - 1], recall[recall.length - 1], 0);
			plot.setColor(Color.black);
			plot.addLabel(0, 0, label);
			Utils.display(title, plot, 0, windowOrganiser);

			windowOrganiser.tile();

			// Create Rois for TP and FP
			if (save)
			{
				double[] matchScore = RankedScoreCalculator.getMatchScore(calc.getScoredAssignments(), nPredicted);
				int matches = 0;
				for (int i = 0; i < matchScore.length; i++)
					if (matchScore[i] != 0)
						matches++;
				if (showTP)
				{
					float[] x = new float[matches];
					float[] y = new float[x.length];
					int n = 0;
					for (int i = 0; i < matchScore.length; i++)
					{
						if (matchScore[i] != 0)
						{
							BasePoint p = (BasePoint) predicted[i];
							x[n] = p.getX() + 0.5f;
							y[n] = p.getY() + 0.5f;
							n++;
						}
					}
					addRoi(0, o, x, y, n, Color.green);
				}
				if (showFP)
				{
					float[] x = new float[nPredicted - matches];
					float[] y = new float[x.length];
					int n = 0;
					for (int i = 0; i < matchScore.length; i++)
					{
						if (matchScore[i] == 0)
						{
							BasePoint p = (BasePoint) predicted[i];
							x[n] = p.getX() + 0.5f;
							y[n] = p.getY() + 0.5f;
							n++;
						}
					}
					addRoi(0, o, x, y, n, Color.red);
				}
			}
		}
		else
		{
			//float[] x = new float[size];
			//float[] y = new float[x.length];
			//for (int i = 0; i < size; i++)
			//{
			//	x[i] = spots[i].x + bounds.x + 0.5f;
			//	y[i] = spots[i].y + bounds.y + 0.5f;
			//}
			//PointRoi roi = new PointRoi(x, y);
			//// Add options to configure colour and labels
			//o.add(roi);

			WindowOrganiser wo = new WindowOrganiser();

			// Option to show the number of neighbours within a set pixel box radius
			int[] count = spotFilterHelper.countNeighbours(spots, width, height, neighbourRadius);

			// Show as histogram the totals...
			int id = Utils.showHistogram(TITLE, new StoredData(count), "Neighbours", 1, 0, 0,
					"Radius = " + neighbourRadius);
			if (Utils.isNewWindow())
				wo.add(id);

			// TODO - Draw n=0, n=1 on the image overlay

			final LUT lut = LUTHelper.createLUT(LutColour.FIRE_LIGHT);
			// These are copied by the ROI
			float[] x = new float[1];
			float[] y = new float[1];
			// Plot the intensity
			double[] intensity = new double[size];
			double[] rank = SimpleArrayUtils.newArray(size, 1, 1.0);
			int top = (topN > 0) ? topN : size;
			int size_1 = size - 1;
			for (int i = 0; i < size; i++)
			{
				intensity[i] = spots[i].intensity;
				if (i < top)
				{
					x[0] = spots[i].x + bounds.x + 0.5f;
					y[0] = spots[i].y + bounds.y + 0.5f;
					Color c = LUTHelper.getColour(lut, size_1 - i, size);
					addRoi(0, o, x, y, 1, c, 2, 1);
				}
			}

			String title = TITLE + " Intensity";
			Plot plot = new Plot(title, "Rank", "Intensity");
			plot.setColor(Color.blue);
			plot.addPoints(rank, intensity, Plot.LINE);
			if (topN > 0 && topN < size)
			{
				plot.setColor(Color.magenta);
				plot.drawLine(topN, 0, topN, intensity[topN - 1]);
			}
			if (select > 0 && select < size)
			{
				plot.setColor(Color.yellow);
				double in = intensity[select - 1];
				plot.drawLine(select, 0, select, in);
				x[0] = spots[select].x + bounds.x + 0.5f;
				y[0] = spots[select].y + bounds.y + 0.5f;
				Color c = LUTHelper.getColour(lut, size_1 - select, size);
				addRoi(0, o, x, y, 1, c, 3, 3);
				plot.setColor(Color.black);
				plot.addLabel(0, 0, "Selected spot intensity = " + Utils.rounded(in));
			}

			Utils.display(title, plot, 0, wo);

			wo.tile();
		}

		imp.setOverlay(o);
	}

	public static void addRoi(int frame, Overlay o, float[] x, float[] y, int n, Color colour)
	{
		// Add as a circle
		addRoi(frame, o, x, y, n, colour, 3);
	}

	public static void addRoi(int frame, Overlay o, float[] x, float[] y, int n, Color colour, int pointType)
	{
		addRoi(frame, o, x, y, n, colour, pointType, 0);
	}

	public static void addRoi(int frame, Overlay o, float[] x, float[] y, int n, Color colour, int pointType, int size)
	{
		if (n == 0)
			return;
		PointRoi roi = new PointRoi(x, y, n);
		roi.setPointType(pointType);
		roi.setFillColor(colour);
		roi.setStrokeColor(colour);
		if (frame != 0)
			roi.setPosition(frame);
		if (size != 0)
			roi.setSize(size);
		o.add(roi);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	@Override
	public void setNPasses(int nPasses)
	{
		// Nothing to do		
	}

	@Override
	public void imageOpened(ImagePlus imp)
	{

	}

	@Override
	public void imageClosed(ImagePlus imp)
	{
	}

	@Override
	public void imageUpdated(ImagePlus imp)
	{
		if (this.imp.getID() == imp.getID() && preview)
		{
			if (imp.getCurrentSlice() != currentSlice && filter != null)
			{
				run(imp.getProcessor(), filter);
			}
		}
	}

	@Override
	public void itemStateChanged(ItemEvent e)
	{
		if (e.getSource() instanceof Choice)
		{
			// Update the settings from the template
			Choice choice = (Choice) e.getSource();
			String templateName = choice.getSelectedItem();
			//System.out.println("Update to " + templateName);

			// Get the configuration template
			TemplateSettings template = ConfigurationTemplate.getTemplate(templateName);

			if (template != null)
			{
				refreshing = true;

				IJ.log("Applying template: " + templateName);

				for (String note : template.getNotesList())
					IJ.log(note);

				boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
				if (template.hasPsf())
				{
					refreshSettings(template.getPsf(), custom);
				}
				if (template.hasFitEngineSettings())
				{
					refreshSettings(template.getFitEngineSettings(), custom);
				}

				refreshing = false;
				//dialogItemChanged(gd, null);
			}
		}
	}

	private void refreshSettings(PSF psf, boolean isCustomTemplate)
	{
		if (!isCustomTemplate || psf == null)
			return;

		// Do not use set() as we support merging a partial PSF
		fitConfig.mergePSF(psf);

		textPSF.select(PSFProtosHelper.getName(fitConfig.getPSFType()));
	}

	/**
	 * Refresh settings.
	 * <p>
	 * If this is a custom template then use all the settings. If a default template then leave some existing spot
	 * settings untouched as the user may have updated them (e.g. PSF width).
	 *
	 * @param fitEngineSettings
	 *            the config
	 * @param isCustomTemplate
	 *            True if a custom template.
	 */
	private void refreshSettings(FitEngineSettings fitEngineSettings, boolean isCustomTemplate)
	{
		// Set the configuration
		// This will clear everything and merge the configuration so
		// remove the fit settings (as we do not care about those).

		this.config.setFitEngineSettings(fitEngineSettings.toBuilder().clearFitSettings().build());
		fitConfig = this.config.getFitConfiguration();

		textCameraModelName.select(fitConfig.getCameraModelName());
		textDataFilterType.select(SettingsManager.getDataFilterTypeNames()[config.getDataFilterType().ordinal()]);
		textDataFilterMethod
				.select(SettingsManager.getDataFilterMethodNames()[config.getDataFilterMethod(0).ordinal()]);
		textSmooth.setText("" + config.getDataFilterParameterValue(0));
		if (config.getDataFiltersCount() > 1)
		{
			textDataFilterMethod2.select(SettingsManager.getDataFilterMethodNames()[config
					.getDataFilterMethod(1, defaultDataFilterMethod).ordinal()]);
			textSmooth2.setText("" + config.getDataFilterParameterValue(1, defaultSmooth));
			// XXX - What about the Abolute/Relative flag?
		}
		textSearch.setText("" + config.getSearch());
		textBorder.setText("" + config.getBorder());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.plugins.PeakFit.FitConfigurationProvider#getFitConfiguration()
	 */
	@Override
	public FitConfiguration getFitConfiguration()
	{
		return fitConfig;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.ExtendedGenericDialog.OptionCollectedListener#optionCollected(ij.gui.ExtendedGenericDialog.
	 * OptionCollectedEvent)
	 */
	@Override
	public void optionCollected(OptionCollectedEvent e)
	{
		// Just run on the current processor
		if (preview)
			run(imp.getProcessor());
	}
}
