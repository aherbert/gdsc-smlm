package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Label;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import gdsc.core.ij.Utils;
import gdsc.core.match.BasePoint;
import gdsc.core.match.Coordinate;
import gdsc.core.match.MatchCalculator;
import gdsc.core.match.MatchResult;
import gdsc.core.match.PointPair;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.engine.DataFilter;
import gdsc.smlm.engine.DataFilterType;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Runs the candidate maxima identification of the image and provides a preview using an overlay
 */
public class SpotFinderPreview implements ExtendedPlugInFilter, DialogListener, ImageListener
{
	private final static String TITLE = "Spot Finder Preview";

	private int flags = DOES_16 | DOES_8G | DOES_32;
	private FitEngineConfiguration config = null;
	private FitConfiguration fitConfig = null;
	private Overlay o = null;
	private ImagePlus imp = null;
	private boolean preview = false;
	private Label label = null;
	private HashMap<Integer, ArrayList<Coordinate>> actualCoordinates = null;
	private static double distance = 1.5;
	private static boolean showTP = true;
	private static boolean showFP = true;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
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
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		this.o = imp.getOverlay();
		this.imp = imp;
		String filename = SettingsManager.getSettingsFilename();
		GlobalSettings settings = SettingsManager.loadSettings(filename);
		config = settings.getFitEngineConfiguration();
		fitConfig = config.getFitConfiguration();

		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Preview candidate maxima");

		String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);

		gd.addStringField("Config_file", filename, 40);

		gd.addNumericField("Initial_StdDev0", fitConfig.getInitialPeakStdDev0(), 3);
		String[] filterTypes = SettingsManager.getNames((Object[]) DataFilterType.values());
		gd.addChoice("Spot_filter_type", filterTypes, filterTypes[config.getDataFilterType().ordinal()]);
		String[] filterNames = SettingsManager.getNames((Object[]) DataFilter.values());
		gd.addChoice("Spot_filter", filterNames, filterNames[config.getDataFilter(0).ordinal()]);
		gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
		gd.addSlider("Search_width", 0.5, 2.5, config.getSearch());
		gd.addSlider("Border", 0.5, 2.5, config.getBorder());

		// Find if this image was created with ground truth data
		if (imp.getID() == CreateData.getImageId())
		{
			MemoryPeakResults results = CreateData.getResults();
			if (results != null)
			{
				gd.addSlider("Match_distance", 0, 2.5, distance);
				gd.addCheckbox("Show TP", showTP);
				gd.addCheckbox("Show FP", showFP);
				gd.addMessage("");
				label = (Label) gd.getMessage();
				// Integer coords
				actualCoordinates = ResultsMatchCalculator.getCoordinates(results.getResults(), true);
			}
		}

		if (!(IJ.isMacro() || java.awt.GraphicsEnvironment.isHeadless()))
		{
			// Listen for changes to an image
			ImagePlus.addImageListener(this);
		}

		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.hideCancelButton();
		gd.setOKLabel("Close");
		gd.showDialog();

		if (!(IJ.isMacro() || java.awt.GraphicsEnvironment.isHeadless()))
			ImagePlus.removeImageListener(this);

		if (!gd.wasCanceled())
		{
			filename = gd.getNextString();
			if (!SettingsManager.saveSettings(settings, filename))
			{
				IJ.error(TITLE, "Failed to save settings to file " + filename);
			}
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
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		gd.getNextChoice();
		gd.getNextString();

		fitConfig.setInitialPeakStdDev0(gd.getNextNumber());
		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		if (label != null)
		{
			distance = gd.getNextNumber();
			showTP = gd.getNextBoolean();
			showFP = gd.getNextBoolean();
		}
		preview = gd.getNextBoolean();
		if (!preview)
		{
			setLabel("");
			this.imp.setOverlay(o);
		}
		return !gd.invalidNumber();
	}

	private void setLabel(String message)
	{
		if (label == null)
			return;
		label.setText(message);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		Rectangle bounds = ip.getRoi();

		MaximaSpotFilter filter = config.createSpotFilter(true);

		// Crop to the ROI
		FloatProcessor fp = ip.crop().toFloat(0, null);

		float[] data = (float[]) fp.getPixels();

		int width = fp.getWidth();
		int height = fp.getHeight();
		Spot[] spots = filter.rank(data, width, height);
		data = filter.getPreprocessedData();

		fp = new FloatProcessor(width, height, data);
		ip = ip.duplicate();
		ip.insert(fp, bounds.x, bounds.y);
		//ip.resetMinAndMax();
		ip.setMinAndMax(fp.getMin(), fp.getMax());

		Overlay o = new Overlay();
		o.add(new ImageRoi(0, 0, ip));

		if (label != null)
		{
			// Get results for frame
			Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, imp.getCurrentSlice());

			Coordinate[] predicted = new Coordinate[spots.length];
			for (int i = 0; i < spots.length; i++)
			{
				predicted[i] = new BasePoint(spots[i].x + bounds.x, spots[i].y + bounds.y);
			}

			// Q. Should this use partial scoring with multi-matches allowed.
			// If so then this needs to be refactored out of the BenchmarkSpotFilter class.

			// TODO - compute AUC and max jaccard and plot			

			// Compute matches
			List<PointPair> matches = new ArrayList<PointPair>(Math.min(actual.length, predicted.length));
			List<Coordinate> FP = new ArrayList<Coordinate>(predicted.length);
			MatchResult result = MatchCalculator.analyseResults2D(actual, predicted,
					distance * fitConfig.getInitialPeakStdDev0(), null, FP, null, matches);

			// Show scores
			setLabel(String.format("P=%s, R=%s, J=%s", Utils.rounded(result.getPrecision()),
					Utils.rounded(result.getRecall()), Utils.rounded(result.getJaccard())));

			// Create Rois for TP and FP
			if (showTP)
			{
				float[] x = new float[matches.size()];
				float[] y = new float[x.length];
				int n = 0;
				for (PointPair pair : matches)
				{
					BasePoint p = (BasePoint) pair.getPoint2();
					x[n] = p.getX() + 0.5f;
					y[n] = p.getY() + 0.5f;
					n++;
				}
				addRoi(0, o, x, y, n, Color.green);
			}
			if (showFP)
			{
				float[] x = new float[predicted.length - matches.size()];
				float[] y = new float[x.length];
				int n = 0;
				for (Coordinate c : FP)
				{
					BasePoint p = (BasePoint) c;
					x[n] = p.getX() + 0.5f;
					y[n] = p.getY() + 0.5f;
					n++;
				}
				addRoi(0, o, x, y, n, Color.red);
			}
		}
		else
		{
			float[] x = new float[spots.length];
			float[] y = new float[x.length];
			for (int i = 0; i < spots.length; i++)
			{
				x[i] = spots[i].x + bounds.x + 0.5f;
				y[i] = spots[i].y + bounds.y + 0.5f;
			}
			PointRoi roi = new PointRoi(x, y);
			// Add options to configure colour and labels
			o.add(roi);
		}
		imp.setOverlay(o);
	}

	public static void addRoi(int frame, Overlay o, float[] x, float[] y, int n, Color colour)
	{
		if (n == 0)
			return;
		PointRoi roi = new PointRoi(x, y, n);
		roi.setFillColor(colour);
		roi.setStrokeColor(colour);
		if (frame != 0)
			roi.setPosition(frame);
		o.add(roi);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	public void setNPasses(int nPasses)
	{
		// Nothing to do		
	}

	public void imageOpened(ImagePlus imp)
	{

	}

	public void imageClosed(ImagePlus imp)
	{
	}

	public void imageUpdated(ImagePlus imp)
	{
		if (this.imp.getID() == imp.getID() && preview)
		{
			run(imp.getProcessor());
		}
	}
}
