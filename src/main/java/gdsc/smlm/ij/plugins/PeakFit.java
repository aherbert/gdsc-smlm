package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
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
import gdsc.smlm.engine.FitEngine;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitJob;
import gdsc.smlm.engine.FitParameters;
import gdsc.smlm.engine.FitParameters.FitTask;
import gdsc.smlm.engine.FitQueue;
import gdsc.smlm.engine.FitWorker;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.filters.SpotFilter;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitCriteria;
import gdsc.smlm.fitting.FitFunction;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import gdsc.smlm.function.CameraNoiseModel;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.SeriesImageSource;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsMode;
import gdsc.smlm.ij.results.ResultsTable;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.PSFCalculatorSettings;
import gdsc.smlm.ij.settings.ResultsSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.IJLogger;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.ij.utils.ImageROIPainter;
import gdsc.smlm.ij.utils.SeriesOpener;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.AggregatedImageSource;
import gdsc.smlm.results.BinaryFilePeakResults;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.FilePeakResults;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.InterlacedImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.utils.NoiseEstimator.Method;
import gdsc.smlm.utils.TextUtils;
import gdsc.smlm.utils.XmlUtils;
import gdsc.smlm.utils.logging.Logger;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JFileChooser;

import org.apache.commons.math3.util.FastMath;

//import ij.io.OpenDialog;

/**
 * Fits local maxima using a 2D Gaussian. Process each frame until a successive number of fits
 * fail to meet the fit criteria.
 */
public class PeakFit implements PlugInFilter, MouseListener, TextListener, ItemListener
{
	private static final String TITLE = "PeakFit";

	private static int FLAGS = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;
	private int plugin_flags;
	private int singleFrame = 0;
	private ImagePlus imp;

	// Used within the run(ImageProcessor) method
	private Rectangle bounds;
	private static boolean optionIgnoreBoundsForNoise = true;
	private boolean ignoreBoundsForNoise = false;
	private Logger logger = null;

	private ImageSource source = null;
	private PeakResultsList results;
	private long time;
	private Calibration calibration;
	private FitEngineConfiguration config = null;
	private FitConfiguration fitConfig;
	private ResultsSettings resultsSettings;
	private boolean silent = false;

	// Flag for extra options in the dialog (shift-key down)
	private boolean extraOptions = false;
	private boolean maximaIdentification = false;
	private boolean fitMaxima = false;
	private boolean simpleFit = false;
	private static int optionIntegrateFrames = 1;
	private int integrateFrames = 1;
	private static boolean optionInterlacedData = false;
	private static int optionDataStart = 1;
	private static int optionDataBlock = 1;
	private static int optionDataSkip = 0;
	private boolean interlacedData = false;
	private int dataStart = 1;
	private int dataBlock = 1;
	private int dataSkip = 0;
	private static boolean optionShowProcessedFrames = false;
	private boolean showProcessedFrames = false;

	private static String inputOption = "";
	private static boolean showTable = true;
	private static boolean showImage = true;
	private static PSFCalculatorSettings calculatorSettings = new PSFCalculatorSettings();

	// Used for the mouse listener
	private TextField textConfigFile;

	// All the fields that will be updated when reloading the configuration file
	private TextField textNmPerPixel;
	private TextField textGain;
	private Checkbox textEMCCD;
	private TextField textExposure;
	private TextField textInitialPeakStdDev0;
	private TextField textInitialPeakStdDev1;
	private TextField textInitialAngleD;
	private Choice textDataFilterType;
	private Choice textDataFilter;
	private TextField textSmooth;
	private TextField textSearch;
	private TextField textBorder;
	private TextField textFitting;
	private Choice textFitSolver;
	private Choice textFitFunction;
	private Checkbox textFitBackground;
	private TextField textFailuresLimit;
	private Checkbox textIncludeNeighbours;
	private TextField textNeighbourHeightThreshold;
	private TextField textResidualsThreshold;
	private TextField textDuplicateDistance;
	private TextField textCoordinateShiftFactor;
	private TextField textSignalStrength;
	private TextField textMinPhotons;
	private TextField textPrecisionThreshold;
	private TextField textNoise;
	private Choice textNoiseMethod;
	private TextField textWidthFactor;
	private Checkbox textLogProgress;
	private Checkbox textShowDeviations;
	private Choice textResultsTable;
	private Choice textResultsImage;
	private Checkbox textWeightedImage;
	private Checkbox textEqualisedImage;
	private TextField textPrecision;
	private TextField textImageScale;
	private TextField textImageRollingWindow;
	private TextField textResultsDirectory;
	private Checkbox textBinaryResults;
	private Checkbox textResultsInMemory;

	public PeakFit()
	{
		init(new FitEngineConfiguration(new FitConfiguration()), new ResultsSettings(), new Calibration());
	}

	public PeakFit(FitEngineConfiguration config)
	{
		init(config, new ResultsSettings(), new Calibration());
	}

	public PeakFit(FitEngineConfiguration config, ResultsSettings resultsSettings)
	{
		init(config, resultsSettings, new Calibration());
	}

	public PeakFit(FitEngineConfiguration config, ResultsSettings resultsSettings, Calibration calibration)
	{
		init(config, resultsSettings, calibration);
	}

	private void init(FitEngineConfiguration config, ResultsSettings resultsSettings, Calibration calibration)
	{
		this.config = config;
		this.resultsSettings = resultsSettings;
		this.calibration = calibration;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		plugin_flags = FLAGS;
		extraOptions = Utils.isExtraOptions();

		maximaIdentification = (arg != null && arg.contains("spot"));
		fitMaxima = (arg != null && arg.contains("maxima"));
		simpleFit = (arg != null && arg.contains("simple"));
		boolean runSeries = (arg != null && arg.contains("series"));

		ImageSource imageSource = null;
		if (fitMaxima)
		{
			imp = null;
			// The maxima will have been identified already. 
			// The image source will be found from the peak results.
			if (!showMaximaDialog())
				return DONE;

			MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
			if (results == null || results.size() == 0)
			{
				IJ.error(TITLE, "No results could be loaded");
				return DONE;
			}

			imageSource = results.getSource();
			plugin_flags |= NO_IMAGE_REQUIRED;
		}
		else if (runSeries)
		{
			imp = null;
			// Select input folder
			String inputDirectory;
			inputDirectory = IJ.getDirectory("Select image series ...");
			//inputDirectory = getInputDirectory("Select image series ...");
			if (inputDirectory == null)
				return DONE;

			// Load input series ...
			SeriesOpener series = new SeriesOpener(inputDirectory, true);
			if (series.getNumberOfImages() == 0)
			{
				IJ.error(TITLE, "No images in the selected directory:\n" + inputDirectory);
				return DONE;
			}

			imageSource = new SeriesImageSource(getName(series.getImageList()), series);
			((SeriesImageSource) imageSource).setLogProgress(true);
			plugin_flags |= NO_IMAGE_REQUIRED;
		}
		else
		{
			if (imp == null)
			{
				IJ.noImage();
				return DONE;
			}

			// Check it is not a previous result
			if (imp.getTitle().endsWith(IJImagePeakResults.IMAGE_SUFFIX))
			{
				IJImageSource tmpImageSource = null;

				// Check the image to see if it has an image source XML structure in the info property
				Object o = imp.getProperty("Info");
				Pattern pattern = Pattern.compile("Source: (<.*IJImageSource>.*<.*IJImageSource>)", Pattern.DOTALL);
				Matcher match = pattern.matcher((o == null) ? "" : o.toString());
				if (match.find())
				{
					ImageSource source = ImageSource.fromXML(match.group(1));
					if (source instanceof IJImageSource)
					{
						tmpImageSource = (IJImageSource) source;
						if (!tmpImageSource.open())
						{
							tmpImageSource = null;
						}
						else
						{
							imp = WindowManager.getImage(tmpImageSource.getName());
						}
					}
				}

				if (tmpImageSource == null)
				{
					// Look for a parent using the title
					String parentTitle = imp.getTitle().substring(0,
							imp.getTitle().length() - IJImagePeakResults.IMAGE_SUFFIX.length() - 1);
					ImagePlus parentImp = WindowManager.getImage(parentTitle);
					if (parentImp != null)
					{
						tmpImageSource = new IJImageSource(parentImp);
						imp = parentImp;
					}
				}
				String message = "The selected image may be a previous fit result";
				if (tmpImageSource != null)
					message += " of: \n \n" + tmpImageSource.getName() + " \n \nFit the parent?";
				else
					message += " \n \nDo you want to continue?";

				YesNoCancelDialog d = new YesNoCancelDialog(null, TITLE, message);
				if (tmpImageSource == null)
				{
					if (!d.yesPressed())
						return DONE;
				}
				else
				{
					if (d.yesPressed())
						imageSource = tmpImageSource;
					if (d.cancelPressed())
						return DONE;
				}
			}

			if (imageSource == null)
				imageSource = new IJImageSource(imp);
		}

		time = -1;

		if (!initialiseImage(imageSource, getBounds(imp), false))
		{
			IJ.error(TITLE, "Failed to initialise the source image: " + imageSource.getName());
			return DONE;
		}

		int flags = showDialog(imp);
		if ((flags & DONE) == 0)
		{
			// Repeat so that we pass in the selected option for ignoring the bounds.
			// This should not be necessary since it is set within the readDialog method.
			//if (ignoreBoundsForNoise)
			//	initialiseImage(imageSource, bounds, ignoreBoundsForNoise);
			initialiseFitting();
		}
		return flags;
	}

	private String getName(String[] imageList)
	{
		String name = imageList[0];
		// Remove directory
		int index = name.lastIndexOf(File.separatorChar);
		if (index > -1)
		{
			name = name.substring(index + 1);
		}
		//// Remove suffix
		//index = name.lastIndexOf('.');
		//if (index > 0)
		//{
		//	name = name.substring(0, index);
		//}
		return "Series " + name;
	}

	private Rectangle getBounds(ImagePlus imp)
	{
		if (imp == null)
			return null;
		Roi roi = imp.getRoi();
		if (roi != null && roi.isArea())
		{
			return roi.getBounds();
		}
		return null;
	}

	/**
	 * @return An input directory containing a series of images
	 */
	@SuppressWarnings("unused")
	private String getInputDirectory(String title)
	{
		final JFileChooser chooser = new JFileChooser()
		{
			private static final long serialVersionUID = 275144634537614122L;

			public void approveSelection()
			{
				if (getSelectedFile().isFile())
				{
					return;
				}
				else
					super.approveSelection();
			}
		};
		if (System.getProperty("os.name").startsWith("Mac OS X"))
		{
			chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		}
		else
		{
			chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		}
		chooser.setDialogTitle(title);
		int returnVal = chooser.showOpenDialog(IJ.getInstance());
		if (returnVal == JFileChooser.APPROVE_OPTION)
		{
			return chooser.getSelectedFile().getPath();
		}
		return null;
	}

	/**
	 * Initialise a new image for fitting and prepare the output results.
	 * <p>
	 * Calls {@link #initialise(ImageSource, Rectangle)} then {@link #initialiseFitting()}.
	 * 
	 * @param imageSource
	 *            The image source
	 * @param bounds
	 *            The region to process from the image
	 * @param ignoreBoundsForNoise
	 *            Set to true if the bounds should be ignored when computing the noise estimate for each frame
	 * @return True if the image was valid and the initialisation was successful
	 */
	public boolean initialise(ImageSource imageSource, Rectangle bounds, boolean ignoreBoundsForNoise)
	{
		if (!initialiseImage(imageSource, bounds, ignoreBoundsForNoise))
			return false;
		return initialiseFitting();
	}

	/**
	 * Initialise a new image.
	 * <p>
	 * Does not set-up for fitting. This can be done using a subsequent call to {@link #initialiseFitting()}.
	 * <p>
	 * This mechanism allows additional result outputs to be added after initialisation using
	 * {@link #addPeakResults(PeakResults)}.
	 * 
	 * @param imageSource
	 *            The image source
	 * @param bounds
	 *            The region to process from the image
	 * @param ignoreBoundsForNoise
	 *            Set to true if the bounds should be ignored when computing the noise estimate for each frame
	 * @return
	 *         True if the image was valid and the initialisation was successful
	 */
	public boolean initialiseImage(ImageSource imageSource, Rectangle bounds, boolean ignoreBoundsForNoise)
	{
		// Initialise for image processing
		if (!setSource(imageSource))
			return false;

		this.ignoreBoundsForNoise = ignoreBoundsForNoise;
		if (bounds == null)
		{
			bounds = new Rectangle(0, 0, source.getWidth(), source.getHeight());
			// No region so no need to ignore the bounds.
			this.ignoreBoundsForNoise = false;
		}
		this.bounds = bounds;
		results = new PeakResultsList();

		time = 0;
		//config = null;

		return true;
	}

	private boolean setSource(ImageSource imageSource)
	{
		// Reset
		this.source = null;
		if (imageSource == null)
			return false;

		// Open the image to ensure it is accessible and the width/height are known
		if (!imageSource.open())
			return false;

		this.source = imageSource;
		return true;
	}

	/**
	 * Set-up the fitting using all the configured properties. Prepare the output results.
	 */
	public boolean initialiseFitting()
	{
		if (source == null)
			return false;

		// Do this to ensure the serialised configuration is correct
		updateFitConfiguration(config);

		results.setSource(source);
		if (maximaIdentification)
			results.setName(source.getName() + " (Maxima)");
		else if (fitMaxima)
			results.setName(source.getName() + " (" + getSolverName() + " Fit Maxima)");
		else
			results.setName(source.getName() + " (" + getSolverName() + ")");
		results.setBounds(bounds);
		Calibration cal = calibration.clone();
		// Account for the frame integration
		// TODO - Should we change this so that if integrate frames is used then the data 
		// are converted to ExtendedPeakResult with a start and end frame
		//cal.exposureTime *= integrateFrames;
		//if (interlacedData)
		//{
		//	cal.exposureTime *= ((double)dataBlock / (dataBlock + dataSkip));
		//}
		results.setCalibration(cal);
		results.setConfiguration(XmlUtils.toXML(config));

		addMemoryResults(results, false);
		addImageResults(results);
		addFileResults(results);
		addTableResults(results);
		addDefaultResults(results);

		results.begin();

		if (simpleFit && showImage)
		{
			for (PeakResults r : results.toArray())
			{
				if (r instanceof IJImagePeakResults)
				{
					ImagePlus i = ((IJImagePeakResults) r).getImagePlus();
					Utils.log("Super-resolution image title = " + i.getTitle());
					WindowManager.toFront(i.getWindow());
				}
			}
		}

		return true;
	}

	private String getSolverName()
	{
		return getSolverName(config.getFitConfiguration());
	}

	public static String getSolverName(FitConfiguration fitConfig)
	{
		FitSolver solver = fitConfig.getFitSolver();
		String name = solver.getShortName();
		if (solver == FitSolver.MLE)
			name += " " + fitConfig.getSearchMethod();
		return name;
	}

	protected void showResults()
	{
		IJ.showProgress(1.0);
		if (time >= 0)
		{
			if (silent)
			{
				results.end();
				return;
			}

			// Check if we are sorting
			IJ.showStatus("Finalising results ...");
			for (PeakResults r : results.toArray())
				if (r instanceof MemoryPeakResults)
				{
					if (((MemoryPeakResults) r).isSortAfterEnd())
						;
					IJ.showStatus("Sorting " + r.size() + " results ...");
					break;
				}

			results.end();

			String textTime = Utils.timeToString(time / 1000000.0);

			int size = getSize();
			String message = String.format("%d Localisation%s. Fitting Time = %s", size, (size == 1) ? "" : "s",
					textTime);
			if (resultsSettings.logProgress)
				IJ.log("-=-=-=-");
			IJ.log(message);
			IJ.showStatus(message);
		}
		else
		{
			IJ.showStatus("");
		}
	}

	private boolean showMaximaDialog()
	{
		int size = MemoryPeakResults.countMemorySize();
		if (size == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return false;
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Select identified maxima for fitting");

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);

		return true;
	}

	@SuppressWarnings("unchecked")
	private int showDialog(ImagePlus imp)
	{
		// Executing as an ImageJ plugin.
		// Override the defaults with those in the configuration file
		String filename = SettingsManager.getSettingsFilename();

		if (simpleFit)
		{
			return showSimpleDialog(filename);
		}

		GlobalSettings settings = SettingsManager.loadSettings(filename);
		calibration = settings.getCalibration();
		config = settings.getFitEngineConfiguration();
		fitConfig = config.getFitConfiguration();
		resultsSettings = settings.getResultsSettings();

		boolean isCrop = (bounds != null && imp != null && (bounds.width < imp.getWidth() || bounds.height < imp
				.getHeight()));

		if (!extraOptions)
		{
			integrateFrames = 1;
			resultsSettings.imageRollingWindow = 0;
			fitConfig.setBackgroundFitting(true);
			fitConfig.setMinIterations(0);
			fitConfig.setNoise(0);
			config.setNoiseMethod(Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);
			showProcessedFrames = false;
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage((maximaIdentification) ? "Identify candidate maxima" : "Fit 2D Gaussian to identified maxima");

		String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);

		gd.addStringField("Config_file", filename, 40);
		gd.addNumericField("Calibration (nm/px)", calibration.nmPerPixel, 2);
		gd.addNumericField("Gain (ADU/photon)", calibration.gain, 2);
		gd.addCheckbox("EM-CCD", calibration.emCCD);
		gd.addNumericField("Exposure_time (ms)", calibration.exposureTime, 2);

		if (isCrop)
			gd.addCheckbox("Ignore_bounds_for_noise", optionIgnoreBoundsForNoise);
		// This is already set to false before the dialog is displayed
		//else
		//	ignoreBoundsForNoise = false;

		gd.addNumericField("Initial_StdDev0", fitConfig.getInitialPeakStdDev0(), 3);
		if (!maximaIdentification)
		{
			gd.addNumericField("Initial_StdDev1", fitConfig.getInitialPeakStdDev1(), 3);
			gd.addNumericField("Initial_Angle", fitConfig.getInitialAngle(), 3);
		}
		String[] filterTypes = SettingsManager.getNames((Object[]) DataFilterType.values());
		gd.addChoice("Spot_filter_type", filterTypes, filterTypes[config.getDataFilterType().ordinal()]);
		String[] filterNames = SettingsManager.getNames((Object[]) DataFilter.values());
		gd.addChoice("Spot_filter", filterNames, filterNames[config.getDataFilter(0).ordinal()]);
		gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
		gd.addSlider("Search_width", 0.5, 2.5, config.getSearch());
		gd.addSlider("Border", 0.5, 2.5, config.getBorder());
		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());
		if (extraOptions && !fitMaxima)
		{
			gd.addCheckbox("Interlaced_data", optionInterlacedData);
			gd.addSlider("Integrate_frames", 1, 5, optionIntegrateFrames);
		}

		Component discardLabel = null;
		if (!maximaIdentification)
		{
			gd.addMessage("--- Gaussian fitting ---");
			String[] solverNames = SettingsManager.getNames((Object[]) FitSolver.values());
			gd.addChoice("Fit_solver", solverNames, solverNames[fitConfig.getFitSolver().ordinal()]);
			String[] functionNames = SettingsManager.getNames((Object[]) FitFunction.values());
			gd.addChoice("Fit_function", functionNames, functionNames[fitConfig.getFitFunction().ordinal()]);
			if (extraOptions)
				gd.addCheckbox("Fit_background", fitConfig.isBackgroundFitting());

			// Parameters specific to each Fit solver are collected in a second dialog 

			gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
			gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
			gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
			gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

			gd.addSlider("Duplicate_distance", 0, 1.5, fitConfig.getDuplicateDistance());

			gd.addMessage("--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");
			discardLabel = gd.getMessage();

			gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
			gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
			gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
			if (extraOptions)
			{
				gd.addNumericField("Noise", fitConfig.getNoise(), 2);
				String[] noiseMethodNames = SettingsManager.getNames((Object[]) Method.values());
				gd.addChoice("Noise_method", noiseMethodNames, noiseMethodNames[config.getNoiseMethod().ordinal()]);
			}
			gd.addSlider("Width_factor", 0.01, 5, fitConfig.getWidthFactor());
			gd.addNumericField("Precision", fitConfig.getPrecisionThreshold(), 2);
		}

		gd.addMessage("--- Results ---");
		gd.addCheckbox("Log_progress", resultsSettings.logProgress);
		if (!maximaIdentification)
		{
			gd.addCheckbox("Show_deviations", resultsSettings.showDeviations);
		}
		String[] tableNames = SettingsManager.getNames((Object[]) ResultsTable.values());
		gd.addChoice("Results_table", tableNames, tableNames[resultsSettings.getResultsTable().ordinal()]);
		String[] imageNames = SettingsManager.getNames((Object[]) ResultsImage.values());
		gd.addMessage("--- Image output ---");
		gd.addChoice("Image", imageNames, imageNames[resultsSettings.getResultsImage().ordinal()]);
		gd.addCheckbox("Weighted", resultsSettings.weightedImage);
		gd.addCheckbox("Equalised", resultsSettings.equalisedImage);
		gd.addSlider("Image_Precision (nm)", 5, 30, resultsSettings.precision);
		gd.addSlider("Image_Scale", 1, 15, resultsSettings.imageScale);
		if (extraOptions)
		{
			gd.addNumericField("Image_window", resultsSettings.imageRollingWindow, 0);
			gd.addCheckbox("Show_processed_frames", optionShowProcessedFrames);
		}
		gd.addMessage("--- File output ---");
		gd.addStringField("Results_dir", resultsSettings.resultsDirectory);
		gd.addCheckbox("Binary_results", resultsSettings.binaryResults);
		gd.addMessage(" ");
		gd.addCheckbox("Results_in_memory", resultsSettings.resultsInMemory);

		// Re-arrange the standard layout which has a GridBagLayout with 2 columns (label,field)
		// to 4 columns: (label,field) x 2

		if (gd.getLayout() != null)
		{
			GridBagLayout grid = (GridBagLayout) gd.getLayout();

			int xOffset = 0, yOffset = 0;
			int lastY = -1, rowCount = 0;
			for (Component comp : gd.getComponents())
			{
				// Check if this should be the second major column
				if (comp == discardLabel)
				{
					xOffset += 2;
					yOffset -= rowCount;
				}
				// Reposition the field
				GridBagConstraints c = grid.getConstraints(comp);
				if (lastY != c.gridy)
					rowCount++;
				lastY = c.gridy;
				c.gridx = c.gridx + xOffset;
				c.gridy = c.gridy + yOffset;
				c.insets.left = c.insets.left + 10 * xOffset;
				c.insets.top = 0;
				c.insets.bottom = 0;
				grid.setConstraints(comp, c);
			}

			if (IJ.isLinux())
				gd.setBackground(new Color(238, 238, 238));
		}

		// Add a mouse listener to the config file field
		if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro()))
		{
			Vector<TextField> texts = (Vector<TextField>) gd.getStringFields();
			Vector<TextField> numerics = (Vector<TextField>) gd.getNumericFields();
			Vector<Checkbox> checkboxes = (Vector<Checkbox>) gd.getCheckboxes();
			Vector<Choice> choices = (Vector<Choice>) gd.getChoices();

			int n = 0;
			int t = 0;
			int b = 0;
			int ch = 0;

			Choice textTemplate = choices.get(ch++);
			textTemplate.addItemListener(this);

			textConfigFile = texts.get(t++);
			textConfigFile.addMouseListener(this);
			textConfigFile.addTextListener(this);

			// TODO: add a value changed listener to detect when typing a new file

			textNmPerPixel = numerics.get(n++);
			textGain = numerics.get(n++);
			textEMCCD = checkboxes.get(b++);
			textExposure = numerics.get(n++);
			textInitialPeakStdDev0 = numerics.get(n++);
			if (!maximaIdentification)
			{
				textInitialPeakStdDev1 = numerics.get(n++);
				textInitialAngleD = numerics.get(n++);
			}
			textDataFilterType = choices.get(ch++);
			textDataFilter = choices.get(ch++);
			textSmooth = numerics.get(n++);
			textSearch = numerics.get(n++);
			textBorder = numerics.get(n++);
			textFitting = numerics.get(n++);
			if (extraOptions && !fitMaxima)
			{
				b++; // Skip over the interlaced data option
				n++; // Skip over the integrate frames option
			}
			if (!maximaIdentification)
			{
				textFitSolver = choices.get(ch++);
				textFitFunction = choices.get(ch++);
				if (extraOptions)
					textFitBackground = checkboxes.get(b++);
				textFailuresLimit = numerics.get(n++);
				textIncludeNeighbours = checkboxes.get(b++);
				textNeighbourHeightThreshold = numerics.get(n++);
				textResidualsThreshold = numerics.get(n++);
				textDuplicateDistance = numerics.get(n++);
				textCoordinateShiftFactor = numerics.get(n++);
				textSignalStrength = numerics.get(n++);
				textMinPhotons = numerics.get(n++);
				textWidthFactor = numerics.get(n++);
				textPrecisionThreshold = numerics.get(n++);
				if (extraOptions)
				{
					textNoise = numerics.get(n++);
					textNoiseMethod = choices.get(ch++);
				}
			}
			textLogProgress = checkboxes.get(b++);
			if (!maximaIdentification)
				textShowDeviations = checkboxes.get(b++);
			textResultsTable = choices.get(ch++);
			textResultsImage = choices.get(ch++);
			textWeightedImage = checkboxes.get(b++);
			textEqualisedImage = checkboxes.get(b++);
			textPrecision = numerics.get(n++);
			textImageScale = numerics.get(n++);
			if (extraOptions)
			{
				textImageRollingWindow = numerics.get(n++);
				b++; // Skip over show processed frames option
			}
			textResultsDirectory = texts.get(t++);
			textResultsDirectory.addMouseListener(this);

			textBinaryResults = checkboxes.get(b++);
			textResultsInMemory = checkboxes.get(b++);
		}

		gd.showDialog();
		
		// The refreshSettings method can be called by the dialog listener.
		// This updates the Calibration, FitEngineConfiguration, and ResultsSettings so set these
		// back in the GlobalSettings object.
		settings.setCalibration(this.calibration);
		settings.setFitEngineConfiguration(this.config);
		settings.setResultsSettings(this.resultsSettings);

		if (gd.wasCanceled() || !readDialog(settings, gd, isCrop))
			return DONE;

		if (imp != null)
		{
			// Store whether the user selected to process all the images.
			int flags = IJ.setupDialog(imp, plugin_flags);

			// Check if cancelled
			if ((flags & DONE) != 0)
				return DONE;

			if ((flags & DOES_STACKS) == 0)
			{
				// Save the slice number for the overlay
				singleFrame = imp.getCurrentSlice();

				// Account for interlaced data
				if (interlacedData)
				{
					int start = singleFrame;

					// Calculate the first frame that is not skipped
					while (ignoreFrame(start) && start > dataStart)
						start--;
					if (start < dataStart)
					{
						log("The current frame (%d) is before the start of the interlaced data", singleFrame);
						return DONE;
					}
					if (start != singleFrame)
						log("Updated the current frame (%d) to a valid interlaced data frame (%d)", singleFrame, start);
					singleFrame = start;
				}

				// Account for integrated frames
				int endFrame = singleFrame;
				if (integrateFrames > 1)
				{
					int totalFrames = 1;
					while (totalFrames < integrateFrames)
					{
						endFrame++;
						if (!ignoreFrame(endFrame))
							totalFrames++;
					}
					log("Updated the image end frame (%d) to %d allow %d integrated frames", singleFrame, endFrame,
							integrateFrames);
				}

				// Create a new image source with the correct frames
				setSource(new IJImageSource(imp, singleFrame, endFrame - singleFrame));

				// Store the image so the results can be added as an overlay
				this.imp = imp;
				this.imp.setOverlay(null);
			}
		}

		// Allow interlaced data by wrapping the image source
		if (interlacedData)
		{
			setSource(new InterlacedImageSource(this.source, dataStart, dataBlock, dataSkip));
		}

		// Allow frame aggregation by wrapping the image source
		if (integrateFrames > 1)
		{
			setSource(new AggregatedImageSource(this.source, integrateFrames));
		}

		// Ask if the user wants to log progress on multiple frame images
		if (resultsSettings.logProgress && source.getFrames() > 1)
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Warning: Log progress on multiple-frame image will be slow");
			gd.addCheckbox("Log_progress", resultsSettings.logProgress);
			gd.showDialog();
			if (gd.wasCanceled())
				return DONE;
			resultsSettings.logProgress = gd.getNextBoolean();
			if (!resultsSettings.logProgress)
				SettingsManager.saveSettings(settings, filename);
		}

		// Get a bias if required
		if (resultsSettings.getResultsTable() == ResultsTable.CALIBRATED && calibration.bias == 0)
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Calibrated results requires a camera bias");
			gd.addNumericField("Camera_bias (ADUs)", calibration.bias, 2);
			gd.showDialog();
			if (!gd.wasCanceled())
			{
				calibration.bias = Math.abs(gd.getNextNumber());
				if (calibration.bias > 0)
					SettingsManager.saveSettings(settings, filename);
			}
		}

		// Return the plugin flags (without the DOES_STACKS flag).
		// The call to run(ImageProcessor) will process the image in 'this.imp' so we only want a 
		// single call to be made.
		return plugin_flags;
	}

	private void log(String format, Object... args)
	{
		if (!silent)
			Utils.log(format, args);
	}

	private int showSimpleDialog(String filename)
	{
		GlobalSettings settings = SettingsManager.loadSettings(filename);
		// Initialise the fit config so that it can be used in the calibration wizard 
		fitConfig = settings.getFitEngineConfiguration().getFitConfiguration();

		boolean requireCalibration = requireCalibration(settings, filename);
		if (requireCalibration)
		{
			if (!showCalibrationWizard(settings, true))
				return DONE;
		}

		// Present dialog with simple output options: Image, Table
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Fit single-molecule localisations");

		if (!requireCalibration)
			gd.addCheckbox("Use_current_calibration", true);
		gd.addCheckbox("Show_table", showTable);
		gd.addCheckbox("Show_image", showImage);

		gd.showDialog();
		if (gd.wasCanceled())
			return DONE;

		boolean useCurrentCalibration = true;
		if (!requireCalibration)
			useCurrentCalibration = gd.getNextBoolean();
		showTable = gd.getNextBoolean();
		showImage = gd.getNextBoolean();

		if (!useCurrentCalibration)
		{
			if (!showCalibrationWizard(settings, false))
				return DONE;
		}

		// Restore fitting to default settings but maintain the calibrated width
		final double sd = fitConfig.getInitialPeakStdDev0();
		config = new FitEngineConfiguration(new FitConfiguration());
		fitConfig = config.getFitConfiguration();
		fitConfig.setInitialPeakStdDev(sd);
		fitConfig.setCoordinateShiftFactor(1);
		resultsSettings = new ResultsSettings();

		// Do simple results output
		resultsSettings.resultsInMemory = true;
		resultsSettings.setResultsTable((showTable) ? ResultsTable.UNCALIBRATED : ResultsTable.NONE);
		if (showImage)
		{
			resultsSettings.setResultsImage(ResultsImage.SIGNAL_INTENSITY);
			resultsSettings.imageScale = Math.ceil(1024 / (FastMath.max(bounds.width, bounds.height)));
			resultsSettings.weightedImage = true;
			resultsSettings.equalisedImage = true;
		}
		else
		{
			resultsSettings.setResultsImage(ResultsImage.NONE);
		}

		// Log the settings we care about:
		calibration = settings.getCalibration();
		IJ.log("-=-=-=-");
		IJ.log("Peak Fit");
		IJ.log("-=-=-=-");
		Utils.log("Pixel pitch = %s", Utils.rounded(calibration.nmPerPixel, 4));
		Utils.log("Exposure Time = %s", Utils.rounded(calibration.exposureTime, 4));
		Utils.log("Gain = %s", Utils.rounded(calibration.gain, 4));
		Utils.log("PSF width = %s", Utils.rounded(fitConfig.getInitialPeakStdDev0(), 4));

		// Save
		settings.setFitEngineConfiguration(config);
		settings.setResultsSettings(resultsSettings);
		SettingsManager.saveSettings(settings, filename);

		return FLAGS;
	}

	/**
	 * Check if the configuration file exists. If not then assume this is the first run of the plugin.
	 * Check the calibration is valid for fitting.
	 * 
	 * @param settings
	 * @param filename
	 * @return
	 */
	private boolean requireCalibration(GlobalSettings settings, String filename)
	{
		if (!new File(filename).exists())
			return true;

		calibration = settings.getCalibration();

		// Check if the calibration contains: Pixel pitch, Gain (can be 1), Exposure time
		if (calibration.nmPerPixel <= 0)
			return true;
		if (calibration.gain <= 0)
			return true;
		if (calibration.exposureTime <= 0)
			return true;

		// Check for a PSF width
		if (fitConfig.getInitialPeakStdDev0() <= 0)
			return true;

		return false;
	}

	private boolean showCalibrationWizard(GlobalSettings settings, boolean showIntroduction)
	{
		if (showIntroduction)
		{
			GenericDialog gd = newWizardDialog("No configuration file could be loaded.",
					"Please follow the configuration wizard to calibrate.");
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
		}

		//Calibration defaultCalibration = new Calibration();
		//if (calibration.nmPerPixel <= 0 || calibration.nmPerPixel == defaultCalibration.nmPerPixel)
		if (!getPixelPitch())
			return false;
		//if (calibration.gain <= 0 || calibration.gain == defaultCalibration.gain)
		if (!getGain())
			return false;
		//if (calibration.exposureTime <= 0 || calibration.exposureTime == defaultCalibration.exposureTime)
		if (!getExposureTime())
			return false;
		// Check for a PSF width other than the default
		//if (fitConfig.getInitialPeakWidth0() == new FitConfiguration().getInitialPeakWidth0())
		if (!getPeakWidth())
			return false;

		// Check parameters
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.nmPerPixel);
			Parameters.isAboveZero("Gain", calibration.gain);
			Parameters.isAboveZero("Exposure time", calibration.exposureTime);
			Parameters.isAboveZero("Initial SD", fitConfig.getInitialPeakStdDev0());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private GenericDialog newWizardDialog(String... messages)
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		final String header = "-=-";
		gd.addMessage(header + " " + TITLE + " Configuration Wizard " + header);
		for (String message : messages)
			gd.addMessage(TextUtils.wrap(message, 80));
		return gd;
	}

	private boolean getPixelPitch()
	{
		GenericDialog gd = newWizardDialog(
				"Enter the size of each pixel. This is required to ensure the dimensions of the image are calibrated.",
				"E.g. a camera with a 6.45um pixel size and a 60x objective will have a pitch of 6450/60 = 107.5nm.");
		// TODO - Add a pop-up calculator...
		gd.addNumericField("Calibration (nm/px)", calibration.nmPerPixel, 2);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.nmPerPixel = Math.abs(gd.getNextNumber());
		return true;
	}

	private boolean getGain()
	{
		GenericDialog gd = newWizardDialog(
				"Enter the total gain.",
				"This is usually supplied with your camera certificate. The gain indicates how many Analogue-to-Digital-Units (ADUs) are recorded at the pixel for each photon registered on the sensor.",
				"The gain is usually expressed using the product of the EM-gain (if applicable), the camera gain and the sensor quantum efficiency.",
				"A value of 1 means no conversion to photons will occur.");
		// TODO - Add a wizard to allow calculation of total gain from EM-gain, camera gain and QE
		gd.addNumericField("Gain (ADU/photon)", calibration.gain, 2);
		gd.addCheckbox("EM-CCD", calibration.emCCD);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.gain = Math.abs(gd.getNextNumber());
		calibration.emCCD = gd.getNextBoolean();
		return true;
	}

	private boolean getExposureTime()
	{
		GenericDialog gd = newWizardDialog(
				"Enter the exposure time. Calibration of the exposure time allows correct reporting of on and off times.",
				"This is the length of time for each frame in the image.");
		gd.addNumericField("Exposure_time (ms)", calibration.exposureTime, 2);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.exposureTime = Math.abs(gd.getNextNumber());
		return true;
	}

	private boolean getPeakWidth()
	{
		GenericDialog gd = newWizardDialog(
				"Enter the expected peak width in pixels.",
				"A point source of light will not be focussed perfectly by the microscope but will appear as a spread out peak. This Point Spread Function (PSF) can be modelled using a 2D Gaussian curve.",
				"An optimised optical system (lens and camera sensor) should have a peak standard deviation of approximately 1 pixel when in focus. This allows the fitting routine to have enough data to identify the centre of the peak without spreading the light over too many pixels (which increases noise).",
				"The peak width can be estimated using the wavelength of light emitted by the single molecules and the parameters of the microscope. Use a PSF calculator by clicking the checkbox below:");
		// Add ability to run the PSF Calculator to get the width
		gd.addCheckbox("Run_PSF_calculator", false);
		gd.addNumericField("Gaussian_SD", fitConfig.getInitialPeakStdDev0(), 3);
		if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro()))
		{
			Checkbox cb = (Checkbox) gd.getCheckboxes().get(0);
			cb.addItemListener(this);
			textInitialPeakStdDev0 = (TextField) gd.getNumericFields().get(0);
		}
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		fitConfig.setInitialPeakStdDev(Math.abs(gd.getNextNumber()));
		return true;
	}

	public void itemStateChanged(ItemEvent e)
	{
		if (e.getSource() instanceof Choice)
		{
			// Update the settings from the template
			Choice choice = (Choice) e.getSource();
			String templateName = choice.getSelectedItem();
			//System.out.println("Update to " + templateName);

			// Get the configuration template
			GlobalSettings template = ConfigurationTemplate.getTemplate(templateName);

			if (template != null)
			{
				boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
				if (template.isFitEngineConfiguration())
				{
					FitEngineConfiguration templateConfig = template.getFitEngineConfiguration().clone();

					// If this is a custom template then use all the settings.
					// If a default template then copy the existing spot settings.
					if (!custom)
					{
						FitConfiguration templateFitConfig = templateConfig.getFitConfiguration();
						templateFitConfig.setInitialPeakStdDev0(fitConfig.getInitialPeakStdDev0());
						templateFitConfig.setInitialPeakStdDev1(fitConfig.getInitialPeakStdDev1());
						templateFitConfig.setInitialAngle(fitConfig.getInitialAngle());
						templateFitConfig.setFitFunction(fitConfig.getFitFunction());
					}

					refreshSettings(templateConfig);
				}
				if (template.isCalibration())
				{
					refreshSettings(template.getCalibration().clone());
				}
				if (template.isResultsSettings())
				{
					refreshSettings(template.getResultsSettings().clone());
				}
			}
		}
		else if (e.getSource() instanceof Checkbox)
		{
			// Run the PSF Calculator
			Checkbox cb = (Checkbox) e.getSource();
			if (cb.getState())
			{
				cb.setState(false);
				PSFCalculator calculator = new PSFCalculator();
				calculatorSettings.pixelPitch = calibration.nmPerPixel / 1000.0;
				calculatorSettings.magnification = 1;
				calculatorSettings.beamExpander = 1;
				double sd = calculator.calculate(calculatorSettings, true);
				if (sd > 0)
					textInitialPeakStdDev0.setText(Double.toString(sd));
			}
		}
	}

	private boolean readDialog(GlobalSettings settings, GenericDialog gd, boolean isCrop)
	{
		// Ignore the template
		gd.getNextChoice();

		String filename = gd.getNextString();

		calibration.nmPerPixel = Math.abs(gd.getNextNumber());
		calibration.gain = Math.abs(gd.getNextNumber());
		calibration.emCCD = gd.getNextBoolean();
		calibration.exposureTime = Math.abs(gd.getNextNumber());
		// Note: The bias and read noise will just end up being what was in the configuration file
		// One fix for this is to save/load only the settings that are required from the configuration file
		// (the others will remain unchanged). This will require a big refactor of the settings save/load.
		// The simple fix is to create a plugin to allow the configuration to be changed for results.
		if (isCrop)
			ignoreBoundsForNoise = optionIgnoreBoundsForNoise = gd.getNextBoolean();

		fitConfig.setInitialPeakStdDev0(gd.getNextNumber());
		if (!maximaIdentification)
		{
			fitConfig.setInitialPeakStdDev1(gd.getNextNumber());
			fitConfig.setInitialAngleD(gd.getNextNumber());
		}
		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		config.setFitting(gd.getNextNumber());
		if (extraOptions && !fitMaxima)
		{
			interlacedData = gd.getNextBoolean();
			integrateFrames = optionIntegrateFrames = (int) gd.getNextNumber();
		}

		if (!maximaIdentification)
		{
			fitConfig.setFitSolver(gd.getNextChoiceIndex());
			fitConfig.setFitFunction(gd.getNextChoiceIndex());
			if (extraOptions)
				fitConfig.setBackgroundFitting(gd.getNextBoolean());
			config.setFailuresLimit((int) gd.getNextNumber());
			config.setIncludeNeighbours(gd.getNextBoolean());
			config.setNeighbourHeightThreshold(gd.getNextNumber());
			config.setResidualsThreshold(gd.getNextNumber());

			fitConfig.setDuplicateDistance(gd.getNextNumber());

			fitConfig.setCoordinateShiftFactor(gd.getNextNumber());
			fitConfig.setSignalStrength(gd.getNextNumber());
			fitConfig.setMinPhotons(gd.getNextNumber());
			if (extraOptions)
			{
				fitConfig.setNoise(gd.getNextNumber());
				config.setNoiseMethod(gd.getNextChoiceIndex());
			}
			fitConfig.setWidthFactor(gd.getNextNumber());
			fitConfig.setPrecisionThreshold(gd.getNextNumber());
		}

		resultsSettings.logProgress = gd.getNextBoolean();
		if (!maximaIdentification)
			resultsSettings.showDeviations = gd.getNextBoolean();

		resultsSettings.setResultsTable(gd.getNextChoiceIndex());
		resultsSettings.setResultsImage(gd.getNextChoiceIndex());
		resultsSettings.weightedImage = gd.getNextBoolean();
		resultsSettings.equalisedImage = gd.getNextBoolean();
		resultsSettings.precision = gd.getNextNumber();
		resultsSettings.imageScale = gd.getNextNumber();
		if (extraOptions)
		{
			resultsSettings.imageRollingWindow = (int) gd.getNextNumber();
			showProcessedFrames = optionShowProcessedFrames = gd.getNextBoolean();
		}
		resultsSettings.resultsDirectory = gd.getNextString();
		resultsSettings.binaryResults = gd.getNextBoolean();
		resultsSettings.resultsInMemory = gd.getNextBoolean();

		// Save to allow dialog state to be maintained even with invalid parameters
		SettingsManager.saveSettings(settings, filename);

		if (gd.invalidNumber())
			return false;

		// Check arguments
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.nmPerPixel);
			Parameters.isAboveZero("Gain", calibration.gain);
			Parameters.isAboveZero("Exposure time", calibration.exposureTime);
			Parameters.isAboveZero("Initial SD0", fitConfig.getInitialPeakStdDev0());
			if (!maximaIdentification)
			{
				Parameters.isAboveZero("Initial SD1", fitConfig.getInitialPeakStdDev1());
				Parameters.isPositive("Initial angle", fitConfig.getInitialAngleD());
			}
			Parameters.isAboveZero("Search_width", config.getSearch());
			Parameters.isAboveZero("Fitting_width", config.getFitting());
			Parameters.isPositive("Integrate frames", integrateFrames);
			if (!maximaIdentification)
			{
				Parameters.isPositive("Failures limit", config.getFailuresLimit());
				Parameters.isPositive("Neighbour height threshold", config.getNeighbourHeightThreshold());
				Parameters.isPositive("Residuals threshold", config.getResidualsThreshold());
				Parameters.isPositive("Duplicate distance", fitConfig.getDuplicateDistance());
				Parameters.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
				Parameters.isPositive("Signal strength", fitConfig.getSignalStrength());
				Parameters.isPositive("Min photons", fitConfig.getMinPhotons());
				if (extraOptions)
					Parameters.isPositive("Noise", fitConfig.getNoise());
				Parameters.isPositive("Width factor", fitConfig.getWidthFactor());
				Parameters.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
			}
			if (resultsSettings.getResultsImage() == ResultsImage.SIGNAL_AV_PRECISION ||
					resultsSettings.getResultsImage() == ResultsImage.LOCALISATIONS_AV_PRECISION)
			{
				Parameters.isAboveZero("Image precision", resultsSettings.precision);
			}
			Parameters.isAboveZero("Image scale", resultsSettings.imageScale);
			if (extraOptions)
				Parameters.isPositive("Image rolling window", resultsSettings.imageRollingWindow);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		// If precision filtering then we need the camera bias
		if (fitConfig.getPrecisionThreshold() > 0)
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Precision filtering can use global noise estimate or local background level.\n \nLocal background requires the camera bias:");
			gd.addCheckbox("Local_background", fitConfig.isPrecisionUsingBackground());
			gd.addNumericField("Camera_bias (ADUs)", calibration.bias, 2);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			fitConfig.setPrecisionUsingBackground(gd.getNextBoolean());
			calibration.bias = Math.abs(gd.getNextNumber());
		}

		if (!configureDataFilter(settings, filename, extraOptions))
			return false;

		// Second dialog for solver dependent parameters
		if (!maximaIdentification)
		{
			if (!configureFitSolver(settings, filename, extraOptions))
				return false;
		}

		// Extra parameters are needed for interlaced data
		if (interlacedData)
		{
			gd = new GenericDialog(TITLE);
			gd.addMessage("Interlaced data requires a repeating pattern of frames to process.\n"
					+ "Describe the regular repeat of the data:\n \n" + "Start = The first frame that contains data\n"
					+ "Block = The number of continuous frames containing data\n"
					+ "Skip = The number of continuous frames to ignore before the next data\n \n"
					+ "E.G. 2:9:1 = Data was imaged from frame 2 for 9 frames, 1 frame to ignore, then repeat.");
			gd.addNumericField("Start", optionDataStart, 0);
			gd.addNumericField("Block", optionDataBlock, 0);
			gd.addNumericField("Skip", optionDataSkip, 0);

			gd.showDialog();
			if (gd.wasCanceled())
				return false;

			if (!gd.wasCanceled())
			{
				dataStart = (int) gd.getNextNumber();
				dataBlock = (int) gd.getNextNumber();
				dataSkip = (int) gd.getNextNumber();

				if (dataStart > 0 && dataBlock > 0 && dataSkip > 0)
				{
					// Store options for next time
					optionInterlacedData = true;
					optionDataStart = dataStart;
					optionDataBlock = dataBlock;
					optionDataSkip = dataSkip;
				}
			}
			else
			{
				interlacedData = false;
			}
		}

		boolean result = SettingsManager.saveSettings(settings, filename);
		if (!result)
		{
			IJ.error(TITLE, "Failed to save settings to file " + filename);
		}
		return result;
	}

	/**
	 * Show a dialog to configure the data filter. The updated settings are saved to the settings file. An error
	 * message is shown if the dialog is cancelled or the configuration is invalid.
	 * 
	 * @param settings
	 * @param filename
	 * @param extraOptions
	 *            True if extra configuration options should be allowed
	 * @return True if the configuration succeeded
	 */
	public static boolean configureDataFilter(GlobalSettings settings, String filename, boolean extraOptions)
	{
		FitEngineConfiguration config = settings.getFitEngineConfiguration();

		int numberOfFilters = 1;
		final int n;
		switch (config.getDataFilterType())
		{
			case JURY:
				n = Integer.MAX_VALUE;
				break;

			case DIFFERENCE:
				n = 2;
				break;

			case SINGLE:
			default:
				n = 1;
		}

		String[] filterNames = SettingsManager.getNames((Object[]) DataFilter.values());

		for (int i = 1; i < n; i++)
		{
			int filter = i + 1;
			GenericDialog gd = new GenericDialog(TITLE);
			gd.enableYesNoCancel("Add", "Continue");
			gd.addMessage(String.format(
					"Configure the %s filter.\nClick continue to proceed with the current set of %d.", config
							.getDataFilterType().toString(), i));
			String fieldName = "Spot_filter" + filter;
			if (IJ.isMacro())
				// Use blank default value so bad macro parameters return nothing
				gd.addStringField(fieldName, "");
			else
				gd.addChoice(fieldName, filterNames, filterNames[config.getDataFilter(i).ordinal()]);
			gd.addSlider("Smoothing" + filter, 0, 4.5, config.getSmooth(i));
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			if (gd.wasOKed())
			{
				int filterIndex = -1;
				if (IJ.isMacro())
				{
					String filterName = gd.getNextString();
					for (int j = 0; j < filterNames.length; j++)
						if (filterNames[j].equals(filterName))
						{
							filterIndex = j;
							break;
						}

					if (filterIndex < 0)
						break;
				}
				else
					filterIndex = gd.getNextChoiceIndex();
				config.setDataFilter(filterIndex, Math.abs(gd.getNextNumber()), i);
				numberOfFilters++;
			}
			else
			{
				break;
			}
		}
		config.setNumberOfFilters(numberOfFilters);
		if (filename != null)
			SettingsManager.saveSettings(settings, filename);
		return true;
	}

	/**
	 * Show a dialog to configure the fit solver. The updated settings are saved to the settings file. An error
	 * message is shown if the dialog is cancelled or the configuration is invalid.
	 * 
	 * @param settings
	 * @param filename
	 * @param extraOptions
	 *            True if extra configuration options should be allowed
	 * @return True if the configuration succeeded
	 */
	/**
	 * @param settings
	 * @param filename
	 * @param extraOptions
	 * @return
	 */
	public static boolean configureFitSolver(GlobalSettings settings, String filename, boolean extraOptions)
	{
		return configureFitSolver(settings, filename, extraOptions, false);
	}

	/**
	 * Show a dialog to configure the fit solver. The updated settings are saved to the settings file. An error
	 * message is shown if the dialog is cancelled or the configuration is invalid.
	 * 
	 * @param settings
	 * @param filename
	 * @param extraOptions
	 *            True if extra configuration options should be allowed
	 * @param ignoreCalibration
	 *            True if the calibration should not be configured
	 * @return True if the configuration succeeded
	 */
	public static boolean configureFitSolver(GlobalSettings settings, String filename, boolean extraOptions,
			boolean ignoreCalibration)
	{
		FitConfiguration fitConfig = settings.getFitEngineConfiguration().getFitConfiguration();
		Calibration calibration = settings.getCalibration();

		if (fitConfig.getFitSolver() == FitSolver.MLE)
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Maximum Likelihood Estimation requires additional parameters");
			if (!ignoreCalibration)
			{
				gd.addNumericField("Camera_bias (ADUs)", calibration.bias, 2);
				gd.addCheckbox("Model_camera_noise", fitConfig.isModelCamera());
				gd.addNumericField("Read_noise (ADUs)", calibration.readNoise, 2);
				gd.addNumericField("Gain (ADU/photon)", calibration.gain, 2);
				gd.addCheckbox("EM-CCD", calibration.emCCD);
			}
			String[] searchNames = SettingsManager.getNames((Object[]) MaximumLikelihoodFitter.SearchMethod.values());
			gd.addChoice("Search_method", searchNames, searchNames[fitConfig.getSearchMethod().ordinal()]);
			gd.addStringField("Relative_threshold", "" + fitConfig.getRelativeThreshold());
			gd.addStringField("Absolute_threshold", "" + fitConfig.getAbsoluteThreshold());
			gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);
			gd.addNumericField("Max_function_evaluations", fitConfig.getMaxFunctionEvaluations(), 0);
			if (extraOptions)
				gd.addCheckbox("Gradient_line_minimisation", fitConfig.isGradientLineMinimisation());
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			if (!ignoreCalibration)
			{
				calibration.bias = Math.abs(gd.getNextNumber());
				fitConfig.setModelCamera(gd.getNextBoolean());
				calibration.readNoise = Math.abs(gd.getNextNumber());
				calibration.gain = Math.abs(gd.getNextNumber());
				calibration.emCCD = gd.getNextBoolean();
				fitConfig.setBias(calibration.bias);
				fitConfig.setReadNoise(calibration.readNoise);
				fitConfig.setGain(calibration.gain);
				fitConfig.setEmCCD(calibration.emCCD);
			}
			fitConfig.setSearchMethod(gd.getNextChoiceIndex());
			try
			{
				fitConfig.setRelativeThreshold(Math.abs(Double.parseDouble(gd.getNextString())));
				fitConfig.setAbsoluteThreshold(Math.abs(Double.parseDouble(gd.getNextString())));
			}
			catch (NumberFormatException e)
			{
				fitConfig.setRelativeThreshold(0);
				fitConfig.setAbsoluteThreshold(0);
			}
			fitConfig.setMaxIterations((int) gd.getNextNumber());
			fitConfig.setMaxFunctionEvaluations((int) gd.getNextNumber());
			if (extraOptions)
				fitConfig.setGradientLineMinimisation(gd.getNextBoolean());
			else
				// This option is for the Conjugate Gradient optimiser and makes it less stable
				fitConfig.setGradientLineMinimisation(false);

			if (filename != null)
				SettingsManager.saveSettings(settings, filename);

			try
			{
				Parameters.isAboveZero("Relative threshold", fitConfig.getRelativeThreshold());
				Parameters.isAboveZero("Absolute threshold", fitConfig.getAbsoluteThreshold());
				Parameters.isAboveZero("Max iterations", fitConfig.getMaxIterations());
				Parameters.isAboveZero("Max function evaluations", fitConfig.getMaxFunctionEvaluations());
				fitConfig.getFunctionSolver();
			}
			catch (IllegalArgumentException e)
			{
				IJ.error(TITLE, e.getMessage());
				return false;
			}
		}
		else if (fitConfig.getFitSolver() == FitSolver.LVM || fitConfig.getFitSolver() == FitSolver.LVM_WEIGHTED)
		{
			// Collect options for LVM fitting
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("LVM requires additional parameters");
			String[] criteriaNames = SettingsManager.getNames((Object[]) FitCriteria.values());
			gd.addChoice("Fit_criteria", criteriaNames, criteriaNames[fitConfig.getFitCriteria().ordinal()]);
			gd.addNumericField("Significant_digits", fitConfig.getSignificantDigits(), 0);
			gd.addNumericField("Coord_delta", fitConfig.getDelta(), 4);
			gd.addNumericField("Lambda", fitConfig.getLambda(), 4);
			if (extraOptions)
				gd.addNumericField("Min_iterations", fitConfig.getMinIterations(), 0);
			gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);

			// Extra parameters are needed for the weighted LVM
			if (fitConfig.getFitSolver() == FitSolver.LVM_WEIGHTED && !ignoreCalibration)
			{
				gd.addMessage("Weighted LVM fitting requires a CCD camera noise model");
				gd.addNumericField("Read_noise (ADUs)", calibration.readNoise, 2);
				gd.addNumericField("Camera_bias (ADUs)", calibration.bias, 2);
			}
			gd.showDialog();
			if (gd.wasCanceled())
				return false;

			fitConfig.setFitCriteria(gd.getNextChoiceIndex());

			fitConfig.setSignificantDigits((int) gd.getNextNumber());
			fitConfig.setDelta(gd.getNextNumber());
			fitConfig.setLambda(gd.getNextNumber());
			if (extraOptions)
				fitConfig.setMinIterations((int) gd.getNextNumber());
			fitConfig.setMaxIterations((int) gd.getNextNumber());

			if (fitConfig.getFitSolver() == FitSolver.LVM_WEIGHTED)
			{
				calibration.readNoise = Math.abs(gd.getNextNumber());
				calibration.bias = Math.abs(gd.getNextNumber());
				fitConfig.setNoiseModel(CameraNoiseModel.createNoiseModel(calibration.readNoise, calibration.bias,
						calibration.emCCD));
			}

			if (filename != null)
				SettingsManager.saveSettings(settings, filename);

			try
			{
				Parameters.isAboveZero("Significant digits", fitConfig.getSignificantDigits());
				Parameters.isAboveZero("Delta", fitConfig.getDelta());
				Parameters.isAboveZero("Lambda", fitConfig.getLambda());
				Parameters.isAboveZero("Max iterations", fitConfig.getMaxIterations());
			}
			catch (IllegalArgumentException e)
			{
				IJ.error(TITLE, e.getMessage());
				return false;
			}
		}
		else if (fitConfig.getFitSolver() == FitSolver.LVM_QUASI_NEWTON)
		{
			// No options yet for Apache LVM fitting. Save options for consistency
			if (filename != null)
				SettingsManager.saveSettings(settings, filename);
		}
		return true;
	}

	/**
	 * Add a result output.
	 * <p>
	 * This can be called after {@link #initialiseImage(ImagePlus)} and before {@link #initialiseFitting()} to add to
	 * the configured result outputs.
	 * 
	 * @param peakResults
	 */
	public void addPeakResults(PeakResults peakResults)
	{
		if (results != null)
			results.addOutput(peakResults);
	}

	private void addImageResults(PeakResultsList resultsList)
	{
		if (resultsSettings.getResultsImage() != ResultsImage.NONE)
		{
			IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(
					resultsSettings.getResultsImage(), resultsSettings.weightedImage, resultsSettings.equalisedImage,
					resultsList.getName(), bounds, calibration.nmPerPixel, calibration.gain,
					resultsSettings.imageScale, resultsSettings.precision, ResultsMode.ADD);
			if (extraOptions)
				image.setRollingWindowSize(resultsSettings.imageRollingWindow);
			resultsList.addOutput(image);
		}
	}

	private void addFileResults(PeakResultsList resultsList)
	{
		String filename = null;
		if (resultsSettings.resultsDirectory != null && new File(resultsSettings.resultsDirectory).exists())
		{
			filename = resultsSettings.resultsDirectory + File.separatorChar + source.getName() +
					((resultsSettings.binaryResults) ? ".results.bin" : ".results.xls");
		}
		else if (resultsSettings.resultsFilename != null && resultsSettings.resultsFilename.length() > 0)
		{
			filename = resultsSettings.resultsFilename;
		}
		if (filename != null)
		{
			FilePeakResults r;
			if (resultsSettings.binaryResults)
				r = new BinaryFilePeakResults(filename, resultsSettings.showDeviations);
			else
				r = new FilePeakResults(filename, resultsSettings.showDeviations);
			r.setSortAfterEnd(Prefs.getThreads() > 1);
			r.setPeakIdColumnName("Frame");
			resultsList.addOutput(r);
		}
	}

	private void addMemoryResults(PeakResultsList resultsList, boolean force)
	{
		if (resultsSettings.resultsInMemory || force)
		{
			MemoryPeakResults results = new MemoryPeakResults();
			results.setSortAfterEnd(Prefs.getThreads() > 1);
			resultsList.addOutput(results);
			MemoryPeakResults.addResults(results);
		}
	}

	private void addTableResults(PeakResultsList resultsList)
	{
		if (resultsSettings.getResultsTable() != ResultsTable.NONE)
		{
			String title = null; // imp.getTitle()
			IJTablePeakResults r = new IJTablePeakResults(resultsSettings.showDeviations, title);
			r.setPeakIdColumnName("Frame");
			r.setShowCalibratedValues(resultsSettings.getResultsTable() == ResultsTable.CALIBRATED);
			r.setClearAtStart(simpleFit);
			r.setShowEndFrame(integrateFrames > 1);
			resultsList.addOutput(r);
		}
	}

	private void addDefaultResults(PeakResultsList resultsList)
	{
		if (resultsList.numberOfOutputs() == 0)
		{
			if (logger != null)
				logger.info("No results output configured. Defaulting to memory");
			addMemoryResults(resultsList, true);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		if (source == null)
		{
			IJ.error(TITLE, "No valid image source configured");
			return;
		}
		if (fitMaxima)
		{
			runMaximaFitting();
		}
		else
		{
			// All setup is already done so simply run
			run();

			addSingleFrameOverlay();
		}
	}

	private void addSingleFrameOverlay()
	{
		// If a single frame was processed add the peaks as an overlay if they are in memory
		if (singleFrame > 0 && imp != null)
		{
			MemoryPeakResults results = null;
			for (PeakResults r : this.results.toArray())
				if (r instanceof MemoryPeakResults)
				{
					results = (MemoryPeakResults) r;
					break;
				}
			if (results == null || results.size() == 0)
				return;

			GenericDialog gd = new GenericDialog(TITLE);
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			gd.addMessage("Add the fitted localisations as an overlay?");
			gd.showDialog();
			if (!gd.wasOKed())
				return;

			float[] ox = new float[results.size()];
			float[] oy = new float[results.size()];
			int points = 0;
			for (PeakResult r : results.getResults())
			{
				ox[points] = r.getXPosition();
				oy[points++] = r.getYPosition();
			}

			ImageROIPainter.addRoi(imp, singleFrame, new PointRoi(ox, oy, points));
		}
	}

	/**
	 * Process the image. The current ROI will be used to define the region processed. The noise can be estimated using
	 * the entire frame or the ROI region.
	 * 
	 * @param imp
	 * @param ignoreBoundsForNoise
	 *            If true estimate the noise from the entire frame, otherwise use only the ROI bounds
	 */
	public void run(ImagePlus imp, boolean ignoreBoundsForNoise)
	{
		run(new IJImageSource(imp), getBounds(imp), ignoreBoundsForNoise);
	}

	/**
	 * Process the image
	 * 
	 * @param imageSource
	 * @param bounds
	 * @param ignoreBoundsForNoise
	 */
	public void run(ImageSource imageSource, Rectangle bounds, boolean ignoreBoundsForNoise)
	{
		if (initialise(imageSource, bounds, ignoreBoundsForNoise))
			run();
	}

	/**
	 * Locate the peaks in the configured image source. Results are saved to the configured output.
	 * <p>
	 * This must be called after initialisation with an image source. Note that each call to this method must be
	 * preceded with initialisation to prepare the image and output options.
	 */
	public void run()
	{
		if (source == null)
			return;

		int totalFrames = source.getFrames();

		ImageStack stack = null;
		if (showProcessedFrames)
			stack = new ImageStack(bounds.width, bounds.height);

		// Do not crop the region from the source if the bounds match the source dimensions
		Rectangle cropBounds = (bounds.x == 0 && bounds.y == 0 && bounds.width == source.getWidth() && bounds.height == source
				.getHeight()) ? null : bounds;

		// Use the FitEngine to allow multi-threading.
		FitEngine engine = createFitEngine(FastMath.min(totalFrames, Prefs.getThreads()));

		final int step = Utils.getProgressInterval(totalFrames);

		boolean shutdown = false;
		int slice = 0;
		while (!shutdown)
		{
			// Noise can optionally be estimated from the entire frame
			float[] data = (ignoreBoundsForNoise) ? source.next() : source.next(cropBounds);
			if (data == null)
				break;

			if (++slice % step == 0)
			{
				if (Utils.showStatus("Slice: " + slice + " / " + totalFrames))
					IJ.showProgress(slice, totalFrames);
			}

			float noise = Float.NaN;
			if (ignoreBoundsForNoise)
			{
				noise = FitWorker.estimateNoise(data, source.getWidth(), source.getHeight(), config.getNoiseMethod());

				// Crop the data to the region
				data = ImageConverter.getData(data, source.getWidth(), source.getHeight(), bounds, null);
			}

			if (showProcessedFrames)
			{
				stack.addSlice(
						String.format("Frame %d - %d", source.getStartFrameNumber(), source.getEndFrameNumber()), data);
			}

			// Get the frame number from the source to allow for interlaced and aggregated data
			engine.run(createJob(source.getStartFrameNumber(), source.getEndFrameNumber(), data, bounds, noise));

			if (escapePressed())
				shutdown = true;
		}

		engine.end(shutdown);
		time = engine.getTime();

		if (showProcessedFrames)
			Utils.display("Processed frames", stack);

		showResults();

		source.close();
	}

	/**
	 * Check if the frame should be ignored (relevant when using interlaced data)
	 * 
	 * @param frame
	 * @return True if the frame should be ignored
	 */
	private boolean ignoreFrame(int frame)
	{
		if (!interlacedData)
			return false;

		// Check if the frame is allowed:
		//    Start
		//      |
		// |----|Block|Skip|Block|Skip|Block|Skip
		// Note the source data is still read so that the source is incremented.
		if (frame < dataStart)
		{
			//System.out.printf("Skipping %d\n", frame);
			return true;
		}
		int frameInBlock = (frame - dataStart) % (dataBlock + dataSkip);
		if (frameInBlock >= dataBlock)
		{
			//System.out.printf("Skipping %d\n", frame);
			return true;
		}
		return false;
	}

	private FitJob createJob(int slice, int endFrame, float[] data, Rectangle bounds2, float noise)
	{
		FitParameters fitParams = null;
		if (slice != endFrame)
		{
			fitParams = new FitParameters();
			fitParams.endT = endFrame;
		}

		if (maximaIdentification)
		{
			if (fitParams == null)
				fitParams = new FitParameters();
			fitParams.fitTask = FitTask.MAXIMA_IDENITIFICATION;
			fitParams.noise = noise;
		}
		else if (!Float.isNaN(noise))
		{
			if (fitParams == null)
				fitParams = new FitParameters();
			fitParams.fitTask = FitTask.PSF_FITTING;
			fitParams.noise = noise;
		}

		if (fitParams != null)
			return new ParameterisedFitJob(fitParams, slice, data, bounds);
		else
			return new FitJob(slice, data, bounds);
	}

	private boolean escapePressed()
	{
		if (IJ.escapePressed())
		{
			IJ.log(TITLE + " stopping ...");
			IJ.beep();
			return true;
		}
		return false;
	}

	/**
	 * Creates a fitting engine using the current configuration.
	 * 
	 * @return The fitting engine
	 */
	public FitEngine createFitEngine()
	{
		return createFitEngine(Prefs.getThreads());
	}

	/**
	 * Creates a fitting engine using the current configuration.
	 * 
	 * @param numberOfThreads
	 * @return The fitting engine
	 */
	public FitEngine createFitEngine(int numberOfThreads)
	{
		// Use a blocking queue to enable progress tracking on the IJ progress bar.
		return createFitEngine(numberOfThreads, FitQueue.BLOCKING);
	}

	/**
	 * Creates a fitting engine using the current configuration.
	 * 
	 * @param numberOfThreads
	 * @param queue
	 * @return The fiting engine
	 */
	public FitEngine createFitEngine(int numberOfThreads, FitQueue queue)
	{
		PeakResults r = results;
		if (results.numberOfOutputs() == 1)
			// Reduce to single object for speed
			r = results.toArray()[0];

		// Update the configuration
		updateFitConfiguration(config);

		FitEngine engine = new FitEngine(config, r, numberOfThreads, queue);

		// Write settings out to the IJ log
		if (resultsSettings.logProgress)
		{
			IJ.log("-=-=-=-");
			IJ.log("Peak Fit");
			IJ.log("-=-=-=-");
			Utils.log("Initial Peak SD = %s,%s", Utils.rounded(fitConfig.getInitialPeakStdDev0()),
					Utils.rounded(fitConfig.getInitialPeakStdDev1()));
			SpotFilter spotFilter = engine.getSpotFilter();
			IJ.log("Spot Filter = " + spotFilter.getDescription());
			int w = 2 * engine.getFitting() + 1;
			Utils.log("Fit window = %d x %d", w, w);
			IJ.log("Coordinate shift = " + Utils.rounded(config.getFitConfiguration().getCoordinateShift()));
			IJ.log("Signal strength = " + Utils.rounded(fitConfig.getSignalStrength()));
			if (extraOptions)
				IJ.log("Noise = " + Utils.rounded(fitConfig.getNoise()));
			IJ.log("Width factor = " + Utils.rounded(fitConfig.getWidthFactor()));
			IJ.log("-=-=-=-");
		}

		return engine;
	}

	/**
	 * Updates the configuration for peak fitting. Configures the calculation of residuals, logging and peak validation.
	 * 
	 * @return
	 */
	private void updateFitConfiguration(FitEngineConfiguration config)
	{
		FitConfiguration fitConfig = config.getFitConfiguration();

		// Adjust the settings that are relevant within the fitting configuration. 
		fitConfig.setComputeResiduals(config.getResidualsThreshold() < 1);
		logger = (resultsSettings.logProgress) ? new IJLogger() : null;
		fitConfig.setLog(logger);

		fitConfig.setComputeDeviations(resultsSettings.showDeviations);

		// Setup peak filtering
		fitConfig.setFitValidation(true);

		// Add the calibration for precision filtering
		fitConfig.setNmPerPixel(calibration.nmPerPixel);
		fitConfig.setGain(calibration.gain);
		fitConfig.setBias(calibration.bias);
		fitConfig.setEmCCD(calibration.emCCD);
	}

	/**
	 * Load the selected results from memory. All multiple frame results are added directly to the results. All single
	 * frame
	 * results are added to a list of candidate maxima per frame and fitted using the configured parameters.
	 */
	private void runMaximaFitting()
	{
		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			log("No results for maxima fitting");
			return;
		}
		// No check for imageSource since this has been check in the calling method

		int totalFrames = source.getFrames();

		// Store the indices of each new time-frame
		results.sort();
		List<PeakResult> candidateMaxima = results.getResults();

		// Use the FitEngine to allow multi-threading.
		FitEngine engine = createFitEngine(FastMath.min(totalFrames, Prefs.getThreads()));

		final int step = Utils.getProgressInterval(totalFrames);

		boolean shutdown = false;
		int slice = candidateMaxima.get(0).peak;
		ArrayList<PeakResult> sliceCandidates = new ArrayList<PeakResult>();
		Iterator<PeakResult> iter = candidateMaxima.iterator();
		while (iter.hasNext())
		{
			PeakResult r = iter.next();
			if (slice != r.peak)
			{
				if (escapePressed())
				{
					shutdown = true;
					break;
				}
				if (slice % step == 0)
				{
					if (Utils.showStatus("Slice: " + slice + " / " + totalFrames))
						IJ.showProgress(slice, totalFrames);
				}

				// Process results
				if (!processResults(engine, sliceCandidates, slice))
					break;

				sliceCandidates.clear();
			}
			slice = r.peak;
			sliceCandidates.add(r);
		}

		// Process final results
		if (!shutdown)
			processResults(engine, sliceCandidates, slice);

		engine.end(shutdown);
		time = engine.getTime();

		showResults();

		source.close();
	}

	private boolean processResults(FitEngine engine, ArrayList<PeakResult> sliceCandidates, int slice)
	{
		// Process results
		int[] maxIndices = new int[sliceCandidates.size()];
		int count = 0;
		ArrayList<PeakResult> processedResults = new ArrayList<PeakResult>(sliceCandidates.size());
		for (PeakResult result : sliceCandidates)
		{
			// Add ExtendedPeakResults to the results if they span multiple frames (they are the result of previous fitting).
			if (result instanceof ExtendedPeakResult && result.peak != result.getEndFrame())
			{
				processedResults.add(result);
			}
			else
			{
				// Fit single frame results.
				maxIndices[count++] = result.origX + bounds.width * result.origY;
			}
		}

		if (!processedResults.isEmpty())
			this.results.addAll(processedResults);

		if (count != 0)
		{
			float[] data = source.get(slice);
			if (data == null)
				return false;

			FitParameters fitParams = new FitParameters();
			fitParams.maxIndices = Arrays.copyOf(maxIndices, count);
			FitJob job = new ParameterisedFitJob(fitParams, slice, data, bounds);

			engine.run(job);
		}

		return true;
	}

	/**
	 * @return The total fitting time in nanoseconds
	 */
	public long getTime()
	{
		return time;
	}

	/**
	 * @return The total number of localisations
	 */
	public int getSize()
	{
		// If only one output in the list it was extracted for the FitEngine to prevent 
		// passing data through all methods in the PeakResultsList
		return (results.numberOfOutputs() == 1) ? results.toArray()[0].size() : results.size();
	}

	/**
	 * @return If true, do not output any log messages (e.g. the final time and count)
	 */
	public boolean isSilent()
	{
		return silent;
	}

	/**
	 * @param silent
	 *            If true, do not output any log messages (e.g. the final time and count)
	 */
	public void setSilent(boolean silent)
	{
		this.silent = silent;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseClicked(java.awt.event.MouseEvent)
	 */
	public void mouseClicked(MouseEvent e)
	{
		if (e.getClickCount() > 1) // Double-click
		{
			if (e.getSource() == textConfigFile)
			{
				String newFilename = Utils.getFilename("Config_File", textConfigFile.getText());
				if (newFilename != null)
				{
					textConfigFile.setText(newFilename);
				}
			}
			else if (e.getSource() == textResultsDirectory)
			{
				String directory = Utils.getDirectory("Results_dir", textResultsDirectory.getText());
				if (directory != null)
					textResultsDirectory.setText(directory);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mousePressed(java.awt.event.MouseEvent)
	 */
	public void mousePressed(MouseEvent e)
	{

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseReleased(java.awt.event.MouseEvent)
	 */
	public void mouseReleased(MouseEvent e)
	{

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseEntered(java.awt.event.MouseEvent)
	 */
	public void mouseEntered(MouseEvent e)
	{

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseExited(java.awt.event.MouseEvent)
	 */
	public void mouseExited(MouseEvent e)
	{

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.TextListener#textValueChanged(java.awt.event.TextEvent)
	 */
	public void textValueChanged(TextEvent e)
	{
		if (e.getSource() == textConfigFile)
		{
			refreshSettings(textConfigFile.getText());
		}
	}

	private void refreshSettings(String newFilename)
	{
		if (newFilename != null && new File(newFilename).exists())
		{
			GenericDialog gd = new GenericDialog(TITLE);
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			gd.addMessage("Reload settings from file");
			gd.showDialog();
			if (gd.wasOKed())
			{
				// Reload the settings and update the GUI
				GlobalSettings settings = SettingsManager.unsafeLoadSettings(newFilename);
				if (settings == null)
					return;

				Calibration calibration = settings.getCalibration();
				refreshSettings(calibration);

				FitEngineConfiguration config = settings.getFitEngineConfiguration();
				refreshSettings(config);

				ResultsSettings resultsSettings = settings.getResultsSettings();
				refreshSettings(resultsSettings);
			}
		}
	}

	private void refreshSettings(Calibration calibration)
	{
		this.calibration = calibration;

		textNmPerPixel.setText("" + calibration.nmPerPixel);
		textGain.setText("" + calibration.gain);
		textEMCCD.setState(calibration.emCCD);
		textExposure.setText("" + calibration.exposureTime);
	}

	private void refreshSettings(FitEngineConfiguration config)
	{
		// Set the configuration
		this.config = config;
		this.fitConfig = config.getFitConfiguration();

		textInitialPeakStdDev0.setText("" + fitConfig.getInitialPeakStdDev0());
		if (!maximaIdentification)
		{
			textInitialPeakStdDev1.setText("" + fitConfig.getInitialPeakStdDev1());
			textInitialAngleD.setText("" + fitConfig.getInitialAngle());
		}
		textDataFilterType.select(config.getDataFilterType().ordinal());
		textDataFilter.select(config.getDataFilter(0).ordinal());
		textSmooth.setText("" + config.getSmooth(0));
		textSearch.setText("" + config.getSearch());
		textBorder.setText("" + config.getBorder());
		textFitting.setText("" + config.getFitting());
		if (!maximaIdentification)
		{
			textFitSolver.select(fitConfig.getFitSolver().ordinal());
			textFitFunction.select(fitConfig.getFitFunction().ordinal());
			if (extraOptions)
				textFitBackground.setState(fitConfig.isBackgroundFitting());
			textFailuresLimit.setText("" + config.getFailuresLimit());
			textIncludeNeighbours.setState(config.isIncludeNeighbours());
			textNeighbourHeightThreshold.setText("" + config.getNeighbourHeightThreshold());
			textResidualsThreshold.setText("" + config.getResidualsThreshold());
			textDuplicateDistance.setText("" + fitConfig.getDuplicateDistance());
			textCoordinateShiftFactor.setText("" + fitConfig.getCoordinateShiftFactor());
			textSignalStrength.setText("" + fitConfig.getSignalStrength());
			textMinPhotons.setText("" + fitConfig.getMinPhotons());
			textWidthFactor.setText("" + fitConfig.getWidthFactor());
			textPrecisionThreshold.setText("" + fitConfig.getPrecisionThreshold());
			if (extraOptions)
			{
				textNoise.setText("" + fitConfig.getNoise());
				textNoiseMethod.select(config.getNoiseMethod().ordinal());
			}
		}
	}

	private void refreshSettings(ResultsSettings resultsSettings)
	{
		this.resultsSettings = resultsSettings;
		textLogProgress.setState(resultsSettings.logProgress);
		if (!maximaIdentification)
			textShowDeviations.setState(resultsSettings.showDeviations);
		textResultsTable.select(resultsSettings.getResultsTable().ordinal());
		textResultsImage.select(resultsSettings.getResultsImage().ordinal());
		textWeightedImage.setState(resultsSettings.weightedImage);
		textEqualisedImage.setState(resultsSettings.equalisedImage);
		textPrecision.setText("" + resultsSettings.precision);
		textImageScale.setText("" + resultsSettings.imageScale);
		if (extraOptions)
			textImageRollingWindow.setText("" + resultsSettings.imageRollingWindow);
		textResultsDirectory.setText("" + resultsSettings.resultsDirectory);
		textBinaryResults.setState(resultsSettings.binaryResults);
		textResultsInMemory.setState(resultsSettings.resultsInMemory);
	}
}
