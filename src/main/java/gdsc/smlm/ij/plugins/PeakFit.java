package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.SystemColor;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JFileChooser;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.IJLogger;
import gdsc.core.ij.SeriesOpener;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Logger;
import gdsc.core.utils.BitFlags;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationConfig.Calibration;
import gdsc.smlm.data.config.CalibrationConfig.CameraType;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitConfig.DataFilterMethod;
import gdsc.smlm.data.config.FitConfig.FitEngineSettings;
import gdsc.smlm.data.config.FitConfig.FitSettings;
import gdsc.smlm.data.config.FitConfig.FitSolver;
import gdsc.smlm.data.config.FitConfig.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitConfigHelper;
import gdsc.smlm.data.config.PSFConfig.PSF;
import gdsc.smlm.data.config.PSFConfig.PSFParameter;
import gdsc.smlm.data.config.PSFConfig.PSFType;
import gdsc.smlm.data.config.PSFConfigHelper;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.ResultsConfig.ResultsFileSettings;
import gdsc.smlm.data.config.ResultsConfig.ResultsImageSettings;
import gdsc.smlm.data.config.ResultsConfig.ResultsImageType;
import gdsc.smlm.data.config.ResultsConfig.ResultsSettings;
import gdsc.smlm.data.config.ResultsConfig.ResultsTableSettings;
import gdsc.smlm.data.config.ResultsConfigHelper;
import gdsc.smlm.data.config.UnitConfig.AngleUnit;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngine;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitJob;
import gdsc.smlm.engine.FitParameters;
import gdsc.smlm.engine.FitParameters.FitTask;
import gdsc.smlm.engine.FitQueue;
import gdsc.smlm.engine.FitWorker;
import gdsc.smlm.engine.ParameterisedFitJob;
import gdsc.smlm.filters.SpotFilter;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.SeriesImageSource;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.PSFCalculatorSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.AggregatedImageSource;
import gdsc.smlm.results.Counter;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.FilePeakResults;
import gdsc.smlm.results.FrameCounter;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.InterlacedImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.filter.DirectFilter;
import gdsc.smlm.results.filter.Filter;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionListener;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;

//import ij.io.OpenDialog;

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

/**
 * Fits local maxima using a 2D Gaussian. Process each frame until a successive number of fits
 * fail to meet the fit criteria.
 */
public class PeakFit implements PlugInFilter, ItemListener
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
	private long time, runTime;
	private FitEngineConfiguration config = null;
	private FitConfiguration fitConfig;
	private ResultsSettings.Builder resultsSettings;
	private boolean silent = false;

	private static int numberOfThreads = 1;

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
	// Testing has shown that we should use ~85% of the total number of cores the system has available.
	// Extra cores then do not make a difference.
	private static double fractionOfThreads = 0.85;

	private static String inputOption = "";
	private static boolean showTable = true;
	private static boolean showImage = true;
	private static PSFCalculatorSettings calculatorSettings = new PSFCalculatorSettings();

	// All the fields that will be updated when reloading the configuration file
	private TextField textNmPerPixel;
	private TextField textGain;
	private Checkbox textEMCCD;
	private TextField textExposure;
	private Choice textPSF;
	private Choice textDataFilterType;
	private Choice textDataFilterMethod;
	private TextField textSmooth;
	private TextField textSearch;
	private TextField textBorder;
	private TextField textFitting;
	private Choice textFitSolver;
	private Checkbox textFitBackground;
	private TextField textFailuresLimit;
	private Checkbox textIncludeNeighbours;
	private TextField textNeighbourHeightThreshold;
	private TextField textResidualsThreshold;
	private TextField textDuplicateDistance;
	private Scrollbar sliderCoordinateShiftFactor;
	private TextField textCoordinateShiftFactor;
	private TextField textSignalStrength;
	private TextField textMinPhotons;
	private TextField textPrecisionThreshold;
	private Checkbox textSmartFilter;
	private Checkbox textDisableSimpleFilter;
	private TextField textNoise;
	private Choice textNoiseMethod;
	private TextField textMinWidthFactor;
	private TextField textWidthFactor;
	private Checkbox textLogProgress;
	private Checkbox textShowDeviations;
	private Checkbox textResultsTable;
	private Choice textResultsImage;
	private TextField textResultsDirectory;
	private Choice textFileFormat;
	private Checkbox textResultsInMemory;

	public PeakFit()
	{
		init(null, null);
	}

	public PeakFit(FitEngineConfiguration config)
	{
		init(config, null);
	}

	public PeakFit(FitEngineConfiguration config, ResultsSettings resultsSettings)
	{
		init(config, resultsSettings);
	}

	private void init(FitEngineConfiguration config, ResultsSettings resultsSettings)
	{
		this.config = (config != null) ? config : new FitEngineConfiguration();
		this.resultsSettings = (resultsSettings != null) ? ResultsSettings.newBuilder(resultsSettings)
				: ResultsSettings.newBuilder();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

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

			// Check for single frame
			singleFrame = results.getFirstFrame();
			final FrameCounter counter = new FrameCounter(singleFrame);
			results.forEach(new PeakResultProcedureX()
			{
				public boolean execute(PeakResult peakResult)
				{
					// The counter will return true (stop execution) if a new frame
					return counter.advance(peakResult.getFrame());
				}
			});
			if (counter.currentFrame() != counter.previousFrame())
				singleFrame = 0;

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
			SeriesOpener series;
			if (extraOptions)
				series = new SeriesOpener(inputDirectory, true, numberOfThreads);
			else
				series = new SeriesOpener(inputDirectory);
			if (series.getNumberOfImages() == 0)
			{
				IJ.error(TITLE, "No images in the selected directory:\n" + inputDirectory);
				return DONE;
			}

			SeriesImageSource seriesImageSource = new SeriesImageSource(getName(series.getImageList()), series);
			seriesImageSource.setLogProgress(true);
			if (extraOptions)
			{
				numberOfThreads = Math.max(1, series.getNumberOfThreads());
				seriesImageSource.setNumberOfThreads(numberOfThreads);
			}
			imageSource = seriesImageSource;

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
				{
					// TODO - Find out why XStream does not serialise the inherited name field in 
					// all ImageSource subclasses. Thjs means that some of the ImageSource details
					// are missing
					if (!Utils.isNullOrEmpty(tmpImageSource.getName()))
						message += " of: \n \n" + tmpImageSource.getName();
					message += " \n \nFit the parent?";
				}
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
	 * 		True if the image was valid and the initialisation was successful
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

		//Calibration cal = calibration.clone();
		// Account for the frame integration
		// TODO - Should we change this so that if integrate frames is used then the data 
		// are converted to ExtendedPeakResult with a start and end frame
		//cal.exposureTime *= integrateFrames;
		//if (interlacedData)
		//{
		//	cal.exposureTime *= ((double)dataBlock / (dataBlock + dataSkip));
		//}

		results.setCalibration(fitConfig.getCalibration());
		results.setPSF(fitConfig.getPSF());
		results.setConfiguration(XmlUtils.toXML(config));

		addTableResults(results);
		ResultsManager.addImageResults(results, resultsSettings.getResultsImageSettings(), bounds,
				(extraOptions) ? ResultsManager.FLAG_EXTRA_OPTIONS : 0);
		addFileResults(results);
		addMemoryResults(results, false);
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
		String name = FitConfigHelper.getName(solver);
		if (solver == FitSolver.MLE)
			name += " " + FitConfigHelper.getName(fitConfig.getSearchMethod());
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
						IJ.showStatus("Sorting " + r.size() + " results ...");
					break;
				}

			results.end();

			String textTime = Utils.timeToString(time / 1000000.0);
			String textRunTime = Utils.timeToString(runTime / 1000000.0);

			int size = getSize();
			String message = String.format("%s. Fitting Time = %s. Run time = %s", Utils.pleural(size, "localisation"),
					textTime, textRunTime);
			if (resultsSettings.getLogProgress())
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

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Select identified maxima for fitting");

		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);

		return true;
	}

	private static PSFType[] _PSFTypeValues;

	public static PSFType[] getPSFTypeValues()
	{
		if (_PSFTypeValues == null)
			initPSFType();
		return _PSFTypeValues;
	}

	private static String[] _PSFTypeNames;

	public static String[] getPSFTypeNames()
	{
		if (_PSFTypeNames == null)
			initPSFType();
		return _PSFTypeNames;
	}

	private static void initPSFType()
	{
		//@formatter:off
		EnumSet<PSFType> d = EnumSet.of(
				PSFType.ONE_AXIS_GAUSSIAN_2D, 
				PSFType.TWO_AXIS_GAUSSIAN_2D, 
				PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
		//@formatter:on
		_PSFTypeValues = d.toArray(new PSFType[d.size()]);
		_PSFTypeNames = new String[_PSFTypeValues.length];
		for (int i = 0; i < _PSFTypeValues.length; i++)
		{
			_PSFTypeNames[i] = PSFConfigHelper.getName(_PSFTypeValues[i]);
		}
	}

	@SuppressWarnings("unchecked")
	private int showDialog(ImagePlus imp)
	{
		// Executing as an ImageJ plugin.

		// Load the settings
		resultsSettings = SettingsManager.readResultsSettings().toBuilder();
		// Settings are within the FitEngineSettings
		config = new FitEngineConfiguration(SettingsManager.readFitEngineSettings());
		fitConfig = config.getFitConfiguration();

		if (simpleFit)
		{
			return showSimpleDialog();
		}

		boolean isCrop = (bounds != null && imp != null &&
				(bounds.width < imp.getWidth() || bounds.height < imp.getHeight()));

		if (!extraOptions)
		{
			integrateFrames = 1;
			resultsSettings.getResultsImageSettingsBuilder().setRollingWindowSize(0);
			fitConfig.setBackgroundFitting(true);
			fitConfig.setNoise(0);
			config.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);
			showProcessedFrames = false;
		}

		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage((maximaIdentification) ? "Identify candidate maxima" : "Fit 2D Gaussian to identified maxima");

		String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);

		// TODO - change this to support camera type
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();
		gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
		gd.addNumericField("Gain (ADU/photon)", calibration.getGain(), 2);
		gd.addCheckbox("EM-CCD", calibration.isEMCCD());
		gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);

		if (isCrop)
			gd.addCheckbox("Ignore_bounds_for_noise", optionIgnoreBoundsForNoise);

		addPSFOptions(gd, fitConfig);

		gd.addChoice("Spot_filter_type", SettingsManager.getDataFilterTypeNames(),
				config.getDataFilterType().ordinal());
		gd.addChoice("Spot_filter", SettingsManager.getDataFilterMethodNames(),
				config.getDataFilterMethod(0).ordinal());
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
			gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());
			if (extraOptions)
				gd.addCheckbox("Fit_background", fitConfig.isBackgroundFitting());

			// Parameters specific to each Fit solver are collected in a second dialog 

			gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
			gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
			gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
			gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

			gd.addSlider("Duplicate_distance", 0, 1.5, config.getDuplicateDistance());

			gd.addMessage("--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");
			discardLabel = gd.getMessage();

			gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
			gd.addCheckbox("Disable_simple_filter", fitConfig.isDisableSimpleFilter());
			gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
			sliderCoordinateShiftFactor = gd.getLastScrollbar();
			gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
			gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
			if (extraOptions)
			{
				gd.addNumericField("Noise", fitConfig.getNoise(), 2);
				gd.addChoice("Noise_method", SettingsManager.getNoiseEstimatorMethodNames(),
						config.getNoiseMethod().ordinal());
			}
			gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
			gd.addSlider("Width_factor", 1.01, 5, fitConfig.getMaxWidthFactor());
			gd.addNumericField("Precision", fitConfig.getPrecisionThreshold(), 2);
		}

		gd.addMessage("--- Results ---");
		gd.addCheckbox("Log_progress", resultsSettings.getLogProgress());
		if (!maximaIdentification)
			gd.addCheckbox("Show_deviations", resultsSettings.getShowDeviations());
		ResultsManager.addTableResultsOptions(gd, resultsSettings);
		ResultsManager.addImageResultsOptions(gd, resultsSettings,
				(extraOptions) ? ResultsManager.FLAG_EXTRA_OPTIONS : 0);
		if (extraOptions)
			gd.addCheckbox("Show_processed_frames", optionShowProcessedFrames);
		ResultsManager.addFileResultsOptions(gd, resultsSettings, ResultsManager.FLAG_RESULTS_DIRECTORY);
		ResultsManager.addInMemoryResultsOptions(gd, resultsSettings);

		if (extraOptions)
		{
			gd.addMessage("--- Misc ---");
			gd.addSlider("Fraction_of_threads", 0.1, 1, fractionOfThreads);
		}

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
		if (Utils.isShowGenericDialog())
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

			textNmPerPixel = numerics.get(n++);
			textGain = numerics.get(n++);
			textEMCCD = checkboxes.get(b++);
			textExposure = numerics.get(n++);
			if (isCrop)
				b++;
			textPSF = choices.get(ch++);
			textDataFilterType = choices.get(ch++);
			textDataFilterMethod = choices.get(ch++);
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
				if (extraOptions)
					textFitBackground = checkboxes.get(b++);
				textFailuresLimit = numerics.get(n++);
				textIncludeNeighbours = checkboxes.get(b++);
				textNeighbourHeightThreshold = numerics.get(n++);
				textResidualsThreshold = numerics.get(n++);
				textDuplicateDistance = numerics.get(n++);
				textSmartFilter = checkboxes.get(b++);
				textDisableSimpleFilter = checkboxes.get(b++);
				textCoordinateShiftFactor = numerics.get(n++);
				textSignalStrength = numerics.get(n++);
				textMinPhotons = numerics.get(n++);
				if (extraOptions)
				{
					textNoise = numerics.get(n++);
					textNoiseMethod = choices.get(ch++);
				}
				textMinWidthFactor = numerics.get(n++);
				textWidthFactor = numerics.get(n++);
				textPrecisionThreshold = numerics.get(n++);

				updateFilterInput();
				textSmartFilter.addItemListener(this);
				textDisableSimpleFilter.addItemListener(this);
			}
			textLogProgress = checkboxes.get(b++);
			if (!maximaIdentification)
				textShowDeviations = checkboxes.get(b++);
			textResultsTable = checkboxes.get(b++);
			textResultsImage = choices.get(ch++);
			if (extraOptions)
			{
				b++; // Skip over show processed frames option
			}
			textResultsDirectory = texts.get(t++);

			textFileFormat = choices.get(ch++);
			textResultsInMemory = checkboxes.get(b++);
		}

		gd.showDialog();

		if (gd.wasCanceled() || !readDialog(gd, isCrop))
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
		if (resultsSettings.getLogProgress() && source.getFrames() > 1)
		{
			ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
			egd.addMessage("Warning: Log progress on multiple-frame image will be slow");
			egd.addCheckbox("Log_progress", true);
			egd.showDialog();
			if (egd.wasCanceled())
				return DONE;
			if (!egd.getNextBoolean())
			{
				resultsSettings.setLogProgress(false);
				SettingsManager.writeSettings(resultsSettings.build());
			}
		}

		// Return the plugin flags (without the DOES_STACKS flag).
		// The call to run(ImageProcessor) will process the image in 'this.imp' so we only want a 
		// single call to be made.
		return plugin_flags;
	}

	public static void addPSFOptions(final ExtendedGenericDialog gd, final FitConfiguration fitConfig)
	{
		gd.addChoice("PSF", getPSFTypeNames(), PSFConfigHelper.getName(fitConfig.getPSFType()),
				new OptionListener<Choice>()
				{
					public void collectOptions(Choice field)
					{
						fitConfig.setPSFType(PeakFit.getPSFTypeValues()[field.getSelectedIndex()]);
						collectOptions();
					}

					public void collectOptions()
					{
						PSFType psfType = fitConfig.getPSFType();
						ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
						PSF psf = fitConfig.getPSF();
						for (PSFParameter p : psf.getParameterList())
							egd.addNumericField(p.getName(), p.getValue(), 3);
						if (psfType == PSFType.ONE_AXIS_GAUSSIAN_2D)
							egd.addCheckbox("Fixed", fitConfig.isFixedPSF());
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return;
						PSF.Builder b = psf.toBuilder();
						int n = b.getParameterCount();
						for (int i = 0; i < n; i++)
							b.getParameterBuilder(i).setValue(egd.getNextNumber());
						fitConfig.setPSF(b.build());
						if (psfType == PSFType.ONE_AXIS_GAUSSIAN_2D)
							fitConfig.setFixedPSF(egd.getNextBoolean());
					}
				});
	}

	private void log(String format, Object... args)
	{
		if (!silent)
			Utils.log(format, args);
	}

	private int showSimpleDialog()
	{
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();
		boolean requireCalibration = requireCalibration(calibration);
		if (requireCalibration)
		{
			if (!showCalibrationWizard(calibration, true))
				return DONE;
		}

		// Present dialog with simple output options: Image, Table
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
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
			if (!showCalibrationWizard(calibration, false))
				return DONE;
		}

		// Restore fitting to default settings but maintain the calibrated width
		final double sd = fitConfig.getInitialXSD();
		config = new FitEngineConfiguration();
		fitConfig = config.getFitConfiguration();
		fitConfig.setInitialPeakStdDev(sd);
		// Allow to move 1 SD
		fitConfig.setCoordinateShiftFactor(1);
		resultsSettings = ResultsSettings.newBuilder();

		// Do simple results output. We only need set non-default values.
		resultsSettings.getResultsInMemorySettingsBuilder().setInMemory(true);
		if (showTable)
		{
			ResultsTableSettings.Builder tableSettings = resultsSettings.getResultsTableSettingsBuilder();
			tableSettings.setShowTable(true);
		}
		if (showImage)
		{
			ResultsImageSettings.Builder imageSettings = resultsSettings.getResultsImageSettingsBuilder();
			imageSettings.setImageType(ResultsImageType.DRAW_INTENSITY);
			imageSettings.setScale(Math.ceil(1024 / (FastMath.max(bounds.width, bounds.height))));
			imageSettings.setWeighted(true);
			imageSettings.setEqualised(true);
		}

		// Log the settings we care about:
		IJ.log("-=-=-=-");
		IJ.log("Peak Fit");
		IJ.log("-=-=-=-");
		Utils.log("Pixel pitch = %s", Utils.rounded(calibration.getNmPerPixel(), 4));
		Utils.log("Exposure Time = %s", Utils.rounded(calibration.getExposureTime(), 4));
		Utils.log("Gain = %s", Utils.rounded(calibration.getGain(), 4));
		Utils.log("PSF width = %s", Utils.rounded(fitConfig.getInitialXSD(), 4));

		// Save
		saveFitEngineSettings();
		SettingsManager.writeSettings(resultsSettings.build());

		return FLAGS;
	}

	private boolean saveFitEngineSettings()
	{
		return saveFitEngineSettings(config);
	}

	private static boolean saveFitEngineSettings(FitEngineConfiguration config)
	{
		FitEngineSettings settings = config.getFitEngineSettings();
		boolean ok = SettingsManager.writeSettings(settings);
		// Write calibration separately as some plugins may just want to load that		
		if (settings.hasFitSettings() && settings.getFitSettings().hasCalibration())
		{
			Calibration c = settings.getFitSettings().getCalibration();
			ok &= SettingsManager.writeSettings(c);
		}
		return ok;
	}

	/**
	 * Check if additional calibration information is required
	 * <p>
	 * Check the calibration is valid for fitting.
	 *
	 * @param calibration
	 *            the calibration
	 * @return True if additional calibration information is required, false if the system is calibrated
	 */
	private boolean requireCalibration(CalibrationWriter calibration)
	{
		// Check if the calibration contains: Pixel pitch, Gain (can be 1), Exposure time
		if (!calibration.hasNmPerPixel())
			return true;
		if (!calibration.hasGain())
			return true;
		if (!calibration.hasExposureTime())
			return true;

		// Check for a PSF width
		if (fitConfig.getInitialXSD() <= 0)
			return true;

		return false;
	}

	private boolean showCalibrationWizard(CalibrationWriter calibration, boolean showIntroduction)
	{
		if (showIntroduction)
		{
			ExtendedGenericDialog gd = newWizardDialog("No configuration file could be loaded.",
					"Please follow the configuration wizard to calibrate.");
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
		}

		//Calibration defaultCalibration = new Calibration();
		//if (calibration.nmPerPixel <= 0 || calibration.nmPerPixel == defaultCalibration.nmPerPixel)
		if (!getPixelPitch(calibration))
			return false;
		//if (calibration.gain <= 0 || calibration.gain == defaultCalibration.gain)
		if (!getGain(calibration))
			return false;
		//if (calibration.exposureTime <= 0 || calibration.exposureTime == defaultCalibration.exposureTime)
		if (!getExposureTime(calibration))
			return false;
		// Check for a PSF width other than the default
		//if (fitConfig.getInitialPeakWidth0() == new FitConfiguration().getInitialPeakWidth0())
		if (!getPeakWidth(calibration))
			return false;

		// Check parameters
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.getNmPerPixel());
			Parameters.isAboveZero("Gain", calibration.getGain());
			Parameters.isAboveZero("Exposure time", calibration.getExposureTime());
			Parameters.isAboveZero("Initial SD", fitConfig.getInitialXSD());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private ExtendedGenericDialog newWizardDialog(String... messages)
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		final String header = "-=-";
		gd.addMessage(header + " " + TITLE + " Configuration Wizard " + header);
		for (String message : messages)
			gd.addMessage(TextUtils.wrap(message, 80));
		return gd;
	}

	private boolean getPixelPitch(CalibrationWriter calibration)
	{
		ExtendedGenericDialog gd = newWizardDialog(
				"Enter the size of each pixel. This is required to ensure the dimensions of the image are calibrated.",
				"E.g. a camera with a 6.45um pixel size and a 60x objective will have a pitch of 6450/60 = 107.5nm.");
		// TODO - Add a pop-up calculator...
		gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
		return true;
	}

	private boolean getGain(CalibrationWriter calibration)
	{
		ExtendedGenericDialog gd = newWizardDialog("Enter the total gain.",
				"This is usually supplied with your camera certificate. The gain indicates how many Analogue-to-Digital-Units (ADUs) are recorded at the pixel for each photon registered on the sensor.",
				"The gain is usually expressed using the product of the EM-gain (if applicable), the camera gain and the sensor quantum efficiency.",
				"A value of 1 means no conversion to photons will occur.");
		// TODO - Add a wizard to allow calculation of total gain from EM-gain, camera gain and QE
		gd.addNumericField("Gain (ADU/photon)", calibration.getGain(), 2);
		gd.addCheckbox("EM-CCD", calibration.isEMCCD());
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.setGain(Math.abs(gd.getNextNumber()));
		calibration.setCameraType((gd.getNextBoolean()) ? CameraType.EMCCD : CameraType.CCD);
		return true;
	}

	private boolean getExposureTime(CalibrationWriter calibration)
	{
		ExtendedGenericDialog gd = newWizardDialog(
				"Enter the exposure time. Calibration of the exposure time allows correct reporting of on and off times.",
				"This is the length of time for each frame in the image.");
		gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.setExposureTime(Math.abs(gd.getNextNumber()));
		return true;
	}

	private boolean getPeakWidth(final CalibrationWriter calibration)
	{
		ExtendedGenericDialog gd = newWizardDialog("Enter the expected peak width in pixels.",
				"A point source of light will not be focussed perfectly by the microscope but will appear as a spread out peak. This Point Spread Function (PSF) can be modelled using a 2D Gaussian curve.",
				"An optimised optical system (lens and camera sensor) should have a peak standard deviation of approximately 1 pixel when in focus. This allows the fitting routine to have enough data to identify the centre of the peak without spreading the light over too many pixels (which increases noise).",
				"The peak width can be estimated using the wavelength of light emitted by the single molecules and the parameters of the microscope. Use a PSF calculator by clicking the checkbox below:");
		// Add ability to run the PSF Calculator to get the width
		gd.addCheckbox("Run_PSF_calculator", false);
		gd.addNumericField("Gaussian_SD", fitConfig.getInitialXSD(), 3);
		if (Utils.isShowGenericDialog())
		{
			final TextField textInitialPeakStdDev0 = (TextField) gd.getNumericFields().get(0);
			gd.addAndGetButton("Run PSF calculator", new ActionListener()
			{
				public void actionPerformed(ActionEvent e)
				{
					// Run the PSF Calculator
					PSFCalculator calculator = new PSFCalculator();
					calculatorSettings.pixelPitch = calibration.getNmPerPixel() / 1000.0;
					calculatorSettings.magnification = 1;
					calculatorSettings.beamExpander = 1;
					double sd = calculator.calculate(calculatorSettings, true);
					if (sd > 0)
						textInitialPeakStdDev0.setText(Double.toString(sd));
				}
			});
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
				String notes = template.getNotes();
				IJ.log("Applying template: " + templateName);
				if (!Utils.isNullOrEmpty(notes))
					IJ.log(notes);

				boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
				if (template.isFitEngineConfiguration())
				{
					//System.out.println(template.getFitEngineConfiguration().getFitEngineSettings());
					refreshSettings(template.getFitEngineConfiguration(), custom);
				}
				// TODO - Fix this
				//				if (template.isCalibration())
				//				{
				//					refreshSettings(template.getCalibration());
				//				}
				//				if (template.isResultsSettings())
				//				{
				//					refreshSettings(template.getResultsSettings());
				//				}
			}
		}
		else if (e.getSource() instanceof Checkbox)
		{
			if (e.getSource() == textSmartFilter)
			{
				// Prevent both filters being enabled
				textDisableSimpleFilter.setState(textSmartFilter.getState());
				updateFilterInput();
			}
			else if (e.getSource() == textDisableSimpleFilter)
			{
				updateFilterInput();
			}
		}
	}

	private void updateFilterInput()
	{
		if (textDisableSimpleFilter.getState())
		{
			sliderCoordinateShiftFactor.setEnabled(false);
			disableEditing(textCoordinateShiftFactor);
			disableEditing(textSignalStrength);
			disableEditing(textMinPhotons);
			// These are used to set bounds
			//disableEditing(textMinWidthFactor);
			//disableEditing(textWidthFactor);
			disableEditing(textPrecisionThreshold);
		}
		else
		{
			sliderCoordinateShiftFactor.setEnabled(true);
			enableEditing(textCoordinateShiftFactor);
			enableEditing(textSignalStrength);
			enableEditing(textMinPhotons);
			//enableEditing(textMinWidthFactor);
			//enableEditing(textWidthFactor);
			enableEditing(textPrecisionThreshold);
		}
	}

	private void disableEditing(TextField textField)
	{
		textField.setEditable(false);
		textField.setBackground(SystemColor.control);
	}

	private void enableEditing(TextField textField)
	{
		textField.setEditable(true);
		textField.setBackground(SystemColor.white);
	}

	private boolean readDialog(ExtendedGenericDialog gd, boolean isCrop)
	{
		// Ignore the template
		gd.getNextChoice();

		CalibrationWriter calibration = fitConfig.getCalibrationWriter();
		calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
		calibration.setGain(Math.abs(gd.getNextNumber()));
		calibration.setCameraType((gd.getNextBoolean()) ? CameraType.EMCCD : CameraType.CCD);
		calibration.setExposureTime(Math.abs(gd.getNextNumber()));

		// Note: The bias and read noise will just end up being what was in the configuration file
		// One fix for this is to save/load only the settings that are required from the configuration file
		// (the others will remain unchanged). This will require a big refactor of the settings save/load.
		// The simple fix is to create a plugin to allow the configuration to be changed for results.
		if (isCrop)
			ignoreBoundsForNoise = optionIgnoreBoundsForNoise = gd.getNextBoolean();

		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), false, 0);
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
			if (extraOptions)
				fitConfig.setBackgroundFitting(gd.getNextBoolean());
			config.setFailuresLimit((int) gd.getNextNumber());
			config.setIncludeNeighbours(gd.getNextBoolean());
			config.setNeighbourHeightThreshold(gd.getNextNumber());
			config.setResidualsThreshold(gd.getNextNumber());

			config.setDuplicateDistance(gd.getNextNumber());

			fitConfig.setSmartFilter(gd.getNextBoolean());
			fitConfig.setDisableSimpleFilter(gd.getNextBoolean());
			fitConfig.setCoordinateShiftFactor(gd.getNextNumber());
			fitConfig.setSignalStrength(gd.getNextNumber());
			fitConfig.setMinPhotons(gd.getNextNumber());
			if (extraOptions)
			{
				fitConfig.setNoise(gd.getNextNumber());
				config.setNoiseMethod(gd.getNextChoiceIndex());
			}
			fitConfig.setMinWidthFactor(gd.getNextNumber());
			fitConfig.setWidthFactor(gd.getNextNumber());
			fitConfig.setPrecisionThreshold(gd.getNextNumber());
		}

		resultsSettings.setLogProgress(gd.getNextBoolean());
		if (!maximaIdentification)
			resultsSettings.setShowDeviations(gd.getNextBoolean());

		resultsSettings.getResultsTableSettingsBuilder().setShowTable(gd.getNextBoolean());
		resultsSettings.getResultsImageSettingsBuilder().setImageTypeValue(gd.getNextChoiceIndex());
		if (extraOptions)
			showProcessedFrames = optionShowProcessedFrames = gd.getNextBoolean();
		resultsSettings.getResultsFileSettingsBuilder().setResultsDirectory(gd.getNextString());
		resultsSettings.getResultsFileSettingsBuilder().setFileFormatValue(gd.getNextChoiceIndex());
		resultsSettings.getResultsInMemorySettingsBuilder().setInMemory(gd.getNextBoolean());
		if (extraOptions)
			fractionOfThreads = Math.abs(gd.getNextNumber());

		gd.collectOptions();

		// Save to allow dialog state to be maintained even with invalid parameters
		saveFitEngineSettings();
		SettingsManager.writeSettings(resultsSettings.build());

		if (gd.invalidNumber())
			return false;

		// Check arguments
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.getNmPerPixel());
			Parameters.isAboveZero("Gain", calibration.getGain());
			Parameters.isAboveZero("Exposure time", calibration.getExposureTime());
			Parameters.isAboveZero("Initial SD0", fitConfig.getInitialXSD());
			if (fitConfig.getPSF().getParameterCount() > 1)
			{
				Parameters.isAboveZero("Initial SD1", fitConfig.getInitialYSD());
			}
			Parameters.isAboveZero("Search_width", config.getSearch());
			Parameters.isAboveZero("Fitting_width", config.getFitting());
			Parameters.isPositive("Integrate frames", integrateFrames);
			if (!maximaIdentification)
			{
				Parameters.isPositive("Failures limit", config.getFailuresLimit());
				Parameters.isPositive("Neighbour height threshold", config.getNeighbourHeightThreshold());
				Parameters.isPositive("Residuals threshold", config.getResidualsThreshold());
				Parameters.isPositive("Duplicate distance", config.getDuplicateDistance());

				if (!fitConfig.isSmartFilter())
				{
					Parameters.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
					Parameters.isPositive("Signal strength", fitConfig.getSignalStrength());
					Parameters.isPositive("Min photons", fitConfig.getMinPhotons());
				}
				if (extraOptions)
					Parameters.isPositive("Noise", fitConfig.getNoise());
				if (!fitConfig.isSmartFilter())
				{
					Parameters.isPositive("Min width factor", fitConfig.getMinWidthFactor());
					Parameters.isPositive("Width factor", fitConfig.getMaxWidthFactor());
					Parameters.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
				}
			}
			final ResultsImageSettings.Builder imageSettings = resultsSettings.getResultsImageSettingsBuilder();
			if (imageSettings.getImageType() == ResultsImageType.DRAW_INTENSITY_AVERAGE_PRECISION ||
					imageSettings.getImageType() == ResultsImageType.DRAW_LOCALISATIONS_AVERAGE_PRECISION)
			{
				Parameters.isAboveZero("Image precision", imageSettings.getAveragePrecision());
			}
			Parameters.isAboveZero("Image scale", imageSettings.getScale());
			if (extraOptions)
				Parameters.isPositive("Image rolling window", imageSettings.getRollingWindowSize());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		int flags = (extraOptions) ? FLAG_EXTRA_OPTIONS : 0;

		// If precision filtering then we need the camera bias
		if (!maximaIdentification)
		{
			if (!configureSmartFilter(config, flags))
				return false;

			if (!fitConfig.isSmartFilter() && fitConfig.getPrecisionThreshold() > 0)
			{
				gd = new ExtendedGenericDialog(TITLE);
				gd.addMessage(
						"Precision filtering can use global noise estimate or local background level.\n \nLocal background requires the camera bias:");
				gd.addCheckbox("Local_background", fitConfig.isPrecisionUsingBackground());
				gd.addNumericField("Camera_bias (ADUs)", calibration.getBias(), 2);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				fitConfig.setPrecisionUsingBackground(gd.getNextBoolean());
				calibration.setBias(Math.abs(gd.getNextNumber()));
			}
		}

		if (!configureDataFilter(config, flags))
			return false;

		// Second dialog for solver dependent parameters
		if (!maximaIdentification)
		{
			if (!configureFitSolver(config, flags))
				return false;
		}

		// Extra parameters are needed for interlaced data
		if (interlacedData)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Interlaced data requires a repeating pattern of frames to process.\n" +
					"Describe the regular repeat of the data:\n \n" + "Start = The first frame that contains data\n" +
					"Block = The number of continuous frames containing data\n" +
					"Skip = The number of continuous frames to ignore before the next data\n \n" +
					"E.G. 2:9:1 = Data was imaged from frame 2 for 9 frames, 1 frame to ignore, then repeat.");
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

		;

		boolean result = saveFitEngineSettings();
		if (!result)
			IJ.error(TITLE, "Failed to save settings");

		return result;
	}

	/**
	 * Show a dialog to configure the smart filter. The updated settings are saved to the settings file.
	 * <p>
	 * If the fit configuration isSmartFilter is not enabled then this method returns true. If it is enabled then a
	 * dialog is shown to input the configuration for a smart filter. If no valid filter can be created from the input
	 * then the method returns false.
	 * <p>
	 * Note: If the smart filter is successfully configured then the user may want to disable the standard fit
	 * validation.
	 *
	 * @param config
	 *            the config
	 * @param flags
	 *            the flags
	 * @return true, if successful
	 */
	public static boolean configureSmartFilter(FitEngineConfiguration config, int flags)
	{
		FitConfiguration fitConfig = config.getFitConfiguration();
		CalibrationWriter calibration = new CalibrationWriter(fitConfig.getCalibration());
		if (!fitConfig.isSmartFilter())
			return true;

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

		String xml = fitConfig.getSmartFilterString();
		if (Utils.isNullOrEmpty(xml))
			xml = fitConfig.getDefaultSmartFilterXML();

		gd.addMessage("Smart filter (used to pick optimum results during fitting)");
		gd.addTextAreas(XmlUtils.convertQuotes(xml), null, 8, 60);
		// Currently we just collect it here even if not needed
		gd.addMessage(
				"Smart filters using precision filtering may require a local background level.\n \nLocal background requires the camera bias:");
		gd.addNumericField("Camera_bias (ADUs)", calibration.getBias(), 2);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		xml = gd.getNextText();
		Filter f = DirectFilter.fromXML(xml);
		if (f == null || !(f instanceof DirectFilter))
			return false;
		fitConfig.setDirectFilter((DirectFilter) f);

		calibration.setBias(Math.abs(gd.getNextNumber()));

		if (BitFlags.anyNotSet(flags, FLAG_NO_SAVE))
		{
			fitConfig.setCalibration(calibration.getCalibration());
			SettingsManager.writeSettings(config.getFitEngineSettings());
		}
		return true;
	}

	/**
	 * Show a dialog to configure the data filter. The updated settings are saved to the settings file. An error
	 * message is shown if the dialog is cancelled or the configuration is invalid.
	 *
	 * @param config
	 *            the config
	 * @param flags
	 *            the flags
	 * @return True if the configuration succeeded
	 */
	public static boolean configureDataFilter(FitEngineConfiguration config, int flags)
	{
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

		String[] filterNames = SettingsManager.getDataFilterMethodNames();
		DataFilterMethod[] filterValues = SettingsManager.getDataFilterMethodValues();

		for (int i = 1; i < n; i++)
		{
			int filter = i + 1;
			ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.enableYesNoCancel("Add", "Continue");
			gd.addMessage(
					String.format("Configure the %s filter.\nClick continue to proceed with the current set of %d.",
							FitConfigHelper.getName(config.getDataFilterType()), i));
			String fieldName = "Spot_filter" + filter;
			if (IJ.isMacro())
				// Use blank default value so bad macro parameters return nothing
				gd.addStringField(fieldName, "");
			else
				gd.addChoice(fieldName, filterNames, filterNames[config.getDataFilterMethod(i).ordinal()]);
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
				config.setDataFilter(filterValues[filterIndex], Math.abs(gd.getNextNumber()), false, i);
				numberOfFilters++;
			}
			else
			{
				break;
			}
		}
		config.setNumberOfFilters(numberOfFilters);
		if (BitFlags.anyNotSet(flags, FLAG_NO_SAVE))
			saveFitEngineSettings(config);
		return true;
	}

	/** Flag to indicate that additional options can be configured. */
	public static final int FLAG_EXTRA_OPTIONS = 0x00000001;
	/** Flag to indicate that the calibration should not be configured. */
	public static final int FLAG_IGNORE_CALIBRATION = 0x00000002;
	/** Flag to indicate that configuration should not be saved. */
	public static final int FLAG_NO_SAVE = 0x00000004;

	/**
	 * Show a dialog to configure the fit solver. The updated settings are saved to the settings file. An error
	 * message is shown if the dialog is cancelled or the configuration is invalid.
	 *
	 * @param config
	 *            the config
	 * @param flags
	 *            the flags
	 * @return True if the configuration succeeded
	 */
	public static boolean configureFitSolver(FitEngineConfiguration config, int flags)
	{
		boolean extraOptions = BitFlags.anySet(flags, FLAG_EXTRA_OPTIONS);
		boolean ignoreCalibration = BitFlags.anySet(flags, FLAG_IGNORE_CALIBRATION);
		boolean saveSettings = BitFlags.anyNotSet(flags, FLAG_NO_SAVE);

		FitConfiguration fitConfig = config.getFitConfiguration();
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();

		FitSolver fitSolver = fitConfig.getFitSolver();

		boolean isLVM = fitSolver == FitSolver.LVM_LSE || fitSolver == FitSolver.LVM_WLSE ||
				fitSolver == FitSolver.LVM_MLE;
		boolean isFastMLE = fitSolver == FitSolver.FAST_MLE || fitSolver == FitSolver.BACKTRACKING_FAST_MLE;
		boolean isSteppingFunctionSolver = isLVM || isFastMLE;

		if (fitSolver == FitSolver.MLE)
		{
			ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Maximum Likelihood Estimation requires additional parameters");
			if (!ignoreCalibration)
			{
				gd.addNumericField("Camera_bias (ADUs)", calibration.getBias(), 2);
				gd.addCheckbox("Model_camera_noise", fitConfig.isModelCamera());
				gd.addNumericField("Read_noise (ADUs)", calibration.getReadNoise(), 2);
				gd.addNumericField("Amplification (ADU/electron)", calibration.getAmplification(), 2);
				gd.addCheckbox("EM-CCD", calibration.isEMCCD());
			}
			String[] searchNames = SettingsManager.getNames((Object[]) MaximumLikelihoodFitter.SearchMethod.values());
			gd.addChoice("Search_method", searchNames, searchNames[fitConfig.getSearchMethod().ordinal()]);
			gd.addStringField("Relative_threshold", Utils.rounded(fitConfig.getRelativeThreshold()));
			gd.addStringField("Absolute_threshold", Utils.rounded(fitConfig.getAbsoluteThreshold()));
			gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);
			gd.addNumericField("Max_function_evaluations", fitConfig.getMaxFunctionEvaluations(), 0);
			if (extraOptions)
				gd.addCheckbox("Gradient_line_minimisation", fitConfig.isGradientLineMinimisation());
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			if (!ignoreCalibration)
			{
				calibration.setBias(Math.abs(gd.getNextNumber()));
				fitConfig.setModelCamera(gd.getNextBoolean());
				calibration.setReadNoise(Math.abs(gd.getNextNumber()));
				calibration.setAmplification(Math.abs(gd.getNextNumber()));
				calibration.setCameraType((gd.getNextBoolean()) ? CameraType.EMCCD : CameraType.CCD);
			}
			fitConfig.setSearchMethod(gd.getNextChoiceIndex());
			fitConfig.setRelativeThreshold(getNumber(gd));
			fitConfig.setAbsoluteThreshold(getNumber(gd));
			fitConfig.setMaxIterations((int) gd.getNextNumber());
			fitConfig.setMaxFunctionEvaluations((int) gd.getNextNumber());
			if (extraOptions)
				fitConfig.setGradientLineMinimisation(gd.getNextBoolean());
			else
				// This option is for the Conjugate Gradient optimiser and makes it less stable
				fitConfig.setGradientLineMinimisation(false);

			if (saveSettings)
				saveFitEngineSettings(config);

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
		else if (isSteppingFunctionSolver)
		{
			boolean requireCalibration = fitSolver != FitSolver.LVM_LSE;

			// Collect options for LVM fitting
			ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage(FitConfigHelper.getName(fitSolver) + " requires additional parameters");
			gd.addStringField("Relative_threshold", Utils.rounded(fitConfig.getRelativeThreshold()));
			gd.addStringField("Absolute_threshold", Utils.rounded(fitConfig.getAbsoluteThreshold()));
			gd.addStringField("Parameter_relative_threshold", Utils.rounded(fitConfig.getParameterRelativeThreshold()));
			gd.addStringField("Parameter_absolute_threshold", Utils.rounded(fitConfig.getParameterAbsoluteThreshold()));
			gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);
			if (isLVM)
				gd.addNumericField("Lambda", fitConfig.getLambda(), 4);
			if (isFastMLE)
				gd.addCheckbox("Fixed_iterations", fitConfig.isFixedIterations());

			// Extra parameters are needed for calibrated fit solvers
			if (requireCalibration)
			{
				// TODO - this must be refactored to get a camera calibration
				if (calibration.getCalibrationOrBuilder().hasCameraCalibration())
				{

				}
				gd.addNumericField("Camera_bias (ADUs)", calibration.getBias(), 2);
				gd.addNumericField("Gain (ADU/photon)", calibration.getGain(), 2);
			}

			gd.addCheckbox("Use_clamping", fitConfig.isUseClamping());
			gd.addCheckbox("Dynamic_clamping", fitConfig.isUseDynamicClamping());
			PSF psf = fitConfig.getPSF();
			boolean isAstigmatism = psf.getPsfType() == PSFType.ASTIGMATIC_GAUSSIAN_2D;
			int nParams = PSFHelper.getParameterCount(psf);
			if (extraOptions)
			{
				gd.addNumericField("Clamp_background", fitConfig.getClampBackground(), 2);
				gd.addNumericField("Clamp_signal", fitConfig.getClampSignal(), 2);
				gd.addNumericField("Clamp_x", fitConfig.getClampX(), 2);
				gd.addNumericField("Clamp_y", fitConfig.getClampY(), 2);
				if (isAstigmatism)
					gd.addNumericField("Clamp_z", fitConfig.getClampZ(), 2);
				else
				{
					if (nParams > 1 || !fitConfig.isFixedPSF())
						gd.addNumericField("Clamp_sx", fitConfig.getClampXSD(), 2);
					if (nParams > 1)
						gd.addNumericField("Clamp_sy", fitConfig.getClampYSD(), 2);
					if (nParams > 2)
						gd.addNumericField("Clamp_angle", fitConfig.getClampAngle(), 2);
				}
			}

			gd.showDialog();
			if (gd.wasCanceled())
				return false;

			fitConfig.setRelativeThreshold(getNumber(gd));
			fitConfig.setAbsoluteThreshold(getNumber(gd));
			fitConfig.setParameterRelativeThreshold(getNumber(gd));
			fitConfig.setParameterAbsoluteThreshold(getNumber(gd));
			fitConfig.setMaxIterations((int) gd.getNextNumber());
			if (isLVM)
				fitConfig.setLambda(gd.getNextNumber());
			if (isFastMLE)
				fitConfig.setFixedIterations(gd.getNextBoolean());

			if (requireCalibration)
			{
				calibration.setBias(Math.abs(gd.getNextNumber()));
				calibration.setGain(Math.abs(gd.getNextNumber()));
			}

			fitConfig.setUseClamping(gd.getNextBoolean());
			fitConfig.setUseDynamicClamping(gd.getNextBoolean());
			if (extraOptions)
			{
				fitConfig.setClampBackground(Math.abs(gd.getNextNumber()));
				fitConfig.setClampSignal(Math.abs(gd.getNextNumber()));
				fitConfig.setClampX(Math.abs(gd.getNextNumber()));
				fitConfig.setClampY(Math.abs(gd.getNextNumber()));
				if (isAstigmatism)
					fitConfig.setClampZ(Math.abs(gd.getNextNumber()));
				else
				{
					if (nParams > 1 || !fitConfig.isFixedPSF())
						fitConfig.setClampXSD(Math.abs(gd.getNextNumber()));
					if (nParams > 1)
						fitConfig.setClampYSD(Math.abs(gd.getNextNumber()));
					if (nParams > 2)
						fitConfig.setClampAngle(Math.abs(gd.getNextNumber()));
				}
			}

			if (saveSettings)
			{
				saveFitEngineSettings(config);
			}

			try
			{
				if (isLVM)
					Parameters.isAboveZero("Lambda", fitConfig.getLambda());
				// This call will check if the convergence criteria are OK
				fitConfig.getFunctionSolver();
			}
			catch (IllegalArgumentException e)
			{
				IJ.error(TITLE, e.getMessage());
				return false;
			}
		}
		else
		{
			IJ.error(TITLE, "Unknown fit solver: " + fitSolver);
			return false;
		}

		if (config.isIncludeNeighbours())
		{
			if (!fitConfig.getFunctionSolver().isBounded())
			{
				IJ.error(TITLE, "Including neighbours requires a bounded fit solver");
				return false;
			}
		}

		return true;
	}

	private static double getNumber(ExtendedGenericDialog gd)
	{
		try
		{
			return Math.abs(Double.parseDouble(gd.getNextString()));
		}
		catch (NumberFormatException e)
		{
			return 0;
		}
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

	private void addTableResults(PeakResultsList resultsList)
	{
		IJTablePeakResults r = ResultsManager.addTableResults(resultsList, resultsSettings.getResultsTableSettings(),
				resultsSettings.getShowDeviations(), false, false);
		if (r != null)
		{
			r.setShowFittingData(true);
			r.setShowNoiseData(true);
			r.setShowZ(PSFHelper.is3D(resultsList.getPSF()));
			r.setClearAtStart(simpleFit);
			r.setShowEndFrame(integrateFrames > 1);
		}
	}

	private void addFileResults(PeakResultsList resultsList)
	{
		ResultsFileSettings resultsSettings = this.resultsSettings.getResultsFileSettings();
		String resultsFilename = null;
		if (resultsSettings.getResultsDirectory() != null && new File(resultsSettings.getResultsDirectory()).exists())
		{
			resultsFilename = resultsSettings.getResultsDirectory() + File.separatorChar + source.getName() +
					".results." + ResultsConfigHelper.getExtension(resultsSettings.getFileFormat());
		}
		else
		{
			resultsFilename = resultsSettings.getResultsFilename();
		}
		PeakResults r = ResultsManager.addFileResults(resultsList, resultsSettings, resultsFilename,
				this.resultsSettings.getShowDeviations(), integrateFrames > 1, false);
		if (r instanceof FilePeakResults)
		{
			FilePeakResults fr = (FilePeakResults) r;
			fr.setSortAfterEnd(Prefs.getThreads() > 1);
		}
	}

	private void addMemoryResults(PeakResultsList resultsList, boolean force)
	{
		if (resultsSettings.getResultsInMemorySettings().getInMemory() || force)
		{
			MemoryPeakResults results = new MemoryPeakResults();
			results.setSortAfterEnd(Prefs.getThreads() > 1);
			resultsList.addOutput(results);
			MemoryPeakResults.addResults(results);
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
		}

		addSingleFrameOverlay();
	}

	private void addSingleFrameOverlay()
	{
		// If a single frame was processed add the peaks as an overlay if they are in memory
		ImagePlus imp = this.imp;

		if (fitMaxima && singleFrame > 0)
		{
			if (source instanceof IJImageSource)
			{
				String title = source.getName();
				imp = WindowManager.getImage(title);
			}
		}

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

			ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			gd.enableYesNoCancel();
			gd.hideCancelButton();
			gd.addMessage("Add the fitted localisations as an overlay?");
			gd.showDialog();
			if (!gd.wasOKed())
				return;

			final LUT lut = LUTHelper.createLUT(LutColour.ICE);
			final Overlay o = new Overlay();
			final int size = results.size();
			final Counter j = new Counter(size);
			final ImagePlus finalImp = imp;
			results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure()
			{
				public void executeXYR(float x, float y, PeakResult r)
				{
					PointRoi roi = new PointRoi(x, y);
					Color c = LUTHelper.getColour(lut, j.decrementAndGet(), size);
					roi.setStrokeColor(c);
					roi.setFillColor(c);
					if (finalImp.getStackSize() > 1)
						roi.setPosition(singleFrame);
					o.add(roi);
				}
			});
			imp.setOverlay(o);
			imp.getWindow().toFront();
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
		Rectangle cropBounds = (bounds.x == 0 && bounds.y == 0 && bounds.width == source.getWidth() &&
				bounds.height == source.getHeight()) ? null : bounds;

		// Use the FitEngine to allow multi-threading.
		FitEngine engine = createFitEngine(getNumberOfThreads(totalFrames));

		final int step = Utils.getProgressInterval(totalFrames);

		runTime = System.nanoTime();
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
				stack.addSlice(String.format("Frame %d - %d", source.getStartFrameNumber(), source.getEndFrameNumber()),
						data);
			}

			// Get the frame number from the source to allow for interlaced and aggregated data
			engine.run(createJob(source.getStartFrameNumber(), source.getEndFrameNumber(), data, bounds, noise));

			if (escapePressed())
				shutdown = true;
		}

		engine.end(shutdown);
		time = engine.getTime();
		runTime = System.nanoTime() - runTime;

		if (showProcessedFrames)
			Utils.display("Processed frames", stack);

		showResults();

		source.close();
	}

	private int getNumberOfThreads(int totalFrames)
	{
		final int t = Math.max(1, (int) (fractionOfThreads * Prefs.getThreads()));
		return FastMath.min(totalFrames, t);
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
		// Use a large queue size to allow images read from disk to be pre-loaded.
		return createFitEngine(numberOfThreads, FitQueue.BLOCKING, numberOfThreads * 10);
	}

	/**
	 * Creates a fitting engine using the current configuration.
	 * 
	 * @param numberOfThreads
	 * @param queue
	 * @param queueSize
	 * @return The fiting engine
	 */
	public FitEngine createFitEngine(int numberOfThreads, FitQueue queue, int queueSize)
	{
		// Ensure thread safety
		PeakResults r = (numberOfThreads > 1) ? results.getThreadSafeList() : results;
		if (results.numberOfOutputs() == 1)
			// Reduce to single object for speed
			r = results.toArray()[0];

		// Update the configuration
		updateFitConfiguration(config);

		FitEngine engine = new FitEngine(config, r, numberOfThreads, queue, queueSize);

		// Write settings out to the IJ log
		if (resultsSettings.getLogProgress())
		{
			IJ.log("-=-=-=-");
			IJ.log("Peak Fit");
			IJ.log("-=-=-=-");
			Utils.log("Initial Peak SD = %s,%s", Utils.rounded(fitConfig.getInitialXSD()),
					Utils.rounded(fitConfig.getInitialYSD()));
			SpotFilter spotFilter = engine.getSpotFilter();
			IJ.log("Spot Filter = " + spotFilter.getDescription());
			int w = 2 * engine.getFitting() + 1;
			Utils.log("Fit window = %d x %d", w, w);
			if (!fitConfig.isDisableSimpleFilter())
			{
				IJ.log("Coordinate shift = " + Utils.rounded(config.getFitConfiguration().getCoordinateShift()));
				IJ.log("Signal strength = " + Utils.rounded(fitConfig.getSignalStrength()));
			}
			if (fitConfig.isDirectFilter())
				IJ.log("Smart filter = " + fitConfig.getSmartFilter().getDescription());
			if (extraOptions)
				IJ.log("Noise = " + Utils.rounded(fitConfig.getNoise()));
			IJ.log("Width factor = " + Utils.rounded(fitConfig.getMaxWidthFactor()));
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
		logger = (resultsSettings.getLogProgress()) ? new IJLogger() : null;
		fitConfig.setLog(logger);

		fitConfig.setComputeDeviations(resultsSettings.getShowDeviations());

		config.configureOutputUnits();
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

		final int totalFrames = source.getFrames();

		// Store the indices of each new time-frame
		results.sort();

		// Use the FitEngine to allow multi-threading.
		final FitEngine engine = createFitEngine(getNumberOfThreads(totalFrames));

		final int step = Utils.getProgressInterval(totalFrames);

		runTime = System.nanoTime();
		final ArrayList<PeakResult> sliceCandidates = new ArrayList<PeakResult>();
		final FrameCounter counter = new FrameCounter(results.getFirstFrame());
		results.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult r)
			{
				if (counter.advance(r.getFrame()))
				{
					if (escapePressed())
						return true;
					int slice = counter.previousFrame();
					if (slice % step == 0)
					{
						if (Utils.showStatus("Slice: " + slice + " / " + totalFrames))
							IJ.showProgress(slice, totalFrames);
					}

					// Process results
					if (!processResults(engine, sliceCandidates, slice))
						return true;

					sliceCandidates.clear();
				}
				sliceCandidates.add(r);

				return false;
			}
		});

		// Process final results
		boolean shutdown = escapePressed();
		if (!shutdown)
			processResults(engine, sliceCandidates, counter.currentFrame());

		engine.end(shutdown);
		time = engine.getTime();
		runTime = System.nanoTime() - runTime;

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
			if (result instanceof ExtendedPeakResult && result.getFrame() != result.getEndFrame())
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
		// passing data through all methods in the PeakResultsList. However the list
		// returns the size using the first entry so this is OK.
		return results.size();
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

	private void refreshSettings(Calibration cal)
	{
		if (cal == null)
			return;

		CalibrationWriter calibration = fitConfig.getCalibrationWriter();
		// Do not use set() as we support merging a partial calibration
		calibration.mergeCalibration(cal);

		if (calibration.hasNmPerPixel())
			textNmPerPixel.setText("" + calibration.getNmPerPixel());
		if (calibration.hasGain())
			textGain.setText("" + calibration.getGain());
		textEMCCD.setState(calibration.isEMCCD());
		if (calibration.hasExposureTime())
			textExposure.setText("" + calibration.getExposureTime());
	}

	/**
	 * Refresh settings.
	 * <p>
	 * If this is a custom template then use all the settings. If a default template then leave some existing spot
	 * settings untouched as the user may have updated them (e.g. PSF width).
	 *
	 * @param pConfig
	 *            the config
	 * @param isCustomTemplate
	 *            True if a custom template.
	 */
	private void refreshSettings(FitEngineConfiguration pConfig, boolean isCustomTemplate)
	{
		// Set the configuration
		// Preserve the current settings for things we do not want to update
		FitSettings current = fitConfig.getFitSettings();
		// This will clear everything and merge the configuration
		this.config.setFitEngineSettings(pConfig.getFitEngineSettings());
		fitConfig = this.config.getFitConfiguration();
		// Reset what we want to preserve
		fitConfig.setCalibration(current.getCalibration());
		fitConfig.setPSF(current.getPsf());

		if (isCustomTemplate)
			textPSF.select(getPSFTypeNames()[fitConfig.getPSFType().ordinal()]);

		textDataFilterType.select(SettingsManager.getDataFilterTypeNames()[config.getDataFilterType().ordinal()]);
		textDataFilterMethod
				.select(SettingsManager.getDataFilterMethodNames()[config.getDataFilterMethod(0).ordinal()]);
		textSmooth.setText("" + config.getSmooth(0));
		textSearch.setText("" + config.getSearch());
		textBorder.setText("" + config.getBorder());
		textFitting.setText("" + config.getFitting());
		if (!maximaIdentification)
		{
			textFitSolver.select(SettingsManager.getFitSolverNames()[fitConfig.getFitSolver().ordinal()]);
			if (extraOptions)
				textFitBackground.setState(fitConfig.isBackgroundFitting());
			textFailuresLimit.setText("" + config.getFailuresLimit());
			textIncludeNeighbours.setState(config.isIncludeNeighbours());
			textNeighbourHeightThreshold.setText("" + config.getNeighbourHeightThreshold());
			textResidualsThreshold.setText("" + config.getResidualsThreshold());
			textDuplicateDistance.setText("" + config.getDuplicateDistance());

			// Filtering
			textSmartFilter.setState(fitConfig.isSmartFilter());
			textDisableSimpleFilter.setState(fitConfig.isDisableSimpleFilter());
			if (!fitConfig.isDisableSimpleFilter())
			{
				textCoordinateShiftFactor.setText("" + fitConfig.getCoordinateShiftFactor());
				textSignalStrength.setText("" + fitConfig.getSignalStrength());
				textWidthFactor.setText("" + fitConfig.getMaxWidthFactor());
				textPrecisionThreshold.setText("" + fitConfig.getPrecisionThreshold());
			}
			// These are used for settings the bounds so they are included
			textMinPhotons.setText("" + fitConfig.getMinPhotons());
			textMinWidthFactor.setText("" + fitConfig.getMinWidthFactor());
			updateFilterInput();

			if (extraOptions)
			{
				textNoise.setText("" + fitConfig.getNoise());
				textNoiseMethod.select(config.getNoiseMethod().ordinal());
			}
		}
	}

	private void refreshSettings(ResultsSettings resultsSettings)
	{
		this.resultsSettings = resultsSettings.toBuilder();
		textLogProgress.setState(resultsSettings.getLogProgress());
		if (!maximaIdentification)
			textShowDeviations.setState(resultsSettings.getShowDeviations());
		textResultsTable.setState(resultsSettings.getResultsTableSettings().getShowTable());
		textResultsImage.select(resultsSettings.getResultsImageSettings().getImageTypeValue());
		textResultsDirectory.setText("" + resultsSettings.getResultsFileSettings().getResultsDirectory());
		textFileFormat.select(resultsSettings.getResultsFileSettings().getFileFormatValue());
		textResultsInMemory.setState(resultsSettings.getResultsInMemorySettings().getInMemory());
	}
}
