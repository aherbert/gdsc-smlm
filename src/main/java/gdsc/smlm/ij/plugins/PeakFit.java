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

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
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

import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.IJLogger;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.SeriesOpener;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Logger;
import gdsc.core.utils.BitFlags;
import gdsc.core.utils.Maths;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import gdsc.smlm.data.config.FitProtos.FitSolver;
import gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import gdsc.smlm.data.config.FitProtosHelper;
import gdsc.smlm.data.config.GUIProtos.PSFCalculatorSettings;
import gdsc.smlm.data.config.GUIProtosHelper;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFParameter;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.ResultsProtos.ResultsFileSettings;
import gdsc.smlm.data.config.ResultsProtos.ResultsImageSettings;
import gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import gdsc.smlm.data.config.ResultsProtosHelper;
import gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
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
import gdsc.smlm.fitting.nonlinear.FastMLESteppingFunctionSolver;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.SeriesImageSource;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.IJImageConverter;
import gdsc.smlm.model.camera.CameraModel;
import gdsc.smlm.model.camera.PerPixelCameraModel;
import gdsc.smlm.results.AggregatedImageSource;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.FilePeakResults;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.InterlacedImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.count.Counter;
import gdsc.smlm.results.count.FrameCounter;
import gdsc.smlm.results.filter.DirectFilter;
import gdsc.smlm.results.filter.Filter;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionCollectedEvent;
import ij.gui.ExtendedGenericDialog.OptionCollectedListener;
import ij.gui.ExtendedGenericDialog.OptionListener;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LutColour;

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
	private String resultsSuffix = null;
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
	private static PSFCalculatorSettings.Builder calculatorSettings = GUIProtosHelper.defaultPSFCalculatorSettings
			.toBuilder();

	// All the fields that will be updated when reloading the configuration file
	private Choice textCameraType;
	private TextField textNmPerPixel;
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
	private TextField textPassRate;
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

	/**
	 * Instantiates a new peak fit.
	 */
	public PeakFit()
	{
		init(null, null);
	}

	/**
	 * Instantiates a new peak fit.
	 *
	 * @param config
	 *            the config
	 */
	public PeakFit(FitEngineConfiguration config)
	{
		init(config, null);
	}

	/**
	 * Instantiates a new peak fit.
	 *
	 * @param config
	 *            the config
	 * @param resultsSettings
	 *            the results settings
	 */
	public PeakFit(FitEngineConfiguration config, ResultsSettings resultsSettings)
	{
		init(config, resultsSettings);
	}

	private void init(FitEngineConfiguration config, ResultsSettings resultsSettings)
	{
		this.config = (config != null) ? config : new FitEngineConfiguration();
		fitConfig = this.config.getFitConfiguration();
		this.resultsSettings = (resultsSettings != null) ? ResultsSettings.newBuilder(resultsSettings)
				: ResultsSettings.newBuilder();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		plugin_flags = FLAGS;
		extraOptions = Utils.isExtraOptions();

		maximaIdentification = (arg != null && arg.contains("spot"));
		fitMaxima = (arg != null && arg.contains("maxima"));
		simpleFit = (arg != null && arg.contains("simple"));
		final boolean runSeries = (arg != null && arg.contains("series"));

		ImageSource imageSource = null;
		if (fitMaxima)
		{
			imp = null;
			// The maxima will have been identified already.
			// The image source will be found from the peak results.
			if (!showMaximaDialog())
				return DONE;

			final MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL);
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
				@Override
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

			final SeriesImageSource seriesImageSource = new SeriesImageSource(getName(series.getImageList()), series);
			seriesImageSource.setTrackProgress(new IJTrackProgress());
			//if (extraOptions)
			//{
			//	numberOfThreads = Math.max(1, series.getNumberOfThreads());
			//	seriesImageSource.setNumberOfThreads(numberOfThreads);
			//}
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
				final Object o = imp.getProperty("Info");
				final Pattern pattern = Pattern.compile("Source: (<.*IJImageSource>.*<.*IJImageSource>)", Pattern.DOTALL);
				final Matcher match = pattern.matcher((o == null) ? "" : o.toString());
				if (match.find())
				{
					final ImageSource source = ImageSource.fromXML(match.group(1));
					if (source instanceof IJImageSource)
					{
						tmpImageSource = (IJImageSource) source;
						if (!tmpImageSource.open())
							tmpImageSource = null;
						else
							imp = WindowManager.getImage(tmpImageSource.getName());
					}
				}

				if (tmpImageSource == null)
				{
					// Look for a parent using the title
					final String parentTitle = imp.getTitle().substring(0,
							imp.getTitle().length() - IJImagePeakResults.IMAGE_SUFFIX.length() - 1);
					final ImagePlus parentImp = WindowManager.getImage(parentTitle);
					if (parentImp != null)
					{
						tmpImageSource = new IJImageSource(parentImp);
						imp = parentImp;
					}
				}
				String message = "The selected image may be a previous fit result";
				if (tmpImageSource != null)
				{
					if (!TextUtils.isNullOrEmpty(tmpImageSource.getName()))
						message += " of: \n \n" + tmpImageSource.getName();
					message += " \n \nFit the parent?";
				}
				else
					message += " \n \nDo you want to continue?";

				final YesNoCancelDialog d = new YesNoCancelDialog(null, TITLE, message);
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
				try
				{
					imageSource = new IJImageSource(imp);
				}
				catch (final IllegalArgumentException e)
				{
					// This can happen if the image has an origin not in integer pixels
					// e.g. the plugin is run on a plot
					IJ.error(TITLE, "Error using image: " + imp.getTitle() + "\n \n" + e.getMessage());
					return DONE;
				}
		}

		time = -1;

		if (!initialiseImage(imageSource, getBounds(imp), false))
		{
			IJ.error(TITLE, "Failed to initialise the source image: " + imageSource.getName());
			return DONE;
		}

		final int flags = showDialog(imp);
		if ((flags & DONE) == 0)
			// Repeat so that we pass in the selected option for ignoring the bounds.
			// This should not be necessary since it is set within the readDialog method.
			//if (ignoreBoundsForNoise)
			//	initialiseImage(imageSource, bounds, ignoreBoundsForNoise);
			initialiseFitting();
		return flags;
	}

	/**
	 * Gets the name.
	 *
	 * @param imageList
	 *            the image list
	 * @return the name
	 */
	static String getName(String[] imageList)
	{
		String name = imageList[0];
		// Remove directory
		final int index = name.lastIndexOf(File.separatorChar);
		if (index > -1)
			name = name.substring(index + 1);
		//// Remove suffix
		//index = name.lastIndexOf('.');
		//if (index > 0)
		//{
		//	name = name.substring(0, index);
		//}
		return "Series " + name;
	}

	private static Rectangle getBounds(ImagePlus imp)
	{
		if (imp == null)
			return null;
		final Roi roi = imp.getRoi();
		if (roi != null && roi.isArea())
			return roi.getBounds();
		return null;
	}

	/**
	 * @return An input directory containing a series of images
	 */
	@SuppressWarnings("unused")
	private static String getInputDirectory(String title)
	{
		final JFileChooser chooser = new JFileChooser()
		{
			private static final long serialVersionUID = 275144634537614122L;

			@Override
			public void approveSelection()
			{
				if (getSelectedFile().isFile())
					return;
				super.approveSelection();
			}
		};
		if (System.getProperty("os.name").startsWith("Mac OS X"))
			chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		else
			chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		chooser.setDialogTitle(title);
		final int returnVal = chooser.showOpenDialog(IJ.getInstance());
		if (returnVal == JFileChooser.APPROVE_OPTION)
			return chooser.getSelectedFile().getPath();
		return null;
	}

	/**
	 * Initialise a new image for fitting and prepare the output results.
	 * <p>
	 * Calls {@link #initialise(ImageSource, Rectangle, boolean)} then {@link #initialiseFitting()}.
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

		// The bounds must fit in the image
		try
		{
			imageSource.checkBounds(bounds);
		}
		catch (final RuntimeException e)
		{
			return false;
		}

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
	 * Sets the results suffix.
	 *
	 * @param resultsSuffix
	 *            the new results suffix
	 */
	public void setResultsSuffix(String resultsSuffix)
	{
		this.resultsSuffix = resultsSuffix;
	}

	/**
	 * Set-up the fitting using all the configured properties. Prepare the output results.
	 *
	 * @return true, if successful
	 */
	public boolean initialiseFitting()
	{
		if (source == null)
			return false;

		// Do this to ensure the serialised configuration is correct
		updateFitConfiguration(config);

		results.setSource(source);
		String name = source.getName();
		if (!TextUtils.isNullOrEmpty(resultsSuffix))
			name += " " + resultsSuffix;
		if (maximaIdentification)
			name += " (Maxima)";
		else if (fitMaxima)
			name += " (" + getSolverName() + " Fit Maxima)";
		else
			name += " (" + getSolverName() + ")";
		results.setName(name);
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
		results.setConfiguration(SettingsManager.toJSON(config.getFitEngineSettings()));

		// This is added first as it cannot be closed. If the table is closed then the
		// number of results at the end is reported incorrectly.
		addMemoryResults(results, false);

		addTableResults(results);
		ResultsManager.addImageResults(results, resultsSettings.getResultsImageSettings(), bounds,
				(extraOptions) ? ResultsManager.FLAG_EXTRA_OPTIONS : 0);
		addFileResults(results);
		addDefaultResults(results);

		results.begin();

		if (simpleFit && showImage)
			for (final PeakResults r : results.toArray())
				if (r instanceof IJImagePeakResults)
				{
					final ImagePlus i = ((IJImagePeakResults) r).getImagePlus();
					Utils.log("Super-resolution image title = " + i.getTitle());
					WindowManager.toFront(i.getWindow());
				}

		return true;
	}

	private String getSolverName()
	{
		return getSolverName(config.getFitConfiguration());
	}

	/**
	 * Gets the solver name.
	 *
	 * @param fitConfig
	 *            the fit config
	 * @return the solver name
	 */
	public static String getSolverName(FitConfiguration fitConfig)
	{
		final FitSolver solver = fitConfig.getFitSolver();
		String name = FitProtosHelper.getName(solver);
		if (solver == FitSolver.MLE)
			name += " " + FitProtosHelper.getName(fitConfig.getSearchMethod());
		return name;
	}

	/**
	 * Show results.
	 */
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
			for (final PeakResults r : results.toArray())
				if (r instanceof MemoryPeakResults)
				{
					if (((MemoryPeakResults) r).isSortAfterEnd())
						IJ.showStatus("Sorting " + r.size() + " results ...");
					break;
				}

			results.end();

			final String textTime = Utils.timeToString(time / 1000000.0);
			final String textRunTime = Utils.timeToString(runTime / 1000000.0);

			final int size = getSize();
			final String message = String.format("%s. Fitting Time = %s. Run time = %s",
					TextUtils.pleural(size, "localisation"), textTime, textRunTime);
			if (resultsSettings.getLogProgress())
				IJ.log("-=-=-=-");
			IJ.log(message);
			IJ.showStatus(message);
		}
		else
			IJ.showStatus("");
	}

	private static boolean showMaximaDialog()
	{
		final int size = MemoryPeakResults.countMemorySize();
		if (size == 0)
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return false;
		}

		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
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

	/**
	 * Gets the PSF type values.
	 *
	 * @return the PSF type values
	 */
	public static PSFType[] getPSFTypeValues()
	{
		if (_PSFTypeValues == null)
			initPSFType();
		return _PSFTypeValues;
	}

	private static String[] _PSFTypeNames;

	/**
	 * Gets the PSF type names.
	 *
	 * @return the PSF type names
	 */
	public static String[] getPSFTypeNames()
	{
		if (_PSFTypeNames == null)
			initPSFType();
		return _PSFTypeNames;
	}

	private static void initPSFType()
	{
		//@formatter:off
		final EnumSet<PSFType> d = EnumSet.of(
				PSFType.ONE_AXIS_GAUSSIAN_2D,
				PSFType.TWO_AXIS_GAUSSIAN_2D,
				PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D,
				PSFType.ASTIGMATIC_GAUSSIAN_2D);
		//@formatter:on
		_PSFTypeValues = d.toArray(new PSFType[d.size()]);
		_PSFTypeNames = new String[_PSFTypeValues.length];
		for (int i = 0; i < _PSFTypeValues.length; i++)
			_PSFTypeNames[i] = PSFProtosHelper.getName(_PSFTypeValues[i]);
	}

	private int showDialog(ImagePlus imp)
	{
		// Executing as an ImageJ plugin.

		// Load the settings
		resultsSettings = SettingsManager.readResultsSettings(0).toBuilder();
		// Settings are within the FitEngineSettings
		config = SettingsManager.readFitEngineConfiguration(0);
		fitConfig = config.getFitConfiguration();

		if (simpleFit)
			return showSimpleDialog();

		final boolean isCrop = (bounds != null && imp != null &&
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

		final String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);

		final CalibrationReader calibration = fitConfig.getCalibrationReader();
		addCameraOptions(gd, 0, fitConfig);
		gd.addNumericField("Calibration", calibration.getNmPerPixel(), 2, 6, "nm/px");
		gd.addNumericField("Exposure_time", calibration.getExposureTime(), 2, 6, "ms");

		if (isCrop)
			gd.addCheckbox("Ignore_bounds_for_noise", optionIgnoreBoundsForNoise);

		final FitConfigurationProvider fitConfigurationProvider = new FitConfigurationProvider()
		{
			@Override
			public FitConfiguration getFitConfiguration()
			{
				return fitConfig;
			}
		};

		final FitEngineConfigurationProvider fitEngineConfigurationProvider = new FitEngineConfigurationProvider()
		{
			@Override
			public FitEngineConfiguration getFitEngineConfiguration()
			{
				return config;
			}
		};

		addPSFOptions(gd, fitConfigurationProvider);
		addDataFilterOptions(gd, fitEngineConfigurationProvider);
		addSearchOptions(gd, fitEngineConfigurationProvider);
		addBorderOptions(gd, fitEngineConfigurationProvider);
		addFittingOptions(gd, fitEngineConfigurationProvider);
		if (extraOptions && !fitMaxima)
		{
			gd.addCheckbox("Interlaced_data", optionInterlacedData);
			gd.addSlider("Integrate_frames", 1, 5, optionIntegrateFrames);
		}

		if (!maximaIdentification)
		{
			gd.addMessage("--- Gaussian fitting ---");
			gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());
			if (extraOptions)
				gd.addCheckbox("Fit_background", fitConfig.isBackgroundFitting());

			// Parameters specific to each Fit solver are collected in a second dialog

			gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
			gd.addNumericField("Pass_rate", config.getPassRate(), 2);
			gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
			gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
			gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

			addDuplicateDistanceOptions(gd, fitEngineConfigurationProvider);

			gd.addMessage("--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");

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
			gd.addSlider("Width_factor", 1, 4.5, fitConfig.getMaxWidthFactor());
			addPrecisionOptions(gd, fitConfigurationProvider);
			// Q. Add dynamically displayed options for z-filtering here?
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

		// Add a mouse listener to the config file field
		if (Utils.isShowGenericDialog())
		{
			final Vector<TextField> texts = gd.getStringFields();
			final Vector<TextField> numerics = gd.getNumericFields();
			final Vector<Checkbox> checkboxes = gd.getCheckboxes();
			final Vector<Choice> choices = gd.getChoices();

			int n = 0;
			int t = 0;
			int b = 0;
			int ch = 0;

			final Choice textTemplate = choices.get(ch++);
			textTemplate.addItemListener(this);

			textCameraType = choices.get(ch++);
			textNmPerPixel = numerics.get(n++);
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
				textPassRate = numerics.get(n++);
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
				b++; // Skip over show processed frames option
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
			final int flags = IJ.setupDialog(imp, plugin_flags);

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
			setSource(new InterlacedImageSource(this.source, dataStart, dataBlock, dataSkip));

		// Allow frame aggregation by wrapping the image source
		if (integrateFrames > 1)
			setSource(new AggregatedImageSource(this.source, integrateFrames));

		// Ask if the user wants to log progress on multiple frame images
		if (resultsSettings.getLogProgress() && source.getFrames() > 1)
		{
			final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
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

	/**
	 * Allow the latest calibration to be provided for update
	 */
	public interface CalibrationProvider
	{
		/**
		 * Gets the calibration.
		 *
		 * @return the calibration
		 */
		public Calibration getCalibration();

		/**
		 * Save the calibration. This is used to save changes to the calibration back to the provider.
		 *
		 * @param calibration
		 *            the calibration
		 */
		public void saveCalibration(Calibration calibration);
	}

	/**
	 * Adds the camera options.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitConfig
	 *            the fit config
	 */
	public static void addCameraOptions(final ExtendedGenericDialog gd, final FitConfiguration fitConfig)
	{
		addCameraOptions(gd, 0, fitConfig);
	}

	/**
	 * Adds the camera options.
	 *
	 * @param gd
	 *            the dialog
	 * @param options
	 *            the options
	 * @param fitConfig
	 *            the fit config
	 */
	public static void addCameraOptions(final ExtendedGenericDialog gd, final int options,
			final FitConfiguration fitConfig)
	{
		addCameraOptions(gd, options, new CalibrationProvider()
		{
			@Override
			public Calibration getCalibration()
			{
				return fitConfig.getCalibration();
			}

			@Override
			public void saveCalibration(Calibration calibration)
			{
				fitConfig.setCalibration(calibration);
			}
		});
	}

	/**
	 * Adds the camera options.
	 *
	 * @param gd
	 *            the dialog
	 * @param calibrationWriter
	 *            the calibration writer
	 */
	public static void addCameraOptions(final ExtendedGenericDialog gd, final CalibrationWriter calibrationWriter)
	{
		addCameraOptions(gd, 0, calibrationWriter);
	}

	/**
	 * Adds the camera options.
	 *
	 * @param gd
	 *            the dialog
	 * @param options
	 *            the options
	 * @param calibrationWriter
	 *            the calibration writer
	 */
	public static void addCameraOptions(final ExtendedGenericDialog gd, final int options,
			final CalibrationWriter calibrationWriter)
	{
		addCameraOptions(gd, options, new PeakFit.CalibrationProvider()
		{
			@Override
			public Calibration getCalibration()
			{
				return calibrationWriter.getCalibration();
			}

			@Override
			public void saveCalibration(Calibration calibration)
			{
				calibrationWriter.setCalibration(calibration);
			}
		});
	}

	/** Flag to indicate that read noise should not be configured. */
	public static final int FLAG_NO_READ_NOISE = 0x00000001;
	/** Flag to indicate that quantum efficiency should be configured. */
	public static final int FLAG_QUANTUM_EFFICIENCY = 0x00000002;
	/** Flag to indicate that gain should not be configured. */
	public static final int FLAG_NO_GAIN = 0x00000004;

	/**
	 * Adds the camera options.
	 *
	 * @param gd
	 *            the dialog
	 * @param options
	 *            the options
	 * @param calibrationProvider
	 *            the calibration provider
	 */
	public static void addCameraOptions(final ExtendedGenericDialog gd, final int options,
			final CalibrationProvider calibrationProvider)
	{
		final CalibrationReader calibration = new CalibrationReader(calibrationProvider.getCalibration());

		gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
				CalibrationProtosHelper.getName(calibration.getCameraType()), new OptionListener<Integer>()
				{
					@Override
					public boolean collectOptions(Integer field)
					{
						final CalibrationWriter calibration = new CalibrationWriter(calibrationProvider.getCalibration());
						final CameraType t = SettingsManager.getCameraTypeValues()[field];
						if (calibration.getCameraType() != t)
						{
							calibration.setCameraType(t);
							calibrationProvider.saveCalibration(calibration.getCalibration());
						}
						final boolean result = collectOptions(false);
						return result;
					}

					@Override
					public boolean collectOptions()
					{
						return collectOptions(true);
					}

					private boolean collectOptions(boolean silent)
					{
						final CalibrationWriter calibration = new CalibrationWriter(calibrationProvider.getCalibration());
						final ExtendedGenericDialog egd = new ExtendedGenericDialog("Camera type options", null);
						if (calibration.isCCDCamera())
						{
							egd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "Count");
							if (BitFlags.anyNotSet(options, FLAG_NO_GAIN))
								egd.addNumericField("Gain", calibration.getCountPerPhoton(), 4, 6, "Count/photon");
							if (BitFlags.anyNotSet(options, FLAG_NO_READ_NOISE))
								egd.addNumericField("Read_noise", calibration.getReadNoise(), 4, 6, "Count");
							if (BitFlags.areSet(options, FLAG_QUANTUM_EFFICIENCY))
								egd.addNumericField("Quantum_efficiency", calibration.getQuantumEfficiency(), 4, 6,
										"electron/photon");
						}
						else if (calibration.isSCMOS())
						{
							final String[] models = CameraModelManager.listCameraModels(true);
							egd.addChoice("Camera_model_name", models, calibration.getCameraModelName());
							if (BitFlags.areSet(options, FLAG_QUANTUM_EFFICIENCY))
								egd.addNumericField("Quantum_efficiency", calibration.getQuantumEfficiency(), 4, 6,
										"electron/photon");
						}
						else
						{
							IJ.error("Unsupported camera type " +
									CalibrationProtosHelper.getName(calibration.getCameraType()));
							return false;
						}
						egd.setSilent(silent);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return false;
						final Calibration old = calibration.getCalibration();
						if (calibration.isCCDCamera())
						{
							calibration.setBias(Math.abs(egd.getNextNumber()));
							if (BitFlags.anyNotSet(options, FLAG_NO_GAIN))
								calibration.setCountPerPhoton(Math.abs(egd.getNextNumber()));
							if (BitFlags.anyNotSet(options, FLAG_NO_READ_NOISE))
								calibration.setReadNoise(Math.abs(egd.getNextNumber()));
							if (BitFlags.areSet(options, FLAG_QUANTUM_EFFICIENCY))
								calibration.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
						}
						else if (calibration.isSCMOS())
						{
							// Note: Since this does not go through the FitConfiguration object the
							// camera model is not invalidated. However any code using this function
							// should later call configureFitSolver(...) which will set the camera model
							// using the camera model name.
							calibration.setCameraModelName(egd.getNextChoice());
							if (BitFlags.areSet(options, FLAG_QUANTUM_EFFICIENCY))
								calibration.setQuantumEfficiency(Math.abs(egd.getNextNumber()));
						}
						final Calibration current = calibration.getCalibration();
						final boolean changed = !old.equals(current);
						if (changed)
							calibrationProvider.saveCalibration(current);
						return changed;
					}
				});
	}

	/**
	 * Allow the latest fitEngineConfiguration to be provided for update
	 */
	public interface FitEngineConfigurationProvider
	{
		/**
		 * Gets the fitEngineConfiguration.
		 *
		 * @return the fitEngineConfiguration
		 */
		public FitEngineConfiguration getFitEngineConfiguration();
	}

	/**
	 * Allow the latest fitConfiguration to be provided for update
	 */
	public interface FitConfigurationProvider
	{
		/**
		 * Gets the fitConfiguration.
		 *
		 * @return the fitConfiguration
		 */
		public FitConfiguration getFitConfiguration();
	}

	/**
	 * Simple implementation of {@link FitEngineConfigurationProvider}.
	 */
	public static class SimpleFitEngineConfigurationProvider implements FitEngineConfigurationProvider
	{
		/** The fit engine configuration. */
		FitEngineConfiguration c;

		/**
		 * Instantiates a new simple fit engine configuration provider.
		 *
		 * @param c
		 *            the configuration
		 */
		SimpleFitEngineConfigurationProvider(FitEngineConfiguration c)
		{
			this.c = c;
		}

		@Override
		public FitEngineConfiguration getFitEngineConfiguration()
		{
			return c;
		}
	}

	/**
	 * Simple implementation of {@link FitConfigurationProvider}.
	 */
	public static class SimpleFitConfigurationProvider implements FitConfigurationProvider
	{
		/** The configuration. */
		FitConfiguration c;

		/**
		 * Instantiates a new simple fit configuration provider.
		 *
		 * @param c
		 *            the configuration
		 */
		SimpleFitConfigurationProvider(FitConfiguration c)
		{
			this.c = c;
		}

		@Override
		public FitConfiguration getFitConfiguration()
		{
			return c;
		}
	}

	/**
	 * Adds the PSF options.
	 * <p>
	 * Note that if an astigmatic PSF is selected then the model must be created with
	 * {@link #configurePSFModel(FitEngineConfiguration, int)}.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitConfiguration
	 *            the fit configuration
	 */
	public static void addPSFOptions(final ExtendedGenericDialog gd, final FitConfiguration fitConfiguration)
	{
		addPSFOptions(gd, new SimpleFitConfigurationProvider(fitConfiguration));
	}

	/**
	 * Adds the PSF options.
	 * <p>
	 * Note that if an astigmatic PSF is selected then the model must be created with
	 * {@link #configurePSFModel(FitEngineConfiguration, int)}.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitConfigurationProvider
	 *            the fit configuration provider
	 */
	public static void addPSFOptions(final ExtendedGenericDialog gd,
			final FitConfigurationProvider fitConfigurationProvider)
	{
		final FitConfiguration fitConfig = fitConfigurationProvider.getFitConfiguration();
		gd.addChoice("PSF", getPSFTypeNames(), PSFProtosHelper.getName(fitConfig.getPSFType()),
				new OptionListener<Integer>()
				{
					@Override
					public boolean collectOptions(Integer field)
					{
						final FitConfiguration fitConfig = fitConfigurationProvider.getFitConfiguration();
						fitConfig.setPSFType(PeakFit.getPSFTypeValues()[field]);
						final boolean result = collectOptions(false);
						return result;
					}

					@Override
					public boolean collectOptions()
					{
						return collectOptions(true);
					}

					private boolean collectOptions(boolean silent)
					{
						final FitConfiguration fitConfig = fitConfigurationProvider.getFitConfiguration();
						final PSFType psfType = fitConfig.getPSFType();
						final ExtendedGenericDialog egd = new ExtendedGenericDialog("PSF Options", null);
						PSF oldPsf = null;
						if (psfType == PSFType.ASTIGMATIC_GAUSSIAN_2D)
						{
							// The PSF is entirely defined in the model
							String[] list = AstigmatismModelManager.listAstigmatismModels(false,
									fitConfig.getCalibrationReader().getNmPerPixel(), 0.1);
							// In case the calibration has not been updated
							if (list.length == 0)
								list = AstigmatismModelManager.listAstigmatismModels(false, true);
							egd.addChoice("Z-model", list, fitConfig.getPSFModelName());
						}
						else
						{
							// Collect the PSF parameters
							oldPsf = fitConfig.getPSF();
							for (int i = 0; i < oldPsf.getParametersCount(); i++)
							{
								final PSFParameter p = oldPsf.getParameters(i);
								egd.addNumericField(String.format("PSF_parameter_%d (%s)", i + 1, p.getName()),
										p.getValue(), 3);
							}
							if (psfType == PSFType.ONE_AXIS_GAUSSIAN_2D)
								egd.addCheckbox("Fixed", fitConfig.isFixedPSF());
						}

						egd.setSilent(silent);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return false;
						if (psfType == PSFType.ASTIGMATIC_GAUSSIAN_2D)
						{
							// The PSF is entirely defined in the model
							fitConfig.setPSFModelName(egd.getNextChoice());
							return true;
						}
						@SuppressWarnings("null")
						final
						PSF.Builder b = oldPsf.toBuilder();
						final int n = b.getParametersCount();
						for (int i = 0; i < n; i++)
							b.getParametersBuilder(i).setValue(egd.getNextNumber());
						final PSF newPsf = b.build();
						fitConfig.setPSF(newPsf);
						boolean changed = !oldPsf.equals(newPsf);
						if (psfType == PSFType.ONE_AXIS_GAUSSIAN_2D)
						{
							final boolean newFixed = egd.getNextBoolean();
							changed = changed || (newFixed != fitConfig.isFixedPSF());
							fitConfig.setFixedPSF(newFixed);
						}
						return changed;
					}
				});
	}

	/**
	 * Used for relative parameters
	 */
	static abstract class RelativeParameterProvider
	{
		/** The min. */
		double min;
		/** The max. */
		double max;
		/** The name. */
		String name;
		/** The fit engine configuration provider. */
		FitEngineConfigurationProvider fitEngineConfigurationProvider;
		/** Set to true to include the value in {@link #getMin()} and {@link #getMax()}. */
		boolean includeValue;

		/**
		 * Instantiates a new relative parameter provider.
		 *
		 * @param min
		 *            the min
		 * @param max
		 *            the max
		 * @param name
		 *            the name
		 * @param fitEngineConfigurationProvider
		 *            the fit engine configuration provider
		 */
		public RelativeParameterProvider(double min, double max, String name,
				FitEngineConfigurationProvider fitEngineConfigurationProvider)
		{
			this(min, max, name, fitEngineConfigurationProvider, false);
		}

		/**
		 * Instantiates a new relative parameter provider.
		 *
		 * @param min
		 *            the min
		 * @param max
		 *            the max
		 * @param name
		 *            the name
		 * @param fitEngineConfigurationProvider
		 *            the fit engine configuration provider
		 * @param includeValue
		 *            the include value
		 */
		public RelativeParameterProvider(double min, double max, String name,
				FitEngineConfigurationProvider fitEngineConfigurationProvider, boolean includeValue)
		{
			this.min = min;
			this.max = max;
			this.name = name;
			this.fitEngineConfigurationProvider = fitEngineConfigurationProvider;
			this.includeValue = includeValue;
		}

		/**
		 * Gets the min.
		 *
		 * @return the min
		 */
		double getMin()
		{
			return (includeValue) ? Math.min(min, getValue()) : min;
		}

		/**
		 * Gets the max.
		 *
		 * @return the max
		 */
		double getMax()
		{
			return (includeValue) ? Math.max(max, getValue()) : max;
		}

		/**
		 * Gets the value.
		 *
		 * @return the value
		 */
		abstract double getValue();

		/**
		 * Checks if is absolute.
		 *
		 * @return true, if is absolute
		 */
		abstract boolean isAbsolute();

		/**
		 * Sets the absolute.
		 *
		 * @param absolute
		 *            the new absolute
		 */
		abstract void setAbsolute(boolean absolute);

		/**
		 * Gets the dialog name.
		 *
		 * @return the dialog name
		 */
		String getDialogName()
		{
			return name.replace(' ', '_');
		}
	}

	/**
	 * Adds the relative parameter options.
	 *
	 * @param gd
	 *            the dialog
	 * @param rp
	 *            the relative parameter
	 */
	static void addRelativeParameterOptions(final ExtendedGenericDialog gd, final RelativeParameterProvider rp)
	{
		final String label = rp.getDialogName();
		gd.addSlider(label, rp.getMin(), rp.getMax(), rp.getValue(), new OptionListener<Double>()
		{
			@Override
			public boolean collectOptions(Double value)
			{
				// Nothing depends on the input double value so just collect the options
				return collectOptions(false);
			}

			@Override
			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				final ExtendedGenericDialog egd = new ExtendedGenericDialog(rp.name + " Options", null);
				final boolean oldValue = rp.isAbsolute();
				egd.addCheckbox(rp.getDialogName() + "_absolute", oldValue);
				egd.setSilent(silent);
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				final boolean newValue = egd.getNextBoolean();
				rp.setAbsolute(newValue);
				return oldValue != newValue;
			}
		});

		// Add a label after the button.
		// The button is added in a panel with a GridBagLayout.
		final Panel p = gd.getLastPanel();
		final GridBagConstraints pc = new GridBagConstraints();
		pc.gridy = 0;
		pc.gridx = 3;
		pc.insets = new Insets(5, 5, 0, 0);
		pc.anchor = GridBagConstraints.EAST;
		final Label flagLabel = new Label();
		updateFlag(flagLabel, rp.isAbsolute());
		p.add(flagLabel, pc);

		gd.addOptionCollectedListener(new OptionCollectedListener()
		{
			@Override
			public void optionCollected(OptionCollectedEvent e)
			{
				if (label.equals(e.getLabel()))
					updateFlag(flagLabel, rp.isAbsolute());
			}
		});
	}

	private static void updateFlag(Label label, boolean absolute)
	{
		// Note: Add spaces so the label has extra width for font width differences
		if (absolute)
			label.setText("Absolute   ");
		else
			label.setText("Relative   ");
	}

	/**
	 * Adds the data filter options for the first filter. Adds to the dialog:
	 * <ul>
	 * <li>a choice of filter type (e.g. single, difference, etc)</li>
	 * <li>a choice of primary filter (e.g. mean, Gaussian, etc)</li>
	 * <li>a single slider for the primary filter parameter</li>
	 * </ul>
	 *
	 * @param gd
	 *            the dialog
	 * @param fitEngineConfigurationProvider
	 *            the fit engine configuration provider
	 */
	public static void addDataFilterOptions(final ExtendedGenericDialog gd,
			final FitEngineConfigurationProvider fitEngineConfigurationProvider)
	{
		final int n = 0;
		final FitEngineConfiguration config = fitEngineConfigurationProvider.getFitEngineConfiguration();
		gd.addChoice("Spot_filter_type", SettingsManager.getDataFilterTypeNames(),
				config.getDataFilterType().ordinal());
		gd.addChoice("Spot_filter", SettingsManager.getDataFilterMethodNames(),
				config.getDataFilterMethod(n).ordinal());
		addRelativeParameterOptions(gd,
				new RelativeParameterProvider(0, 2.5, "Smoothing", fitEngineConfigurationProvider, true)
				{
					@Override
					void setAbsolute(boolean absolute)
					{
						final FitEngineConfiguration c = fitEngineConfigurationProvider.getFitEngineConfiguration();
						final DataFilterMethod m = c.getDataFilterMethod(n);
						final double smooth = c.getDataFilterParameter(n).getValue();
						c.setDataFilter(m, smooth, absolute, n);
					}

					@Override
					boolean isAbsolute()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration()
								.getDataFilterParameterAbsolute(n);
					}

					@Override
					double getValue()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getDataFilterParameter(n)
								.getValue();
					}
				});
	}

	/**
	 * Adds the search options. A single slider for the search parameter is added to the dialog.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitEngineConfigurationProvider
	 *            the fit engine configuration provider
	 */
	public static void addSearchOptions(final ExtendedGenericDialog gd,
			final FitEngineConfigurationProvider fitEngineConfigurationProvider)
	{
		addRelativeParameterOptions(gd,
				new RelativeParameterProvider(0.5, 2.5, "Search Width", fitEngineConfigurationProvider, true)
				{
					@Override
					void setAbsolute(boolean absolute)
					{
						fitEngineConfigurationProvider.getFitEngineConfiguration().setSearchAbsolute(absolute);
					}

					@Override
					boolean isAbsolute()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getSearchAbsolute();
					}

					@Override
					double getValue()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getSearch();
					}
				});
	}

	/**
	 * Adds the border options. A single slider for the border parameter is added to the dialog.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitEngineConfigurationProvider
	 *            the fit engine configuration provider
	 */
	public static void addBorderOptions(final ExtendedGenericDialog gd,
			final FitEngineConfigurationProvider fitEngineConfigurationProvider)
	{
		addRelativeParameterOptions(gd,
				new RelativeParameterProvider(0.5, 2.5, "Border Width", fitEngineConfigurationProvider, true)
				{
					@Override
					void setAbsolute(boolean absolute)
					{
						fitEngineConfigurationProvider.getFitEngineConfiguration().setBorderAbsolute(absolute);
					}

					@Override
					boolean isAbsolute()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getBorderAbsolute();
					}

					@Override
					double getValue()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getBorder();
					}
				});
	}

	/**
	 * Adds the fitting options. A single slider for the fitting parameter is added to the dialog.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitEngineConfigurationProvider
	 *            the fit engine configuration provider
	 */
	public static void addFittingOptions(final ExtendedGenericDialog gd,
			final FitEngineConfigurationProvider fitEngineConfigurationProvider)
	{
		// For this we allow the slider range to increase as the user may have a large fit width
		addRelativeParameterOptions(gd,
				new RelativeParameterProvider(2, 4.5, "Fitting Width", fitEngineConfigurationProvider, true)
				{
					@Override
					void setAbsolute(boolean absolute)
					{
						fitEngineConfigurationProvider.getFitEngineConfiguration().setFittingAbsolute(absolute);
					}

					@Override
					boolean isAbsolute()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getFittingAbsolute();
					}

					@Override
					double getValue()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getFitting();
					}
				});
	}

	/**
	 * Adds the duplicate distance options. A single slider for the duplicate distance parameter is added to the dialog.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitEngineConfigurationProvider
	 *            the fit engine configuration provider
	 */
	public static void addDuplicateDistanceOptions(final ExtendedGenericDialog gd,
			final FitEngineConfigurationProvider fitEngineConfigurationProvider)
	{
		addRelativeParameterOptions(gd,
				new RelativeParameterProvider(0, 1.5, "Duplicate Distance", fitEngineConfigurationProvider)
				{
					@Override
					void setAbsolute(boolean absolute)
					{
						fitEngineConfigurationProvider.getFitEngineConfiguration()
								.setDuplicateDistanceAbsolute(absolute);
					}

					@Override
					boolean isAbsolute()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration()
								.getDuplicateDistanceAbsolute();
					}

					@Override
					double getValue()
					{
						return fitEngineConfigurationProvider.getFitEngineConfiguration().getDuplicateDistance();
					}
				});
	}

	/**
	 * Adds the precision options. A single numeric field for the precision is added. A pop-up is added to allow the
	 * precision method to be configured.
	 *
	 * @param gd
	 *            the dialog
	 * @param fitConfigurationProvider
	 *            the fit configuration provider
	 */
	public static void addPrecisionOptions(final ExtendedGenericDialog gd,
			final FitConfigurationProvider fitConfigurationProvider)
	{
		gd.addNumericField("Precision", fitConfigurationProvider.getFitConfiguration().getPrecisionThreshold(), 2,
				new OptionListener<Double>()
				{
					@Override
					public boolean collectOptions(Double field)
					{
						final FitConfiguration fitConfig = fitConfigurationProvider.getFitConfiguration();
						fitConfig.setPrecisionThreshold(field);
						final boolean result = collectOptions(false);
						return result;
					}

					@Override
					public boolean collectOptions()
					{
						return collectOptions(true);
					}

					private boolean collectOptions(boolean silent)
					{
						final FitConfiguration fitConfig = fitConfigurationProvider.getFitConfiguration();
						final ExtendedGenericDialog egd = new ExtendedGenericDialog("Precision Options", null);
						final int oldIndex = fitConfig.getPrecisionMethod().ordinal();
						egd.addChoice("Precision_method", SettingsManager.getPrecisionMethodNames(), oldIndex);
						egd.setSilent(silent);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return false;
						final int newIndex = egd.getNextChoiceIndex();
						fitConfig.setPrecisionMethod(newIndex);
						final boolean changed = oldIndex != newIndex;
						return changed;
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
		// Just support circular fitting
		fitConfig.setPSF(PSFProtosHelper.defaultOneAxisGaussian2DPSF);
		fitConfig.setFixedPSF(false);

		// TODO - Support sCMOS camera. This may be 'too difficult' as the
		// user will need to have created a per-pixel calibration image

		final CalibrationWriter calibration = fitConfig.getCalibrationWriter();
		final boolean requireCalibration = requireCalibration(calibration);
		if (requireCalibration)
			if (!showCalibrationWizard(calibration, true))
				return DONE;

		// Present dialog with simple output options: Image, Table
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
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
			if (!showCalibrationWizard(calibration, false))
				return DONE;

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
			final ResultsTableSettings.Builder tableSettings = resultsSettings.getResultsTableSettingsBuilder();
			tableSettings.setShowTable(true);
		}
		if (showImage)
		{
			final ResultsImageSettings.Builder imageSettings = resultsSettings.getResultsImageSettingsBuilder();
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
		Utils.log("PSF width = %s", Utils.rounded(fitConfig.getInitialXSD(), 4));

		// Save
		fitConfig.setCalibration(calibration.getCalibration());
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
		return SettingsManager.writeSettings(config, 0);
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
		// Check for a supported camera
		if (!calibration.isCCDCamera())
			return true;
		// Check if the calibration contains: Pixel pitch, Gain (can be 1), Exposure time
		if (!calibration.hasNmPerPixel())
			return true;
		// Bias can be zero (but this is unlikely)
		if (!calibration.hasCountPerPhoton() || !(calibration.getBias() > 0))
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
			final ExtendedGenericDialog gd = newWizardDialog("No configuration file could be loaded.",
					"Please follow the configuration wizard to calibrate.");
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
		}

		if (!getCameraType(calibration))
			return false;
		if (!getPixelPitch(calibration))
			return false;
		if (!getGain(calibration))
			return false;
		if (!getExposureTime(calibration))
			return false;
		if (!getPeakWidth(calibration))
			return false;

		// Check parameters
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.getNmPerPixel());
			// We allow bias to be zero
			Parameters.isAboveZero("Gain", calibration.getCountPerPhoton());
			Parameters.isAboveZero("Exposure time", calibration.getExposureTime());
			Parameters.isAboveZero("Initial SD", fitConfig.getInitialXSD());
		}
		catch (final IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		return true;
	}

	private static ExtendedGenericDialog newWizardDialog(String... messages)
	{
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		final String header = "-=-";
		gd.addMessage(header + " " + TITLE + " Configuration Wizard " + header);
		for (final String message : messages)
			gd.addMessage(TextUtils.wrap(message, 80));
		return gd;
	}

	private static boolean getCameraType(CalibrationWriter calibration)
	{
		final ExtendedGenericDialog gd = newWizardDialog("Enter the type of camera.");
		gd.addChoice("Camera_type", SettingsManager.getCameraTypeNames(),
				CalibrationProtosHelper.getName(calibration.getCameraType()));
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
		if (!calibration.isCCDCamera())
		{
			// TODO - Support sCMOS camera

			IJ.error("Unsupported camera type " + CalibrationProtosHelper.getName(calibration.getCameraType()));
			return false;
		}
		return true;
	}

	private static boolean getPixelPitch(CalibrationWriter calibration)
	{
		final ExtendedGenericDialog gd = newWizardDialog(
				"Enter the size of each pixel. This is required to ensure the dimensions of the image are calibrated.",
				"E.g. a camera with a 6.45um pixel size and a 60x objective will have a pitch of 6450/60 = 107.5nm.");
		// TODO - Add a pop-up calculator...
		gd.addNumericField("Calibration", calibration.getNmPerPixel(), 2, 6, "nm/px");
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
		return true;
	}

	private static boolean getGain(CalibrationWriter calibration)
	{
		final ExtendedGenericDialog gd = newWizardDialog("Enter the bias and total gain.",
				"This is usually supplied with your camera certificate. The bias is a fixed offset added to the camera counts. The gain indicates how many Analogue-to-Digital-Units (Count) are recorded at the pixel for each photon registered on the sensor.",
				"The gain is usually expressed using the product of the EM-gain (if applicable), the camera gain and the sensor quantum efficiency.",
				"A value of 1 means no conversion to photons will occur.");
		// TODO - Add a wizard to allow calculation of total gain from EM-gain, camera gain and QE
		gd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "Count");
		gd.addNumericField("Gain", calibration.getCountPerPhoton(), 2, 6, "Count/photon");
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.setBias(Math.abs(gd.getNextNumber()));
		calibration.setCountPerPhoton(Math.abs(gd.getNextNumber()));
		return true;
	}

	private static boolean getExposureTime(CalibrationWriter calibration)
	{
		final ExtendedGenericDialog gd = newWizardDialog(
				"Enter the exposure time. Calibration of the exposure time allows correct reporting of on and off times.",
				"This is the length of time for each frame in the image.");
		gd.addNumericField("Exposure_time", calibration.getExposureTime(), 2, 6, "ms");
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		calibration.setExposureTime(Math.abs(gd.getNextNumber()));
		return true;
	}

	private boolean getPeakWidth(final CalibrationWriter calibration)
	{
		final ExtendedGenericDialog gd = newWizardDialog("Enter the expected peak width in pixels.",
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
				@Override
				public void actionPerformed(ActionEvent e)
				{
					// Run the PSF Calculator
					final PSFCalculator calculator = new PSFCalculator();
					calculatorSettings.setPixelPitch(calibration.getNmPerPixel() / 1000.0);
					calculatorSettings.setMagnification(1);
					calculatorSettings.setBeamExpander(1);
					final double sd = calculator.calculate(calculatorSettings.build(), true);
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

	@Override
	public void itemStateChanged(ItemEvent e)
	{
		if (e.getSource() instanceof Choice)
		{
			// Update the settings from the template
			final Choice choice = (Choice) e.getSource();
			final String templateName = choice.getSelectedItem();
			//System.out.println("Update to " + templateName);

			// Get the configuration template
			final TemplateSettings template = ConfigurationTemplate.getTemplate(templateName);

			if (template != null)
			{
				IJ.log("Applying template: " + templateName);

				for (final String note : template.getNotesList())
					IJ.log(note);

				final boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
				if (template.hasCalibration())
					refreshSettings(template.getCalibration());
				if (template.hasPsf())
					refreshSettings(template.getPsf(), custom);
				if (template.hasFitEngineSettings())
					refreshSettings(template.getFitEngineSettings(), custom);
				if (template.hasResultsSettings())
					refreshSettings(template.getResultsSettings());
			}
		}
		else if (e.getSource() instanceof Checkbox)
			if (e.getSource() == textSmartFilter)
			{
				// Prevent both filters being enabled
				textDisableSimpleFilter.setState(textSmartFilter.getState());
				updateFilterInput();
			}
			else if (e.getSource() == textDisableSimpleFilter)
				updateFilterInput();
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

	private static void disableEditing(TextField textField)
	{
		textField.setEditable(false);
		textField.setBackground(SystemColor.control);
	}

	private static void enableEditing(TextField textField)
	{
		textField.setEditable(true);
		textField.setBackground(Color.white);
	}

	private boolean readDialog(ExtendedGenericDialog gd, boolean isCrop)
	{
		// Ignore the template
		gd.getNextChoice();

		final CalibrationWriter calibration = fitConfig.getCalibrationWriter();
		calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
		calibration.setNmPerPixel(Math.abs(gd.getNextNumber()));
		calibration.setExposureTime(Math.abs(gd.getNextNumber()));
		fitConfig.setCalibration(calibration.getCalibration());

		// Note: The bias and read noise will just end up being what was in the configuration file
		// One fix for this is to save/load only the settings that are required from the configuration file
		// (the others will remain unchanged). This will require a big refactor of the settings save/load.
		// The simple fix is to create a plugin to allow the configuration to be changed for results.
		if (isCrop)
			ignoreBoundsForNoise = optionIgnoreBoundsForNoise = gd.getNextBoolean();

		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		config.setDataFilterType(gd.getNextChoiceIndex());
		// Note: The absolute flag is set in extra options
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
			if (extraOptions)
				fitConfig.setBackgroundFitting(gd.getNextBoolean());
			config.setFailuresLimit((int) gd.getNextNumber());
			config.setPassRate(gd.getNextNumber());
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
		resultsSettings.getResultsFileSettingsBuilder().setFileFormatValue(gd.getNextChoiceIndex());
		resultsSettings.getResultsFileSettingsBuilder().setResultsDirectory(gd.getNextString());
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
			// No check on camera calibration. This is left to the FitConfiguration to
			// error if the settings are incorrect

			Parameters.isAboveZero("nm per pixel", calibration.getNmPerPixel());
			Parameters.isAboveZero("Exposure time", calibration.getExposureTime());
			if (fitConfig.getPSFTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE)
			{
				Parameters.isAboveZero("Initial SD0", fitConfig.getInitialXSD());
				if (fitConfig.getPSF().getParametersCount() > 1)
					Parameters.isAboveZero("Initial SD1", fitConfig.getInitialYSD());
			}
			Parameters.isAboveZero("Search_width", config.getSearch());
			Parameters.isAboveZero("Fitting_width", config.getFitting());
			Parameters.isPositive("Integrate frames", integrateFrames);
			if (!maximaIdentification)
			{
				// This can be negative to disable, i.e. fit everything
				//Parameters.isPositive("Failures limit", config.getFailuresLimit());
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
					if (fitConfig.getPrecisionThreshold() > 0)
					{
						if (fitConfig.getPrecisionMethod() == PrecisionMethod.PRECISION_METHOD_NA)
							throw new IllegalArgumentException("Precision filter requires a precision method");
						if (fitConfig.isPrecisionUsingBackground() && calibration.isCCDCamera() &&
								(calibration.getBias() == 0 || !calibration.hasCountPerPhoton()))
							throw new IllegalArgumentException(
									"Precision using the local background requires the camera bias");
					}
				}
			}
			final ResultsImageSettings.Builder imageSettings = resultsSettings.getResultsImageSettingsBuilder();
			if (imageSettings.getImageType() == ResultsImageType.DRAW_INTENSITY_AVERAGE_PRECISION ||
					imageSettings.getImageType() == ResultsImageType.DRAW_LOCALISATIONS_AVERAGE_PRECISION)
				Parameters.isAboveZero("Image precision", imageSettings.getAveragePrecision());
			Parameters.isAboveZero("Image scale", imageSettings.getScale());
			if (extraOptions)
				Parameters.isPositive("Image rolling window", imageSettings.getRollingWindowSize());
		}
		catch (final IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		final int flags = (extraOptions) ? FLAG_EXTRA_OPTIONS : 0;

		// If precision filtering then we need the camera bias
		if (!maximaIdentification)
		{
			if (!configurePSFModel(config, flags))
				return false;

			if (!configureResultsFilter(config, flags))
				return false;
		}

		if (!configureDataFilter(config, flags))
			return false;

		// Second dialog for solver dependent parameters
		if (!maximaIdentification)
			if (!configureFitSolver(config, source.getBounds(), bounds, flags))
				return false;

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
				interlacedData = false;
		}

		final boolean result = saveFitEngineSettings();
		if (!result)
			IJ.error(TITLE, "Failed to save settings");

		return result;
	}

	/**
	 * Show a dialog to configure the PSF model. The updated settings are saved to the settings file.
	 * <p>
	 * If the configuration is for a 3D PSF then a dialog to configure the z model is shown.
	 *
	 * @param config
	 *            the config
	 * @return true, if successful
	 */
	public static boolean configurePSFModel(FitEngineConfiguration config)
	{
		return configurePSFModel(config, FLAG_NO_SAVE);
	}

	/**
	 * Show a dialog to configure the PSF model. The updated settings are saved to the settings file.
	 * <p>
	 * If the configuration is for a 3D PSF then a dialog to configure the z model is shown.
	 *
	 * @param config
	 *            the config
	 * @param flags
	 *            the flags
	 * @return true, if successful
	 */
	public static boolean configurePSFModel(FitEngineConfiguration config, int flags)
	{
		final FitConfiguration fitConfig = config.getFitConfiguration();
		if (fitConfig.getPSFTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE)
			return true;

		// Get the astigmatism z-model
		final AstigmatismModel model = AstigmatismModelManager.getModel(fitConfig.getPSFModelName());
		if (model == null)
		{
			IJ.error(TITLE, "Failed to load the model: " + fitConfig.getPSFModelName());
			return false;
		}

		// Conversion to the correct units in pixels is done within the FitConfiguration object.
		fitConfig.setAstigmatismModel(model);

		if (BitFlags.anyNotSet(flags, FLAG_NO_SAVE))
			SettingsManager.writeSettings(config, 0);

		return true;
	}

	/**
	 * Show a dialog to configure the results filter. The updated settings are saved to the settings file.
	 * <p>
	 * If the configuration is for a 3D PSF then a dialog to configure the z range for the results is shown (see
	 * {@link #configureZFilter(FitEngineConfiguration, int)}).
	 * <p>
	 * If the configuration is for a smart filter then a dialog to configure the smart filter is shown (see
	 * {@link #configureSmartFilter(FitEngineConfiguration, int)}).
	 *
	 * @param config
	 *            the config
	 * @param flags
	 *            the flags
	 * @return true, if successful
	 */
	public static boolean configureResultsFilter(FitEngineConfiguration config, int flags)
	{
		boolean result = configureZFilter(config, flags);
		result = result && configureSmartFilter(config, flags);
		return result;
	}

	/**
	 * Show a dialog to configure the results z filter. The updated settings are saved to the settings file.
	 * <p>
	 * If the fit configuration PSF is not 3D or the simple filter is disabled then this method returns true. If it is
	 * enabled then a dialog is shown to input the configuration for the z filter.
	 * <p>
	 * Note: The PSF and any z-model must be correctly configured for fitting in pixel units.
	 *
	 * @param config
	 *            the config
	 * @param flags
	 *            the flags
	 * @return true, if successful
	 */
	public static boolean configureZFilter(FitEngineConfiguration config, int flags)
	{
		final FitConfiguration fitConfig = config.getFitConfiguration();
		if (fitConfig.isDisableSimpleFilter() || !fitConfig.is3D())
			return true;

		// Create a converter to map the model units in pixels to nm for the dialog.
		// Note the output units of pixels may not yet be set in the calibration so we assume it is pixels.
		//TypeConverter<DistanceUnit> c = fitConfig.getCalibrationReader().getDistanceConverter(DistanceUnit.NM);
		final TypeConverter<DistanceUnit> c = UnitConverterFactory.createConverter(DistanceUnit.PIXEL, DistanceUnit.NM,
				fitConfig.getCalibrationReader().getNmPerPixel());

		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

		gd.addMessage("3D filter");
		gd.addNumericField("Min_z", c.convert(fitConfig.getMinZ()), 0, 6, "nm");
		gd.addNumericField("Max_z", c.convert(fitConfig.getMaxZ()), 0, 6, "nm");

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		final double minZ = gd.getNextNumber();
		final double maxZ = gd.getNextNumber();

		if (gd.invalidNumber() || minZ > maxZ)
		{
			IJ.error(TITLE, "Min Z must be equal or below the max Z");
			return false;
		}

		// Map back
		fitConfig.setMinZ(c.convertBack(minZ));
		fitConfig.setMaxZ(c.convertBack(maxZ));

		if (BitFlags.anyNotSet(flags, FLAG_NO_SAVE))
			SettingsManager.writeSettings(config, 0);

		return true;
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
		final FitConfiguration fitConfig = config.getFitConfiguration();
		if (!fitConfig.isSmartFilter())
			return true;

		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

		String xml = fitConfig.getSmartFilterString();
		if (TextUtils.isNullOrEmpty(xml))
			xml = fitConfig.getDefaultSmartFilterXML();

		gd.addMessage("Smart filter (used to pick optimum results during fitting)");
		gd.addTextAreas(gdsc.core.utils.XmlUtils.convertQuotes(xml), null, 8, 60);
		// Add message about precision filtering
		gd.addMessage(TextUtils.wrap("Note: Smart filters using precision may require a local background level. " +
				"Ensure the camera calibration is correct including any bias.", 80));

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		xml = gd.getNextText();
		final Filter f = Filter.fromXML(xml);
		if (f == null || !(f instanceof DirectFilter))
			return false;

		fitConfig.setDirectFilter((DirectFilter) f);

		if (BitFlags.anyNotSet(flags, FLAG_NO_SAVE))
			SettingsManager.writeSettings(config, 0);
		return true;
	}

	/**
	 * Show a dialog to configure the data filter. The data filter type and the first data filter must ALREADY be set in
	 * the configuration. The subsequent filters are then configured, e.g. for difference and jury filters.
	 * <p>
	 * The updated settings are saved to the settings file. An error message is shown if the dialog is cancelled or the
	 * configuration is invalid.
	 * <p>
	 * If the configuration is for a per-pixel camera type (e.g. sCMOS) then the camera model will be loaded using the
	 * configured camera model name. This will be used to validate the filter to check the filter supports the per-pixel
	 * camera type.
	 *
	 * @param config
	 *            the config
	 * @param flags
	 *            the flags
	 * @return True if the configuration succeeded
	 */
	public static boolean configureDataFilter(final FitEngineConfiguration config, int flags)
	{
		int numberOfFilters = 1;
		int n;
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

		final String[] filterNames = SettingsManager.getDataFilterMethodNames();
		final DataFilterMethod[] filterValues = SettingsManager.getDataFilterMethodValues();

		// We use the previous value in the event the configuration does not have any current values.
		// Check we have at least the first filter.
		if (config.getDataFiltersCount() == 0)
			throw new IllegalStateException("No primary filter is configured");

		final FitEngineConfigurationProvider fitEngineConfigurationProvider = new FitEngineConfigurationProvider()
		{
			@Override
			public FitEngineConfiguration getFitEngineConfiguration()
			{
				return config;
			}
		};

		for (int i = 1; i < n; i++)
		{
			final int filter = i + 1;
			final int ii = i;

			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			if (filter == n)
				// This is maximum filter count so no continue option
				gd.addMessage(
						String.format("Configure the %s filter.", FitProtosHelper.getName(config.getDataFilterType())));
			else
			{
				gd.enableYesNoCancel("Add", "Continue");
				gd.addMessage(
						String.format("Configure the %s filter.\nClick continue to proceed with the current set of %d.",
								FitProtosHelper.getName(config.getDataFilterType()), i));
			}

			final String fieldName = "Spot_filter" + filter;
			if (IJ.isMacro())
				// Use blank default value so bad macro parameters return nothing
				gd.addStringField(fieldName, "");
			else
				gd.addChoice(fieldName, filterNames,
						filterNames[config.getDataFilterMethod(ii, config.getDataFilterMethod(ii - 1)).ordinal()]);
			addRelativeParameterOptions(gd,
					new RelativeParameterProvider(0, 4.5, "Smoothing" + filter, fitEngineConfigurationProvider, true)
					{
						@Override
						void setAbsolute(boolean absolute)
						{
							// Get the current settings
							final FitEngineConfiguration c = fitEngineConfigurationProvider.getFitEngineConfiguration();
							final DataFilterMethod m = c.getDataFilterMethod(ii);
							final double smooth = c.getDataFilterParameter(ii).getValue();
							// Reset with the new absolute value
							c.setDataFilter(m, smooth, absolute, ii);
						}

						@Override
						boolean isAbsolute()
						{
							final FitEngineConfiguration c = fitEngineConfigurationProvider.getFitEngineConfiguration();
							return c.getDataFilterParameterAbsolute(ii, c.getDataFilterParameterAbsolute(ii - 1));
						}

						@Override
						double getValue()
						{
							final FitEngineConfiguration c = fitEngineConfigurationProvider.getFitEngineConfiguration();
							return c.getDataFilterParameterValue(ii, c.getDataFilterParameterValue(ii - 1));
						}
					});
			//gd.addSlider("Smoothing" + filter, 0, 4.5, config.getDataFilterParameterValue(i, defaultSmooth));
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			if (gd.wasOKed())
			{
				int filterIndex = -1;
				if (IJ.isMacro())
				{
					final String filterName = gd.getNextString();
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
				// Note: The absolute flag is set in extra options
				config.setDataFilter(filterValues[filterIndex], Math.abs(gd.getNextNumber()), i);
				gd.collectOptions();
				numberOfFilters++;
			}
			else
				break;
		}
		config.setNumberOfFilters(numberOfFilters);

		if (BitFlags.anyNotSet(flags, FLAG_NO_SAVE))
			saveFitEngineSettings(config);

		final FitConfiguration fitConfig = config.getFitConfiguration();
		final CalibrationReader calibration = fitConfig.getCalibrationReader();
		if (calibration.isSCMOS())
			fitConfig.setCameraModel(CameraModelManager.load(fitConfig.getCameraModelName()));

		try
		{
			config.createSpotFilter();
		}
		catch (final Exception e) // IllegalStateException
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

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
	 * <p>
	 * The bounds are used to validate the camera model. The camera model must be large enough to cover the source
	 * bounds. If larger then it will be cropped. Optionally an internal region of the input image can be specifed. This
	 * is relative to the width and height of the input image. If no camera model is present then the bounds can be
	 * null.
	 *
	 * @param config
	 *            the config
	 * @param sourceBounds
	 *            the source image bounds (used to validate the camera model dimensions)
	 * @param bounds
	 *            the crop bounds (relative to the input image, used to validate the camera model dimensions)
	 * @param flags
	 *            the flags
	 * @return True if the configuration succeeded
	 */
	public static boolean configureFitSolver(FitEngineConfiguration config, Rectangle sourceBounds, Rectangle bounds,
			int flags)
	{
		final boolean extraOptions = BitFlags.anySet(flags, FLAG_EXTRA_OPTIONS);
		final boolean ignoreCalibration = BitFlags.anySet(flags, FLAG_IGNORE_CALIBRATION);
		final boolean saveSettings = BitFlags.anyNotSet(flags, FLAG_NO_SAVE);

		final FitConfiguration fitConfig = config.getFitConfiguration();
		final CalibrationWriter calibration = fitConfig.getCalibrationWriter();

		final FitSolver fitSolver = fitConfig.getFitSolver();

		final boolean isLVM = fitSolver == FitSolver.LVM_LSE || fitSolver == FitSolver.LVM_WLSE ||
				fitSolver == FitSolver.LVM_MLE;
		final boolean isFastMLE = fitSolver == FitSolver.FAST_MLE || fitSolver == FitSolver.BACKTRACKING_FAST_MLE;
		final boolean isSteppingFunctionSolver = isLVM || isFastMLE;

		if (fitSolver == FitSolver.MLE)
		{
			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			if (!ignoreCalibration)
			{
				gd.addMessage("Maximum Likelihood Estimation requires CCD-type camera parameters");
				gd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "count");
				gd.addCheckbox("Model_camera_noise", fitConfig.isModelCamera());
				gd.addNumericField("Read_noise", calibration.getReadNoise(), 2, 6, "count");
				gd.addNumericField("Quantum_efficiency", calibration.getQuantumEfficiency(), 2, 6, "electron/photon");
				gd.addCheckbox("EM-CCD", calibration.isEMCCD());
			}
			else
				gd.addMessage("Maximum Likelihood Estimation requires additional parameters");
			// This works because the proto configuration enum matches the named enum
			final String[] searchNames = SettingsManager.getNames((Object[]) MaximumLikelihoodFitter.SearchMethod.values());
			gd.addChoice("Search_method", searchNames, searchNames[fitConfig.getSearchMethod().getNumber()]);
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
				calibration.setQuantumEfficiency(Math.abs(gd.getNextNumber()));
				calibration.setCameraType((gd.getNextBoolean()) ? CameraType.EMCCD : CameraType.CCD);
				fitConfig.setCalibration(calibration.getCalibration());
			}
			fitConfig.setSearchMethod(gd.getNextChoiceIndex());
			fitConfig.setRelativeThreshold(getThresholdNumber(gd));
			fitConfig.setAbsoluteThreshold(getThresholdNumber(gd));
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
			catch (final Exception e) // IllegalArgumentException, IllegalStateException
			{
				IJ.error(TITLE, e.getMessage());
				return false;
			}
		}
		else if (isSteppingFunctionSolver)
		{
			final boolean requireCalibration = !ignoreCalibration && fitSolver != FitSolver.LVM_LSE;

			// Collect options for LVM fitting
			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
			final String fitSolverName = FitProtosHelper.getName(fitSolver);
			gd.addMessage(fitSolverName + " requires additional parameters");
			gd.addStringField("Relative_threshold", Utils.rounded(fitConfig.getRelativeThreshold()));
			gd.addStringField("Absolute_threshold", Utils.rounded(fitConfig.getAbsoluteThreshold()));
			gd.addStringField("Parameter_relative_threshold", Utils.rounded(fitConfig.getParameterRelativeThreshold()));
			gd.addStringField("Parameter_absolute_threshold", Utils.rounded(fitConfig.getParameterAbsoluteThreshold()));
			gd.addNumericField("Max_iterations", fitConfig.getMaxIterations(), 0);
			if (isLVM)
				gd.addNumericField("Lambda", fitConfig.getLambda(), 4);
			if (isFastMLE)
			{
				gd.addCheckbox("Fixed_iterations", fitConfig.isFixedIterations());
				// This works because the proto configuration enum matches the named enum
				final String[] lineSearchNames = SettingsManager
						.getNames((Object[]) FastMLESteppingFunctionSolver.LineSearchMethod.values());
				gd.addChoice("Line_search_method", lineSearchNames,
						lineSearchNames[fitConfig.getLineSearchMethod().getNumber()]);
			}

			gd.addCheckbox("Use_clamping", fitConfig.isUseClamping());
			gd.addCheckbox("Dynamic_clamping", fitConfig.isUseDynamicClamping());
			final PSF psf = fitConfig.getPSF();
			final boolean isAstigmatism = psf.getPsfType() == PSFType.ASTIGMATIC_GAUSSIAN_2D;
			final int nParams = PSFHelper.getParameterCount(psf);
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

			// Extra parameters are needed for calibrated fit solvers
			if (requireCalibration)
			{
				switch (calibration.getCameraType())
				{
					case CCD:
					case EMCCD:
					case SCMOS:
						break;
					default:
						IJ.error(TITLE, fitSolverName + " requires camera calibration");
						return false;
				}

				gd.addMessage(fitSolverName + " requires calibration for camera: " +
						CalibrationProtosHelper.getName(calibration.getCameraType()));
				if (calibration.isSCMOS())
				{
					final String[] models = CameraModelManager.listCameraModels(true);
					gd.addChoice("Camera_model_name", models, fitConfig.getCameraModelName());
				}
				else
				{
					gd.addNumericField("Camera_bias", calibration.getBias(), 2, 6, "Count");
					gd.addNumericField("Gain", calibration.getCountPerPhoton(), 2, 6, "Count/photon");
				}
			}

			gd.showDialog();
			if (gd.wasCanceled())
				return false;

			fitConfig.setRelativeThreshold(getThresholdNumber(gd));
			fitConfig.setAbsoluteThreshold(getThresholdNumber(gd));
			fitConfig.setParameterRelativeThreshold(getThresholdNumber(gd));
			fitConfig.setParameterAbsoluteThreshold(getThresholdNumber(gd));
			fitConfig.setMaxIterations((int) gd.getNextNumber());
			if (isLVM)
				fitConfig.setLambda(gd.getNextNumber());
			if (isFastMLE)
			{
				fitConfig.setFixedIterations(gd.getNextBoolean());
				fitConfig.setLineSearchMethod(gd.getNextChoiceIndex());
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

			if (requireCalibration)
				if (calibration.isSCMOS())
					fitConfig.setCameraModelName(gd.getNextChoice());
				else
				{
					calibration.setBias(Math.abs(gd.getNextNumber()));
					calibration.setCountPerPhoton(Math.abs(gd.getNextNumber()));
					fitConfig.setCalibration(calibration.getCalibration());
				}

			// Do this even if collection of calibration settings was ignored. This ensures the
			// camera model is set.
			if (calibration.isSCMOS())
			{
				fitConfig.setCameraModel(CameraModelManager.load(fitConfig.getCameraModelName()));
				if (!checkCameraModel(fitConfig, sourceBounds, bounds, true))
					return false;
			}

			if (saveSettings)
				saveFitEngineSettings(config);

			try
			{
				if (isLVM)
					Parameters.isAboveZero("Lambda", fitConfig.getLambda());
				// This call will check if the configuration is OK (including convergence criteria)
				fitConfig.getFunctionSolver();
			}
			catch (final Exception e) // IllegalArgumentException, IllegalStateException
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
			if (!fitConfig.getFunctionSolver().isBounded())
			{
				IJ.error(TITLE, "Including neighbours requires a bounded fit solver");
				return false;
			}

		return true;
	}

	/**
	 * Gets the threshold number from the next string field. If invalid then -1 is returned (which deactivates the
	 * threshold).
	 *
	 * @param gd
	 *            the generic dialog
	 * @return the number
	 */
	private static double getThresholdNumber(ExtendedGenericDialog gd)
	{
		try
		{
			return Double.parseDouble(gd.getNextString());
		}
		catch (final NumberFormatException e)
		{
			return -1;
		}
	}

	/**
	 * Check the camera model covers the region of the source.
	 *
	 * @param fitConfig
	 *            the fit config
	 * @param sourceBounds
	 *            the source bounds of the input image
	 * @param cropBounds
	 *            the crop bounds (relative to the input image)
	 * @param initialise
	 *            the initialise flag
	 * @return true, if successful
	 */
	private static boolean checkCameraModel(FitConfiguration fitConfig, Rectangle sourceBounds, Rectangle cropBounds,
			boolean initialise)
	{
		final CalibrationReader calibration = fitConfig.getCalibrationReader();
		if (calibration.isSCMOS() && sourceBounds != null)
		{
			CameraModel cameraModel = fitConfig.getCameraModel();
			if (cameraModel == null)
				throw new IllegalStateException("No camera model for camera type: " + calibration.getCameraType());

			// The camera model origin must be reset to be relative to the source bounds origin
			cameraModel = cropCameraModel(cameraModel, sourceBounds, cropBounds, true);
			if (cameraModel == null)
				return false;

			if (initialise && cameraModel instanceof PerPixelCameraModel)
				((PerPixelCameraModel) cameraModel).initialise();
			fitConfig.setCameraModel(cameraModel);
		}
		return true;
	}

	/**
	 * Combine the source bounds with an internal crop bounds to create the region required from the camera model.
	 *
	 * @param sourceBounds
	 *            the source bounds
	 * @param cropBounds
	 *            the crop bounds
	 * @return the rectangle
	 */
	public static Rectangle combineBounds(Rectangle sourceBounds, Rectangle cropBounds)
	{
		if (sourceBounds == null)
			throw new NullPointerException("No source bounds");
		if (cropBounds == null)
			return sourceBounds;

		// Check the crop is inside the source
		if (!new Rectangle(sourceBounds.width, sourceBounds.height).contains(cropBounds))
			throw new IllegalArgumentException("Crop bounds does not fit within the source width x height");

		// Make relative
		if (sourceBounds.x != 0 || sourceBounds.y != 0)
		{
			cropBounds = (Rectangle) cropBounds.clone();
			cropBounds.x += sourceBounds.x;
			cropBounds.y += sourceBounds.y;
		}

		return cropBounds;
	}

	/**
	 * Crop a camera model for processing data from a cropped image frame of the given bounds. The target bounds are
	 * created by combining the crop with the source bounds. The camera model bounds will be checked to verify that the
	 * target fits within the model.
	 * <p>
	 * If the model is smaller then an error is thrown.
	 * <p>
	 * If the model is larger then a crop will be made using a dialog to select the crop.
	 * <p>
	 * If the model is the same size then no crop is made, even if the origin is incorrect.
	 * <p>
	 * Optionally the model can be updated so that the origin is relative to the source bounds. If no crop is used the
	 * origin will be 0,0. Otherwise it will be equal to the crop origin.
	 * <p>
	 * This method can be used to prepare a camera model for processing images frames of crop width x height.
	 *
	 * @param cameraModel
	 *            the camera model
	 * @param sourceBounds
	 *            the source bounds
	 * @param cropBounds
	 *            the crop bounds (relative to the input image). If null then the full width x height of the source is
	 *            used.
	 * @param resetOrigin
	 *            the reset origin flag (to set the output camera model origin to match the source bounds)
	 * @return the camera model (or null if the dialog was cancelled)
	 * @throws IllegalArgumentException
	 *             If the model is null or the crop cannot be done
	 */
	public static CameraModel cropCameraModel(CameraModel cameraModel, Rectangle sourceBounds, Rectangle cropBounds,
			boolean resetOrigin) throws IllegalArgumentException
	{
		if (cameraModel == null)
			throw new IllegalArgumentException("No camera model");
		Rectangle modelBounds = cameraModel.getBounds();
		if (modelBounds == null || sourceBounds == null)
			return cameraModel;

		// Combine the source bounds with the crop
		sourceBounds = combineBounds(sourceBounds, cropBounds);

		final int width = sourceBounds.width;
		final int height = sourceBounds.height;
		if (modelBounds.width < width || modelBounds.height < height)
			throw new IllegalArgumentException(String.format(
					"Camera model bounds [x=%d,y=%d,width=%d,height=%d] is smaller than image size [%dx%d]",
					modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height, width, height));
		else if (modelBounds.width > width || modelBounds.height > height)
		{
			final GenericDialog gd2 = new GenericDialog("Crop Camera Model");
			//@formatter:off
			gd2.addMessage(String.format(
					"WARNING:\n \nCamera model bounds\n[x=%d,y=%d,width=%d,height=%d]\nare larger than the image size [%dx%d].\n \nCrop the model?",
					modelBounds.x, modelBounds.y, modelBounds.width, modelBounds.height,
					width, height
					));
			//@formatter:on
			final int upperx = modelBounds.x + modelBounds.width - width;
			final int uppery = modelBounds.y + modelBounds.height - height;
			int ox = sourceBounds.x;
			int oy = sourceBounds.y;
			gd2.addSlider("Origin_x", modelBounds.x, upperx, Maths.clip(modelBounds.x, upperx, ox));
			gd2.addSlider("Origin_y", modelBounds.y, uppery, Maths.clip(modelBounds.y, uppery, oy));
			gd2.showDialog();
			if (gd2.wasCanceled())
				return null;
				//throw new IllegalArgumentException("Unknown camera model crop");
			ox = (int) gd2.getNextNumber();
			oy = (int) gd2.getNextNumber();

			final Rectangle bounds = new Rectangle(ox, oy, width, height);
			cameraModel = cameraModel.crop(bounds, false);
			modelBounds = cameraModel.getBounds();
			if (modelBounds.width != bounds.width || modelBounds.height != bounds.height)
				throw new IllegalArgumentException("Failed to crop camera model using bounds: " + bounds);
		}

		if (resetOrigin)
		{
			// Reset origin to the source origin
			cameraModel = cameraModel.copy();
			if (cropBounds == null)
				cameraModel.setOrigin(0, 0);
			else
				cameraModel.setOrigin(cropBounds.x, cropBounds.y);
		}

		return cameraModel;
	}

	/**
	 * Add a result output.
	 * <p>
	 * This can be called after {@link #initialiseImage(ImageSource, Rectangle, boolean)} and before
	 * {@link #initialiseFitting()} to add to
	 * the configured result outputs.
	 *
	 * @param peakResults
	 *            the peak results
	 */
	public void addPeakResults(PeakResults peakResults)
	{
		if (results != null)
			results.addOutput(peakResults);
	}

	private void addTableResults(PeakResultsList resultsList)
	{
		final IJTablePeakResults r = ResultsManager.addTableResults(resultsList, resultsSettings.getResultsTableSettings(),
				resultsSettings.getShowDeviations(), false, false, false);
		if (r != null)
		{
			r.setShowZ(PSFHelper.is3D(resultsList.getPSF()));
			r.setClearAtStart(simpleFit);
			r.setShowEndFrame(integrateFrames > 1);
		}
	}

	private void addFileResults(PeakResultsList resultsList)
	{
		final ResultsFileSettings resultsSettings = this.resultsSettings.getResultsFileSettings();
		if (resultsSettings.getFileFormat().getNumber() > 0)
		{
			String resultsFilename = null;
			if (resultsSettings.getResultsDirectory() != null &&
					new File(resultsSettings.getResultsDirectory()).exists())
				resultsFilename = resultsSettings.getResultsDirectory() + File.separatorChar + source.getName() +
						".results." + ResultsProtosHelper.getExtension(resultsSettings.getFileFormat());
			else // This is used for running via other code calling PeakFit methods,
			// i.e. not as an ImageJ plugin.
			if (plugin_flags == 0)
				resultsFilename = resultsSettings.getResultsFilename();
			final PeakResults r = ResultsManager.addFileResults(resultsList, resultsSettings, resultsFilename,
					this.resultsSettings.getShowDeviations(), integrateFrames > 1, false);
			if (r instanceof FilePeakResults)
			{
				final FilePeakResults fr = (FilePeakResults) r;
				fr.setSortAfterEnd(Prefs.getThreads() > 1);
			}
		}
	}

	private void addMemoryResults(PeakResultsList resultsList, boolean force)
	{
		if (resultsSettings.getResultsInMemorySettings().getInMemory() || force)
		{
			final MemoryPeakResults results = new MemoryPeakResults();
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
	@Override
	public void run(ImageProcessor ip)
	{
		if (source == null)
		{
			IJ.error(TITLE, "No valid image source configured");
			return;
		}
		if (fitMaxima)
			runMaximaFitting();
		else
			// All setup is already done so simply run
			run();

		addSingleFrameOverlay();
	}

	private void addSingleFrameOverlay()
	{
		// If a single frame was processed add the peaks as an overlay if they are in memory
		ImagePlus imp = this.imp;

		if (fitMaxima && singleFrame > 0)
			if (source instanceof IJImageSource)
			{
				final String title = source.getName();
				imp = WindowManager.getImage(title);
			}

		if (singleFrame > 0 && imp != null)
		{
			MemoryPeakResults results = null;
			for (final PeakResults r : this.results.toArray())
				if (r instanceof MemoryPeakResults)
				{
					results = (MemoryPeakResults) r;
					break;
				}
			if (results == null || results.size() == 0)
				return;

			final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
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
				@Override
				public void executeXYR(float x, float y, PeakResult r)
				{
					final PointRoi roi = new PointRoi(x, y);
					final Color c = LUTHelper.getColour(lut, j.decrementAndGet(), size);
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
	 *            the imp
	 * @param ignoreBoundsForNoise
	 *            If true estimate the noise from the entire frame, otherwise use only the ROI bounds
	 */
	public void run(ImagePlus imp, boolean ignoreBoundsForNoise)
	{
		run(new IJImageSource(imp), getBounds(imp), ignoreBoundsForNoise);
	}

	/**
	 * Process the image.
	 *
	 * @param imageSource
	 *            the image source
	 * @param bounds
	 *            the bounds
	 * @param ignoreBoundsForNoise
	 *            the ignore bounds for noise
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
	@SuppressWarnings("null")
	public void run()
	{
		if (source == null)
			return;

		final int totalFrames = source.getFrames();

		ImageStack stack = null;
		if (showProcessedFrames)
			stack = new ImageStack(bounds.width, bounds.height);

		// Do not crop the region from the source if the bounds match the source dimensions
		final Rectangle cropBounds = (bounds.x == 0 && bounds.y == 0 && bounds.width == source.getWidth() &&
				bounds.height == source.getHeight()) ? null : bounds;

		// Use the FitEngine to allow multi-threading.
		final FitEngine engine = createFitEngine(getNumberOfThreads(totalFrames));
		if (engine == null)
			return;

		final int step = Utils.getProgressInterval(totalFrames);

		// To pre-process data for noise estimation
		boolean isFitCameraCounts = false;
		CameraModel cameraModel = null;
		if (ignoreBoundsForNoise)
		{
			isFitCameraCounts = fitConfig.isFitCameraCounts();
			cameraModel = fitConfig.getCameraModel();
		}

		runTime = System.nanoTime();
		boolean shutdown = false;
		int slice = 0;
		final String format = String.format("Slice: %%d / %d (Results=%%d)", totalFrames);
		while (!shutdown)
		{
			// Noise can optionally be estimated from the entire frame
			float[] data = (ignoreBoundsForNoise) ? source.next() : source.next(cropBounds);
			if (data == null)
				break;

			if (++slice % step == 0)
				if (Utils.showStatus(String.format(format, slice, results.size())))
					IJ.showProgress(slice, totalFrames);

			float noise = Float.NaN;
			if (ignoreBoundsForNoise)
			{
				// We must pre-process the data before noise estimation
				final float[] data2 = data.clone();
				if (isFitCameraCounts)
					cameraModel.removeBias(data2);
				else
					cameraModel.removeBiasAndGain(data2);

				noise = FitWorker.estimateNoise(data2, source.getWidth(), source.getHeight(), config.getNoiseMethod());

				// Crop the data to the region
				data = IJImageConverter.getData(data, source.getWidth(), source.getHeight(), bounds, null);
			}

			if (showProcessedFrames)
				stack.addSlice(String.format("Frame %d - %d", source.getStartFrameNumber(), source.getEndFrameNumber()),
						data);

			// Get the frame number from the source to allow for interlaced and aggregated data
			engine.run(createJob(source.getStartFrameNumber(), source.getEndFrameNumber(), data, bounds, noise));

			shutdown = escapePressed();
		}

		engine.end(shutdown);
		time = engine.getTime();
		runTime = System.nanoTime() - runTime;

		if (showProcessedFrames)
			Utils.display("Processed frames", stack);

		showResults();

		source.close();
	}

	/**
	 * Gets the number of threads.
	 *
	 * @param totalFrames
	 *            the total frames
	 * @return the number of threads
	 */
	private static int getNumberOfThreads(int totalFrames)
	{
		final int t = Math.max(1, (int) (fractionOfThreads * Prefs.getThreads()));
		return FastMath.min(totalFrames, t);
	}

	/**
	 * Check if the frame should be ignored (relevant when using interlaced data).
	 *
	 * @param frame
	 *            the frame
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
			//System.out.printf("Skipping %d\n", frame);
			return true;
		final int frameInBlock = (frame - dataStart) % (dataBlock + dataSkip);
		if (frameInBlock >= dataBlock)
			//System.out.printf("Skipping %d\n", frame);
			return true;
		return false;
	}

	private FitJob createJob(int startFrame, int endFrame, float[] data, Rectangle bounds, float noise)
	{
		FitParameters fitParams = null;
		if (startFrame != endFrame)
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
			return new ParameterisedFitJob(fitParams, startFrame, data, bounds);
		return new FitJob(startFrame, data, bounds);
	}

	private static boolean escapePressed()
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
	 *            the number of threads
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
	 *            the number of threads
	 * @param queue
	 *            the queue
	 * @param queueSize
	 *            the queue size
	 * @return The fiting engine
	 */
	public FitEngine createFitEngine(int numberOfThreads, FitQueue queue, int queueSize)
	{
		// Ensure thread safety
		final PeakResultsList list = (numberOfThreads > 1) ? results.getThreadSafeList() : results;

		// Reduce to single object for speed
		final PeakResults r = (results.numberOfOutputs() == 1) ? list.toArray()[0] : list;

		// Update the configuration
		if (!updateFitConfiguration(config))
			return null;

		final FitEngine engine = new FitEngine(config, r, numberOfThreads, queue, queueSize);

		// Write settings out to the IJ log
		if (resultsSettings.getLogProgress())
		{
			IJ.log("-=-=-=-");
			IJ.log("Peak Fit");
			IJ.log("-=-=-=-");
			Utils.log("Initial Peak SD = %s,%s", Utils.rounded(fitConfig.getInitialXSD()),
					Utils.rounded(fitConfig.getInitialYSD()));
			final SpotFilter spotFilter = engine.getSpotFilter();
			IJ.log("Spot Filter = " + spotFilter.getDescription());
			final int w = 2 * engine.getFitting() + 1;
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
	 * @param config
	 *            the config
	 * @return true, if successful
	 */
	private boolean updateFitConfiguration(FitEngineConfiguration config)
	{
		final FitConfiguration fitConfig = config.getFitConfiguration();

		// Adjust the settings that are relevant within the fitting configuration.
		fitConfig.setComputeResiduals(config.getResidualsThreshold() < 1);
		logger = (resultsSettings.getLogProgress()) ? new IJLogger() : null;
		fitConfig.setLog(logger);

		if (resultsSettings.getShowDeviations())
			// Note: This may already by true if the deviations are needed for the smart filter
			fitConfig.setComputeDeviations(resultsSettings.getShowDeviations());

		config.configureOutputUnits();

		if (!checkCameraModel(fitConfig, source.getBounds(), bounds, true))
			return false;
		return true;
	}

	/**
	 * Load the selected results from memory. All multiple frame results are added directly to the results. All single
	 * frame results are added to a list of candidate maxima per frame and fitted using the configured parameters.
	 */
	private void runMaximaFitting()
	{
		final MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false, DistanceUnit.PIXEL);
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
		if (engine == null)
			return;

		final int step = Utils.getProgressInterval(totalFrames);

		runTime = System.nanoTime();
		final ArrayList<PeakResult> sliceCandidates = new ArrayList<>();
		final FrameCounter counter = new FrameCounter(results.getFirstFrame());
		results.forEach(new PeakResultProcedureX()
		{
			@Override
			public boolean execute(PeakResult r)
			{
				if (counter.advance(r.getFrame()))
				{
					if (escapePressed())
						return true;
					final int slice = counter.previousFrame();
					if (slice % step == 0)
						if (Utils.showStatus("Slice: " + slice + " / " + totalFrames))
							IJ.showProgress(slice, totalFrames);

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
		final boolean shutdown = escapePressed();
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
		final int[] maxIndices = new int[sliceCandidates.size()];
		int count = 0;
		final ArrayList<PeakResult> processedResults = new ArrayList<>(sliceCandidates.size());
		for (final PeakResult result : sliceCandidates)
			// Add ExtendedPeakResults to the results if they span multiple frames (they are the result of previous fitting).
			if (result instanceof ExtendedPeakResult && result.getFrame() != result.getEndFrame())
				processedResults.add(result);
			else
				// Fit single frame results.
				maxIndices[count++] = result.getOrigX() + bounds.width * result.getOrigY();

		if (!processedResults.isEmpty())
			this.results.addAll(processedResults);

		if (count != 0)
		{
			final float[] data = source.get(slice);
			if (data == null)
				return false;

			final FitParameters fitParams = new FitParameters();
			fitParams.maxIndices = Arrays.copyOf(maxIndices, count);
			final FitJob job = new ParameterisedFitJob(fitParams, slice, data, bounds);

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

		// Do not use set() as we support merging a partial calibration
		fitConfig.mergeCalibration(cal);
		final CalibrationReader calibration = fitConfig.getCalibrationReader();

		textCameraType.select(CalibrationProtosHelper.getName(calibration.getCameraType()));
		if (calibration.hasNmPerPixel())
			textNmPerPixel.setText("" + calibration.getNmPerPixel());
		if (calibration.hasExposureTime())
			textExposure.setText("" + calibration.getExposureTime());
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
		// This will clear everything and merge the configuration
		this.config.setFitEngineSettings(fitEngineSettings);
		fitConfig = this.config.getFitConfiguration();

		textDataFilterType.select(SettingsManager.getDataFilterTypeNames()[config.getDataFilterType().ordinal()]);
		textDataFilterMethod
				.select(SettingsManager.getDataFilterMethodNames()[config.getDataFilterMethod(0).ordinal()]);
		textSmooth.setText("" + config.getDataFilterParameterValue(0));
		textSearch.setText("" + config.getSearch());
		textBorder.setText("" + config.getBorder());
		textFitting.setText("" + config.getFitting());
		if (!maximaIdentification)
		{
			textFitSolver.select(SettingsManager.getFitSolverNames()[fitConfig.getFitSolver().ordinal()]);
			if (extraOptions)
				textFitBackground.setState(fitConfig.isBackgroundFitting());
			textFailuresLimit.setText("" + config.getFailuresLimit());
			textPassRate.setText("" + config.getPassRate());
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
