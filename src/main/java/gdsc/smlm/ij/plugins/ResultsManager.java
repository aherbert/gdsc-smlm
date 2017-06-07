package gdsc.smlm.ij.plugins;

import java.awt.Choice;
import java.awt.EventQueue;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumSet;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFileChooser;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.utils.TurboList;

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

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.results.ResultsFileFormat;
import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsMode;
import gdsc.smlm.ij.results.ResultsTable;
import gdsc.smlm.ij.settings.Constants;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.ResultsSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.BinaryFilePeakResults;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.FileFormat;
import gdsc.smlm.results.FilePeakResults;
import gdsc.smlm.results.MALKFilePeakResults;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.PeakResultsReader;
import gdsc.smlm.results.ResultOption;
import gdsc.smlm.results.TSFPeakResultsWriter;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionListener;
import ij.gui.GenericDialog;
import ij.gui.YesNoCancelDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.util.Java2;

/**
 * Opens peaks results and displays/converts them
 */
public class ResultsManager implements PlugIn
{
	public enum InputSource
	{
		//@formatter:off
		FILE{ public String getName() { return "File"; }}, 
		MEMORY{ public String getName() { return "Memory"; }}, 
		MEMORY_MULTI_FRAME{ public String getName() { return "Memory (Multi-Frame)"; }},
		MEMORY_SINGLE_FRAME{ public String getName() { return "Memory (Single-Frame)"; }},
		MEMORY_CLUSTERED{ public String getName() { return "Memory (Clustered)"; }}, 
		NONE{ public String getName() { return "None"; }};
		//@formatter:on

		@Override
		public String toString()
		{
			return getName();
		}

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();
	}

	private static String TITLE = "Peak Results Manager";

	static String INPUT_FILE = "File";
	static String INPUT_MEMORY = "Memory";
	static String INPUT_NONE = "[None]";

	private static String inputOption = "";
	private static String inputFilename = Prefs.get(Constants.inputFilename, "");
	private static boolean chooseRoi = false;
	private static String roiImage = "";
	private Rectangle roiBounds;
	private int roiImageWidth, roiImageHeight;

	private ResultsSettings resultsSettings = new ResultsSettings();
	private boolean extraOptions;

	private boolean fileInput = false;

	private static double input_nmPerPixel = Prefs.get(Constants.inputNmPerPixel, 0);
	private static double input_gain = Prefs.get(Constants.inputGain, 1);
	private static double input_exposureTime = Prefs.get(Constants.inputExposureTime, 0);
	private static float input_noise = (float) Prefs.get(Constants.inputNoise, 0);
	private static ArrayList<String> selected;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		extraOptions = Utils.isExtraOptions();
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if ("load".equals(arg))
		{
			batchLoad();
			return;
		}

		if (arg != null && arg.startsWith("clear"))
		{
			Collection<MemoryPeakResults> allResults;
			boolean removeAll = false;
			if (arg.contains("multi"))
			{
				MultiDialog md = new MultiDialog(TITLE, new MultiDialog.MemoryResultsItems());
				md.addSelected(selected);
				md.showDialog();
				if (md.wasCanceled())
					return;
				selected = md.getSelectedResults();
				if (selected.isEmpty())
					return;
				allResults = new ArrayList<MemoryPeakResults>(selected.size());
				for (String name : selected)
				{
					MemoryPeakResults r = MemoryPeakResults.getResults(name);
					if (r != null)
						allResults.add(r);
				}
			}
			else
			{
				removeAll = true;
				allResults = MemoryPeakResults.getAllResults();
			}
			if (allResults.isEmpty())
				return;

			long memorySize = 0;
			int size = 0;
			for (MemoryPeakResults results : allResults)
			{
				memorySize += MemoryPeakResults.estimateMemorySize(results.getResults());
				size += results.size();
			}
			String memory = MemoryPeakResults.memorySizeString(memorySize);
			String count = Utils.pleural(size, "result");
			String sets = Utils.pleural(allResults.size(), "set");

			GenericDialog gd = new GenericDialog(TITLE);

			gd.addMessage(String.format("Do you want to remove %s from memory (%s, %s)?", count, sets, memory));
			gd.enableYesNoCancel();
			gd.showDialog();
			if (!gd.wasOKed())
				return;

			if (removeAll)
				MemoryPeakResults.clearMemory();
			else
			{
				for (MemoryPeakResults results : allResults)
					MemoryPeakResults.removeResults(results.getName());
			}

			SummariseResults.clearSummaryTable();
			IJ.log(String.format("Cleared %s (%s, %s)", count, sets, memory));
			return;
		}

		if (!showDialog())
			return;

		MemoryPeakResults results = loadResults(inputOption);

		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		results = cropToRoi(results);
		if (results.size() == 0)
		{
			IJ.error(TITLE, "No results within the crop region");
			return;
		}

		if (resultsSettings.resultsInMemory && fileInput)
			MemoryPeakResults.addResults(results);

		Rectangle bounds = results.getBounds(true);
		boolean showDeviations = resultsSettings.showDeviations && canShowDeviations(results);
		boolean showEndFrame = canShowEndFrame(results);
		boolean showId = canShowId(results);

		// Display the configured output
		PeakResultsList outputList = new PeakResultsList();

		outputList.copySettings(results);
		//String title = results.getSource();
		//if (title == null || title.length() == 0)
		//	output.setSource(TITLE);

		addTableResults(results, outputList, showDeviations, showEndFrame);
		addImageResults(outputList, results.getName(), bounds, results.getNmPerPixel(), results.getGain());
		addFileResults(outputList, showDeviations, showEndFrame, showId);

		if (outputList.numberOfOutputs() == 0)
		{
			// Error if results were loaded since we have no outputs
			if (!(resultsSettings.resultsInMemory && fileInput))
			{
				IJ.error(TITLE, "No output selected");
				return;
			}
		}
		else
		{
			IJ.showStatus("Processing outputs ...");

			// Reduce to single object for speed
			PeakResults output = (outputList.numberOfOutputs() == 1) ? outputList.toArray()[0] : outputList;

			output.begin();

			// Process in batches to provide progress
			List<PeakResult> list = results.getResults();
			int progress = 0;
			int totalProgress = list.size();
			int stepProgress = Utils.getProgressInterval(totalProgress);
			TurboList<PeakResult> batch = new TurboList<PeakResult>(stepProgress);
			for (PeakResult result : list)
			{
				if (progress % stepProgress == 0)
				{
					IJ.showProgress(progress, totalProgress);
				}
				progress++;
				batch.addf(result);

				if (batch.size() == stepProgress)
				{
					output.addAll(batch);
					batch.clearf();
					if (isInterrupted())
						break;
				}
			}
			IJ.showProgress(1);
			output.end();
		}

		IJ.showStatus(String.format("Processed %d result%s", results.size(), (results.size() > 1) ? "s" : ""));
	}

	/**
	 * Check if the escape key has been pressed. Show a status aborted message if true.
	 * 
	 * @return True if aborted
	 */
	public static boolean isInterrupted()
	{
		if (IJ.escapePressed())
		{
			IJ.beep();
			IJ.showStatus("Aborted");
			return true;
		}
		return false;
	}

	private boolean canShowDeviations(MemoryPeakResults results)
	{
		for (PeakResult r : results)
			if (r.paramsStdDev != null)
				return true;
		return false;
	}

	private boolean canShowEndFrame(MemoryPeakResults results)
	{
		for (PeakResult r : results.getResults())
			if (r.getFrame() != r.getEndFrame())
				return true;
		return false;
	}

	private boolean canShowId(MemoryPeakResults results)
	{
		final int id = results.getHead().getId();
		for (PeakResult r : results.getResults())
			if (id != r.getId())
				return true;
		return false;
	}

	private void addTableResults(MemoryPeakResults results, PeakResultsList resultsList, boolean showDeviations,
			boolean showEndFrame)
	{
		if (resultsSettings.getResultsTable() != ResultsTable.NONE)
		{
			IJTablePeakResults r = new IJTablePeakResults(showDeviations);
			r.setPeakIdColumnName("Frame");
			r.setShowCalibratedValues(resultsSettings.getResultsTable() == ResultsTable.CALIBRATED);
			// Get a bias if required
			Calibration calibration = results.getCalibration();
			if (r.isShowCalibratedValues() && calibration.getBias() == 0)
			{
				GenericDialog gd = new GenericDialog(TITLE);
				gd.addMessage("Calibrated results requires a camera bias");
				gd.addNumericField("Camera_bias (ADUs)", calibration.getBias(), 2);
				gd.showDialog();
				if (!gd.wasCanceled())
				{
					calibration.setBias(Math.abs(gd.getNextNumber()));
				}
			}
			r.setShowEndFrame(showEndFrame);
			resultsList.addOutput(r);
		}
	}

	private void addImageResults(PeakResultsList resultsList, String title, Rectangle bounds, double nmPerPixel,
			double gain)
	{
		if (resultsSettings.getResultsImage() != ResultsImage.NONE)
		{
			IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(resultsSettings.getResultsImage(),
					resultsSettings.weightedImage, resultsSettings.equalisedImage, title, bounds, nmPerPixel, gain,
					resultsSettings.imageScale, resultsSettings.precision, ResultsMode.ADD);
			if (extraOptions)
				image.setRollingWindowSize(resultsSettings.imageRollingWindow);
			image.setRepaintDelay(2000);
			resultsList.addOutput(image);
		}
	}

	private void addFileResults(PeakResultsList resultsList, boolean showDeviations, boolean showEndFrame,
			boolean showId)
	{
		if (resultsSettings.resultsFilename != null && resultsSettings.resultsFilename.length() > 0)
		{
			// Remove extension
			resultsSettings.resultsFilename = Utils.replaceExtension(resultsSettings.resultsFilename,
					resultsSettings.getResultsFileFormat().getExtension());

			if (fileInput && inputFilename.equals(resultsSettings.resultsFilename))
			{
				IJ.log(TITLE + ": Input and output files are the same, skipping output ...");
				return;
			}

			// Check if file exists
			File file = new File(resultsSettings.resultsFilename);
			if (file.exists())
			{
				YesNoCancelDialog d = new YesNoCancelDialog(IJ.getInstance(), TITLE,
						"Overwrite existing file?\n" + resultsSettings.resultsFilename);
				if (!d.yesPressed())
					return;
			}

			File parent = file.getParentFile();
			if (parent != null && parent.exists())
			{
				PeakResults r;
				switch (resultsSettings.getResultsFileFormat())
				{
					case GDSC_BINARY:
						r = new BinaryFilePeakResults(resultsSettings.resultsFilename, showDeviations, showEndFrame,
								showId);
						break;
					case GDSC_TEXT:
						r = new FilePeakResults(resultsSettings.resultsFilename, showDeviations, showEndFrame, showId);
						break;
					case MALK:
						r = new MALKFilePeakResults(resultsSettings.resultsFilename);
						break;
					case TSF:
						r = new TSFPeakResultsWriter(resultsSettings.resultsFilename);
						break;
					default:
						throw new RuntimeException(
								"Unsupported file format: " + resultsSettings.getResultsFileFormat());
				}
				if (r instanceof FilePeakResults)
				{
					FilePeakResults fr = (FilePeakResults) r;
					fr.setPeakIdColumnName("Frame");
				}
				resultsList.addOutput(r);
			}
		}
	}

	private boolean showDialog()
	{
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		// Build a list of all images with a region ROI
		List<String> titles = new LinkedList<String>();
		if (WindowManager.getWindowCount() > 0)
		{
			for (int imageID : WindowManager.getIDList())
			{
				ImagePlus imp = WindowManager.getImage(imageID);
				if (imp != null && imp.getRoi() != null && imp.getRoi().isArea())
					titles.add(imp.getTitle());
			}
		}

		GlobalSettings settings = SettingsManager.loadSettings();
		resultsSettings = settings.getResultsSettings();

		gd.addMessage("Read the Peak Results and output to a new format");

		gd.addMessage("Select the Peak Results");
		addInput(gd, inputOption, InputSource.MEMORY, InputSource.FILE);

		if (!titles.isEmpty())
			gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", chooseRoi);

		gd.addMessage("--- Table output ---");
		String[] tableNames = SettingsManager.getNames((Object[]) ResultsTable.values());
		gd.addChoice("Results_table", tableNames, tableNames[resultsSettings.getResultsTable().ordinal()],
				new OptionListener<Choice>()
				{
					public void collectOptions(Choice field)
					{
						resultsSettings.setResultsTable(field.getSelectedIndex());
						collectOptions();
					}

					public void collectOptions()
					{
						ResultsTable resultsTable = resultsSettings.getResultsTable();
						if (resultsTable == ResultsTable.NONE)
						{
							return;
						}
						ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
						egd.addCheckbox("Show_deviations", resultsSettings.showDeviations);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return;
						resultsSettings.showDeviations = egd.getNextBoolean();
					}
				});

		gd.addMessage("--- Image output ---");
		String[] imageNames = SettingsManager.getNames((Object[]) ResultsImage.values());
		final EnumSet<ResultsImage> requirePrecision = EnumSet.of(ResultsImage.LOCALISATIONS_AV_PRECISION,
				ResultsImage.SIGNAL_AV_PRECISION);
		final EnumSet<ResultsImage> requireWeighted = EnumSet.of(ResultsImage.LOCALISATIONS,
				ResultsImage.SIGNAL_INTENSITY, ResultsImage.FRAME_NUMBER, ResultsImage.ERROR);
		gd.addChoice("Image", imageNames, imageNames[resultsSettings.getResultsImage().ordinal()],
				new OptionListener<Choice>()
				{
					public void collectOptions(Choice field)
					{
						resultsSettings.setResultsImage(field.getSelectedIndex());
						collectOptions();
					}

					public void collectOptions()
					{
						ResultsImage resultsImage = resultsSettings.getResultsImage();
						if (resultsImage == ResultsImage.NONE)
						{
							return;
						}
						ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
						if (requireWeighted.contains(resultsImage))
							egd.addCheckbox("Weighted", resultsSettings.weightedImage);
						egd.addCheckbox("Equalised", resultsSettings.equalisedImage);
						if (requirePrecision.contains(resultsImage))
							egd.addSlider("Image_Precision (nm)", 5, 30, resultsSettings.precision);
						egd.addSlider("Image_Scale", 1, 15, resultsSettings.imageScale);
						if (extraOptions)
							egd.addNumericField("Image_Window", resultsSettings.imageRollingWindow, 0);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return;
						if (requireWeighted.contains(resultsImage))
							resultsSettings.weightedImage = egd.getNextBoolean();
						resultsSettings.equalisedImage = egd.getNextBoolean();
						if (requirePrecision.contains(resultsImage))
							resultsSettings.precision = egd.getNextNumber();
						resultsSettings.imageScale = egd.getNextNumber();
						if (extraOptions)
							resultsSettings.imageRollingWindow = (int) egd.getNextNumber();
					}
				});

		gd.addMessage("--- File output ---");
		// Do not add a results file to prevent constant overwrite messages
		gd.addFilenameField("Results_file", "");
		// TODO - move these to an options pop-uo
		String[] formatNames = SettingsManager.getNames((Object[]) ResultsFileFormat.values());
		gd.addChoice("Results_format", formatNames, formatNames[resultsSettings.getResultsFileFormat().ordinal()]);
		gd.addMessage(" ");
		gd.addCheckbox("Save_to_memory (file input only)", resultsSettings.resultsInMemory);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		inputFilename = gd.getNextString();
		if (!titles.isEmpty())
			chooseRoi = gd.getNextBoolean();
		resultsSettings.setResultsTable(gd.getNextChoiceIndex());
		resultsSettings.setResultsImage(gd.getNextChoiceIndex());
		resultsSettings.resultsFilename = gd.getNextString();
		resultsSettings.setResultsFileFormat(gd.getNextChoiceIndex());
		resultsSettings.resultsInMemory = gd.getNextBoolean();

		// Check arguments
		try
		{
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

		Prefs.set(Constants.inputFilename, inputFilename);

		if (!titles.isEmpty() && chooseRoi)
		{
			if (titles.size() == 1)
			{
				roiImage = titles.get(0);
				Recorder.recordOption("Image", roiImage);
			}
			else
			{
				String[] items = titles.toArray(new String[titles.size()]);
				ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
				egd.addMessage("Select the source image for the ROI");
				egd.addChoice("Image", items, roiImage);
				egd.showDialog();
				if (egd.wasCanceled())
					return false;
				roiImage = egd.getNextChoice();
			}
			ImagePlus imp = WindowManager.getImage(roiImage);

			roiBounds = imp.getRoi().getBounds();
			roiImageWidth = imp.getWidth();
			roiImageHeight = imp.getHeight();
		}
		else
		{
			roiBounds = null;
		}

		SettingsManager.saveSettings(settings);

		return true;
	}

	/**
	 * Add a list of input sources to the generic dialog. The choice field will be named 'input'. If a file input option
	 * is selected then a field will be added name 'Input_file'.
	 * <p>
	 * If the source is a memory source then it will not be added if it is empty. If not empty then a summary of the
	 * number of localisation is added as a message to the dialog.
	 * 
	 * @param gd
	 * @param inputOption
	 * @param inputs
	 */
	public static void addInput(ExtendedGenericDialog gd, String inputOption, InputSource... inputs)
	{
		addInput(gd, "Input", inputOption, inputs);
	}

	/**
	 * Add a list of input sources to the generic dialog. The choice field will be named inputName. If a file input
	 * option
	 * is selected then a field will be added name 'Input_file'.
	 * <p>
	 * If the source is a memory source then it will not be added if it is empty. If not empty then a summary of the
	 * number of localisation is added as a message to the dialog.
	 * 
	 * @param gd
	 * @param inputName
	 * @param inputOption
	 * @param inputs
	 */
	public static void addInput(ExtendedGenericDialog gd, String inputName, String inputOption, InputSource... inputs)
	{
		ArrayList<String> source = new ArrayList<String>(3);
		boolean fileInput = false;
		for (InputSource input : inputs)
		{
			ResultsManager.addInputSource(source, input);
			if (input == InputSource.FILE)
				fileInput = true;
		}
		if (source.isEmpty())
			addInputSource(source, InputSource.NONE);

		addInputSourceToDialog(gd, inputName, inputOption, source, fileInput);
	}

	/**
	 * Add a list of input sources to the generic dialog. The choice field will be named inputName. If the file input
	 * option
	 * is true then a field will be added name 'Input_file'.
	 * 
	 * @param gd
	 * @param inputName
	 * @param inputOption
	 *            The option to select by default
	 * @param source
	 * @param fileInput
	 */
	@SuppressWarnings("unused")
	public static void addInputSourceToDialog(final ExtendedGenericDialog gd, String inputName, String inputOption,
			ArrayList<String> source, boolean fileInput)
	{
		String[] options = source.toArray(new String[source.size()]);
		// Find the option
		inputOption = removeFormatting(inputOption);

		int optionIndex = 0;
		for (int i = 0; i < options.length; i++)
		{
			String name = removeFormatting(options[i]);
			if (name.equals(inputOption))
			{
				optionIndex = i;
				break;
			}
		}
		final Choice c = gd.addAndGetChoice(inputName, options, options[optionIndex]);
		if (fileInput)
		{
			final TextField tf = gd.addFilenameField("Input_file", inputFilename);
			final JButton b = gd.getLastOptionButton();

			// Add a listener to the choice to enable the file input field.
			// Currently we hide the filename field and pack the dialog. 
			// We may wish to just disable the fields and leave them there.
			// This could be a user configured option in a global GDSC settings class.
			if (Utils.isShowGenericDialog())
			{
				final Label l = gd.getLastLabel();
				final Panel p = gd.getLastPanel();
				ItemListener listener = new ItemListener()
				{
					public void itemStateChanged(ItemEvent e)
					{
						boolean enable = INPUT_FILE.equals(c.getSelectedItem());
						if (enable != l.isVisible())
						{
							l.setVisible(enable);
							p.setVisible(enable);
							//tf.setVisible(enable);
							//b.setVisible(enable);
							gd.pack();
						}
					}
				};

				// Run once to set up
				listener.itemStateChanged(null);

				c.addItemListener(listener);
			}
		}
	}

	/**
	 * Remove the extra information added to a name for use in dialogs
	 * 
	 * @param name
	 *            the formatted name
	 * @return The name
	 */
	public static String removeFormatting(String name)
	{
		int index = name.lastIndexOf('[');
		if (index > 0)
			name = name.substring(0, index - 1);
		return name;
	}

	/**
	 * Add an input source the list. If the source is a memory source then it will not be added if it is
	 * empty. If not empty then a summary of the number of localisation is added as a message to the dialog.
	 * 
	 * @param source
	 * @param input
	 */
	public static void addInputSource(ArrayList<String> source, InputSource input)
	{
		switch (input)
		{
			case NONE:
				source.add(INPUT_NONE);
				break;

			case FILE:
				source.add(INPUT_FILE);
				break;

			case MEMORY:
			case MEMORY_MULTI_FRAME:
			case MEMORY_SINGLE_FRAME:
			case MEMORY_CLUSTERED:
				for (String name : MemoryPeakResults.getResultNames())
				{
					addInputSource(source, MemoryPeakResults.getResults(name), input);
				}
				break;
		}
	}

	/**
	 * Add a memory input source to the list
	 * 
	 * @param source
	 * @param memoryResults
	 * @param input
	 *            MEMORY_MULTI_FRAME : Select only those results with at least one result spanning frames,
	 *            MEMORY_CLUSTERED : Select only those results which have at least some IDs above zero (allowing zero to
	 *            be a valid cluster Id for no cluster)
	 */
	public static void addInputSource(ArrayList<String> source, MemoryPeakResults memoryResults, InputSource input)
	{
		if (memoryResults.size() > 0)
		{
			switch (input)
			{
				case MEMORY_MULTI_FRAME:
					if (!isMultiFrame(memoryResults))
						return;
					break;
				case MEMORY_SINGLE_FRAME:
					if (isMultiFrame(memoryResults))
						return;
					break;
				case MEMORY_CLUSTERED:
					if (!hasID(memoryResults))
						return;
					break;
				default:
			}
			source.add(getName(memoryResults));
		}
	}

	/**
	 * Get the name of the results for use in dialogs
	 * 
	 * @param memoryResults
	 * @return The name
	 */
	public static String getName(MemoryPeakResults memoryResults)
	{
		return memoryResults.getName() + " [" + memoryResults.size() + "]";
	}

	/**
	 * Load results from memory using a name from a dialog
	 * 
	 * @param name
	 * @return The results
	 */
	public static MemoryPeakResults loadMemoryResults(String name)
	{
		return MemoryPeakResults.getResults(removeFormatting(name));
	}

	/**
	 * Check for multi-frame results.
	 * 
	 * @param memoryResults
	 * @return True if at least one result spanning frames
	 */
	public static boolean isMultiFrame(MemoryPeakResults memoryResults)
	{
		for (PeakResult r : memoryResults.getResults())
			if (r.getFrame() < r.getEndFrame())
				return true;
		return false;
	}

	/**
	 * Check for any IDs above zero.
	 * 
	 * @param memoryResults
	 * @return True if any results have IDs above zero
	 */
	public static boolean hasID(MemoryPeakResults memoryResults)
	{
		for (PeakResult r : memoryResults.getResults())
			if (r.getId() > 0)
				return true;
		return false;
	}

	/**
	 * Check for all IDs above zero.
	 * 
	 * @param memoryResults
	 * @return True if all results have IDs above zero
	 */
	public static boolean isID(MemoryPeakResults memoryResults)
	{
		for (PeakResult r : memoryResults.getResults())
			if (r.getId() <= 0)
				return false;
		return true;
	}

	/**
	 * All results must be an ExtendedPeakResult.
	 * 
	 * @param memoryResults
	 * @return True if all ExtendedPeakResult
	 */
	public static boolean isExtended(MemoryPeakResults memoryResults)
	{
		for (PeakResult r : memoryResults.getResults())
			if (!(r instanceof ExtendedPeakResult))
				return false;
		return true;
	}

	/**
	 * Gets the name of the next input source from the dialog
	 * 
	 * @param gd
	 */
	public static String getInputSource(GenericDialog gd)
	{
		String source = gd.getNextChoice();
		return removeFormatting(source);
	}

	/**
	 * Load the results from the named input option
	 * 
	 * @param inputOption
	 * @param checkCalibration
	 *            Set to true to ensure the results have a valid calibration
	 * @return
	 */
	public static MemoryPeakResults loadInputResults(String inputOption, boolean checkCalibration)
	{
		MemoryPeakResults results = null;
		PeakResultsReader reader = null;
		if (inputOption.equals(INPUT_NONE))
		{

		}
		else if (inputOption.equals(INPUT_FILE))
		{
			IJ.showStatus("Reading results file ...");
			reader = new PeakResultsReader(inputFilename);
			IJ.showStatus("Reading " + reader.getFormat() + " results file ...");
			ResultOption[] options = reader.getOptions();
			if (options != null)
				collectOptions(reader, options);
			reader.setTracker(new IJTrackProgress());
			results = reader.getResults();
			reader.getTracker().progress(1.0);

			if (results != null && results.size() > 0)
			{
				// If the name contains a .tif suffix then create an image source
				if (results.getName() != null && results.getName().contains(".tif") && results.getSource() == null)
				{
					int index = results.getName().indexOf(".tif");
					results.setSource(new IJImageSource(results.getName().substring(0, index)));
				}
			}
		}
		else
		{
			results = loadMemoryResults(inputOption);
		}

		if (results != null && results.size() > 0 && checkCalibration)
		{
			if (!checkCalibration(results, reader))
				results = null;
		}
		IJ.showStatus("");
		return results;
	}

	private static void collectOptions(PeakResultsReader reader, ResultOption[] options)
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Options required for file format: " + reader.getFormat().getName());
		for (ResultOption option : options)
		{
			if (option.hasValues())
			{
				String[] items = new String[option.values.length];
				for (int i = 0; i < items.length; i++)
					items[i] = option.values[i].toString();
				gd.addChoice(getName(option), items, option.getValue().toString());
			}
			else if (option.getValue() instanceof Number)
			{
				Number n = (Number) option.getValue();
				if (n.doubleValue() == n.intValue())
				{
					gd.addNumericField(getName(option), n.intValue(), 0);
				}
				else
				{
					String value = n.toString();
					int sig = 0;
					int index = value.indexOf('.');
					if (index != -1)
					{
						// There is a decimal point. Count the digits after it
						while (++index < value.length())
						{
							if (!Character.isDigit(value.charAt(index)))
							{
								// A non-digit character after the decimal point is for scientific notation
								sig = -sig;
								break;
							}
							sig++;
						}
					}
					gd.addNumericField(getName(option), n.doubleValue(), sig);
				}
			}
			else if (option.getValue() instanceof String)
			{
				gd.addStringField(getName(option), (String) option.getValue());
			}
			else if (option.getValue() instanceof Boolean)
			{
				gd.addCheckbox(getName(option), (Boolean) option.getValue());
			}
			else
			{
				IJ.log(TITLE + ": Unsupported reader option: " + option.name + "=" + option.getValue().toString());
			}
		}
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		try
		{
			for (ResultOption option : options)
			{
				if (option.hasValues())
				{
					option.setValue(option.values[gd.getNextChoiceIndex()]);
				}
				else if (option.getValue() instanceof Number)
				{
					double d = gd.getNextNumber();
					// Convert to the correct type using the String value constructor for the number
					option.setValue(
							option.getValue().getClass().getConstructor(String.class).newInstance(Double.toString(d)));
				}
				else if (option.getValue() instanceof String)
				{
					option.setValue(gd.getNextString());
				}
				else if (option.getValue() instanceof Boolean)
				{
					option.setValue(gd.getNextBoolean());
				}
			}
			reader.setOptions(options);
		}
		catch (Exception e)
		{
			// This can occur if the options are not valid
			IJ.log(TITLE + ": Failed to configure reader options: " + e.getMessage());
		}
	}

	private static String getName(ResultOption option)
	{
		return option.name.replace(' ', '_');
	}

	/**
	 * Check the calibration of the results exists, if not then prompt for it with a dialog
	 * 
	 * @param results
	 *            The results
	 * @return True if OK; false if calibration dialog cancelled
	 */
	public static boolean checkCalibration(MemoryPeakResults results)
	{
		return checkCalibration(results, null);
	}

	/**
	 * Check the calibration of the results exists, if not then prompt for it with a dialog
	 * 
	 * @param results
	 *            The results
	 * @param reader
	 *            Used to determine the file type
	 * @return True if OK; false if calibration dialog cancelled
	 */
	private static boolean checkCalibration(MemoryPeakResults results, PeakResultsReader reader)
	{
		// Check for Calibration
		Calibration calibration = results.getCalibration();
		String msg = "partially calibrated";
		if (calibration == null)
		{
			// Make sure the user knows all the values have not been set
			calibration = new Calibration();
			msg = "uncalibrated";
		}
		else
		{
			// Validate to set the valid flags
			calibration.validate();
		}

		final float noise = getNoise(results);
		// Only check for essential calibration settings (i.e. not readNoise, bias, emCCD, amplification)
		if (!calibration.hasNmPerPixel() || !calibration.hasGain() || !calibration.hasExposureTime() || noise <= 0)
		{
			boolean convert = false;

			// We may have results that are within configured bounds. If so we do not need the conversion
			boolean showConvert = true;
			if (results.getBounds(false) != null)
			{
				final Rectangle bounds = results.getBounds(false);
				showConvert = false;
				for (PeakResult r : results.getResults())
				{
					if (!bounds.contains(r.getXPosition(), r.getYPosition()))
					{
						showConvert = true;
						break;
					}
				}
			}

			if (!calibration.hasNmPerPixel())
				calibration.setNmPerPixel(input_nmPerPixel);
			if (!calibration.hasGain())
				calibration.setGain(input_gain);
			if (!calibration.hasExposureTime())
				calibration.setExposureTime(input_exposureTime);

			Rectangle2D.Float dataBounds = results.getDataBounds();

			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage(
					String.format("Results are %s.\nData bounds = (%s,%s) to (%s,%s)", msg, Utils.rounded(dataBounds.x),
							Utils.rounded(dataBounds.y), Utils.rounded(dataBounds.y + dataBounds.getWidth()),
							Utils.rounded(dataBounds.x + dataBounds.getHeight())));
			gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
			gd.addNumericField("Gain (ADU/photon)", calibration.getGain(), 2);
			gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);
			if (noise <= 0)
				gd.addNumericField("Noise (ADU)", input_noise, 2);
			if (showConvert)
				gd.addCheckbox("Convert_nm_to_pixels", convert);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			input_nmPerPixel = Math.abs(gd.getNextNumber());
			input_gain = Math.abs(gd.getNextNumber());
			input_exposureTime = Math.abs(gd.getNextNumber());
			if (noise == 0)
				input_noise = Math.abs((float) gd.getNextNumber());
			if (showConvert)
				convert = gd.getNextBoolean();

			Prefs.set(Constants.inputNmPerPixel, input_nmPerPixel);
			Prefs.set(Constants.inputGain, input_gain);
			Prefs.set(Constants.inputExposureTime, input_exposureTime);
			Prefs.set(Constants.inputNoise, input_noise);

			results.setCalibration(new Calibration(input_nmPerPixel, input_gain, input_exposureTime));

			if (convert && input_nmPerPixel > 0)
			{
				// Note: NSTORM stores 2xSD
				final double widthConversion = (reader != null && reader.getFormat() == FileFormat.NSTORM)
						? 1.0 / (2 * input_nmPerPixel) : 1.0 / input_nmPerPixel;
				for (PeakResult p : results.getResults())
				{
					p.params[Gaussian2DFunction.X_POSITION] /= input_nmPerPixel;
					p.params[Gaussian2DFunction.Y_POSITION] /= input_nmPerPixel;
					p.params[Gaussian2DFunction.X_SD] *= widthConversion;
					p.params[Gaussian2DFunction.Y_SD] *= widthConversion;
				}
			}
			if (noise == 0)
			{
				for (PeakResult p : results.getResults())
				{
					p.noise = input_noise;
				}
			}
		}
		return true;
	}

	/**
	 * get the first non-zero noise value
	 * 
	 * @param results
	 * @return The noise (zero if no results have a noise value)
	 */
	private static float getNoise(MemoryPeakResults results)
	{
		for (PeakResult r : results.getResults())
		{
			if (r.noise != 0)
				return r.noise;
		}
		return 0;
	}

	/**
	 * Load the results from the named input option
	 * 
	 * @param inputOption
	 * @return
	 */
	private MemoryPeakResults loadResults(String inputOption)
	{
		if (inputOption.equals(INPUT_FILE))
		{
			fileInput = true;
		}
		return loadInputResults(inputOption, true);
	}

	private MemoryPeakResults cropToRoi(MemoryPeakResults results)
	{
		if (roiBounds == null)
			return results;

		// Adjust bounds relative to input results image
		double xscale = roiImageWidth / results.getBounds().width;
		double yscale = roiImageHeight / results.getBounds().height;
		roiBounds.x /= xscale;
		roiBounds.width /= xscale;
		roiBounds.y /= yscale;
		roiBounds.height /= yscale;

		float minX = (int) (roiBounds.x);
		float maxX = (int) Math.ceil(roiBounds.x + roiBounds.width);
		float minY = (int) (roiBounds.y);
		float maxY = (int) Math.ceil(roiBounds.y + roiBounds.height);

		// Create a new set of results within the bounds
		MemoryPeakResults newResults = new MemoryPeakResults();
		newResults.begin();
		for (PeakResult peakResult : results.getResults())
		{
			float x = peakResult.params[Gaussian2DFunction.X_POSITION];
			float y = peakResult.params[Gaussian2DFunction.Y_POSITION];
			if (x < minX || x > maxX || y < minY || y > maxY)
				continue;
			newResults.add(peakResult);
		}
		newResults.end();
		newResults.copySettings(results);
		newResults.setBounds(new Rectangle((int) minX, (int) minY, (int) (maxX - minX), (int) (maxY - minY)));
		return newResults;
	}

	/**
	 * Gets the input filename that will be used in {@link ResultsManager#loadInputResults(String, boolean)}.
	 *
	 * @return the input filename
	 */
	static String getInputFilename()
	{
		return inputFilename;
	}

	/**
	 * Sets the input filename that will be used in {@link ResultsManager#loadInputResults(String, boolean)}.
	 *
	 * @param inputFilename
	 *            the new input filename
	 */
	static void setInputFilename(String inputFilename)
	{
		ResultsManager.inputFilename = inputFilename;
	}

	private String omDirectory;
	private File[] omFiles;

	/**
	 * Batch load a set of results files.
	 */
	private void batchLoad()
	{
		// Adapted from ij.io.Opener.openMultiple

		String resetInputFilename = inputFilename;

		Java2.setSystemLookAndFeel();
		// run JFileChooser in a separate thread to avoid possible thread deadlocks
		try
		{
			EventQueue.invokeAndWait(new Runnable()
			{
				public void run()
				{
					JFileChooser fc = new JFileChooser();
					fc.setMultiSelectionEnabled(true);
					File dir = null;
					String sdir = OpenDialog.getDefaultDirectory();
					if (sdir != null)
						dir = new File(sdir);
					if (dir != null)
						fc.setCurrentDirectory(dir);
					int returnVal = fc.showOpenDialog(IJ.getInstance());
					if (returnVal != JFileChooser.APPROVE_OPTION)
						return;
					omFiles = fc.getSelectedFiles();
					if (omFiles.length == 0)
					{ // getSelectedFiles does not work on some JVMs
						omFiles = new File[1];
						omFiles[0] = fc.getSelectedFile();
					}
					omDirectory = fc.getCurrentDirectory().getPath() + File.separator;
				}
			});
		}
		catch (Exception e)
		{
		}
		if (omDirectory == null)
			return;
		OpenDialog.setDefaultDirectory(omDirectory);
		for (int i = 0; i < omFiles.length; i++)
		{
			String path = omDirectory + omFiles[i].getName();
			load(path);
		}

		inputFilename = resetInputFilename;
	}

	private void load(String path)
	{
		inputFilename = path;
		// Record this as a single load of the results manager.
		// This should support any dialogs that are presented in loadInputResults(...)
		// to get the calibration.
		if (Recorder.record)
		{
			Recorder.setCommand("Results Manager");
			Recorder.recordOption("input", INPUT_FILE);
			Recorder.recordOption("input_file", path);
			Recorder.recordOption("results_table", ResultsTable.NONE.getName());
			Recorder.recordOption("image", ResultsImage.NONE.getName());
			Recorder.recordOption("results_file", "[]");
			Recorder.recordOption("save_to_memory");
		}
		MemoryPeakResults results = loadInputResults(INPUT_FILE, true);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded from " + path);
		}
		else
		{
			if (Recorder.record)
				Recorder.saveCommand();
			MemoryPeakResults.addResults(results);
		}
	}
}
