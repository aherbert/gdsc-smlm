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

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.IJTrackProgress;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.results.ResultsImage;
import gdsc.smlm.ij.results.ResultsMode;
import gdsc.smlm.ij.results.ResultsTable;
import gdsc.smlm.ij.settings.Constants;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.ResultsSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.BinaryFilePeakResults;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.FileFormat;
import gdsc.smlm.results.FilePeakResults;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.PeakResultsReader;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.YesNoCancelDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;

import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

/**
 * Opens peaks results and displays/converts them
 */
public class ResultsManager implements PlugIn, MouseListener
{
	public enum InputSource
	{
		FILE("File"), MEMORY("Memory"), MEMORY_TRACED("Memory (Traced)"), NONE("None");

		private String name;

		private InputSource(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
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

	private boolean fileInput = false;

	private GenericDialog gd;
	private TextField text1;
	private TextField text2;

	private static double input_nmPerPixel = Prefs.get(Constants.inputNmPerPixel, 0);
	private static double input_gain = Prefs.get(Constants.inputGain, 1);
	private static double input_exposureTime = Prefs.get(Constants.inputExposureTime, 0);
	private static float input_noise = (float) Prefs.get(Constants.inputNoise, 0);

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (arg != null && arg.equals("clear"))
		{
			long memorySize = MemoryPeakResults.estimateMemorySize();
			if (memorySize > 0)
			{
				String memory = MemoryPeakResults.memorySizeString(memorySize);
				int size = MemoryPeakResults.countMemorySize();

				GenericDialog gd = new GenericDialog(TITLE);
				gd.addMessage(String.format("Do you want to remove all %d result%s from memory (%s)", size,
						(size > 1) ? "s" : "", memory));
				gd.enableYesNoCancel();
				gd.showDialog();
				if (!gd.wasOKed())
					return;

				IJ.log(String.format("Cleared %d result%s (%s)", size, (size > 1) ? "s" : "", memory));
				MemoryPeakResults.clearMemory();
				SummariseResults.clearSummaryTable();
			}
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

		IJ.showStatus("Processing outputs ...");

		Rectangle bounds = results.getBounds(true);
		boolean showDeviations = resultsSettings.showDeviations && canShowDeviations(results);
		boolean showEndFrame = canShowEndFrame(results);
		boolean showId = canShowId(results);

		// Display the configured output
		PeakResultsList output = new PeakResultsList();

		output.copySettings(results);
		//String title = results.getSource();
		//if (title == null || title.length() == 0)
		//	output.setSource(TITLE);

		addTableResults(output, showDeviations, showEndFrame);
		addImageResults(output, results.getName(), bounds, results.getNmPerPixel(), results.getGain());
		addFileResults(output, showDeviations, showEndFrame, showId);

		output.begin();

		// Process in batches to provide progress
		int batchSize = 100;
		ArrayList<PeakResult> batch = new ArrayList<PeakResult>(batchSize);
		List<PeakResult> list = results.getResults();
		int size = 0;
		for (; size < list.size(); size += batchSize)
		{
			IJ.showProgress(size, list.size());
			for (int j = size; j < FastMath.min(list.size(), size + batchSize); j++)
				batch.add(list.get(j));
			output.addAll(batch);
			batch.clear();
			if (isInterrupted())
				break;
		}
		IJ.showProgress(1);
		output.end();

		if (size > results.size())
			size = results.size();

		IJ.showStatus(String.format("Processed %d result%s", size, (size > 1) ? "s" : ""));
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
			if (r.peak != r.getEndFrame())
				return true;
		return false;
	}

	private boolean canShowId(MemoryPeakResults results)
	{
		final int id = results.getResults().get(0).getId();
		for (PeakResult r : results.getResults())
			if (id != r.getId())
				return true;
		return false;
	}

	private void addTableResults(PeakResultsList resultsList, boolean showDeviations, boolean showEndFrame)
	{
		if (resultsSettings.getResultsTable() != ResultsTable.NONE)
		{
			IJTablePeakResults r = new IJTablePeakResults(showDeviations);
			r.setPeakIdColumnName("Frame");
			r.setShowCalibratedValues(resultsSettings.getResultsTable() == ResultsTable.CALIBRATED);
			r.setShowEndFrame(showEndFrame);
			resultsList.addOutput(r);
		}
	}

	private void addImageResults(PeakResultsList resultsList, String title, Rectangle bounds, double nmPerPixel,
			double gain)
	{
		if (resultsSettings.getResultsImage() != ResultsImage.NONE)
		{
			IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(
					resultsSettings.getResultsImage(), resultsSettings.weightedImage, resultsSettings.equalisedImage,
					title, bounds, nmPerPixel, gain, resultsSettings.imageScale, resultsSettings.precision,
					ResultsMode.ADD);
			image.setRollingWindowSize(resultsSettings.imageRollingWindow);
			resultsList.addOutput(image);
		}
	}

	private void addFileResults(PeakResultsList resultsList, boolean showDeviations, boolean showEndFrame,
			boolean showId)
	{
		if (resultsSettings.resultsFilename != null && resultsSettings.resultsFilename.length() > 0)
		{
			// Remove extension
			File file = new File(resultsSettings.resultsFilename);
			int index = file.getName().lastIndexOf('.');
			if (index > 0)
			{
				resultsSettings.resultsFilename = file.getParent() + File.separatorChar +
						file.getName().substring(0, index);
			}

			// Add used extension
			resultsSettings.resultsFilename += (resultsSettings.binaryResults) ? ".bin" : ".xls";

			if (fileInput && inputFilename.equals(resultsSettings.resultsFilename))
			{
				IJ.log(TITLE + ": Input and output files are the same, skipping output ...");
				return;
			}

			// Check if file exists
			file = new File(resultsSettings.resultsFilename);
			if (file.exists())
			{
				YesNoCancelDialog d = new YesNoCancelDialog(IJ.getInstance(), TITLE, "Overwrite existing file?\n" +
						resultsSettings.resultsFilename);
				if (!d.yesPressed())
					return;
			}

			File parent = file.getParentFile();
			if (parent != null && parent.exists())
			{
				FilePeakResults r;
				if (resultsSettings.binaryResults)
					r = new BinaryFilePeakResults(resultsSettings.resultsFilename, showDeviations, showEndFrame, showId);
				else
					r = new FilePeakResults(resultsSettings.resultsFilename, showDeviations, showEndFrame, showId);
				r.setPeakIdColumnName("Frame");
				resultsList.addOutput(r);
			}
		}
	}

	private boolean showDialog()
	{
		gd = new GenericDialog(TITLE);
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
		gd.addChoice("Results_table", tableNames, tableNames[resultsSettings.getResultsTable().ordinal()]);
		gd.addCheckbox("Show_deviations", resultsSettings.showDeviations);

		gd.addMessage("--- Image output ---");
		String[] imageNames = SettingsManager.getNames((Object[]) ResultsImage.values());
		gd.addChoice("Image", imageNames, imageNames[resultsSettings.getResultsImage().ordinal()]);
		gd.addCheckbox("Weighted", resultsSettings.weightedImage);
		gd.addCheckbox("Equalised", resultsSettings.equalisedImage);
		gd.addSlider("Image_Precision (nm)", 0, 1.5, resultsSettings.precision);
		gd.addSlider("Image_Scale", 1, 15, resultsSettings.imageScale);
		gd.addNumericField("Image_Window", resultsSettings.imageRollingWindow, 0);

		gd.addMessage("--- File output ---");
		gd.addStringField("Results_file", resultsSettings.resultsFilename);
		gd.addCheckbox("Binary_results (file output only)", resultsSettings.binaryResults);
		gd.addMessage(" ");
		gd.addCheckbox("Results_in_memory (file input only)", resultsSettings.resultsInMemory);

		// Dialog to allow double click to select files using a file chooser
		if (!java.awt.GraphicsEnvironment.isHeadless())
		{
			text1 = (TextField) gd.getStringFields().get(0);
			text2 = (TextField) gd.getStringFields().get(1);
			text1.addMouseListener(this);
			text2.addMouseListener(this);
			text1.setColumns(30);
			text2.setColumns(30);
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		inputFilename = gd.getNextString();
		if (!titles.isEmpty())
			chooseRoi = gd.getNextBoolean();
		resultsSettings.setResultsTable(gd.getNextChoiceIndex());
		resultsSettings.showDeviations = gd.getNextBoolean();
		resultsSettings.setResultsImage(gd.getNextChoiceIndex());
		resultsSettings.weightedImage = gd.getNextBoolean();
		resultsSettings.equalisedImage = gd.getNextBoolean();
		resultsSettings.precision = gd.getNextNumber();
		resultsSettings.imageScale = gd.getNextNumber();
		resultsSettings.imageRollingWindow = (int) gd.getNextNumber();
		resultsSettings.resultsFilename = gd.getNextString();
		resultsSettings.binaryResults = gd.getNextBoolean();
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
				gd = new GenericDialog(TITLE);
				gd.addMessage("Select the source image for the ROI");
				gd.addChoice("Image", items, roiImage);
				gd.showDialog();
				if (gd.wasCanceled())
					return false;
				roiImage = gd.getNextChoice();
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

		return SettingsManager.saveSettings(settings);
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
	public static void addInput(GenericDialog gd, String inputOption, InputSource... inputs)
	{
		addInput(gd, "Input", inputOption, inputs);
	}

	/**
	 * Add a list of input sources to the generic dialog. The choice field will be named 'input'. If a file input option
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
	public static void addInput(GenericDialog gd, String inputName, String inputOption, InputSource... inputs)
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

		String[] options = source.toArray(new String[source.size()]);
		// Find the option
		inputOption = removeSizeSuffix(inputOption);

		int optionIndex = 0;
		for (int i = 0; i < options.length; i++)
		{
			String name = removeSizeSuffix(options[i]);
			if (name.equals(inputOption))
			{
				optionIndex = i;
				break;
			}
		}
		gd.addChoice(inputName, options, options[optionIndex]);
		if (fileInput)
			gd.addStringField("Input_file", inputFilename, 30);
	}

	private static String removeSizeSuffix(String name)
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
				for (String name : MemoryPeakResults.getResultNames())
				{
					addInputSource(source, MemoryPeakResults.getResults(name), false);
				}
				break;

			case MEMORY_TRACED:
				for (String name : MemoryPeakResults.getResultNames())
				{
					addInputSource(source, MemoryPeakResults.getResults(name), true);
				}
				break;
		}
	}

	/**
	 * Add a memory input source to the list
	 * 
	 * @param source
	 * @param memoryResults
	 * @param traced
	 *            Select only those results with ExtendedPeakResult objects
	 */
	public static void addInputSource(ArrayList<String> source, MemoryPeakResults memoryResults, boolean traced)
	{
		if (memoryResults.size() > 0)
		{
			if (traced)
			{
				// All results must be an ExtendedPeakResult
				for (PeakResult r : memoryResults.getResults())
					if (!(r instanceof ExtendedPeakResult))
						return;
			}
			source.add(memoryResults.getName() + " [" + memoryResults.size() + "]");
		}
	}

	/**
	 * Gets the name of the next input source from the dialog
	 * 
	 * @param gd
	 */
	public static String getInputSource(GenericDialog gd)
	{
		String source = gd.getNextChoice();
		return removeSizeSuffix(source);
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
			reader.setTracker(new IJTrackProgress());
			results = reader.getResults();
			reader.getTracker().progress(1.0);

			if (results != null && results.size() > 0)
			{
				// If the name contains a .tif suffix then create an image source
				if (results.getName() != null && results.getName().contains(".tif"))
				{
					results.setSource(new IJImageSource(results.getName(), results.getName()));
				}
			}
		}
		else
		{
			inputOption = removeSizeSuffix(inputOption);
			results = MemoryPeakResults.getResults(inputOption);
		}

		if (results != null && results.size() > 0 && checkCalibration)
		{
			// Check for Calibration
			Calibration calibration = results.getCalibration();
			final float noise = getNoise(results);
			// Do not check for calibration readNoise, bias, emCCD			
			if (calibration == null || calibration.nmPerPixel <= 0 || calibration.gain <= 0 ||
					calibration.exposureTime <= 0 || noise <= 0)
			{
				String msg = "partially calibrated";
				boolean convert = false;
				if (calibration == null)
				{
					// Make sure the user knows all the values have not been set
					calibration = new Calibration(0, 0, 0);
					msg = "uncalibrated";
				}

				if (calibration.nmPerPixel <= 0)
					calibration.nmPerPixel = input_nmPerPixel;
				if (calibration.gain <= 0)
					calibration.gain = input_gain;
				if (calibration.exposureTime <= 0)
					calibration.exposureTime = input_exposureTime;

				GenericDialog gd = new GenericDialog(TITLE);
				gd.addMessage("Results are " + msg);
				gd.addNumericField("Calibration (nm/px)", calibration.nmPerPixel, 2);
				gd.addNumericField("Gain (ADU/photon)", calibration.gain, 2);
				gd.addNumericField("Exposure_time (ms)", calibration.exposureTime, 2);
				if (noise <= 0)
					gd.addNumericField("Noise (ADU)", input_noise, 2);
				gd.addCheckbox("Convert_nm_to_pixels", convert);
				gd.showDialog();
				if (gd.wasCanceled())
					return null;
				input_nmPerPixel = Math.abs(gd.getNextNumber());
				input_gain = Math.abs(gd.getNextNumber());
				input_exposureTime = Math.abs(gd.getNextNumber());
				if (noise == 0)
					input_noise = Math.abs((float) gd.getNextNumber());
				convert = gd.getNextBoolean();

				Prefs.set(Constants.inputNmPerPixel, input_nmPerPixel);
				Prefs.set(Constants.inputGain, input_gain);
				Prefs.set(Constants.inputExposureTime, input_exposureTime);
				Prefs.set(Constants.inputNoise, input_noise);

				results.setCalibration(new Calibration(input_nmPerPixel, input_gain, input_exposureTime));

				if (convert && input_nmPerPixel > 0)
				{
					// Note: NSTORM stores 2xSD
					final double widthConversion = (reader != null && reader.getFormat() == FileFormat.NSTORM) ? 1.0 / (2 * input_nmPerPixel)
							: 1.0 / input_nmPerPixel;
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
		}
		return results;
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

	public void mouseClicked(MouseEvent e)
	{
		if (e.getClickCount() > 1) // Double-click
		{
			TextField text = (e.getSource() == text1) ? text1 : text2;
			String[] path = Utils.decodePath(text.getText());
			OpenDialog chooser = new OpenDialog("Coordinate file", path[0], path[1]);
			if (chooser.getFileName() != null)
			{
				text.setText(chooser.getDirectory() + chooser.getFileName());
			}
		}
	}

	public void mousePressed(MouseEvent e)
	{

	}

	public void mouseReleased(MouseEvent e)
	{

	}

	public void mouseEntered(MouseEvent e)
	{

	}

	public void mouseExited(MouseEvent e)
	{

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
}
