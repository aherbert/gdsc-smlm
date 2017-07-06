package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
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

import javax.swing.JButton;
import javax.swing.JFileChooser;

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

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.utils.BitFlags;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.ResultsConfig.ResultsFileFormat;
import gdsc.smlm.data.config.ResultsConfig.ResultsFileSettings;
import gdsc.smlm.data.config.ResultsConfig.ResultsImageMode;
import gdsc.smlm.data.config.ResultsConfig.ResultsImageSettings;
import gdsc.smlm.data.config.ResultsConfig.ResultsImageType;
import gdsc.smlm.data.config.ResultsConfig.ResultsInMemorySettings;
import gdsc.smlm.data.config.ResultsConfig.ResultsSettings;
import gdsc.smlm.data.config.ResultsConfig.ResultsSettings.Builder;
import gdsc.smlm.data.config.ResultsConfig.ResultsTableSettings;
import gdsc.smlm.data.config.ResultsConfigHelper;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;
import gdsc.smlm.data.config.UnitConfig.IntensityUnit;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.ij.results.ImagePeakResultsFactory;
import gdsc.smlm.ij.settings.Constants;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.BinaryFilePeakResults;
import gdsc.smlm.results.Counter;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.FixedPeakResultList;
import gdsc.smlm.results.MALKFilePeakResults;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResults;
import gdsc.smlm.results.PeakResultsList;
import gdsc.smlm.results.PeakResultsReader;
import gdsc.smlm.results.ResultOption;
import gdsc.smlm.results.TSFPeakResultsWriter;
import gdsc.smlm.results.TextFilePeakResults;
import gdsc.smlm.results.procedures.PeakResultProcedure;
import gdsc.smlm.results.procedures.PeakResultProcedureX;
import ij.IJ;
import ij.Prefs;
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

	private ResultsSettings.Builder resultsSettings = ResultsSettings.newBuilder();
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
				memorySize += MemoryPeakResults.estimateMemorySize(results);
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
			Utils.log("Cleared %s (%s, %s)", count, sets, memory);
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

		IJ.showStatus("Loaded " + Utils.pleural(results.size(), "result"));

		boolean saved = false;
		if (resultsSettings.getResultsInMemorySettings().getInMemory() && fileInput)
		{
			MemoryPeakResults.addResults(results);
			saved = true;
		}

		if (!resultsSettings.getResultsTableSettings().getShowTable() &&
				resultsSettings.getResultsImageSettings().getImageType().getNumber() <= 0 &&
				Utils.isNullOrEmpty(resultsSettings.getResultsFileSettings().getResultsFilename()))
		{
			// No outputs. Error if results were not saved to memory
			if (!saved)
			{
				IJ.error(TITLE, "No output selected");
			}
			return;
		}

		Rectangle bounds = results.getBounds(true);
		boolean showDeviations = resultsSettings.getShowDeviations() && canShowDeviations(results);
		boolean showEndFrame = canShowEndFrame(results);
		boolean showId = canShowId(results);

		// Display the configured output
		PeakResultsList outputList = new PeakResultsList();

		outputList.copySettings(results);
		//String title = results.getSource();
		//if (title == null || title.length() == 0)
		//	output.setSource(TITLE);

		addTableResults(outputList, resultsSettings.getResultsTableSettings(), showDeviations, showEndFrame,
				results.is3D());
		addImageResults(outputList, resultsSettings.getResultsImageSettings(), bounds,
				(extraOptions) ? FLAG_EXTRA_OPTIONS : 0);
		addFileResults(outputList, showDeviations, showEndFrame, showId);

		IJ.showStatus("Processing outputs ...");

		// Reduce to single object for speed
		final PeakResults output = (outputList.numberOfOutputs() == 1) ? outputList.toArray()[0] : outputList;

		output.begin();

		// Note: We could add a batch iterator to the MemoryPeakResults.
		// However the speed increase will be marginal as the main time
		// taken is in processing the outputs.

		// Process in batches to provide progress
		final Counter progress = new Counter();
		final int totalProgress = results.size();
		final int batchSize = Math.max(100, totalProgress / 10);
		final FixedPeakResultList batch = new FixedPeakResultList(batchSize);
		IJ.showProgress(0);
		results.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult result)
			{
				batch.add(result);
				if (batch.size == batchSize)
				{
					if (IJ.escapePressed())
					{
						batch.clear();
						return true;
					}
					output.addAll(batch.results);
					batch.clear();
					IJ.showProgress(progress.incrementAndGet(batchSize), totalProgress);
				}
				return false;
			}
		});

		// Will be empty if interrupted
		if (batch.isNotEmpty())
			output.addAll(batch.toArray());

		IJ.showProgress(1);
		output.end();

		if (output.size() == results.size())
			IJ.showStatus("Processed " + Utils.pleural(results.size(), "result"));
		else
			IJ.showStatus(String.format("A %d/%s", output.size(), Utils.pleural(results.size(), "result")));
	}

	private boolean canShowDeviations(MemoryPeakResults results)
	{
		return results.hasDeviations();
	}

	private boolean canShowEndFrame(MemoryPeakResults results)
	{
		return results.hasEndFrame();
	}

	private boolean canShowId(MemoryPeakResults results)
	{
		return results.hasId();
	}

	public static IJTablePeakResults addTableResults(PeakResultsList resultsList, ResultsTableSettings resultsSettings,
			boolean showDeviations, boolean showEndFrame, boolean showZ)
	{
		if (resultsSettings.getShowTable())
		{
			IJTablePeakResults r = new IJTablePeakResults(showDeviations);
			r.setDistanceUnit(resultsSettings.getDistanceUnit());
			r.setIntensityUnit(resultsSettings.getIntensityUnit());
			r.setAngleUnit(resultsSettings.getAngleUnit());
			r.setComputePrecision(resultsSettings.getComputePrecision());
			r.setShowEndFrame(showEndFrame);
			r.setRoundingPrecision(resultsSettings.getRoundingPrecision());
			r.setShowZ(showZ);
			r.setShowFittingData(resultsSettings.getShowFittingData());
			r.setShowNoiseData(resultsSettings.getShowNoiseData());
			resultsList.addOutput(r);
			return r;
		}
		return null;
	}

	public static void addImageResults(PeakResultsList resultsList, ResultsImageSettings resultsSettings,
			Rectangle bounds, int flags)
	{
		if (resultsSettings.getImageType().getNumber() > 0)
		{
			IJImagePeakResults image = ImagePeakResultsFactory.createPeakResultsImage(resultsSettings.getImageType(),
					resultsSettings.getWeighted(), resultsSettings.getEqualised(), resultsList.getName(), bounds,
					resultsList.getNmPerPixel(), resultsList.getGain(), resultsSettings.getScale(),
					resultsSettings.getAveragePrecision(), ResultsImageMode.IMAGE_ADD);
			if (BitFlags.anySet(flags, FLAG_EXTRA_OPTIONS))
				image.setRollingWindowSize(resultsSettings.getRollingWindowSize());
			image.setRepaintDelay(2000);
			resultsList.addOutput(image);
		}
	}

	private void addFileResults(PeakResultsList resultsList, boolean showDeviations, boolean showEndFrame,
			boolean showId)
	{
		ResultsFileSettings resultsSettings = this.resultsSettings.getResultsFileSettings();
		if (!Utils.isNullOrEmpty(resultsSettings.getResultsFilename()))
		{
			// Remove extension
			String resultsFilename = Utils.replaceExtension(resultsSettings.getResultsFilename(),
					ResultsConfigHelper.getExtension(resultsSettings.getFileFormat()));

			if (fileInput && inputFilename.equals(resultsFilename))
			{
				IJ.log(TITLE + ": Input and output files are the same, skipping output ...");
				return;
			}

			// Update
			this.resultsSettings.getResultsFileSettingsBuilder().setResultsFilename(resultsFilename);
			SettingsManager.writeSettings(this.resultsSettings.build());

			// Check if file exists
			File file = new File(resultsFilename);
			if (file.exists())
			{
				YesNoCancelDialog d = new YesNoCancelDialog(IJ.getInstance(), TITLE,
						"Overwrite existing file?\n" + resultsFilename);
				if (!d.yesPressed())
					return;
			}

			addFileResults(resultsList, resultsSettings, resultsFilename, showDeviations, showEndFrame, showId);
		}
	}

	public static PeakResults addFileResults(PeakResultsList resultsList, ResultsFileSettings resultsSettings,
			String resultsFilename, boolean showDeviations, boolean showEndFrame, boolean showId)
	{
		File file = new File(resultsFilename);
		File parent = file.getParentFile();
		if (parent != null && parent.exists())
		{
			PeakResults r;
			switch (resultsSettings.getFileFormat())
			{
				case BINARY:
					r = new BinaryFilePeakResults(resultsFilename, showDeviations, showEndFrame, showId);
					break;
				case TEXT:
					TextFilePeakResults f = new TextFilePeakResults(resultsFilename, showDeviations, showEndFrame,
							showId);
					f.setDistanceUnit(resultsSettings.getDistanceUnit());
					f.setIntensityUnit(resultsSettings.getIntensityUnit());
					f.setAngleUnit(resultsSettings.getAngleUnit());
					f.setComputePrecision(resultsSettings.getComputePrecision());
					r = f;
					break;
				case MALK:
					r = new MALKFilePeakResults(resultsFilename);
					break;
				case TSF:
					r = new TSFPeakResultsWriter(resultsFilename);
					break;
				default:
					throw new RuntimeException("Unsupported file format: " + resultsSettings.getFileFormat());
			}
			resultsList.addOutput(r);
			return r;
		}
		return null;
	}

	private boolean showDialog()
	{
		final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		resultsSettings = SettingsManager.readResultsSettings(0).toBuilder();

		gd.addMessage("Read the Peak Results and output to a new format");

		gd.addMessage("Select the Peak Results");
		addInput(gd, inputOption, InputSource.MEMORY, InputSource.FILE);

		final Choice inputChoice = gd.getLastChoice();

		addTableResultsOptions(gd, resultsSettings);
		addImageResultsOptions(gd, resultsSettings);
		addFileResultsOptions(gd, resultsSettings, 0);
		addInMemoryResultsOptions(gd, resultsSettings);

		final Label messageLabel = (Label) gd.getMessage();
		final Checkbox saveCheckbox = gd.getLastCheckbox();

		// Hide the in-memory settings if the input is not a file
		if (Utils.isShowGenericDialog())
		{
			final Label saveLabel = gd.getLastLabel();
			ItemListener listener = new ItemListener()
			{
				public void itemStateChanged(ItemEvent e)
				{
					boolean enable = INPUT_FILE.equals(inputChoice.getSelectedItem());
					if (enable != messageLabel.isVisible())
					{
						messageLabel.setVisible(enable);
						saveCheckbox.setVisible(enable);
						saveLabel.setVisible(enable);
						gd.pack();
					}
				}
			};

			// Run once to set up
			listener.itemStateChanged(null);

			inputChoice.addItemListener(listener);
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		inputFilename = gd.getNextString();
		resultsSettings.getResultsTableSettingsBuilder().setShowTable(gd.getNextBoolean());
		resultsSettings.getResultsImageSettingsBuilder().setImageTypeValue(gd.getNextChoiceIndex());
		resultsSettings.getResultsFileSettingsBuilder().setResultsFilename(gd.getNextString());
		resultsSettings.getResultsFileSettingsBuilder().setFileFormatValue(gd.getNextChoiceIndex());
		resultsSettings.getResultsInMemorySettingsBuilder().setInMemory(gd.getNextBoolean());

		gd.collectOptions();
		
		// Check arguments
		try
		{
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

		Prefs.set(Constants.inputFilename, inputFilename);

		SettingsManager.writeSettings(resultsSettings.build());

		return true;
	}

	public static void addTableResultsOptions(final ExtendedGenericDialog gd, final Builder resultsSettings)
	{
		gd.addMessage("--- Table output ---");
		final ResultsTableSettings.Builder tableSettings = resultsSettings.getResultsTableSettingsBuilder();
		gd.addCheckbox("Show_results_table", tableSettings.getShowTable(), new OptionListener<Checkbox>()
		{
			public void collectOptions(Checkbox field)
			{
				tableSettings.setShowTable(field.getState());
				collectOptions();
			}

			public void collectOptions()
			{
				if (!tableSettings.getShowTable())
				{
					return;
				}
				ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
				egd.addChoice("Table_distance_unit", SettingsManager.getDistanceUnitNames(),
						tableSettings.getDistanceUnit().getNumber());
				egd.addChoice("Table_intensity_unit", SettingsManager.getIntensityUnitNames(),
						tableSettings.getIntensityUnit().getNumber());
				egd.addChoice("Table_angle_unit", SettingsManager.getAngleUnitNames(),
						tableSettings.getAngleUnit().getNumber());
				egd.addCheckbox("Table_show_fitting_data", tableSettings.getShowFittingData());
				egd.addCheckbox("Table_show_noise_data", tableSettings.getShowNoiseData());
				egd.addCheckbox("Table_show_precision", tableSettings.getComputePrecision());
				egd.addSlider("Table_precision", 0, 10, tableSettings.getRoundingPrecision());
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return;
				tableSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
				tableSettings.setIntensityUnitValue(egd.getNextChoiceIndex());
				tableSettings.setAngleUnitValue(egd.getNextChoiceIndex());
				tableSettings.setShowFittingData(egd.getNextBoolean());
				tableSettings.setShowNoiseData(egd.getNextBoolean());
				tableSettings.setComputePrecision(egd.getNextBoolean());
				tableSettings.setRoundingPrecision((int) egd.getNextNumber());
			}
		});
	}

	/** Use this to add extra options to the dialog. */
	public static final int FLAG_EXTRA_OPTIONS = 0x00000001;
	/** Use this to add the results directory to the file results dialog. */
	public static final int FLAG_RESULTS_DIRECTORY = 0x00000002;
	/** Use this to add the results file to the file results dialog. */
	public static final int FLAG_RESULTS_FILE = 0x00000004;

	private void addImageResultsOptions(final ExtendedGenericDialog gd, final Builder resultsSettings)
	{
		addImageResultsOptions(gd, resultsSettings, (extraOptions) ? FLAG_EXTRA_OPTIONS : 0);
	}

	public static void addImageResultsOptions(final ExtendedGenericDialog gd, final Builder resultsSettings,
			final int flags)
	{
		gd.addMessage("--- Image output ---");
		final ResultsImageSettings.Builder imageSettings = resultsSettings.getResultsImageSettingsBuilder();
		final EnumSet<ResultsImageType> requirePrecision = EnumSet.of(
				ResultsImageType.DRAW_LOCALISATIONS_AVERAGE_PRECISION,
				ResultsImageType.DRAW_INTENSITY_AVERAGE_PRECISION);
		final EnumSet<ResultsImageType> requireWeighted = EnumSet.of(ResultsImageType.DRAW_LOCALISATIONS,
				ResultsImageType.DRAW_INTENSITY, ResultsImageType.DRAW_FRAME_NUMBER, ResultsImageType.DRAW_FIT_ERROR);
		gd.addChoice("Image", SettingsManager.getResultsImageTypeNames(), imageSettings.getImageType().getNumber(),
				new OptionListener<Choice>()
				{
					public void collectOptions(Choice field)
					{
						imageSettings.setImageTypeValue(field.getSelectedIndex());
						collectOptions();
					}

					public void collectOptions()
					{
						ResultsImageType resultsImage = imageSettings.getImageType();
						if (resultsImage.getNumber() <= 0)
						{
							return;
						}
						boolean extraOptions = BitFlags.anySet(flags, FLAG_EXTRA_OPTIONS);
						ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
						if (requireWeighted.contains(resultsImage))
							egd.addCheckbox("Weighted", imageSettings.getWeighted());
						egd.addCheckbox("Equalised", imageSettings.getEqualised());
						if (requirePrecision.contains(resultsImage))
							egd.addSlider("Image_Precision (nm)", 5, 30, imageSettings.getAveragePrecision());
						egd.addSlider("Image_Scale", 1, 15, imageSettings.getScale());
						if (extraOptions)
							egd.addNumericField("Image_Window", imageSettings.getRollingWindowSize(), 0);
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return;
						if (requireWeighted.contains(resultsImage))
							imageSettings.setWeighted(egd.getNextBoolean());
						imageSettings.setEqualised(egd.getNextBoolean());
						if (requirePrecision.contains(resultsImage))
							imageSettings.setAveragePrecision(egd.getNextNumber());
						imageSettings.setScale(egd.getNextNumber());
						if (extraOptions)
							imageSettings.setRollingWindowSize((int) egd.getNextNumber());
					}
				});
	}

	public static void addFileResultsOptions(final ExtendedGenericDialog gd, final Builder resultsSettings,
			final int flags)
	{
		gd.addMessage("--- File output ---");
		final ResultsFileSettings.Builder fileSettings = resultsSettings.getResultsFileSettingsBuilder();
		if (BitFlags.anySet(flags, FLAG_RESULTS_DIRECTORY))
			gd.addDirectoryField("Results_directory", fileSettings.getResultsDirectory());
		if (BitFlags.anySet(flags, FLAG_RESULTS_FILE))
			gd.addFilenameField("Results_file", fileSettings.getResultsFilename());
		else
			// Do not add a results file to prevent constant overwrite messages
			gd.addFilenameField("Results_file", "");
		gd.addChoice("Results_format", SettingsManager.getResultsFileFormatNames(),
				fileSettings.getFileFormat().getNumber(), new OptionListener<Choice>()
				{
					public void collectOptions(Choice field)
					{
						fileSettings.setFileFormatValue(field.getSelectedIndex());
						collectOptions();
					}

					public void collectOptions()
					{
						ResultsFileFormat resultsFileFormat = fileSettings.getFileFormat();
						if (!ResultsConfigHelper.isGDSC(resultsFileFormat))
						{
							return;
						}
						ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
						if (resultsFileFormat == ResultsFileFormat.TEXT)
						{
							egd.addChoice("File_distance_unit", SettingsManager.getDistanceUnitNames(),
									fileSettings.getDistanceUnit().ordinal());
							egd.addChoice("File_intensity_unit", SettingsManager.getIntensityUnitNames(),
									fileSettings.getIntensityUnit().ordinal());
							egd.addChoice("File_angle_unit", SettingsManager.getAngleUnitNames(),
									fileSettings.getAngleUnit().ordinal());
							egd.addCheckbox("File_show_precision", fileSettings.getComputePrecision());
						}
						egd.addCheckbox("Show_deviations", resultsSettings.getShowDeviations());
						egd.showDialog(true, gd);
						if (egd.wasCanceled())
							return;
						fileSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
						fileSettings.setIntensityUnitValue(egd.getNextChoiceIndex());
						fileSettings.setAngleUnitValue(egd.getNextChoiceIndex());
						fileSettings.setComputePrecision(egd.getNextBoolean());
						resultsSettings.setShowDeviations(egd.getNextBoolean());
					}
				});
	}

	public static void addInMemoryResultsOptions(final ExtendedGenericDialog gd, final Builder resultsSettings)
	{
		gd.addMessage("--- Memory output ---");
		final ResultsInMemorySettings.Builder memorySettings = resultsSettings.getResultsInMemorySettingsBuilder();
		gd.addCheckbox("Save_to_memory", memorySettings.getInMemory());
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
		return memoryResults.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult r)
			{
				return r.getFrame() < r.getEndFrame();
			}
		});
	}

	/**
	 * Check for any IDs above zero.
	 * 
	 * @param memoryResults
	 * @return True if any results have IDs above zero
	 */
	public static boolean hasID(MemoryPeakResults memoryResults)
	{
		return memoryResults.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult r)
			{
				return r.getId() > 0;
			}
		});
	}

	/**
	 * Check for all IDs above zero.
	 * 
	 * @param memoryResults
	 * @return True if all results have IDs above zero
	 */
	public static boolean isID(MemoryPeakResults memoryResults)
	{
		return !memoryResults.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult r)
			{
				return r.getId() <= 0;
			}
		});
	}

	/**
	 * All results must be an ExtendedPeakResult.
	 * 
	 * @param memoryResults
	 * @return True if all are an ExtendedPeakResult
	 */
	public static boolean isExtended(MemoryPeakResults memoryResults)
	{
		return !memoryResults.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult r)
			{
				return !(r instanceof ExtendedPeakResult);
			}
		});
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
	 * Load the results from the named input option. If the results are not empty then a check can be made for
	 * calibration, and data using the legacy standard units (distance in Pixel and intensity in Count).
	 * <p>
	 * If the calibration cannot be obtained or the units are incorrect then the null will be returned.
	 *
	 * @param inputOption
	 *            the input option
	 * @param checkCalibration
	 *            Set to true to ensure the results have a valid calibration
	 * @return the results
	 */
	public static MemoryPeakResults loadInputResults(String inputOption, boolean checkCalibration)
	{
		return loadInputResults(inputOption, checkCalibration, DistanceUnit.PIXEL, IntensityUnit.COUNT);
	}

	/**
	 * Load the results from the named input option. If the results are not empty then a check can be made for
	 * calibration, and data using the specified units. If the calibration cannot be obtained or the units are incorrect
	 * then the null will be returned.
	 *
	 * @param inputOption
	 *            the input option
	 * @param checkCalibration
	 *            Set to true to ensure the results have a valid calibration
	 * @param distanceUnit
	 *            the required distance unit for the results
	 * @param intensityUnit
	 *            the required intensity unit for the results
	 * @return the results
	 */
	public static MemoryPeakResults loadInputResults(String inputOption, boolean checkCalibration,
			DistanceUnit distanceUnit, IntensityUnit intensityUnit)
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

		try
		{
			if (results == null)
				return null;
			if (results.isEmpty())
				return results;

			if (checkCalibration)
			{
				if (!checkCalibration(results, reader))
					return null;
			}
			if (distanceUnit != null && results.getDistanceUnit() != distanceUnit)
				return null;
			if (intensityUnit != null && results.getIntensityUnit() != intensityUnit)
				return null;
		}
		finally
		{
			IJ.showStatus("");
		}
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
		CalibrationWriter calibration = results.getCalibrationWriterSafe();
		String msg = "partially calibrated";
		if (results.hasCalibration())
		{
			// Make sure the user knows all the values have not been set
			msg = "uncalibrated";
		}

		// Only check for essential calibration settings (i.e. not readNoise, bias, emCCD, amplification)
		if (!calibration.hasNmPerPixel() || !calibration.hasGain() || !calibration.hasExposureTime())
		{
			final float noise = getNoise(results);

			if (!calibration.hasNmPerPixel())
				calibration.setNmPerPixel(input_nmPerPixel);
			if (!calibration.hasGain())
				calibration.setGain(input_gain);
			if (!calibration.hasExposureTime())
				calibration.setExposureTime(input_exposureTime);

			Rectangle2D.Float dataBounds = results.getDataBounds(null);

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
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			input_nmPerPixel = Math.abs(gd.getNextNumber());
			input_gain = Math.abs(gd.getNextNumber());
			input_exposureTime = Math.abs(gd.getNextNumber());
			if (noise <= 0)
				input_noise = Math.abs((float) gd.getNextNumber());

			Prefs.set(Constants.inputNmPerPixel, input_nmPerPixel);
			Prefs.set(Constants.inputGain, input_gain);
			Prefs.set(Constants.inputExposureTime, input_exposureTime);
			Prefs.set(Constants.inputNoise, input_noise);

			calibration.setNmPerPixel(input_nmPerPixel);
			calibration.setGain(input_gain);
			calibration.setExposureTime(input_exposureTime);

			results.setCalibration(calibration.getCalibration());

			if (noise == 0)
			{
				results.forEach(new PeakResultProcedure()
				{
					public void execute(PeakResult p)
					{
						p.noise = input_noise;
					}
				});
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
		final float[] noise = new float[1];
		results.forEach(new PeakResultProcedureX()
		{
			public boolean execute(PeakResult r)
			{
				if (r.noise != 0)
				{
					noise[0] = r.noise;
					return true;
				}
				return false;
			}
		});
		return noise[0];
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
		return loadInputResults(inputOption, true, null, null);
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
			Recorder.recordOption("image",
					SettingsManager.getResultsImageTypeNames()[ResultsImageType.DRAW_NONE.getNumber()]);
			Recorder.recordOption("results_file", "[]");
			Recorder.recordOption("save_to_memory");
		}
		MemoryPeakResults results = loadInputResults(INPUT_FILE, true, null, null);
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
