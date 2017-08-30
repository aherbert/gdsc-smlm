package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Color;
import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.SystemColor;
import java.awt.TextField;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.Vector;

import gdsc.core.ij.Utils;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

/**
 * Adjust the configuration used for fitting.
 */
public class Configuration implements PlugIn, ItemListener
{
	private static final String TITLE = "Fit Configuration";
	private static String templateFilename = "";
	private static String notes = "";

	private FitEngineConfiguration config;
	private FitConfiguration fitConfig;

	private boolean configurationChanged = false;

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
	private TextField textFailuresLimit;
	private Checkbox textIncludeNeighbours;
	private TextField textNeighbourHeightThreshold;
	private TextField textResidualsThreshold;
	private TextField textDuplicateDistance;
	private Checkbox textSmartFilter;
	private Checkbox textDisableSimpleFilter;
	private TextField textCoordinateShiftFactor;
	private TextField textSignalStrength;
	private TextField textMinPhotons;
	private TextField textPrecisionThreshold;
	private TextField textMinWidthFactor;
	private TextField textWidthFactor;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		showDialog();
	}

	/**
	 * Show the current properties
	 */
	@SuppressWarnings("unchecked")
	public void showDialog()
	{
		configurationChanged = false;

		config = SettingsManager.readFitEngineConfiguration(0);
		fitConfig = config.getFitConfiguration();
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Configuration settings for the single-molecule localisation microscopy plugins");

		String[] templates = ConfigurationTemplate.getTemplateNames(true);
		gd.addChoice("Template", templates, templates[0]);

		PeakFit.addCameraOptions(gd, calibration);
		gd.addNumericField("Calibration (nm/px)", calibration.getNmPerPixel(), 2);
		gd.addNumericField("Exposure_time (ms)", calibration.getExposureTime(), 2);

		gd.addMessage("--- Gaussian parameters ---");
		PeakFit.addPSFOptions(gd, fitConfig);

		gd.addMessage("--- Maxima identification ---");
		gd.addChoice("Spot_filter_type", SettingsManager.getDataFilterTypeNames(),
				config.getDataFilterType().ordinal());
		gd.addChoice("Spot_filter", SettingsManager.getDataFilterMethodNames(),
				config.getDataFilterMethod(0).ordinal());
		gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
		gd.addSlider("Search_width", 0.5, 2.5, config.getSearch());
		gd.addSlider("Border", 0.5, 2.5, config.getBorder());
		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());

		gd.addMessage("--- Gaussian fitting ---");
		Component splitLabel = gd.getMessage();

		gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());

		// Parameters specific to each Fit solver are collected in a second dialog 

		gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
		gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
		gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
		gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

		gd.addSlider("Duplicate_distance", 0, 1.5, config.getDuplicateDistance());

		gd.addMessage("--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");

		gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
		gd.addCheckbox("Disable_simple_filter", fitConfig.isDisableSimpleFilter());
		gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
		gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
		gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
		gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
		gd.addSlider("Width_factor", 1.01, 5, fitConfig.getMaxWidthFactor());
		gd.addNumericField("Precision_threshold", fitConfig.getPrecisionThreshold(), 2);

		// Add a mouse listener to the config file field
		if (Utils.isShowGenericDialog())
		{
			Vector<TextField> numerics = (Vector<TextField>) gd.getNumericFields();
			Vector<Checkbox> checkboxes = (Vector<Checkbox>) gd.getCheckboxes();
			Vector<Choice> choices = (Vector<Choice>) gd.getChoices();

			int n = 0;
			int b = 0;
			int ch = 0;

			Choice textTemplate = choices.get(ch++);
			textTemplate.addItemListener(this);

			textCameraType = choices.get(ch++);
			textNmPerPixel = numerics.get(n++);
			textExposure = numerics.get(n++);
			textPSF = choices.get(ch++);
			textDataFilterType = choices.get(ch++);
			textDataFilterMethod = choices.get(ch++);
			textSmooth = numerics.get(n++);
			textSearch = numerics.get(n++);
			textBorder = numerics.get(n++);
			textFitting = numerics.get(n++);
			textFitSolver = choices.get(ch++);
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
			textMinWidthFactor = numerics.get(n++);
			textWidthFactor = numerics.get(n++);
			textPrecisionThreshold = numerics.get(n++);

			updateFilterInput();
			textSmartFilter.addItemListener(this);
			textDisableSimpleFilter.addItemListener(this);
		}

		if (gd.getLayout() != null)
		{
			GridBagLayout grid = (GridBagLayout) gd.getLayout();

			int xOffset = 0, yOffset = 0;
			int lastY = -1, rowCount = 0;
			for (Component comp : gd.getComponents())
			{
				// Check if this should be the second major column
				if (comp == splitLabel)
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

		gd.enableYesNoCancel("Save", "Save template");
		//gd.setCancelLabel("Close");

		gd.showDialog();

		if (gd.wasCanceled())
			return;

		// In case a template update the calibration
		calibration = fitConfig.getCalibrationWriter();

		// Ignore the template
		gd.getNextChoice();

		calibration.setCameraType(SettingsManager.getCameraTypeValues()[gd.getNextChoiceIndex()]);
		calibration.setNmPerPixel(gd.getNextNumber());
		calibration.setExposureTime(gd.getNextNumber());
		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), false, 0);
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		config.setFitting(gd.getNextNumber());

		fitConfig.setFitSolver(gd.getNextChoiceIndex());

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
		fitConfig.setMinWidthFactor(gd.getNextNumber());
		fitConfig.setWidthFactor(gd.getNextNumber());
		fitConfig.setPrecisionThreshold(gd.getNextNumber());

		gd.collectOptions();

		// Check arguments
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.getNmPerPixel());
			Parameters.isAboveZero("Gain", calibration.getCountPerPhoton());
			Parameters.isAboveZero("Exposure time", calibration.getExposureTime());
			Parameters.isAboveZero("Initial SD0", fitConfig.getInitialXSD());
			Parameters.isAboveZero("Initial SD1", fitConfig.getInitialYSD());
			Parameters.isAboveZero("Search_width", config.getSearch());
			Parameters.isAboveZero("Fitting_width", config.getFitting());
			Parameters.isPositive("Failures limit", config.getFailuresLimit());
			Parameters.isPositive("Neighbour height threshold", config.getNeighbourHeightThreshold());
			Parameters.isPositive("Residuals threshold", config.getResidualsThreshold());
			Parameters.isPositive("Duplicate distance", config.getDuplicateDistance());
			Parameters.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
			Parameters.isPositive("Signal strength", fitConfig.getSignalStrength());
			Parameters.isPositive("Min photons", fitConfig.getMinPhotons());
			Parameters.isPositive("Min width factor", fitConfig.getMinWidthFactor());
			Parameters.isPositive("Width factor", fitConfig.getMaxWidthFactor());
			Parameters.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return;
		}

		if (gd.invalidNumber())
			return;

		boolean saveToFile = !gd.wasOKed();

		int flags = PeakFit.FLAG_NO_SAVE;
		if (!PeakFit.configureSmartFilter(config, flags))
			return;
		if (!PeakFit.configureDataFilter(config, flags))
			return;
		PeakFit.configureFitSolver(config, 0, 0, flags);

		if (saveToFile)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addFilenameField("Template_filename", templateFilename);
			gd.addMessage("Add notes to the template ...");
			gd.addTextAreas(notes, null, 10, 60);
			gd.showDialog();
			if (gd.wasCanceled())
				return;
			String filename = gd.getNextString();
			notes = gd.getNextText();

			if (filename != null)
			{
				templateFilename = Utils.replaceExtension(filename, ".txt");
				File file = new File(templateFilename);
				String name = Utils.removeExtension(file.getName());
				TemplateSettings.Builder settings = TemplateSettings.newBuilder();
				settings.addNotes(notes);
				settings.setCalibration(fitConfig.getCalibration());
				settings.setFitEngineSettings(config.getFitEngineSettings());
				// Note: No results settings are currently supported
				settings.setPsf(fitConfig.getPSF());
				if (!ConfigurationTemplate.saveTemplate(name, settings.build(), file))
					IJ.error(TITLE, "Failed to save to file: " + templateFilename);
			}
		}
		else
		{
			SettingsManager.writeSettings(config, 0);
			configurationChanged = true;
		}
	}

	public boolean isConfigurationChanged()
	{
		return configurationChanged;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ItemListener#itemStateChanged(java.awt.event.ItemEvent)
	 */
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
				IJ.log("Applying template: " + templateName);
				File file = ConfigurationTemplate.getTemplateFile(templateName);
				if (file != null)
					templateFilename = file.getPath();

				StringBuilder sb = new StringBuilder();
				for (String note : template.getNotesList())
				{
					sb.append(note);
					if (!note.endsWith("\n"))
						sb.append("\n");
					IJ.log(note);
				}
				notes = sb.toString();

				boolean custom = ConfigurationTemplate.isCustomTemplate(templateName);
				if (template.hasCalibration())
				{
					refreshSettings(template.getCalibration());
				}
				if (template.hasPsf())
				{
					refreshSettings(template.getPsf(), custom);
				}
				if (template.hasFitEngineSettings())
				{
					refreshSettings(template.getFitEngineSettings(), custom);
				}
			}
		}
		else if (e.getSource() instanceof Checkbox)
		{
			if (e.getSource() == textSmartFilter)
			{
				textDisableSimpleFilter.setState(textSmartFilter.getState());
			}
			updateFilterInput();
		}
	}

	private void updateFilterInput()
	{
		if (textDisableSimpleFilter.getState())
		{
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

	private void refreshSettings(Calibration cal)
	{
		if (cal == null)
			return;

		// Do not use set() as we support merging a partial calibration
		fitConfig.mergeCalibration(cal);
		CalibrationWriter calibration = fitConfig.getCalibrationWriter();

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
		textSmooth.setText("" + config.getSmooth(0));
		textSearch.setText("" + config.getSearch());
		textBorder.setText("" + config.getBorder());
		textFitting.setText("" + config.getFitting());
		textFitSolver.select(SettingsManager.getFitSolverNames()[fitConfig.getFitSolver().ordinal()]);
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
	}
}
