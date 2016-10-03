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
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitFunction;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.core.ij.Utils;
import gdsc.smlm.results.Calibration;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.YesNoCancelDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

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
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.io.File;
import java.util.Vector;

/**
 * Adjust the configuration used for fitting.
 */
// TODO - This could be incorporated into the PeakFit plugin as another run mode.
public class Configuration implements PlugIn, MouseListener, TextListener, ItemListener
{
	private static final String TITLE = "Fit Configuration";

	private boolean configurationChanged = false;

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
	private TextField textFailuresLimit;
	private Checkbox textIncludeNeighbours;
	private TextField textNeighbourHeightThreshold;
	private TextField textResidualsThreshold;
	private TextField textDuplicateDistance;
	private Checkbox textSmartFilter;
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

		String filename = SettingsManager.getSettingsFilename();

		GlobalSettings settings = SettingsManager.loadSettings(filename);
		FitEngineConfiguration config = settings.getFitEngineConfiguration();
		FitConfiguration fitConfig = config.getFitConfiguration();
		Calibration calibration = settings.getCalibration();

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);
		gd.addMessage("Configuration settings for the single-molecule localisation microscopy plugins");

		gd.addStringField("Config_file", filename, 40);
		gd.addNumericField("Calibration (nm/px)", calibration.nmPerPixel, 2);
		gd.addNumericField("Gain", calibration.gain, 2);
		gd.addCheckbox("EM-CCD", calibration.emCCD);
		gd.addNumericField("Exposure_time (ms)", calibration.exposureTime, 2);

		gd.addMessage("--- Gaussian parameters ---");
		gd.addNumericField("Initial_StdDev0", fitConfig.getInitialPeakStdDev0(), 3);
		gd.addNumericField("Initial_StdDev1", fitConfig.getInitialPeakStdDev1(), 3);
		gd.addNumericField("Initial_Angle", fitConfig.getInitialAngle(), 3);

		gd.addMessage("--- Maxima identification ---");
		String[] filterTypes = SettingsManager.getNames((Object[]) DataFilterType.values());
		gd.addChoice("Spot_filter_type", filterTypes, filterTypes[config.getDataFilterType().ordinal()]);
		String[] filterNames = SettingsManager.getNames((Object[]) DataFilter.values());
		gd.addChoice("Spot_filter", filterNames, filterNames[config.getDataFilter(0).ordinal()]);
		gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
		gd.addSlider("Search_width", 0.5, 2.5, config.getSearch());
		gd.addSlider("Border", 0.5, 2.5, config.getBorder());
		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());

		gd.addMessage("--- Gaussian fitting ---");
		Component splitLabel = gd.getMessage();

		String[] solverNames = SettingsManager.getNames((Object[]) FitSolver.values());
		gd.addChoice("Fit_solver", solverNames, solverNames[fitConfig.getFitSolver().ordinal()]);
		String[] functionNames = SettingsManager.getNames((Object[]) FitFunction.values());
		gd.addChoice("Fit_function", functionNames, functionNames[fitConfig.getFitFunction().ordinal()]);

		// Parameters specific to each Fit solver are collected in a second dialog 

		gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
		gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
		gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
		gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

		gd.addSlider("Duplicate_distance", 0, 1.5, fitConfig.getDuplicateDistance());

		gd.addMessage("--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");

		gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
		gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
		gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
		gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
		gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
		gd.addSlider("Width_factor", 1.01, 5, fitConfig.getWidthFactor());
		gd.addNumericField("Precision_threshold", fitConfig.getPrecisionThreshold(), 2);

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

			textConfigFile = texts.get(t++);
			textConfigFile.addMouseListener(this);
			textConfigFile.addTextListener(this);

			// TODO: add a value changed listener to detect when typing a new file

			textNmPerPixel = numerics.get(n++);
			textGain = numerics.get(n++);
			textEMCCD = checkboxes.get(b++);
			textExposure = numerics.get(n++);
			textInitialPeakStdDev0 = numerics.get(n++);
			textInitialPeakStdDev1 = numerics.get(n++);
			textInitialAngleD = numerics.get(n++);
			textDataFilterType = choices.get(ch++);
			textDataFilter = choices.get(ch++);
			textSmooth = numerics.get(n++);
			textSearch = numerics.get(n++);
			textBorder = numerics.get(n++);
			textFitting = numerics.get(n++);
			textFitSolver = choices.get(ch++);
			textFitFunction = choices.get(ch++);
			textFailuresLimit = numerics.get(n++);
			textIncludeNeighbours = checkboxes.get(b++);
			textNeighbourHeightThreshold = numerics.get(n++);
			textResidualsThreshold = numerics.get(n++);
			textDuplicateDistance = numerics.get(n++);
			textSmartFilter = checkboxes.get(b++);
			textCoordinateShiftFactor = numerics.get(n++);
			textSignalStrength = numerics.get(n++);
			textMinPhotons = numerics.get(n++);
			textMinWidthFactor = numerics.get(n++);
			textWidthFactor = numerics.get(n++);
			textPrecisionThreshold = numerics.get(n++);

			updateFilterInput();
			textSmartFilter.addItemListener(this);
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

		gd.showDialog();

		if (gd.wasCanceled())
			return;

		filename = gd.getNextString();

		calibration.nmPerPixel = gd.getNextNumber();
		calibration.gain = gd.getNextNumber();
		calibration.emCCD = gd.getNextBoolean();
		calibration.exposureTime = gd.getNextNumber();
		fitConfig.setInitialPeakStdDev0(gd.getNextNumber());
		fitConfig.setInitialPeakStdDev1(gd.getNextNumber());
		fitConfig.setInitialAngleD(gd.getNextNumber());
		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		config.setFitting(gd.getNextNumber());

		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		fitConfig.setFitFunction(gd.getNextChoiceIndex());

		config.setFailuresLimit((int) gd.getNextNumber());
		config.setIncludeNeighbours(gd.getNextBoolean());
		config.setNeighbourHeightThreshold(gd.getNextNumber());
		config.setResidualsThreshold(gd.getNextNumber());

		fitConfig.setDuplicateDistance(gd.getNextNumber());

		fitConfig.setSmartFilter(gd.getNextBoolean());
		fitConfig.setCoordinateShiftFactor(gd.getNextNumber());
		fitConfig.setSignalStrength(gd.getNextNumber());
		fitConfig.setMinPhotons(gd.getNextNumber());
		fitConfig.setMinWidthFactor(gd.getNextNumber());
		fitConfig.setWidthFactor(gd.getNextNumber());
		fitConfig.setPrecisionThreshold(gd.getNextNumber());

		// Check arguments
		try
		{
			Parameters.isAboveZero("nm per pixel", calibration.nmPerPixel);
			Parameters.isAboveZero("Gain", calibration.gain);
			Parameters.isAboveZero("Exposure time", calibration.exposureTime);
			Parameters.isAboveZero("Initial SD0", fitConfig.getInitialPeakStdDev0());
			Parameters.isAboveZero("Initial SD1", fitConfig.getInitialPeakStdDev1());
			Parameters.isPositive("Initial angle", fitConfig.getInitialAngleD());
			Parameters.isAboveZero("Search_width", config.getSearch());
			Parameters.isAboveZero("Fitting_width", config.getFitting());
			Parameters.isPositive("Failures limit", config.getFailuresLimit());
			Parameters.isPositive("Neighbour height threshold", config.getNeighbourHeightThreshold());
			Parameters.isPositive("Residuals threshold", config.getResidualsThreshold());
			Parameters.isPositive("Duplicate distance", fitConfig.getDuplicateDistance());
			Parameters.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
			Parameters.isPositive("Signal strength", fitConfig.getSignalStrength());
			Parameters.isPositive("Min photons", fitConfig.getMinPhotons());
			Parameters.isPositive("Min width factor", fitConfig.getMinWidthFactor());
			Parameters.isPositive("Width factor", fitConfig.getWidthFactor());
			Parameters.isPositive("Precision threshold", fitConfig.getPrecisionThreshold());
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return;
		}

		if (gd.invalidNumber())
			return;

		configurationChanged = SettingsManager.saveSettings(settings, filename);
		if (configurationChanged)
			SettingsManager.saveSettingsFilename(filename);

		if (!PeakFit.configureSmartFilter(settings, filename))
			return;
		if (!PeakFit.configureDataFilter(settings, filename, false))
			return;
		PeakFit.configureFitSolver(settings, filename, false);
	}

	public boolean isConfigurationChanged()
	{
		return configurationChanged;
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
				String[] path = Utils.decodePath(textConfigFile.getText());
				OpenDialog chooser = new OpenDialog("Configuration_File", path[0], path[1]);
				if (chooser.getFileName() != null)
				{
					String newFilename = chooser.getDirectory() + chooser.getFileName();
					textConfigFile.setText(newFilename);
				}
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
			YesNoCancelDialog d = new YesNoCancelDialog(IJ.getInstance(), TITLE, "Reload settings from file");
			d.setVisible(true);
			if (d.yesPressed())
			{
				// Reload the settings and update the GUI
				// XXX : This does not deal with loading settings into fields that are not displayed,
				// e.g. for configuring the Fit Solvers. This could be done by writing into
				// a class scope settings instance (loaded in showDialog()). However the user would not 
				// see all the changes that have been written, since the later dialogs are shown depending 
				// on what options are initially configured. 

				GlobalSettings settings = SettingsManager.unsafeLoadSettings(newFilename);
				if (settings == null)
					return;
				FitEngineConfiguration config = settings.getFitEngineConfiguration();
				FitConfiguration fitConfig = config.getFitConfiguration();
				Calibration calibration = settings.getCalibration();

				textNmPerPixel.setText("" + calibration.nmPerPixel);
				textGain.setText("" + calibration.gain);
				textEMCCD.setState(calibration.emCCD);
				textExposure.setText("" + calibration.exposureTime);
				textInitialPeakStdDev0.setText("" + fitConfig.getInitialPeakStdDev0());
				textInitialPeakStdDev1.setText("" + fitConfig.getInitialPeakStdDev1());
				textInitialAngleD.setText("" + fitConfig.getInitialAngle());
				textDataFilterType.select(config.getDataFilterType().ordinal());
				textDataFilter.select(config.getDataFilter(0).ordinal());
				textSmooth.setText("" + config.getSmooth(0));
				textSearch.setText("" + config.getSearch());
				textBorder.setText("" + config.getBorder());
				textFitting.setText("" + config.getFitting());
				textFitSolver.select(fitConfig.getFitSolver().ordinal());
				textFitFunction.select(fitConfig.getFitFunction().ordinal());
				textFailuresLimit.setText("" + config.getFailuresLimit());
				textIncludeNeighbours.setState(config.isIncludeNeighbours());
				textNeighbourHeightThreshold.setText("" + config.getNeighbourHeightThreshold());
				textResidualsThreshold.setText("" + config.getResidualsThreshold());
				textDuplicateDistance.setText("" + fitConfig.getDuplicateDistance());
				textCoordinateShiftFactor.setText("" + fitConfig.getCoordinateShiftFactor());
				textSignalStrength.setText("" + fitConfig.getSignalStrength());
				textMinPhotons.setText("" + fitConfig.getMinPhotons());
				textMinWidthFactor.setText("" + fitConfig.getMinWidthFactor());
				textWidthFactor.setText("" + fitConfig.getWidthFactor());
				textPrecisionThreshold.setText("" + fitConfig.getPrecisionThreshold());
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ItemListener#itemStateChanged(java.awt.event.ItemEvent)
	 */
	public void itemStateChanged(ItemEvent e)
	{
		if (e.getSource() instanceof Checkbox)
			updateFilterInput();
	}

	private void updateFilterInput()
	{
		if (textSmartFilter.getState())
		{
			disableEditing(textCoordinateShiftFactor);
			disableEditing(textSignalStrength);
			disableEditing(textMinPhotons);
			disableEditing(textMinWidthFactor);
			disableEditing(textWidthFactor);
			disableEditing(textPrecisionThreshold);
		}
		else
		{
			enableEditing(textCoordinateShiftFactor);
			enableEditing(textSignalStrength);
			enableEditing(textMinPhotons);
			enableEditing(textMinWidthFactor);
			enableEditing(textWidthFactor);
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
}
