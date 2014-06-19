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

import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.PSFCalculatorSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.awt.AWTEvent;
import java.awt.Label;
import java.awt.SystemColor;
import java.awt.TextField;

/**
 * Calculates the expected PSF width for a Gaussian approximation to the Airy disk.
 */
public class PSFCalculator implements PlugIn, DialogListener
{
	/**
	 * This is the proportionality factor (f) calculated from fitting quantum dots on the imaging microscope at the
	 * University of Sussex. The Gaussian PSF Standard Deviation = f * lambda / (2 * NA).
	 */
	public static final double DEFAULT_PROPORTIONALITY_FACTOR = 1.52;

	private PSFCalculatorSettings settings;
	private GenericDialog gd;
	private Label pixelPitchLabel;
	private TextField widthNmText;
	private TextField widthPixelsText;
	private TextField sdNmText;
	private TextField sdPixelsText;

	public void run(String arg)
	{
		GlobalSettings globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getPsfCalculatorSettings();

		double sd = calculate(settings, false);
		if (sd < 0)
			return;

		globalSettings.getFitEngineConfiguration().getFitConfiguration().setInitialPeakStdDev((float) sd);
		globalSettings.getFitEngineConfiguration().getFitConfiguration().setInitialAngle(0);
		globalSettings.getCalibration().nmPerPixel = getPixelPitch();
		SettingsManager.saveSettings(globalSettings);
	}

	/**
	 * Present an interactive dialog that allows the user to calculate the Gaussian PSF standard deviation using the
	 * provided settings.
	 * 
	 * @param settings
	 * @param simpleMode
	 *            Only present a wavelength, NA and proportionality factor fields.
	 * @return the PSF standard deviation
	 */
	public double calculate(PSFCalculatorSettings settings, boolean simpleMode)
	{
		gd = new GenericDialog("PSF Calculator");
		gd.addHelp(About.HELP_URL);

		this.settings = settings;

		if (!simpleMode)
		{
			gd.addNumericField("Pixel_pitch (um)", settings.pixelPitch, 2);
			gd.addNumericField("Magnification", settings.magnification, 0);
			gd.addNumericField("Beam_Expander", settings.beamExpander, 2);
		}
		gd.addMessage(getPixelPitchLabel(settings.pixelPitch * 1000.0));

		gd.addNumericField("Wavelength (nm)", settings.wavelength, 0);
		gd.addNumericField("Numerical_Aperture (NA)", settings.numericalAperture, 2);
		gd.addNumericField("Proportionality_factor", settings.proportionalityFactor, 2);

		if (!simpleMode)
		{
			gd.addNumericField("Width (nm)",
					calculateWidth(settings.wavelength, settings.numericalAperture, settings.proportionalityFactor), 3);
			gd.addNumericField(
					"Width (pixels)",
					calculateWidth(settings.pixelPitch, settings.magnification * settings.beamExpander,
							settings.wavelength, settings.numericalAperture, settings.proportionalityFactor), 3);
		}
		gd.addNumericField("StdDev (nm)",
				calculateStdDev(settings.wavelength, settings.numericalAperture, settings.proportionalityFactor), 3);
		gd.addNumericField(
				"StdDev (pixels)",
				calculateStdDev(settings.pixelPitch, settings.magnification * settings.beamExpander,
						settings.wavelength, settings.numericalAperture, settings.proportionalityFactor), 3);

		if (!simpleMode)
		{
			pixelPitchLabel = (Label) gd.getMessage();
			pixelPitchLabel.setText(getPixelPitchLabel());
			widthNmText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 4);
			widthPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 3);
		}
		sdNmText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 2);
		sdPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 1);

		//widthNmText
		if (!simpleMode)
		{
			disableEditing(widthNmText);
			disableEditing(widthPixelsText);
		}
		disableEditing(sdNmText);
		disableEditing(sdPixelsText);

		if (!simpleMode)
			gd.addMessage("Save StdDev pixel width to the fitting properties");

		gd.addDialogListener(this);

		gd.showDialog();

		if (gd.wasCanceled())
		{
			return -1;
		}

		return calculateStdDev(settings.pixelPitch, settings.magnification * settings.beamExpander,
				settings.wavelength, settings.numericalAperture, settings.proportionalityFactor);
	}

	private void disableEditing(TextField textField)
	{
		textField.setEditable(false);
		textField.setBackground(SystemColor.control);
	}

	private boolean readDialog()
	{
		if (widthNmText != null)
		{
			settings.pixelPitch = gd.getNextNumber();
			settings.magnification = gd.getNextNumber();
			settings.beamExpander = gd.getNextNumber();
		}
		settings.wavelength = gd.getNextNumber();
		settings.numericalAperture = gd.getNextNumber();
		settings.proportionalityFactor = gd.getNextNumber();

		// Check arguments
		try
		{
			if (widthNmText != null)
			{
				Parameters.isAboveZero("Pixel pitch", settings.pixelPitch);
				Parameters.isAboveZero("Magnification", settings.magnification);
				Parameters.isEqualOrAbove("Beam expander", settings.beamExpander, 1);
			}
			Parameters.isAboveZero("Wavelength", settings.wavelength);
			Parameters.isAboveZero("Numerical aperture", settings.numericalAperture);
			Parameters.isEqualOrAbove("Proportionality factor", settings.proportionalityFactor, 1);
		}
		catch (IllegalArgumentException ex)
		{
			// Q. Is logging the error necessary given that we will be in a live update preview?
			//IJ.log(ex.getMessage());
			return false;
		}

		return !gd.invalidNumber();
	}

	/**
	 * Calculates the expected PSF standard deviation (nm) for a Gaussian approximation to the Airy disk.
	 * <p>
	 * <a href="http://en.wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile">http://en.
	 * wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile</a>
	 * 
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 */
	public static double calculateStdDev(double wavelength, double numericalAperture, double proportionalityFactor)
	{
		double fwhm = calculateWidth(wavelength, numericalAperture, proportionalityFactor);
		return fwhm / (2.0 * Math.sqrt(2.0 * Math.log(2.0)));
	}

	/**
	 * Calculates the expected PSF standard deviation (pixels) for a Gaussian approximation to the Airy disk.
	 * 
	 * @param pixelPitch
	 *            Camera pixel pitch in micrometers (um)
	 * @param magnification
	 *            Objective magnification
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 */
	public static double calculateStdDev(double pixelPitch, double magnification, double wavelength,
			double numericalAperture, double proportionalityFactor)
	{
		double sd = calculateStdDev(wavelength, numericalAperture, proportionalityFactor);
		return convertToPixels(pixelPitch, magnification, sd);
	}

	/**
	 * @param pixelPitch
	 *            Camera pixel pitch in micrometers (um)
	 * @param magnification
	 *            Objective magnification
	 * @param size
	 *            Size in nm
	 * @return
	 */
	private static double convertToPixels(double pixelPitch, double magnification, double size)
	{
		double pixelPitchInNm = pixelPitch * 1000 / magnification;
		return size / pixelPitchInNm;
	}

	/**
	 * Calculates the expected PSF peak width at half maximum for a Gaussian approximation to the Airy disk.
	 * <p>
	 * FWHM is proportional to the lambda / 2NA. This method computes the FWHM as proportionalityFactor * lambda / 2NA.
	 * <p>
	 * The proportionality factor should be determined by fitting the PSF of quantum dots on the microscope.
	 * 
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 */
	public static double calculateWidth(double wavelength, double numericalAperture, double proportionalityFactor)
	{
		// FWHM is proportional to the lambda / 2NA. 
		double fwhm = proportionalityFactor * wavelength / (2.0 * numericalAperture);
		return fwhm;
	}

	/**
	 * Calculates the expected PSF peak width at half maximum for a Gaussian approximation to the Airy disk.
	 * 
	 * @param pixelPitch
	 *            Camera pixel pitch in micrometers (nm)
	 * @param magnification
	 *            Objective magnification
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 * @return
	 */
	public static double calculateWidth(double pixelPitch, double magnification, double wavelength,
			double numericalAperture, double proportionalityFactor)
	{
		double fwhm = calculateWidth(wavelength, numericalAperture, proportionalityFactor);
		return convertToPixels(pixelPitch, magnification, fwhm);
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		if (e == null || e.getSource() == sdNmText || e.getSource() == sdPixelsText || e.getSource() == widthPixelsText)
			return true;

		// Update width
		if (readDialog())
		{
			if (widthNmText != null)
			{
				pixelPitchLabel.setText(getPixelPitchLabel());
				widthNmText.setText(IJ
						.d2s(calculateWidth(settings.wavelength, settings.numericalAperture,
								settings.proportionalityFactor), 3));
				widthPixelsText.setText(IJ.d2s(
						calculateWidth(settings.pixelPitch, settings.magnification * settings.beamExpander,
								settings.wavelength, settings.numericalAperture, settings.proportionalityFactor), 3));
			}
			sdNmText.setText(IJ
					.d2s(calculateStdDev(settings.wavelength, settings.numericalAperture,
							settings.proportionalityFactor), 3));
			sdPixelsText.setText(IJ.d2s(
					calculateStdDev(settings.pixelPitch, settings.magnification * settings.beamExpander,
							settings.wavelength, settings.numericalAperture, settings.proportionalityFactor), 3));
			return true;
		}

		return false;
	}

	private String getPixelPitchLabel(double pixelPitch)
	{
		return "Pixel pitch (nm) = " + IJ.d2s(pixelPitch, 3);
	}

	private String getPixelPitchLabel()
	{
		return getPixelPitchLabel(getPixelPitch());
	}

	private double getPixelPitch()
	{
		return settings.pixelPitch * 1000 / (settings.magnification * settings.beamExpander);
	}
}
