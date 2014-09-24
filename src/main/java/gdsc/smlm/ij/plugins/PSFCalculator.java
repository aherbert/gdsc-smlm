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
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.utils.AiryPattern;
import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Label;
import java.awt.SystemColor;
import java.awt.TextField;

import xal.tools.math.BesselFunction;

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

	/**
	 * The factor for convering a Gaussian standard deviation to Full Width at Half Maxima (FWHM)
	 */
	public static final double SD_TO_FWHM_FACTOR = (2.0 * Math.sqrt(2.0 * Math.log(2.0)));

	private PSFCalculatorSettings settings;
	private GenericDialog gd;
	private Label pixelPitchLabel;
	private TextField widthNmText;
	private TextField widthPixelsText;
	private TextField sdNmText;
	private TextField sdPixelsText;

	// Used for the PSF profile plot
	private double[] x = null, y, y2;

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

		plotProfile(calculateAiryWidth(settings.pixelPitch, settings.magnification * settings.beamExpander,
				settings.wavelength, settings.numericalAperture, settings.proportionalityFactor));

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
	 *            Microscope numerical aperture (NA)
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 * @return the SD in nm
	 */
	public static double calculateStdDev(double wavelength, double numericalAperture, double proportionalityFactor)
	{
		double fwhm = calculateWidth(wavelength, numericalAperture, proportionalityFactor);
		return fwhm / SD_TO_FWHM_FACTOR;
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
	 *            Microscope numerical aperture (NA)
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 * @return the SD in pixels
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
	 * @return the size in pixels
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
	 *            Microscope numerical aperture (NA)
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 * @return the width in nm
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
	 *            Microscope numerical aperture (NA)
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2NA
	 * @return the width in pixels
	 */
	public static double calculateWidth(double pixelPitch, double magnification, double wavelength,
			double numericalAperture, double proportionalityFactor)
	{
		double fwhm = calculateWidth(wavelength, numericalAperture, proportionalityFactor);
		return convertToPixels(pixelPitch, magnification, fwhm);
	}

	/**
	 * Calculates the PSF peak width for the Airy disk.
	 * <p>
	 * Width is proportional to the lambda / (2*pi*NA). This method computes the width as proportionalityFactor * lambda
	 * / (2 * pi * NA).
	 * <p>
	 * The proportionality factor should be determined by fitting the PSF of quantum dots on the microscope.
	 * 
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture (NA)
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2piNA
	 * @return the width in nm
	 */
	public static double calculateAiryWidth(double wavelength, double numericalAperture, double proportionalityFactor)
	{
		double width = proportionalityFactor * wavelength / (2.0 * Math.PI * numericalAperture);
		return width;
	}

	/**
	 * Calculates the PSF peak width for the Airy disk.
	 * 
	 * @param pixelPitch
	 *            Camera pixel pitch in micrometers (nm)
	 * @param magnification
	 *            Objective magnification
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture (NA)
	 * @param proportionalityFactor
	 *            Expresses the proportional relationship between the diffraction limit and lambda/2piNA
	 * @return the width in pixels
	 */
	public static double calculateAiryWidth(double pixelPitch, double magnification, double wavelength,
			double numericalAperture, double proportionalityFactor)
	{
		double width = calculateAiryWidth(wavelength, numericalAperture, proportionalityFactor);
		return convertToPixels(pixelPitch, magnification, width);
	}

	private boolean lock = false;

	private synchronized boolean aquireLock()
	{
		if (lock)
			return false;
		return lock = true;
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		if (e == null)
			return true;
		Object o = e.getSource();
		if (o == sdNmText || o == sdPixelsText)
			return true;
		if (widthNmText != null)
		{
			if (o == pixelPitchLabel || o == widthNmText || o == widthPixelsText)
				return true;
		}

		if (aquireLock())
		{
			try
			{
				// Continue while the parameter is changing
				boolean parametersChanged = true;
				while (parametersChanged)
				{
					if (!readDialog())
						return false;

					// Store the parameters to be processed
					double pixelPitch = settings.pixelPitch;
					double magnification = settings.magnification;
					double beamExpander = settings.beamExpander;
					double wavelength = settings.wavelength;
					double numericalAperture = settings.numericalAperture;
					double proportionalityFactor = settings.proportionalityFactor;

					// Do something with parameters
					if (widthNmText != null)
					{
						pixelPitchLabel.setText(getPixelPitchLabel());
						widthNmText.setText(IJ.d2s(
								calculateWidth(wavelength, numericalAperture, proportionalityFactor), 3));
						widthPixelsText.setText(IJ.d2s(
								calculateWidth(pixelPitch, magnification * beamExpander, wavelength, numericalAperture,
										proportionalityFactor), 3));
					}
					sdNmText.setText(IJ.d2s(calculateStdDev(wavelength, numericalAperture, proportionalityFactor), 3));
					sdPixelsText.setText(IJ.d2s(
							calculateStdDev(pixelPitch, magnification * beamExpander, wavelength, numericalAperture,
									proportionalityFactor), 3));

					plotProfile(calculateAiryWidth(pixelPitch, magnification * beamExpander, wavelength,
							numericalAperture, proportionalityFactor));

					// Check if the parameters have changed again
					parametersChanged = (pixelPitch != settings.pixelPitch) ||
							(magnification != settings.magnification) || (beamExpander != settings.beamExpander) ||
							(wavelength != settings.wavelength) || (numericalAperture != settings.numericalAperture) ||
							(proportionalityFactor != settings.proportionalityFactor);
				}
			}
			finally
			{
				// Ensure the running flag is reset
				lock = false;
			}
		}

		return true;
	}

	private void plotProfile(double scale)
	{
		if (x == null)
		{
			x = Utils.newArray(200, -10, 0.1);
			y = new double[x.length];
			y2 = new double[x.length];
			for (int i = 0; i < x.length; i++)
			{
				y[i] = AiryPattern.intensity(x[i]);
				y2[i] = AiryPattern.intensityGaussian(x[i]);
			}
		}
		double[] x2 = new double[x.length];
		for (int i = 0; i < x2.length; i++)
		{
			x2[i] = x[i] * scale;
		}
		String title = "PSF profile";
		Plot p = new Plot(title, "px", "", x2, y);
		p.setColor(Color.RED);
		p.addPoints(x2, y2, Plot.LINE);
		p.setColor(Color.BLUE);
		Utils.display(title, p);
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
