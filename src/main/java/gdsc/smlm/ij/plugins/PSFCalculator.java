package gdsc.smlm.ij.plugins;

import gdsc.smlm.function.gaussian.GaussianFunction;

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
import gdsc.core.ij.Utils;
import gdsc.smlm.model.AiryPattern;
import ij.IJ;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.plugin.PlugIn;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Label;
import java.awt.SystemColor;
import java.awt.TextField;

import org.apache.commons.math3.util.FastMath;

/**
 * Calculates the expected PSF width for a Gaussian approximation to the Airy disk.
 */
public class PSFCalculator implements PlugIn, DialogListener
{
	private static final String TITLE = "PSF Calculator";
	/**
	 * This is the factor (f) for scaling an Airy width to a Gaussian approximation. The Gaussian PSF
	 * Standard Deviation = f * lambda / (2 * pi * NA).
	 */
	public static final double AIRY_TO_GAUSSIAN = 1.323;

	private PSFCalculatorSettings settings;
	private GenericDialog gd;
	private Label abbeLimitLabel;
	private Label pixelPitchLabel;
	private TextField widthNmText;
	private TextField widthPixelsText;
	private TextField sdNmText;
	private TextField sdPixelsText;
	private TextField fwhmPixelsText;

	// Used for the PSF profile plot
	private double[] x = null, y, y2;

	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
		
		GlobalSettings globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getPsfCalculatorSettings();

		double sd = calculate(settings, false);
		if (sd < 0)
			return;

		globalSettings.getFitEngineConfiguration().getFitConfiguration().setInitialPeakStdDev((float) sd);
		globalSettings.getFitEngineConfiguration().getFitConfiguration().setInitialAngle(0);
		globalSettings.getCalibration().setNmPerPixel(getPixelPitch());
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
		gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		this.settings = settings;

		if (!simpleMode)
		{
			gd.addNumericField("Pixel_pitch (um)", settings.pixelPitch, 2);
			gd.addNumericField("Magnification", settings.magnification, 0);
			gd.addNumericField("Beam_Expander", settings.beamExpander, 2);
			gd.addMessage(getPixelPitchLabel());
			pixelPitchLabel = (Label) gd.getMessage();
			//pixelPitchLabel.setText(getPixelPitchLabel());
		}

		gd.addSlider("Wavelength (nm)", 400, 750, settings.wavelength);
		gd.addSlider("Numerical_Aperture (NA)", 1, 1.5, settings.numericalAperture);
		gd.addMessage(getAbbeLimitLabel());
		abbeLimitLabel = (Label) gd.getMessage();
		//abbe.setText(getAbbeLimitLabel());
		gd.addMessage("*** Account for optical aberations and focus error ***");
		gd.addSlider("Proportionality_factor", 1, 2.5, settings.proportionalityFactor);
		gd.addCheckbox("Adjust_for_square_pixels", settings.adjustForSquarePixels);

		if (!simpleMode)
		{
			gd.addNumericField("Airy Width (nm)", calculateAiryWidth(settings.wavelength, settings.numericalAperture),
					3);
			gd.addNumericField(
					"Airy Width (pixels)",
					calculateAiryWidth(settings.pixelPitch, settings.magnification * settings.beamExpander,
							settings.wavelength, settings.numericalAperture), 3);
		}
		gd.addNumericField("StdDev (nm)",
				calculateStdDev(settings.wavelength, settings.numericalAperture, settings.proportionalityFactor), 3);
		double sd = calculateStdDev(settings.pixelPitch, settings.magnification * settings.beamExpander,
				settings.wavelength, settings.numericalAperture, settings.proportionalityFactor,
				settings.adjustForSquarePixels);
		gd.addNumericField("StdDev (pixels)", sd, 3);
		gd.addNumericField("HWHM (pixels)", sd * GaussianFunction.SD_TO_HWHM_FACTOR, 3);

		if (!simpleMode)
		{
			widthNmText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 5);
			widthPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 4);
		}
		sdNmText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 3);
		sdPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 2);
		fwhmPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 1);

		//widthNmText
		if (!simpleMode)
		{
			disableEditing(widthNmText);
			disableEditing(widthPixelsText);
		}
		disableEditing(sdNmText);
		disableEditing(sdPixelsText);
		disableEditing(fwhmPixelsText);

		if (!simpleMode)
			gd.addMessage("Save StdDev pixel width to the fitting properties");

		gd.addDialogListener(this);

		double s = calculateStdDev(settings.pixelPitch, settings.magnification * settings.beamExpander, settings.wavelength,
				settings.numericalAperture, 1, false);
		plotProfile(
				calculateAiryWidth(settings.pixelPitch, settings.magnification * settings.beamExpander,
						settings.wavelength, settings.numericalAperture), sd / s);

		gd.showDialog();

		if (gd.wasCanceled())
		{
			return -1;
		}

		return calculateStdDev(settings.pixelPitch, settings.magnification * settings.beamExpander,
				settings.wavelength, settings.numericalAperture, settings.proportionalityFactor,
				settings.adjustForSquarePixels);
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
		settings.adjustForSquarePixels = gd.getNextBoolean();

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
			Parameters.isAboveZero("Proportionality factor", settings.proportionalityFactor);
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
	 * @return the SD in nm
	 */
	public static double calculateStdDev(double wavelength, double numericalAperture)
	{
		return AIRY_TO_GAUSSIAN * calculateAiryWidth(wavelength, numericalAperture);
	}

	/**
	 * Calculates the expected PSF standard deviation (nm) for a Gaussian approximation to the Airy disk.
	 * <p>
	 * <a href="http://en.wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile">http://en.
	 * wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile</a>
	 * <p>
	 * Using a scale factor of 1 results in the best match of the Gaussian to the Airy profile. A value above 1 should
	 * be used to increase the width to account for deviation of the optical system from the theoretical limit and
	 * out-of-focus objects.
	 * 
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture (NA)
	 * @param scaleFactor
	 *            Scale factor to account for deviation of the optical system and out-of-focus objects (should be >=1)
	 * @return the SD in nm
	 */
	public static double calculateStdDev(double wavelength, double numericalAperture, double scaleFactor)
	{
		return scaleFactor * AIRY_TO_GAUSSIAN * calculateAiryWidth(wavelength, numericalAperture);
	}

	/**
	 * If the pixel size (a) is provided the standard deviation (s) is adjusted to account for square pixels:
	 * 
	 * <pre>
	 * sa^2 = s^2 + a^2/12.
	 * </pre>
	 * 
	 * @param s
	 * @param a
	 * @return sa
	 */
	public static double squarePixelAdjustment(double s, final double a)
	{
		if (a > 0)
			s = Math.sqrt(s * s + a * a / 12.0);
		return s;
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
	 * @param scaleFactor
	 *            Scale factor to account for deviation of the optical system and out-of-focus objects (should be >=1)
	 * @param adjust
	 *            Adjust for square pixels
	 * @return the SD in pixels
	 */
	private static double calculateStdDev(double pixelPitch, double magnification, double wavelength,
			double numericalAperture, double scaleFactor, boolean adjust)
	{
		double s = calculateStdDev(wavelength, numericalAperture, scaleFactor);
		final double pixelPitchInNm = pixelPitch * 1000 / magnification;
		if (adjust)
			s = squarePixelAdjustment(s, pixelPitchInNm);
		s /= pixelPitchInNm;
		return s;
	}

	/**
	 * Calculates the PSF peak width for the Airy disk at the first dark ring (zero intersection).
	 * <p>
	 * Width is the lambda / (2*pi*NA).
	 * 
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture (NA)
	 * @return the width in nm
	 */
	public static double calculateAiryWidth(double wavelength, double numericalAperture)
	{
		double width = wavelength / (2.0 * Math.PI * numericalAperture);
		return width;
	}

	/**
	 * Calculates the PSF peak width for the Airy disk at the first dark ring (zero intersection).
	 * 
	 * @param pixelPitch
	 *            Camera pixel pitch in micrometers (um)
	 * @param magnification
	 *            Objective magnification
	 * @param wavelength
	 *            Wavelength of light in nanometers (nm)
	 * @param numericalAperture
	 *            Microscope numerical aperture (NA)
	 * @return the width in pixels
	 */
	public static double calculateAiryWidth(double pixelPitch, double magnification, double wavelength,
			double numericalAperture)
	{
		double width = calculateAiryWidth(wavelength, numericalAperture);
		final double pixelPitchInNm = pixelPitch * 1000 / magnification;
		return width / pixelPitchInNm;
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
		if (o == sdNmText || o == sdPixelsText || o == fwhmPixelsText || o == abbeLimitLabel)
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
					boolean adjustForSquarePixels = settings.adjustForSquarePixels;
					double proportionalityFactor = settings.proportionalityFactor;

					// Do something with parameters
					if (widthNmText != null)
					{
						pixelPitchLabel.setText(getPixelPitchLabel());
						widthNmText.setText(IJ.d2s(calculateAiryWidth(wavelength, numericalAperture), 3));
						widthPixelsText.setText(IJ.d2s(
								calculateAiryWidth(pixelPitch, magnification * beamExpander, wavelength,
										numericalAperture), 3));
					}
					abbeLimitLabel.setText(getAbbeLimitLabel());
					sdNmText.setText(IJ.d2s(calculateStdDev(wavelength, numericalAperture, proportionalityFactor), 3));
					double sd = calculateStdDev(pixelPitch, magnification * beamExpander, wavelength,
							numericalAperture, proportionalityFactor, adjustForSquarePixels);
					sdPixelsText.setText(IJ.d2s(sd, 3));
					fwhmPixelsText.setText(IJ.d2s(sd * GaussianFunction.SD_TO_HWHM_FACTOR, 3));

					double s = calculateStdDev(pixelPitch, magnification * beamExpander, wavelength,
							numericalAperture, 1, false);
					plotProfile(
							calculateAiryWidth(pixelPitch, magnification * beamExpander, wavelength, numericalAperture),
							sd / s);

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

	/**
	 * @param airyWidth The Airy width
	 * @param factor Factor used to scale the Airy approximation using the Gaussian 
	 */
	private void plotProfile(double airyWidth, double factor)
	{
		if (x == null)
		{
			x = Utils.newArray(200, -10, 0.1);
			y = new double[x.length];
			y2 = new double[x.length];
			for (int i = 0; i < x.length; i++)
			{
				y[i] = AiryPattern.intensity(x[i]);
			}
		}
		double[] x2 = new double[x.length];
		for (int i = 0; i < x2.length; i++)
		{
			x2[i] = x[i] * airyWidth;
			y2[i] = AiryPattern.intensityGaussian(x[i] / factor);
		}
		String title = "PSF profile";
		Plot2 p = new Plot2(title, "px", "", x2, y);
		p.addLabel(0, 0, "Blue = Airy; Red = Gaussian");
		p.setColor(Color.RED);
		p.addPoints(x2, y2, Plot2.LINE);
		final double sd = airyWidth * AIRY_TO_GAUSSIAN * factor;
		final double sdHeight = 0.606530659; //intensityGaussian(1);
		p.drawLine(-sd, 0, -sd, sdHeight);
		p.drawLine(sd, 0, sd, sdHeight);
		p.setColor(Color.BLUE);
		Utils.display(title, p);
	}

	/**
	 * Calculate the intensity of the Gaussian at distance x from the centre
	 * 
	 * @param x
	 * @return The intensity
	 */
	public static double intensityGaussian(double x)
	{
		if (x == 0)
			return 1;
		return FastMath.exp(-0.5 * (x * x));
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
	
	private String getAbbeLimitLabel(double abbeLimit)
	{
		return "Abbe limit (nm) = " + IJ.d2s(abbeLimit, 3);
	}

	private String getAbbeLimitLabel()
	{
		return getAbbeLimitLabel(getAbbeLimit());
	}

	private double getAbbeLimit()
	{
		return settings.wavelength / (2 * settings.numericalAperture);
	}
}
