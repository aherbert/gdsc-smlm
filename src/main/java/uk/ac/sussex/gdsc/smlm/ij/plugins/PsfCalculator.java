/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

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
import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.PSFCalculatorSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.model.AiryPattern;

/**
 * Calculates the expected PSF width for a Gaussian approximation to the Airy disk.
 */
public class PsfCalculator implements PlugIn, DialogListener {
  private static final String TITLE = "PSF Calculator";
  /**
   * This is the factor (f) for scaling an Airy width to a Gaussian approximation. The Gaussian PSF
   * Standard Deviation = f * lambda / (2 * pi * NA).
   */
  public static final double AIRY_TO_GAUSSIAN = 1.323;

  private PSFCalculatorSettings.Builder settingsBuilder;
  private GenericDialog gd;
  private Label abbeLimitLabel;
  private Label pixelPitchLabel;
  private TextField widthNmText;
  private TextField widthPixelsText;
  private TextField sdNmText;
  private TextField sdPixelsText;
  private TextField fwhmPixelsText;

  // Used for the PSF profile plot
  private double[] x;
  private double[] y;
  private double[] y2;

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    final PSFCalculatorSettings settings = SettingsManager.readPsfCalculatorSettings(0);

    final double sd = calculate(settings, false);
    if (sd < 0) {
      return;
    }

    SettingsManager.writeSettings(this.settingsBuilder);

    final FitEngineConfiguration config = SettingsManager.readFitEngineConfiguration(0);
    final FitConfiguration fitConfig = config.getFitConfiguration();
    fitConfig.setNmPerPixel(getPixelPitch());
    fitConfig.setPsfType(PSFType.ONE_AXIS_GAUSSIAN_2D);
    fitConfig.setInitialPeakStdDev(sd);
    SettingsManager.writeSettings(config, 0);
  }

  /**
   * Present an interactive dialog that allows the user to calculate the Gaussian PSF standard
   * deviation using the provided settings.
   *
   * @param settings the settings
   * @param simpleMode Only present a wavelength, NA and proportionality factor fields.
   * @return the PSF standard deviation
   */
  public double calculate(PSFCalculatorSettings settings, boolean simpleMode) {
    gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("psf-calculator"));

    this.settingsBuilder = settings.toBuilder();

    if (!simpleMode) {
      gd.addNumericField("Pixel_pitch (um)", settings.getPixelPitch(), 2);
      gd.addNumericField("Magnification", settings.getMagnification(), 0);
      gd.addNumericField("Beam_Expander", settings.getBeamExpander(), 2);
      gd.addMessage(getPixelPitchLabel());
      pixelPitchLabel = (Label) gd.getMessage();
    }

    gd.addSlider("Wavelength (nm)", 400, 750, settings.getWavelength());
    gd.addSlider("Numerical_Aperture (NA)", 1, 1.5, settings.getNumericalAperture());
    gd.addMessage(getAbbeLimitLabel());
    abbeLimitLabel = (Label) gd.getMessage();
    gd.addMessage("*** Account for optical aberations and focus error ***");
    gd.addSlider("Proportionality_factor", 1, 2.5, settings.getProportionalityFactor());
    gd.addCheckbox("Adjust_for_square_pixels", settings.getAdjustForSquarePixels());

    if (!simpleMode) {
      gd.addNumericField("Airy Width (nm)",
          calculateAiryWidth(settings.getWavelength(), settings.getNumericalAperture()), 3);
      gd.addNumericField("Airy Width (pixels)",
          calculateAiryWidth(settings.getPixelPitch(),
              settings.getMagnification() * settings.getBeamExpander(), settings.getWavelength(),
              settings.getNumericalAperture()),
          3);
    }
    gd.addNumericField("StdDev (nm)", calculateStdDev(settings.getWavelength(),
        settings.getNumericalAperture(), settings.getProportionalityFactor()), 3);
    final double sd = calculateStdDev(settings.getPixelPitch(),
        settings.getMagnification() * settings.getBeamExpander(), settings.getWavelength(),
        settings.getNumericalAperture(), settings.getProportionalityFactor(),
        settings.getAdjustForSquarePixels());
    gd.addNumericField("StdDev (pixels)", sd, 3);
    gd.addNumericField("HWHM (pixels)", sd * Gaussian2DFunction.SD_TO_HWHM_FACTOR, 3);

    if (!simpleMode) {
      widthNmText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 5);
      widthPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 4);
    }
    sdNmText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 3);
    sdPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 2);
    fwhmPixelsText = (TextField) gd.getNumericFields().get(gd.getNumericFields().size() - 1);

    // widthNmText
    if (!simpleMode) {
      disableEditing(widthNmText);
      disableEditing(widthPixelsText);
    }
    disableEditing(sdNmText);
    disableEditing(sdPixelsText);
    disableEditing(fwhmPixelsText);

    if (!simpleMode) {
      gd.addMessage("Save StdDev pixel width to the fitting properties");
    }

    gd.addDialogListener(this);

    final double s = calculateStdDev(settings.getPixelPitch(),
        settings.getMagnification() * settings.getBeamExpander(), settings.getWavelength(),
        settings.getNumericalAperture(), 1, false);
    plotProfile(calculateAiryWidth(settings.getPixelPitch(),
        settings.getMagnification() * settings.getBeamExpander(), settings.getWavelength(),
        settings.getNumericalAperture()), sd / s);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return -1;
    }

    return calculateStdDev(settings.getPixelPitch(),
        settings.getMagnification() * settings.getBeamExpander(), settings.getWavelength(),
        settings.getNumericalAperture(), settings.getProportionalityFactor(),
        settings.getAdjustForSquarePixels());
  }

  private static void disableEditing(TextField textField) {
    textField.setEditable(false);
    textField.setBackground(SystemColor.control);
  }

  private boolean readDialog() {
    if (widthNmText != null) {
      settingsBuilder.setPixelPitch(gd.getNextNumber());
      settingsBuilder.setMagnification(gd.getNextNumber());
      settingsBuilder.setBeamExpander(gd.getNextNumber());
    }
    settingsBuilder.setWavelength(gd.getNextNumber());
    settingsBuilder.setNumericalAperture(gd.getNextNumber());
    settingsBuilder.setProportionalityFactor(gd.getNextNumber());
    settingsBuilder.setAdjustForSquarePixels(gd.getNextBoolean());

    // Check arguments
    try {
      if (widthNmText != null) {
        ParameterUtils.isAboveZero("Pixel pitch", settingsBuilder.getPixelPitch());
        ParameterUtils.isAboveZero("Magnification", settingsBuilder.getMagnification());
        ParameterUtils.isEqualOrAbove("Beam expander", settingsBuilder.getBeamExpander(), 1);
      }
      ParameterUtils.isAboveZero("Wavelength", settingsBuilder.getWavelength());
      ParameterUtils.isAboveZero("Numerical aperture", settingsBuilder.getNumericalAperture());
      ParameterUtils.isAboveZero("Proportionality factor",
          settingsBuilder.getProportionalityFactor());
    } catch (final IllegalArgumentException ex) {
      // Q. Is logging the error necessary given that we will be in a live update preview?
      return false;
    }

    return !gd.invalidNumber();
  }

  /**
   * Calculates the expected PSF standard deviation (nm) for a Gaussian approximation to the Airy
   * disk.
   *
   * <p><a
   * href="http://en.wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile">http://en.
   * wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile</a>
   *
   * @param wavelength Wavelength of light in nanometers (nm)
   * @param numericalAperture Microscope numerical aperture (NA)
   * @return the SD in nm
   */
  public static double calculateStdDev(double wavelength, double numericalAperture) {
    return AIRY_TO_GAUSSIAN * calculateAiryWidth(wavelength, numericalAperture);
  }

  /**
   * Calculates the expected PSF standard deviation (nm) for a Gaussian approximation to the Airy
   * disk.
   *
   * <p><a
   * href="http://en.wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile">http://en.
   * wikipedia.org/wiki/Airy_disk#Approximation_using_a_Gaussian_profile</a>
   *
   * <p>Using a scale factor of 1 results in the best match of the Gaussian to the Airy profile. A
   * value above 1 should be used to increase the width to account for deviation of the optical
   * system from the theoretical limit and out-of-focus objects.
   *
   * @param wavelength Wavelength of light in nanometers (nm)
   * @param numericalAperture Microscope numerical aperture (NA)
   * @param scaleFactor Scale factor to account for deviation of the optical system and out-of-focus
   *        objects (should be {@code >=1}s)
   * @return the SD in nm
   */
  public static double calculateStdDev(double wavelength, double numericalAperture,
      double scaleFactor) {
    return scaleFactor * AIRY_TO_GAUSSIAN * calculateAiryWidth(wavelength, numericalAperture);
  }

  /**
   * Calculates the expected PSF standard deviation (pixels) for a Gaussian approximation to the
   * Airy disk.
   *
   * @param pixelPitch Camera pixel pitch in micrometers (um)
   * @param magnification Objective magnification
   * @param wavelength Wavelength of light in nanometers (nm)
   * @param numericalAperture Microscope numerical aperture (NA)
   * @param scaleFactor Scale factor to account for deviation of the optical system and out-of-focus
   *        objects (should be >=1)
   * @param adjust Adjust for square pixels
   * @return the SD in pixels
   */
  private static double calculateStdDev(double pixelPitch, double magnification, double wavelength,
      double numericalAperture, double scaleFactor, boolean adjust) {
    double sd = calculateStdDev(wavelength, numericalAperture, scaleFactor);
    final double pixelPitchInNm = pixelPitch * 1000 / magnification;
    if (adjust) {
      sd = squarePixelAdjustment(sd, pixelPitchInNm);
    }
    sd /= pixelPitchInNm;
    return sd;
  }

  /**
   * If the pixel size (a) is provided the standard deviation (s) is adjusted to account for square
   * pixels:
   *
   * <pre>
   * sa^2 = s^2 + a^2/12.
   * </pre>
   *
   * <p>This is relevant if using a single Gaussian evaluated at the centre of the pixel (0.5,0.5)
   * to represent the value over the entire pixel. If using a complete Gaussian function using the
   * integral of the error function (erf) then this is not needed.
   *
   * @param sd Gaussian standard deviation
   * @param pixelSize The pixel size
   * @return sa The adjusted standard deviation
   */
  public static double squarePixelAdjustment(double sd, final double pixelSize) {
    if (pixelSize > 0) {
      sd = Math.sqrt(sd * sd + pixelSize * pixelSize / 12.0);
    }
    return sd;
  }

  /**
   * Calculates the PSF peak width for the Airy disk at the first dark ring (zero intersection).
   *
   * <p>Width is the lambda / (2*pi*NA).
   *
   * @param wavelength Wavelength of light in nanometers (nm)
   * @param numericalAperture Microscope numerical aperture (NA)
   * @return the width in nm
   */
  public static double calculateAiryWidth(double wavelength, double numericalAperture) {
    return wavelength / (2.0 * Math.PI * numericalAperture);
  }

  /**
   * Calculates the PSF peak width for the Airy disk at the first dark ring (zero intersection).
   *
   * @param pixelPitch Camera pixel pitch in micrometers (um)
   * @param magnification Objective magnification
   * @param wavelength Wavelength of light in nanometers (nm)
   * @param numericalAperture Microscope numerical aperture (NA)
   * @return the width in pixels
   */
  public static double calculateAiryWidth(double pixelPitch, double magnification,
      double wavelength, double numericalAperture) {
    final double width = calculateAiryWidth(wavelength, numericalAperture);
    final double pixelPitchInNm = pixelPitch * 1000 / magnification;
    return width / pixelPitchInNm;
  }

  private boolean lock;

  private synchronized boolean aquireLock() {
    if (lock) {
      return false;
    }
    return lock = true;
  }

  @Override
  public boolean dialogItemChanged(GenericDialog gd, AWTEvent event) {
    if (event == null) {
      return true;
    }
    final Object o = event.getSource();
    if (o == sdNmText || o == sdPixelsText || o == fwhmPixelsText || o == abbeLimitLabel) {
      return true;
    }
    if (widthNmText != null && (o == pixelPitchLabel || o == widthNmText || o == widthPixelsText)) {
      return true;
    }

    if (aquireLock()) {
      try {
        // Continue while the parameter is changing
        boolean parametersChanged = true;
        while (parametersChanged) {
          if (!readDialog()) {
            return false;
          }

          // Store the parameters to be processed
          final double pixelPitch = settingsBuilder.getPixelPitch();
          final double magnification = settingsBuilder.getMagnification();
          final double beamExpander = settingsBuilder.getBeamExpander();
          final double wavelength = settingsBuilder.getWavelength();
          final double numericalAperture = settingsBuilder.getNumericalAperture();
          final boolean adjustForSquarePixels = settingsBuilder.getAdjustForSquarePixels();
          final double proportionalityFactor = settingsBuilder.getProportionalityFactor();

          // Do something with parameters
          if (widthNmText != null) {
            pixelPitchLabel.setText(getPixelPitchLabel());
            widthNmText.setText(IJ.d2s(calculateAiryWidth(wavelength, numericalAperture), 3));
            widthPixelsText.setText(IJ.d2s(calculateAiryWidth(pixelPitch,
                magnification * beamExpander, wavelength, numericalAperture), 3));
          }
          abbeLimitLabel.setText(getAbbeLimitLabel());
          sdNmText.setText(
              IJ.d2s(calculateStdDev(wavelength, numericalAperture, proportionalityFactor), 3));
          final double sd = calculateStdDev(pixelPitch, magnification * beamExpander, wavelength,
              numericalAperture, proportionalityFactor, adjustForSquarePixels);
          sdPixelsText.setText(IJ.d2s(sd, 3));
          fwhmPixelsText.setText(IJ.d2s(sd * Gaussian2DFunction.SD_TO_HWHM_FACTOR, 3));

          final double s = calculateStdDev(pixelPitch, magnification * beamExpander, wavelength,
              numericalAperture, 1, false);
          plotProfile(calculateAiryWidth(pixelPitch, magnification * beamExpander, wavelength,
              numericalAperture), sd / s);

          // Check if the parameters have changed again
          parametersChanged = (pixelPitch != settingsBuilder.getPixelPitch())
              || (magnification != settingsBuilder.getMagnification())
              || (beamExpander != settingsBuilder.getBeamExpander())
              || (wavelength != settingsBuilder.getWavelength())
              || (numericalAperture != settingsBuilder.getNumericalAperture())
              || (proportionalityFactor != settingsBuilder.getProportionalityFactor());
        }
      } finally {
        // Ensure the running flag is reset
        lock = false;
      }
    }

    return true;
  }

  /**
   * Plot the profile.
   *
   * @param airyWidth The Airy width
   * @param factor Factor used to scale the Airy approximation using the Gaussian
   */
  private void plotProfile(double airyWidth, double factor) {
    if (x == null) {
      x = SimpleArrayUtils.newArray(200, -10, 0.1);
      y = new double[x.length];
      y2 = new double[x.length];
      for (int i = 0; i < x.length; i++) {
        y[i] = AiryPattern.intensity(x[i]);
      }
    }
    final double[] x2 = new double[x.length];
    for (int i = 0; i < x2.length; i++) {
      x2[i] = x[i] * airyWidth;
      y2[i] = AiryPattern.intensityGaussian(x[i] / factor);
    }
    final String title = "PSF profile";
    final Plot p = new Plot(title, "px", "");
    p.addLabel(0, 0, "Blue = Airy; Red = Gaussian");
    p.setColor(Color.BLUE);
    p.addPoints(x2, y, Plot.LINE);
    p.setColor(Color.RED);
    p.addPoints(x2, y2, Plot.LINE);
    final double sd = airyWidth * AIRY_TO_GAUSSIAN * factor;
    final double sdHeight = 0.606530659; // intensityGaussian(1)
    p.drawLine(-sd, 0, -sd, sdHeight);
    p.drawLine(sd, 0, sd, sdHeight);
    ImageJUtils.display(title, p);
  }

  /**
   * Calculate the intensity of the Gaussian at distance x from the centre.
   *
   * @param x the x
   * @return The intensity
   */
  public static double intensityGaussian(double x) {
    if (x == 0) {
      return 1;
    }
    return FastMath.exp(-0.5 * (x * x));
  }

  private static String getPixelPitchLabel(double pixelPitch) {
    return "Pixel pitch (nm) = " + IJ.d2s(pixelPitch, 3);
  }

  private String getPixelPitchLabel() {
    return getPixelPitchLabel(getPixelPitch());
  }

  private double getPixelPitch() {
    return settingsBuilder.getPixelPitch() * 1000
        / (settingsBuilder.getMagnification() * settingsBuilder.getBeamExpander());
  }

  private static String getAbbeLimitLabel(double abbeLimit) {
    return "Abbe limit (nm) = " + IJ.d2s(abbeLimit, 3);
  }

  private String getAbbeLimitLabel() {
    return getAbbeLimitLabel(getAbbeLimit());
  }

  private double getAbbeLimit() {
    return settingsBuilder.getWavelength() / (2 * settingsBuilder.getNumericalAperture());
  }
}
