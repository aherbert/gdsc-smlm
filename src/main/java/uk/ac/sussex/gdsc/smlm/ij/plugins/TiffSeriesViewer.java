/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SeriesOpener;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.io.ExtendedFileInfo;
import uk.ac.sussex.gdsc.core.logging.TrackProgressAdaptor;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.ij.SeriesImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.Constants;
import uk.ac.sussex.gdsc.smlm.results.ImageSource.ReadHint;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.Prefs;
import ij.VirtualStack;
import ij.io.FileInfo;
import ij.io.TiffEncoder;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;

import java.awt.Choice;
import java.awt.Font;
import java.awt.Label;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Reads a TIFF image using the series image source and presents it using a read-only virtual stack
 * image.
 */
public class TiffSeriesViewer implements PlugIn {
  private static final String TITLE = "Tiff Series Viewer";

  private Label label;
  private Label label2;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    private static final String[] MODE = {"Directory", "File"};
    private static final String[] OUTPUT_MODE = {"Image", "Files"};

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    private int inputMode;
    private String inputDirectory;
    private String inputFile;
    private boolean logProgress;
    private int outputMode;
    private int imageCount;
    private String outputDirectory;

    Settings() {
      // Set defaults
      inputMode = (int) Prefs.get(Constants.tiffSeriesMode, 0);
      inputDirectory = Prefs.get(Constants.tiffSeriesDirectory, "");
      inputFile = Prefs.get(Constants.tiffSeriesFile, "");
      logProgress = Prefs.getBoolean(Constants.tiffSeriesLogProgress, false);
      outputMode = (int) Prefs.get(Constants.tiffSeriesOutputMode, 0);
      imageCount = (int) Prefs.get(Constants.tiffSeriesOutputNImages, 1);
      outputDirectory = Prefs.get(Constants.tiffSeriesOutputDirectory, "");
    }

    Settings(Settings source) {
      inputMode = source.inputMode;
      inputDirectory = source.inputDirectory;
      inputFile = source.inputFile;
      logProgress = source.logProgress;
      outputMode = source.outputMode;
      imageCount = source.imageCount;
      outputDirectory = source.outputDirectory;
    }

    Settings copy() {
      return new Settings(this);
    }

    /**
     * Load a copy of the settings.
     *
     * @return the settings
     */
    static Settings load() {
      return lastSettings.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      lastSettings.set(this);
      Prefs.set(Constants.tiffSeriesMode, inputMode);
      Prefs.set(Constants.tiffSeriesDirectory, inputDirectory);
      Prefs.set(Constants.tiffSeriesFile, inputFile);
      Prefs.set(Constants.tiffSeriesLogProgress, logProgress);
      Prefs.set(Constants.tiffSeriesOutputMode, outputMode);
      Prefs.set(Constants.tiffSeriesOutputNImages, imageCount);
      Prefs.set(Constants.tiffSeriesOutputDirectory, outputDirectory);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    settings = Settings.load();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addChoice("Mode", Settings.MODE, settings.inputMode, new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        settings.inputMode = value;
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        // This has limited silent support to fake running in a macro
        if (settings.inputMode == 0) {
          String dir = null;
          final String title = "Select image series ...";
          if (silent) {
            final String macroOptions = Macro.getOptions();
            if (macroOptions != null) {
              dir = Macro.getValue(macroOptions, title, null);
            }
          } else {
            dir = ImageJUtils.getDirectory(title, settings.inputDirectory);
          }
          if (TextUtils.isNullOrEmpty(dir)) {
            return false;
          }
          settings.inputDirectory = dir;
        } else {
          String file = null;
          final String title = "Select image ...";
          if (silent) {
            final String macroOptions = Macro.getOptions();
            if (macroOptions != null) {
              file = Macro.getValue(macroOptions, title, null);
            }
          } else {
            file = ImageJUtils.getFilename(title, settings.inputFile);
          }
          if (TextUtils.isNullOrEmpty(file)) {
            return false;
          }
          settings.inputFile = file;
        }
        updateLabel();
        return true;
      }
    });
    gd.addMessage("");
    label = gd.getLastLabel();
    if (ImageJUtils.isShowGenericDialog()) {
      final Choice choice = gd.getLastChoice();
      choice.addItemListener(event -> {
        settings.inputMode = choice.getSelectedIndex();
        updateLabel();
      });
      updateLabel();
    }
    gd.addCheckbox("Log_progress", settings.logProgress);
    gd.addChoice("Output_mode", Settings.OUTPUT_MODE, settings.outputMode,
        new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer value) {
            settings.outputMode = value;
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            if (settings.outputMode == 0) {
              // Nothing to do
              return false;
            }
            final ExtendedGenericDialog egd = new ExtendedGenericDialog("Output Options");
            egd.addNumericField("Slices_per_image", settings.imageCount, 0);
            egd.addDirectoryField("Output_directory", settings.outputDirectory);
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            settings.imageCount = (int) egd.getNextNumber();
            settings.outputDirectory = egd.getNextString();
            updateLabel2();
            return true;
          }
        });
    gd.addMessage("");
    label2 = gd.getLastLabel();
    if (ImageJUtils.isShowGenericDialog()) {
      final Choice choice = gd.getLastChoice();
      choice.addItemListener(event -> {
        settings.outputMode = choice.getSelectedIndex();
        updateLabel2();
      });
      updateLabel2();
    }

    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    settings.inputMode = gd.getNextChoiceIndex();
    settings.logProgress = gd.getNextBoolean();
    settings.outputMode = gd.getNextChoiceIndex();

    settings.save();

    SeriesImageSource source;
    if (settings.inputMode == 0) {
      final SeriesOpener series = new SeriesOpener(settings.inputDirectory);
      if (series.getNumberOfImages() == 0) {
        IJ.error(TITLE, "No images in the selected directory:\n" + settings.inputDirectory);
        return;
      }
      source = new SeriesImageSource(PeakFit.getName(series.getImageList()), series);
    } else {
      source = new SeriesImageSource(FileUtils.getName(settings.inputFile),
          new String[] {settings.inputFile});
    }

    source.setBufferLimit(0); // No memory buffer
    source.setReadHint(ReadHint.NONSEQUENTIAL);

    if (!source.isTiffSeries) {
      IJ.error(TITLE, "Not a TIFF image");
      return;
    }
    ImageJUtils.showStatus("Opening TIFF ...");
    final TrackProgressAdaptor progress = new TrackProgressAdaptor() {
      @Override
      public void progress(double fraction) {
        IJ.showProgress(fraction);
      }

      @Override
      public void progress(long position, long total) {
        IJ.showProgress((double) position / total);
      }

      @Override
      public void log(String format, Object... args) {
        if (settings.logProgress) {
          ImageJUtils.log(format, args);
        }
      }

      @Override
      public void status(String format, Object... args) {
        ImageJUtils.showStatus(() -> String.format(format, args));
      }

      @Override
      public boolean isLog() {
        return settings.logProgress;
      }
    };
    source.setTrackProgress(progress);
    if (!source.open()) {
      IJ.error(TITLE, "Cannot open the image");
      return;
    }
    ImageJUtils.showStatus("");

    // Create a virtual stack
    final TiffSeriesVirtualStack stack = new TiffSeriesVirtualStack(source);
    if (settings.outputMode == 0) {
      stack.show();
    } else {
      final int nImages = Math.max(1, settings.imageCount);
      final ImagePlus imp = stack.createImp();
      // The calibration only has the offset so ignore for speed.
      // Calibration cal = imp.getCalibration();
      final int size = stack.getSize();

      // Create the format string
      final int digits = String.format("%d", size).length();
      final String format =
          new File(settings.outputDirectory, imp.getShortTitle() + "%0" + digits + "d.tif")
              .getPath();

      IJ.showStatus("Saving image ...");
      try {
        for (int i = 1; i <= size; i += nImages) {
          if (ImageJUtils.isInterrupted()) {
            break;
          }
          ImageJUtils.showSlowProgress(i, size);
          final String path = String.format(format, i);
          final ImageStack out = new ImageStack(source.getWidth(), source.getHeight());
          for (int j = 0, k = i; j < nImages && k <= size; j++, k++) {
            out.addSlice(null, stack.getPixels(k));
          }
          final ImagePlus outImp = new ImagePlus(path, out);
          // outImp.setCalibration(cal);
          saveAsTiff(outImp, path);
        }
        IJ.showStatus("Saved image");
      } catch (final IOException ex) {
        IJ.log(ExceptionUtils.getStackTrace(ex));
        IJ.error(TITLE, "Failed to save image: " + ex.getMessage());
        IJ.showStatus("Failed to save image");
      } finally {
        ImageJUtils.clearSlowProgress();
      }
    }
  }

  private static void saveAsTiff(ImagePlus imp, String path) throws IOException {
    final FileInfo fi = imp.getFileInfo();
    fi.nImages = imp.getStackSize();
    try (OutputStream out = new BufferedOutputStream(new FileOutputStream(path))) {
      new TiffEncoder(fi).write(out);
    }
  }

  private void updateLabel() {
    if (settings.inputMode == 0) {
      label.setText(settings.inputDirectory);
    } else {
      label.setText(settings.inputFile);
    }
  }

  private void updateLabel2() {
    if (settings.outputMode == 0) {
      label2.setText("");
    } else {
      label2.setText(String.format("Slices per image = %d : %s", settings.imageCount,
          settings.outputDirectory));
    }
  }

  /**
   * Override methods in the ij.VirtualStack class to provide the pixels from a TIFF series. The
   * stack cannot be modified.
   */
  public static class TiffSeriesVirtualStack extends VirtualStack {
    /** The image source. */
    SeriesImageSource source;

    /**
     * Instantiates a new tiff series virtual stack with a source. The source must have been
     * successfully opened.
     *
     * @param source the source
     */
    public TiffSeriesVirtualStack(SeriesImageSource source) {
      super(source.getWidth(), source.getHeight(), null, null);
      if (!source.isValid(1)) {
        throw new IllegalArgumentException("Source has no frames");
      }
      this.source = source;
      final Object pixels = source.getRaw(1);
      if (pixels == null) {
        throw new IllegalArgumentException("Source has no first frame");
      }
      setBitDepth(ImageJUtils.getBitDepth(pixels));
    }

    /**
     * Wrap the series in a ImagePlus object.
     *
     * @return the image plus
     */
    public ImagePlus createImp() {
      final ImagePlus imp = new ImagePlus(source.getName(), this);
      addInfo(imp);
      return imp;
    }

    /**
     * Show the series in a ImagePlus object.
     *
     * @return the image plus
     */
    public ImagePlus show() {
      final ImagePlus imp = ImageJUtils.display(source.getName(), this);
      addInfo(imp);
      return imp;
    }

    private void addInfo(ImagePlus imp) {
      // So the FileSaver can save the stack make sure the FileInfo is not null
      imp.getFileInfo();

      // Get metadata from the source
      final ExtendedFileInfo[] fileInfo = source.getFileInfo(0);
      if (ArrayUtils.getLength(fileInfo) > 0) {
        final ExtendedFileInfo efi = fileInfo[0];
        if (efi.getExtendedMetaData() != null) {
          imp.setProperty("Info", efi.getExtendedMetaData());
        } else if (efi.info != null) {
          imp.setProperty("Info", efi.info);
        }
      }
      if (source.getXOrigin() != 0 || source.getYOrigin() != 0) {
        final Calibration cal = imp.getLocalCalibration();
        cal.xOrigin = -source.getXOrigin();
        cal.yOrigin = -source.getYOrigin();
      }
    }

    /**
     * Does nothing.
     */
    @Override
    public void addSlice(String name) {
      // Do nothing
    }

    /**
     * Does nothing.
     */
    @Override
    public void deleteSlice(int n) {
      // Do nothing
    }

    @Override
    public ImageProcessor getProcessor(int n) {
      final Object pixels = source.getRaw(n);
      ImageProcessor ip = null;
      int depthThisImage = 0;
      if (pixels != null) {
        ip = ImageJUtils.createProcessor(getWidth(), getHeight(), pixels);
      } else {
        ip = new ByteProcessor(getWidth(), getHeight());
        ip.invert();
        int size = getHeight() / 20;
        if (size < 9) {
          size = 9;
        }
        final Font font = new Font("Helvetica", Font.PLAIN, size);
        ip.setFont(font);
        ip.setAntialiasedText(true);
        ip.setColor(0);
        ip.drawString("Error opening frame " + n, size, size * 2);
        depthThisImage = 8;
      }
      // Convert to the correct bit depth
      if (depthThisImage != getBitDepth()) {
        switch (getBitDepth()) {
          case 8:
            ip = ip.convertToByte(true);
            break;
          case 16:
            ip = ip.convertToShort(true);
            break;
          case 24:
            ip = ip.convertToRGB();
            break;
          case 32:
            ip = ip.convertToFloat();
            break;
          default:
            // Unknown bit depth. Leave as is.
            break;
        }
      }
      // Note: width/height checks are not required as the source checks the dimensions
      return ip;
    }

    @Override
    public int saveChanges(int n) {
      return -1; // Not implemented
    }

    @Override
    public int getSize() {
      return source.getFrames();
    }

    @Override
    public String getSliceLabel(int n) {
      return null;
    }

    @Override
    public String getDirectory() {
      return null;
    }

    @Override
    public String getFileName(int n) {
      return null;
    }

    @Override
    public ImageStack sortDicom(String[] strings, String[] info, int maxDigits) {
      // Don't sort
      return this;
    }
  }
}
