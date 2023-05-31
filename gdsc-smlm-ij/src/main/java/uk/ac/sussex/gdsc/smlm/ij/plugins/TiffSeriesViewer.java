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
import java.awt.Choice;
import java.awt.Font;
import java.awt.Label;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SeriesOpener;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.io.ExtendedFileInfo;
import uk.ac.sussex.gdsc.core.logging.TrackProgressAdapter;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.ij.SeriesImageSource;
import uk.ac.sussex.gdsc.smlm.results.ImageSource.ReadHint;

/**
 * Reads a TIFF image using the series image source and presents it using a read-only virtual stack
 * image.
 */
public class TiffSeriesViewer implements PlugIn {
  private static final String TITLE = "Tiff Series Viewer";

  /** The 1st label in the dialog. */
  private Label label;
  /** The 2nd label in the dialog. */
  private Label label2;

  /** The plugin settings. */
  Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String[] MODE = {"Directory", "File"};
    static final String[] OUTPUT_MODE = {"Image", "Files", "Validate"};
    static final int INPUT_DIRECTORY = 0;
    static final int OUTPUT_FILES = 1;
    static final int OUTPUT_VALIDATE = 2;

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    int inputMode;
    String inputDirectory;
    String inputFile;
    boolean logProgress;
    int outputMode;
    int imageCount;
    String outputDirectory;

    Settings() {
      // Set defaults
      inputMode = (int) Prefs.get(PrefsKey.tiffSeriesMode, 0);
      inputDirectory = Prefs.get(PrefsKey.tiffSeriesDirectory, "");
      inputFile = Prefs.get(PrefsKey.tiffSeriesFile, "");
      logProgress = Prefs.getBoolean(PrefsKey.tiffSeriesLogProgress, false);
      outputMode = (int) Prefs.get(PrefsKey.tiffSeriesOutputMode, 0);
      imageCount = (int) Prefs.get(PrefsKey.tiffSeriesOutputNImages, 1);
      outputDirectory = Prefs.get(PrefsKey.tiffSeriesOutputDirectory, "");
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
      return INSTANCE.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      INSTANCE.set(this);
      Prefs.set(PrefsKey.tiffSeriesMode, inputMode);
      Prefs.set(PrefsKey.tiffSeriesDirectory, inputDirectory);
      Prefs.set(PrefsKey.tiffSeriesFile, inputFile);
      Prefs.set(PrefsKey.tiffSeriesLogProgress, logProgress);
      Prefs.set(PrefsKey.tiffSeriesOutputMode, outputMode);
      Prefs.set(PrefsKey.tiffSeriesOutputNImages, imageCount);
      Prefs.set(PrefsKey.tiffSeriesOutputDirectory, outputDirectory);
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
            if (settings.outputMode != Settings.OUTPUT_FILES) {
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

    gd.addHelp(HelpUrls.getUrl("tiff-series-viewer"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }

    settings.inputMode = gd.getNextChoiceIndex();
    settings.logProgress = gd.getNextBoolean();
    settings.outputMode = gd.getNextChoiceIndex();

    settings.save();

    SeriesImageSource source;
    if (settings.inputMode == Settings.INPUT_DIRECTORY) {
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
    final TrackProgressAdapter progress = new TrackProgressAdapter() {
      private final boolean logProgress = settings.logProgress;

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
        if (logProgress) {
          ImageJUtils.log(format, args);
        }
      }

      @Override
      public void status(String format, Object... args) {
        ImageJUtils.showStatus(() -> String.format(format, args));
      }

      @Override
      public boolean isLog() {
        return logProgress;
      }
    };
    source.setTrackProgress(progress);
    if (settings.outputMode == Settings.OUTPUT_VALIDATE) {
      if (validate(source, progress)) {
        ImageJUtils.finished(source.getName() + " valid");
      } else {
        IJ.error(TITLE, "Not a valid image series");
      }
      return;
    }
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

      try {
        for (int i = 1; i <= size; i += nImages) {
          if (ImageJUtils.isInterrupted()) {
            break;
          }
          ImageJUtils.showSlowProgress(i, size);
          final int from = i;
          final int to = Math.min(size, i + nImages - 1);
          ImageJUtils.showStatus(() -> "Creating sub-image ... " + from + " - " + to);
          final String path = String.format(format, i);
          final ImageStack out = new ImageStack(source.getWidth(), source.getHeight());
          for (int k = from; k <= to; k++) {
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
    try (OutputStream out = new BufferedOutputStream(Files.newOutputStream(Paths.get(path)))) {
      new TiffEncoder(fi).write(out);
    }
  }

  /**
   * Update the 1st label in the dialog.
   */
  void updateLabel() {
    if (settings.inputMode == Settings.INPUT_DIRECTORY) {
      label.setText(settings.inputDirectory);
    } else {
      label.setText(settings.inputFile);
    }
  }

  /**
   * Update the 2nd label in the dialog.
   */
  void updateLabel2() {
    if (settings.outputMode == Settings.OUTPUT_FILES) {
      label2.setText(String.format("Slices per image = %d : %s", settings.imageCount,
          settings.outputDirectory));
    } else {
      label2.setText("");
    }
  }

  /**
   * Validate the series. This method will open the series; get the last raw frame from each image;
   * and close the series.
   *
   * @param source the source
   * @param progress the progress
   * @return true if valid
   */
  private static boolean validate(SeriesImageSource source, TrackProgressAdapter progress) {
    if (!source.open()) {
      progress.log("Failed to open image series: " + source.getName());
      return false;
    }
    // All images have non-zero size
    int frame = 0;
    final int n = source.getSeriesSize();
    try {
      for (int i = 0; i < n; i++) {
        ImageJUtils.showSlowProgress(i, n);
        final int size = source.getImageSize(i);
        if (size <= 0) {
          final FileInfo[] fi = source.getFileInfo(i);
          if (fi != null) {
            progress.log("Unknown size for image: [%d] %s", i, fi[0].fileName);
          } else {
            progress.log("Unknown size for image: [%d]", i);
          }
          return false;
        }
        frame += size;
        if (source.getRaw(frame) == null) {
          return false;
        }
      }
    } finally {
      ImageJUtils.clearSlowProgress();
      source.close();
    }
    return true;
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
      ImageProcessor ip;
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
