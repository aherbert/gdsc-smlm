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
import ij.WindowManager;
import ij.gui.ImageWindow;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.text.TextPanel;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Label;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageAdapter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.utils.MemoryUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrentMonoStack;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJTablePeakResults;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultView;

/**
 * Produces a summary table of the results that are stored in memory.
 */
public class OverlayResults implements PlugIn {
  private static final String TITLE = "Overlay Results";

  private String[] names;
  private int[] ids;
  private Choice choice;
  private Checkbox checkbox;
  private Label label;

  private final ConcurrentMonoStack<Integer> inbox = new ConcurrentMonoStack<>();

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String name;
    boolean showTable;

    Settings() {
      // Set defaults
      name = "";
    }

    Settings(Settings source) {
      name = source.name;
      showTable = source.showTable;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  private class Worker extends ImageAdapter implements ItemListener, Runnable {
    private boolean running = true;
    private final boolean[] error = new boolean[ids.length];
    // The results text window (so we can close it)
    private TextWindow tw;
    private Rectangle windowBounds;

    PeakResultView view;
    TypeConverter<DistanceUnit> converter;

    private int currentIndex;
    private int currentSlice = -1;

    @Override
    public void itemStateChanged(ItemEvent event) {
      // Read other options from the dialog
      settings.showTable = checkbox.getState();
      refresh();
    }

    @Override
    public void imageUpdated(ImagePlus imp) {
      if (imp == null) {
        return;
      }
      // If this is a change to the current image then refresh
      if (ids[currentIndex] == imp.getID()) {
        refresh();
      }
    }

    /**
     * Refresh the overlay with the currently selected image from the list.
     */
    void refresh() {
      inbox.insert(choice.getSelectedIndex());
    }

    @Override
    public void run() {
      while (running) {
        try {
          final Integer job = inbox.pop();
          if (job == null || !running) {
            break;
          }
          final int index = job.intValue();
          if (index == 0) {
            // This is selection of no image
            clearOldOverlay();
            continue;
          }

          drawOverlay(index);
        } catch (final InterruptedException ex) {
          running = false;
          Logger.getLogger(OverlayResults.class.getName()).log(Level.WARNING,
              "Unexpected intteruption", ex);
          Thread.currentThread().interrupt();
        }
      }
      clearOldOverlay();
      closeTextWindow();
    }

    private void clearOldOverlay() {
      if (currentIndex != 0) {
        final ImagePlus oldImp = WindowManager.getImage(ids[currentIndex]);
        if (oldImp != null) {
          oldImp.setOverlay(null);
        }
      }
      view = null;
      currentSlice = -1;
      currentIndex = 0;
    }

    private void closeTextWindow() {
      if (tw != null) {
        windowBounds = tw.getBounds();
        tw.close();
        tw = null;
      }
    }

    /**
     * Draw the overlay.
     *
     * <p>This is only called when index > 0.
     *
     * @param newIndex the new index
     */
    private void drawOverlay(int newIndex) {
      // Get new and old image Id
      final int oldIndex = currentIndex;
      final int oldId = ids[oldIndex];
      final int newId = ids[newIndex];

      // Check identity of the image
      final boolean newImage = oldId != newId;
      if (newImage) {
        clearOldOverlay();
      }

      currentIndex = newIndex;

      final ImagePlus imp = WindowManager.getImage(newId);
      final String name = names[newIndex];

      if (imp == null) {
        // Image has been closed.
        logError("Image not available", name);
        return;
      }

      final int newSlice = imp.getCurrentSlice();
      // If same results then check if the slice has changed
      if (newIndex == oldIndex) {
        if (currentSlice == newSlice) {
          final boolean isShowing = tw != null;
          if (settings.showTable == isShowing) {
            // No change from last time
            return;
          }
        }
      } else {
        // Get a new snapshot view
        view = null;
      }
      currentSlice = newSlice;

      final MemoryPeakResults results = MemoryPeakResults.getResults(name);
      if (results == null) {
        // Results have been cleared from memory (or renamed).
        logError("Results not available", name);
        return;
      }
      clearError();

      final ImageJTablePeakResults table;
      IntOpenHashSet selectedId = null;
      if (settings.showTable) {
        final boolean hasId = results.hasId();

        // Old selected item
        boolean is3D = false;
        if (hasId && tw != null) {
          final TextPanel tp = tw.getTextPanel();
          final int idColumn = ImageJUtils.getColumn(tp, "Id");
          final int start = tp.getSelectionStart();
          if (start != -1 && idColumn != -1) {
            selectedId = new IntOpenHashSet();
            final int end = tp.getSelectionEnd();
            for (int index = start; index <= end; index++) {
              final String text = tp.getLine(index).split("\t")[idColumn];
              selectedId.add(Integer.parseInt(text));
            }
          }
          // Keep the z column to avoid table redraw
          is3D = tp.getColumnHeadings().contains("\tZ");
        }

        // New table
        is3D = is3D || results.is3D();
        table = new ImageJTablePeakResults(false);
        table.setTableTitle(TITLE);
        table.copySettings(results);
        table.setClearAtStart(true);
        table.setAddCounter(true);
        table.setHideSourceText(true);
        table.setShowZ(is3D);
        table.setShowId(hasId);
        table.begin();

        tw = table.getResultsWindow();
        if (windowBounds != null) {
          tw.setBounds(windowBounds);
        } else {
          // Position under the window
          final ImageWindow win = imp.getWindow();
          final Point p = win.getLocation();
          p.y += win.getHeight();
          tw.setLocation(p);
        }
      } else {
        table = null;
        closeTextWindow();
      }

      if (view == null) {
        view = results.getSnapshotPeakResultView();
        converter = results.getDistanceConverter(DistanceUnit.PIXEL);
      }
      int select = -1;
      final PeakResult[] frameResults = view.getResultsByFrame(currentSlice);
      int size = 0;
      float[] ox = new float[11];
      float[] oy = new float[11];
      for (int i = 0; i < frameResults.length; i++) {
        final PeakResult r = frameResults[i];
        if (ox.length == size) {
          ox = Arrays.copyOf(ox, MemoryUtils.createNewCapacity(size + 1, size));
          oy = Arrays.copyOf(oy, ox.length);
        }
        ox[size] = converter.convert(r.getXPosition());
        oy[size] = converter.convert(r.getYPosition());
        size++;
        if (table != null) {
          table.add(r);
          if (selectedId != null && selectedId.contains(r.getId())) {
            // For now just preserve the first selected ID.
            // This at least allows tracking a single localisation.
            select = i;
            selectedId = null;
          }
        }
      }

      // Coords are copied in PolygonRoi
      final PointRoi roi = new OffsetPointRoi(ox, oy, size);
      roi.setPointType(3);
      if (newImage) {
        // New windows to the front
        imp.getWindow().toFront();
      }
      imp.setOverlay(new Overlay(roi));

      if (table != null) {
        table.end();
        final TextPanel tp = table.getResultsWindow().getTextPanel();
        tp.scrollToTop();

        // Reselect the same Id
        if (select != -1) {
          table.select(select);
        }
      }
    }

    private void logError(String msg, String name) {
      if (!error[currentIndex]) {
        ImageJUtils.log("%s Error: %s for results '%s'", TITLE, msg, name);
        label.setText("Error: " + msg + ". Restart this plugin to refresh.");
      }
      error[currentIndex] = true;
    }

    private void clearError() {
      error[currentIndex] = false;
      if (!TextUtils.isNullOrEmpty(label.getText())) {
        label.setText("");
      }
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    names = new String[MemoryPeakResults.getResultNames().size() + 1];
    ids = new int[names.length];
    int count = 0;
    names[count++] = "(None)";
    for (final MemoryPeakResults results : MemoryPeakResults.getAllResults()) {
      if (results.getSource() != null
          && results.getSource().getOriginal() instanceof IJImageSource) {
        final IJImageSource source = (IJImageSource) (results.getSource().getOriginal());
        final ImagePlus imp = WindowManager.getImage(source.getName());
        if (imp != null) {
          ids[count] = imp.getID();
          names[count++] = results.getName();
        }
      }
    }
    if (count == 1) {
      IJ.error(TITLE, "There are no result images available");
      return;
    }
    names = Arrays.copyOf(names, count);

    Thread thread = null;
    Worker worker = null;
    final NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
    settings = Settings.load();
    settings.save();
    gd.addMessage("Overlay results on current image frame");
    gd.addChoice("Results", names, (settings.name == null) ? "" : settings.name);
    gd.addCheckbox("Show_table", settings.showTable);
    gd.addMessage("");
    gd.addHelp(HelpUrls.getUrl("overlay-results"));
    gd.hideCancelButton();
    gd.setOKLabel("Close");
    if (!(IJ.isMacro() || java.awt.GraphicsEnvironment.isHeadless())) {
      worker = new Worker();

      choice = (Choice) gd.getChoices().get(0);
      choice.addItemListener(worker);
      checkbox = (Checkbox) gd.getCheckboxes().get(0);
      checkbox.addItemListener(worker);
      label = (Label) gd.getMessage();

      // Initialise
      worker.refresh();

      // Listen for changes to an image
      ImagePlus.addImageListener(worker);

      thread = new Thread(worker);
      thread.setDaemon(true);
      thread.start();
    }
    gd.showDialog();
    if (worker != null) {
      ImagePlus.removeImageListener(worker);
    }
    if (!gd.wasCanceled()) {
      settings.name = gd.getNextChoice();
      settings.showTable = gd.getNextBoolean();
    }
    if (thread != null && worker != null) {
      worker.running = false;
      inbox.close(true);
      try {
        thread.join(0);
      } catch (final InterruptedException ex) {
        Logger.getLogger(getClass().getName()).log(Level.WARNING, "Unexpected interruption", ex);
        Thread.currentThread().interrupt();
      }
    }
  }
}
