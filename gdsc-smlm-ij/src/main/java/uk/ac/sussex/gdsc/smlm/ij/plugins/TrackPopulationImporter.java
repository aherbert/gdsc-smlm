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
import ij.Prefs;
import ij.plugin.PlugIn;
import it.unimi.dsi.fastutil.longs.Long2IntMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.regex.Pattern;
import org.apache.commons.rng.core.util.NumberFactory;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.OpenHashMaps.CustomLong2IntOpenHashMap;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Plugin to import population categories assigned to an existing track dataset. Maps the category
 * to a unique localisation result using the key 'Frame:ID'.
 */
public class TrackPopulationImporter implements PlugIn {
  private static final String TITLE = "Track Population Importer";
  private static final Pattern CSV = Pattern.compile(", *");

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    private static final String KEY_FILENAME = "gdsc.smlm.trackpopulationimporter.filename";
    private static final String KEY_INPUT_OPTION = "gdsc.smlm.trackpopulationimporter.inputOption";
    private static final String KEY_IGNORE_UNMAPPED =
        "gdsc.smlm.trackpopulationimporter.ignoreUnmapped";
    private static final String KEY_NEW_DATASET = "gdsc.smlm.trackpopulationimporter.newDataset";

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String filename;
    String inputOption;
    boolean ignoreUnmapped;
    boolean newDataset;

    Settings() {
      // Set defaults
      filename = Prefs.get(KEY_FILENAME, "");
      inputOption = Prefs.get(KEY_INPUT_OPTION, "");
      ignoreUnmapped = Prefs.get(KEY_IGNORE_UNMAPPED, false);
      newDataset = Prefs.get(KEY_NEW_DATASET, false);
    }

    Settings(Settings source) {
      filename = source.filename;
      inputOption = source.inputOption;
      ignoreUnmapped = source.ignoreUnmapped;
      newDataset = source.newDataset;
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
      Prefs.set(KEY_FILENAME, filename);
      Prefs.set(KEY_INPUT_OPTION, inputOption);
      Prefs.set(KEY_IGNORE_UNMAPPED, ignoreUnmapped);
      Prefs.set(KEY_NEW_DATASET, newDataset);
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    final MemoryResultsList items = new MemoryResultsList(MemoryPeakResults::hasId);

    if (items.isEmpty()) {
      IJ.error(TITLE, "No traced localisations in memory");
      return;
    }

    // Get input options
    if (!showDialog()) {
      return;
    }

    final MemoryPeakResults results = MemoryPeakResults.getResults(settings.inputOption);

    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    // Overwrite or create new dataset
    final MemoryPeakResults out;
    if (settings.newDataset) {
      out = results.copy(true);
      out.setName(out.getName() + " (Categorised)");
    } else {
      out = results;
    }
    // Working data
    final PeakResult[] data = out.toArray();

    // Mapping methods will throw exceptions
    try {
      // Create a map of 'Frame:ID -> result'
      final Long2IntMap resultMap = createResultMap(data);

      // Read the mapping 'Frame:ID -> Category'
      CustomLong2IntOpenHashMap categoryMap;
      try (BufferedReader br = Files.newBufferedReader(Paths.get(settings.filename))) {
        categoryMap = createCategoryMap(br);
      }

      final BitSet assigned = assignCategory(data, resultMap, categoryMap);
      ImageJPluginLoggerHelper.getLogger(TrackPopulationImporter.class).log(Level.INFO,
          () -> String.format("Assigned %d/%d", assigned.cardinality(), data.length));

      // Set all other categories to zero (optional)
      if (!settings.ignoreUnmapped) {
        assigned.flip(0, data.length);
        assigned.stream().forEach(i -> data[i] = changeCategory(data[i], 0));
      }
    } catch (final Exception ex) {
      handleException(ex);
    }

    // Save
    out.begin();
    out.addAll(data);
    out.end();
  }

  private boolean showDialog() {
    settings = Settings.load();
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Import track populations.\n \nInput format:\nFrame,ID,Category\n \n"
        + "Assigns 'Frame:ID' -> Category in the specified results.\n"
        + "ID/Category <= 0 are ignored.");
    gd.addFilenameField("Filename", settings.filename, 30);
    ResultsManager.addInput(gd, "Results", settings.inputOption, InputSource.MEMORY_CLUSTERED);
    gd.addCheckbox("Ignore_unmapped", settings.ignoreUnmapped);
    gd.addCheckbox("New_dataset", settings.newDataset);
    gd.addHelp(HelpUrls.getUrl("track-population-importer"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.filename = gd.getNextString();
    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.ignoreUnmapped = gd.getNextBoolean();
    settings.newDataset = gd.getNextBoolean();
    settings.save();
    return true;
  }

  /**
   * Creates the result map.
   *
   * <pre>
   * Frame:ID -> result index
   * </pre>
   *
   * @param data the data
   * @return the map
   * @throws IllegalArgumentException if the results contain duplicate keys
   */
  static Long2IntMap createResultMap(PeakResult[] data) {
    final Long2IntOpenHashMap map = new Long2IntOpenHashMap(data.length);
    map.defaultReturnValue(-1);
    for (int i = 0; i < data.length; i++) {
      final PeakResult r = data[i];
      final int id = r.getId();
      if (id <= 0) {
        continue;
      }
      final int frame = r.getFrame();
      final long key = getKey(frame, id);
      if (map.put(key, i) != -1) {
        throw new IllegalArgumentException("Duplicate key: " + frame + ":" + id);
      }
    }
    return map;
  }

  /**
   * Creates the category map.
   *
   * <pre>
   * Frame:ID -> category
   * </pre>
   *
   * @param reader the reader
   * @return the map
   * @throws IOException Signals that an I/O exception has occurred.
   * @throws IllegalArgumentException if the category map does not contain 3 integer fields per
   *         line; or the map contains duplicate keys
   */
  static CustomLong2IntOpenHashMap createCategoryMap(BufferedReader reader) throws IOException {
    final CustomLong2IntOpenHashMap map = new CustomLong2IntOpenHashMap(16);
    map.defaultReturnValue(-1);
    for (String line = reader.readLine(); line != null; line = reader.readLine()) {
      final String[] fields = CSV.split(line);
      if (fields.length != 3) {
        throw new IllegalArgumentException("Expected 3 fields: " + line);
      }
      try {
        final int frame = Integer.parseInt(fields[0]);
        final int id = Integer.parseInt(fields[1]);
        final int category = Integer.parseInt(fields[2]);
        if (id <= 0 || category <= 0) {
          continue;
        }
        final long key = getKey(frame, id);
        if (map.put(key, category) != -1) {
          throw new IllegalArgumentException("Duplicate key: " + frame + ":" + id);
        }
      } catch (final NumberFormatException ex) {
        throw new IllegalArgumentException("Expected integer fields: " + line, ex);
      }
    }
    return map;
  }

  /**
   * Assign the category to the result.
   *
   * @param data the data
   * @param resultMap the result map
   * @param categoryMap the category map
   * @return the bit set
   * @throws IllegalArgumentException if the category map contains a key not in the results
   */
  static BitSet assignCategory(PeakResult[] data, Long2IntMap resultMap,
      CustomLong2IntOpenHashMap categoryMap) {
    final BitSet assigned = new BitSet(data.length);
    categoryMap.forEach((long key, int category) -> {
      final int index = resultMap.get(key);
      if (index < 0) {
        throw new IllegalArgumentException(
            "Key not present in the results: " + getFrame(key) + ":" + getId(key));
      }
      // Create a new result (if necessary)
      data[index] = changeCategory(data[index], category);
      assigned.set(index);
    });

    return assigned;
  }

  /**
   * Gets the key.
   *
   * @param frame the frame
   * @param id the id
   * @return the key
   */
  static long getKey(int frame, int id) {
    return NumberFactory.makeLong(frame, id);
  }

  /**
   * Gets the frame.
   *
   * @param key the key
   * @return the frame
   */
  static int getFrame(long key) {
    return NumberFactory.extractHi(key);
  }

  /**
   * Gets the id.
   *
   * @param key the key
   * @return the id
   */
  static int getId(long key) {
    return NumberFactory.extractLo(key);
  }

  private static void handleException(Exception ex) {
    IJ.error(TITLE, "Failed to import: " + ex.getMessage());
    ImageJPluginLoggerHelper.getLogger(TrackPopulationImporter.class).log(Level.SEVERE,
        "Failed to import", ex);
  }

  /**
   * Change the category of the result.
   *
   * @param result the result
   * @param category the category
   * @return the updated result
   */
  static PeakResult changeCategory(PeakResult result, int category) {
    if (result.getCategory() != category) {
      final AttributePeakResult r2 =
          result instanceof AttributePeakResult ? (AttributePeakResult) result
              : new AttributePeakResult(result);
      r2.setCategory(category);
      return r2;
    }
    return result;
  }
}
