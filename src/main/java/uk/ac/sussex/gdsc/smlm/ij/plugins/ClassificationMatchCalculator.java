/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import com.thoughtworks.xstream.converters.ConversionException;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.Prefs;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;
import uk.ac.sussex.gdsc.core.ij.BufferedTextWindow;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.MatchCalculator;
import uk.ac.sussex.gdsc.core.match.PointPair;
import uk.ac.sussex.gdsc.core.match.RandIndex;
import uk.ac.sussex.gdsc.core.match.Resequencer;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultPoint;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyzrResultProcedure;

/**
 * Compares the classification of localisations in two sets of results and computes the match
 * statistics.
 */
public class ClassificationMatchCalculator implements PlugIn {
  private static final String TITLE = "Classification Calculator";

  private static AtomicReference<TextWindow> resultsWindowRef = new AtomicReference<>();

  /** An empty coordinate array. */
  private static final Coordinate[] EMPTY_COORD_ARRAY = new Coordinate[0];

  private static final String[] ANALYSIS_OPTION =
      SettingsManager.getNames((Object[]) ClassAnalysis.values());

  /** The plugin settings. */
  private Settings settings;

  private enum ClassAnalysis implements NamedObject {
    IGNORE("Ignore"), IGNORE_ZERO("Ignore zero"), ALL("All");

    static ClassAnalysis fromNumber(int number) {
      if (number == ClassAnalysis.ALL.ordinal()) {
        return ClassAnalysis.ALL;
      }
      if (number == ClassAnalysis.IGNORE_ZERO.ordinal()) {
        return ClassAnalysis.IGNORE_ZERO;
      }
      return ClassAnalysis.IGNORE;
    }

    final String name;

    ClassAnalysis(String name) {
      this.name = name;
    }

    @Override
    public String getName() {
      return name;
    }
  }

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    private static final String KEY_INPUT_OPTION1 =
        "gdsc.smlm.classificationmatchcalculator.inputOption1";
    private static final String KEY_INPUT_OPTION2 =
        "gdsc.smlm.classificationmatchcalculator.inputOption2";
    private static final String KEY_MATCH_DISTANCE =
        "gdsc.smlm.classificationmatchcalculator.distanceThreshold";
    private static final String KEY_USE_ID = "gdsc.smlm.classificationmatchcalculator.useId";
    private static final String KEY_USE_CATEGORY =
        "gdsc.smlm.classificationmatchcalculator.useCategory";

    String inputOption1;
    String inputOption2;
    double matchDistance;
    ClassAnalysis useId;
    ClassAnalysis useCategory;

    Settings() {
      inputOption1 = Prefs.get(KEY_INPUT_OPTION1, "");
      inputOption2 = Prefs.get(KEY_INPUT_OPTION2, "");
      matchDistance = Prefs.get(KEY_MATCH_DISTANCE, 0.01);
      useId = ClassAnalysis.fromNumber(Prefs.getInt(KEY_USE_ID, ClassAnalysis.IGNORE.ordinal()));
      useCategory = ClassAnalysis
          .fromNumber(Prefs.getInt(KEY_USE_CATEGORY, ClassAnalysis.IGNORE_ZERO.ordinal()));
    }

    Settings(Settings source) {
      this.inputOption1 = source.inputOption1;
      this.inputOption2 = source.inputOption2;
      this.matchDistance = source.matchDistance;
      this.useId = source.useId;
      this.useCategory = source.useCategory;
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
      Prefs.set(KEY_INPUT_OPTION1, inputOption1);
      Prefs.set(KEY_INPUT_OPTION2, inputOption2);
      Prefs.set(KEY_MATCH_DISTANCE, matchDistance);
      Prefs.set(KEY_USE_ID, useId.ordinal());
      Prefs.set(KEY_USE_CATEGORY, useCategory.ordinal());
    }
  }

  private interface Mapper {
    /**
     * The number of unique mappings.
     *
     * @return the int
     */
    int size();

    /**
     * Map the key to a value.
     *
     * @param key the key
     * @return the value
     */
    int map(int key);

    /**
     * Return a mapper that maps all keys to a single value.
     *
     * @return the mapper
     */
    static Mapper single() {
      return new Mapper() {
        @Override
        public int size() {
          return 1;
        }

        @Override
        public int map(int key) {
          return 0;
        }
      };
    }

    /**
     * Return a mapper that maps all keys to an offset value.
     *
     * @param size the size
     * @param offset the offset
     * @return the mapper
     */
    static Mapper offset(int size, int offset) {
      if (offset == 0) {
        // Identity mapper
        return new Mapper() {
          @Override
          public int size() {
            return size;
          }

          @Override
          public int map(int key) {
            return key;
          }
        };
      }
      return new Mapper() {
        @Override
        public int size() {
          return size;
        }

        @Override
        public int map(int key) {
          return key - offset;
        }
      };
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    // Load the results
    final MemoryPeakResults results1 =
        ResultsManager.loadInputResults(settings.inputOption1, false, null, null);
    final MemoryPeakResults results2 =
        ResultsManager.loadInputResults(settings.inputOption2, false, null, null);
    IJ.showStatus("");
    if (results1 == null || results1.size() == 0) {
      IJ.error(TITLE, "No results 1 could be loaded");
      return;
    }
    if (results2 == null || results2.size() == 0) {
      IJ.error(TITLE, "No results 2 could be loaded");
      return;
    }
    // Check the results can be loaded in the pixels
    try {
      results1.getDistanceConverter(DistanceUnit.PIXEL);
      results2.getDistanceConverter(DistanceUnit.PIXEL);
    } catch (ConversionException | ConfigurationException ex) {
      IJ.error(TITLE, "Distances cannot be loaded in pixels");
      return;
    }

    final long start = System.nanoTime();
    runCompareClassifications(results1, results2);
    final long nanoseconds = System.nanoTime() - start;

    IJ.showStatus(String.format("%s = %s", TITLE, TextUtils.nanosToString(nanoseconds)));
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    settings = Settings.load();
    gd.addMessage("Compare the points in two results sets\nand compute the match statistics");
    ResultsManager.addInput(gd, "Results1", settings.inputOption1, InputSource.MEMORY_CLUSTERED,
        InputSource.MEMORY_CATEGORY);
    ResultsManager.addInput(gd, "Results2", settings.inputOption2, InputSource.MEMORY_CLUSTERED,
        InputSource.MEMORY_CATEGORY);
    gd.addNumericField("Match_distance", settings.matchDistance, -2, 6, "px");
    gd.addChoice("Use_id", ANALYSIS_OPTION, settings.useId.ordinal());
    gd.addChoice("Use_category", ANALYSIS_OPTION, settings.useCategory.ordinal());

    gd.addHelp(HelpUrls.getUrl("classification-match-calculator"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption1 = gd.getNextChoice();
    settings.inputOption2 = gd.getNextChoice();
    settings.matchDistance = gd.getNextNumber();
    settings.useId = ClassAnalysis.fromNumber(gd.getNextChoiceIndex());
    settings.useCategory = ClassAnalysis.fromNumber(gd.getNextChoiceIndex());
    settings.save();

    if (settings.useId.ordinal() + settings.useCategory.ordinal() == 0) {
      IJ.error(TITLE, "No classifications specified (id or category)");
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isPositive("Match distance", settings.matchDistance);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private void runCompareClassifications(MemoryPeakResults results1, MemoryPeakResults results2) {
    final List<PointPair> allMatches = new LinkedList<>();

    // Optionally exclude results which do not have an id and/or category
    Predicate<PeakResult> test =
        settings.useId == ClassAnalysis.IGNORE_ZERO ? r -> r.getId() != 0 : null;
    if (settings.useCategory == ClassAnalysis.IGNORE_ZERO) {
      final Predicate<PeakResult> test2 = r -> r.getCategory() != 0;
      test = test == null ? test2 : test.and(test2);
    } else if (test == null) {
      test = r -> true;
    }

    // Divide the results into time points
    final TIntObjectHashMap<List<PeakResultPoint>> coordinates1 = getCoordinates(results1, test);
    final TIntObjectHashMap<List<PeakResultPoint>> coordinates2 = getCoordinates(results2, test);

    // Process each time point
    int n1 = 0;
    int n2 = 0;
    for (final int t : getTimepoints(coordinates1, coordinates2)) {
      final Coordinate[] c1 = getCoordinates(coordinates1, t);
      final Coordinate[] c2 = getCoordinates(coordinates2, t);
      n1 += c1.length;
      n2 += c2.length;

      final List<PointPair> matches = new LinkedList<>();

      MatchCalculator.analyseResults3D(c1, c2, settings.matchDistance, null, null, null, matches);

      allMatches.addAll(matches);
    }

    if (allMatches.isEmpty()) {
      IJ.error(TITLE, "No localisation matches between the two results sets");
      return;
    }

    // Get the unique Ids and Categories in the matches.
    final Mapper ids = getMapper(allMatches, PeakResult::getId, settings.useId);
    final Mapper cats = getMapper(allMatches, PeakResult::getCategory, settings.useCategory);

    // Map id/category to an index = stride * cat + id
    final int stride = ids.size();

    // Any integer is allowed as an index
    if ((long) stride * cats.size() > 1L << 32) {
      IJ.error(TITLE, "Too many combinations of id and category to assigne unique labels");
      return;
    }


    // Extract indices
    final int[] set1 = new int[allMatches.size()];
    final int[] set2 = new int[allMatches.size()];
    int i = 0;

    for (final PointPair r : allMatches) {
      set1[i] = toIndex(stride, ids, cats, ((PeakResultPoint) r.getPoint1()).getPeakResult());
      set2[i] = toIndex(stride, ids, cats, ((PeakResultPoint) r.getPoint2()).getPeakResult());
      i++;
    }

    final Resequencer re = new Resequencer();
    re.setCacheMap(true);
    re.renumber(set1);
    re.renumber(set2);

    // Compare
    final RandIndex r = new RandIndex().compute(set1, set2);

    final TextWindow resultsWindow = ImageJUtils.refresh(resultsWindowRef,
        () -> new TextWindow(TITLE + " Results",
            "Results1\tResults2\tID\tCategory\tn1\tc1\tn2\tc2\tMatched\tRand Index\tAdjusted RI", "",
            900, 300));
    try (BufferedTextWindow bw = new BufferedTextWindow(resultsWindow)) {
      final StringBuilder sb = new StringBuilder(2048);
      sb.append(results1.getName()).append('\t');
      sb.append(results2.getName()).append('\t');
      sb.append(ANALYSIS_OPTION[settings.useId.ordinal()]).append('\t');
      sb.append(ANALYSIS_OPTION[settings.useCategory.ordinal()]).append('\t');
      sb.append(n1).append('\t');
      sb.append(MathUtils.max(set1) + 1).append('\t');
      sb.append(n2).append('\t');
      sb.append(MathUtils.max(set2) + 1).append('\t');
      sb.append(set1.length).append('\t');
      sb.append(MathUtils.rounded(r.getRandIndex())).append('\t');
      sb.append(MathUtils.rounded(r.getAdjustedRandIndex())).append('\t');
      bw.append(sb.toString());
    }
  }

  /**
   * Gets the mapper that can create a value from a natural sequence starting from 0 for each unique
   * key in the results. If the analysis is set to ignore then a single mapping to zero is created.
   *
   * @param allMatches the all matches
   * @param fun the function to get the key value
   * @param analysis the type of analysis
   * @return the mapper
   */
  private static Mapper getMapper(List<PointPair> allMatches, ToIntFunction<PeakResult> fun,
      ClassAnalysis analysis) {
    if (analysis == ClassAnalysis.IGNORE) {
      return Mapper.single();
    }

    // Find the unique values
    final TIntHashSet set = new TIntHashSet();
    for (final PointPair r : allMatches) {
      set.add(fun.applyAsInt(((PeakResultPoint) r.getPoint1()).getPeakResult()));
      set.add(fun.applyAsInt(((PeakResultPoint) r.getPoint2()).getPeakResult()));
    }

    // Edge case of 1 value
    if (set.size() == 1) {
      return Mapper.single();
    }

    // Map to a natural sequence from zero
    final int[] keys = set.toArray();
    Arrays.sort(keys);

    // Check if a discrete sequence already
    if (keys[keys.length - 1] - keys[0] == set.size() - 1) {
      return Mapper.offset(set.size(), keys[0]);
    }

    // Map each key to a value starting from 0
    final TIntIntHashMap map = new TIntIntHashMap(keys.length);
    for (final int k : keys) {
      map.put(k, map.size());
    }

    return new Mapper() {
      @Override
      public int size() {
        return map.size();
      }

      @Override
      public int map(int key) {
        return map.get(key);
      }
    };
  }

  /**
   * Convert the peak to an index representing the id and/or category combination.
   *
   * @param stride the stride
   * @param ids the ids mapper
   * @param cats the categories mapper
   * @param r the peak
   * @return the index
   */
  private static int toIndex(int stride, Mapper ids, Mapper cats, PeakResult r) {
    return stride * cats.map(r.getCategory()) + ids.map(r.getId());
  }

  /**
   * Build a map between the peak id (time point) and a list of coordinates that pass the filter.
   *
   * @param results the results
   * @param test the test
   * @return the coordinates
   */
  public static TIntObjectHashMap<List<PeakResultPoint>> getCoordinates(MemoryPeakResults results,
      Predicate<PeakResult> test) {
    final TIntObjectHashMap<List<PeakResultPoint>> coords = new TIntObjectHashMap<>();
    if (results.size() > 0) {
      // Do not use HashMap directly to build the coords object since there
      // will be many calls to getEntry(). Instead sort the results and use
      // a new list for each time point
      results.sort();

      // Create list
      final LocalList<PeakResultPoint> tmpCoords = new LocalList<>();

      // Add the results for each frame
      final FrameCounter counter = results.newFrameCounter();
      results.forEach(DistanceUnit.PIXEL, (XyzrResultProcedure) (x, y, z, r) -> {
        if (counter.advance(r.getFrame()) && !tmpCoords.isEmpty()) {
          coords.put(counter.previousFrame(), tmpCoords.copy());
          tmpCoords.clear();
        }
        if (test.test(r)) {
          tmpCoords.add(new PeakResultPoint(r.getFrame(), x, y, z, r));
        }
      });
      if (!tmpCoords.isEmpty()) {
        coords.put(counter.currentFrame(), tmpCoords.copy());
      }
    }
    return coords;
  }

  /**
   * Return an array of coordinates for the given time point. Returns an empty array if there are no
   * coordinates.
   *
   * @param coords the coords
   * @param time the time
   * @return the coordinates
   */
  public static Coordinate[] getCoordinates(TIntObjectHashMap<List<PeakResultPoint>> coords,
      int time) {
    final List<PeakResultPoint> tmp = coords.get(time);
    if (tmp != null) {
      return tmp.toArray(EMPTY_COORD_ARRAY);
    }
    return EMPTY_COORD_ARRAY;
  }

  /**
   * Merge the time points from each map into a single sorted list of unique time points.
   *
   * @param actualCoordinates the actual coordinates
   * @param predictedCoordinates the predicted coordinates
   * @return a list of time points
   */
  private static int[] getTimepoints(TIntObjectHashMap<List<PeakResultPoint>> actualCoordinates,
      TIntObjectHashMap<List<PeakResultPoint>> predictedCoordinates) {

    // Do inline to avoid materialising the keys arrays
    final TIntHashSet hashset =
        new TIntHashSet(Math.max(actualCoordinates.size(), predictedCoordinates.size()));
    final TIntProcedure p = value -> {
      hashset.add(value);
      return true;
    };
    actualCoordinates.forEachKey(p);
    predictedCoordinates.forEachKey(p);
    final int[] set = hashset.toArray();

    Arrays.sort(set);
    return set;
  }
}
