/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import java.awt.Point;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.MatchCalculator;
import uk.ac.sussex.gdsc.core.match.MatchResult;
import uk.ac.sussex.gdsc.core.match.PointPair;
import uk.ac.sussex.gdsc.core.match.Pulse;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageRoiPainter;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.utils.CoordinateProvider;

/**
 * Compares the coordinates in sets of traced results and computes the match statistics.
 */
public class TraceMatchCalculator implements PlugIn {
  private static final String TITLE = "Trace Match Calculator";

  private static AtomicBoolean writeHeader = new AtomicBoolean(true);
  private static final AtomicReference<TextWindow> RESULTS_WINDOW_REF = new AtomicReference<>();
  private static final AtomicReference<WindowAndPainter> PAIRS_WINDOW_REF = new AtomicReference<>();
  private static final AtomicReference<WindowAndPainter> TRIPLES_WINDOW_REF =
      new AtomicReference<>();

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    private static final String[] SORT_OPTIONS = {"Score", "Time"};

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    String inputOption1;
    String inputOption2;
    String inputOption3;
    double distanceThreshold;
    double beta;
    boolean showPairs;
    int sortIndex;

    Settings() {
      // Set defaults
      inputOption1 = "";
      inputOption2 = "";
      inputOption3 = "";
      distanceThreshold = 1;
      beta = 4;
      sortIndex = 1;
    }

    Settings(Settings source) {
      inputOption1 = source.inputOption1;
      inputOption2 = source.inputOption2;
      inputOption3 = source.inputOption3;
      distanceThreshold = source.distanceThreshold;
      beta = source.beta;
      showPairs = source.showPairs;
      sortIndex = source.sortIndex;
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
    }
  }

  /**
   * Class to allow atomic update of the text window and the painter.
   */
  private static class WindowAndPainter {
    final TextWindow textWindow;
    final ImageRoiPainter painter;

    WindowAndPainter(TextWindow textWindow, ImageRoiPainter painter) {
      this.textWindow = textWindow;
      this.painter = painter;
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
    IJ.showStatus("");
    if (results1 == null || results1.size() == 0) {
      IJ.error(TITLE, "No results 1 could be loaded");
      return;
    }
    final MemoryPeakResults results2 =
        ResultsManager.loadInputResults(settings.inputOption2, false, null, null);
    if (results2 == null || results2.size() == 0) {
      IJ.error(TITLE, "No results 2 could be loaded");
      return;
    }
    if (results1.getDistanceUnit() != results2.getDistanceUnit()) {
      IJ.error(TITLE, "Distance unit should be the same for the results 1 & 2");
      return;
    }
    final MemoryPeakResults results3 =
        ResultsManager.loadInputResults(settings.inputOption3, false, null, null);
    if (results3 != null && results1.getDistanceUnit() != results3.getDistanceUnit()) {
      IJ.error(TITLE, "Distance unit should be the same for the results 1 & 3");
      return;
    }

    final long start = System.nanoTime();
    compareCoordinates(results1, results2, results3, settings.distanceThreshold);
    final double seconds = (System.nanoTime() - start) / 1000000000.0;

    IJ.showStatus(String.format("%s = %ss", TITLE, MathUtils.rounded(seconds, 4)));
  }

  private boolean showDialog() {
    settings = Settings.load();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addMessage("Compare the traced points in results sets\nand computes the match statistics");
    ResultsManager.addInput(gd, "Results1", settings.inputOption1, InputSource.MEMORY_MULTI_FRAME);
    ResultsManager.addInput(gd, "Results2", settings.inputOption2, InputSource.MEMORY_MULTI_FRAME);
    ResultsManager.addInput(gd, "Results3", settings.inputOption3, InputSource.NONE,
        InputSource.MEMORY_MULTI_FRAME);
    gd.addNumericField("Distance", settings.distanceThreshold, 2, 6, "px");

    gd.addNumericField("Beta", settings.beta, 2);
    gd.addCheckbox("Show_pairs", settings.showPairs);
    gd.addChoice("Sort_pairs", Settings.SORT_OPTIONS, settings.sortIndex);

    gd.addHelp(HelpUrls.getUrl("trace-match-calculator"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption1 = gd.getNextChoice();
    settings.inputOption2 = gd.getNextChoice();
    settings.inputOption3 = gd.getNextChoice();
    settings.distanceThreshold = gd.getNextNumber();
    settings.beta = gd.getNextNumber();
    settings.showPairs = gd.getNextBoolean();
    settings.sortIndex = gd.getNextChoiceIndex();
    settings.save();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Distance threshold", settings.distanceThreshold);
      ParameterUtils.isPositive("Beta", settings.beta);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  @SuppressWarnings("null")
  private void compareCoordinates(MemoryPeakResults results1, MemoryPeakResults results2,
      MemoryPeakResults results3, double distanceThreshold) {
    final Pulse[] p1 = extractPulses(results1);
    final Pulse[] p2 = extractPulses(results2);
    final Pulse[] p3 = extractPulses(results3);

    final List<Pulse> tp = null;
    List<Pulse> fp = null;
    List<Pulse> fn = null;
    List<PointPair> pairs = null;

    final List<Pulse> tp2 = null;
    List<Pulse> fp2 = null;
    List<Pulse> fn2 = null;
    List<PointPair> pairs2 = null;

    if (settings.showPairs) {
      pairs = new LinkedList<>();
      fp = new LinkedList<>();
      fn = new LinkedList<>();
      pairs2 = new LinkedList<>();
      fp2 = new LinkedList<>();
      fn2 = new LinkedList<>();
    }

    final MatchResult result =
        MatchCalculator.analyseResults2D(p1, p2, distanceThreshold, tp, fp, fn, pairs);
    final MatchResult result2 = p3 != null
        ? MatchCalculator.analyseResults2D(p1, p3, distanceThreshold, tp2, fp2, fn2, pairs2)
        : null;

    // Create output
    Consumer<String> resultsOutput;
    if (!java.awt.GraphicsEnvironment.isHeadless()) {
      final TextWindow resultsWindow = ImageJUtils.refresh(RESULTS_WINDOW_REF,
          () -> new TextWindow(TITLE + " Results", createResultsHeader(), "", 900, 300));
      resultsOutput = resultsWindow::append;

      if (settings.showPairs) {
        if (p3 == null) {
          // Produce a pairs output
          final WindowAndPainter wap = refresh(PAIRS_WINDOW_REF, true, resultsWindow, results1);

          // Add the unmatched points
          WindowManager.getIDList();

          for (final Coordinate c : fn) {
            pairs.add(new PointPair(c, null));
          }
          for (final Coordinate c : fp) {
            pairs.add(new PointPair(null, c));
          }

          final List<? extends PointPair> sortedPairs = sortPointPairs(pairs);

          for (final PointPair pair : sortedPairs) {
            addPairResult(wap.textWindow, pair);
          }
        } else {
          // Produce a triple output
          final WindowAndPainter wap = refresh(TRIPLES_WINDOW_REF, false, resultsWindow, results1);

          final HashMap<Pulse, Triple> map = new HashMap<>();
          final ArrayList<Triple> triples = new ArrayList<>(pairs.size());
          for (final PointPair pair : pairs) {
            final Pulse p = (Pulse) pair.getPoint1();
            final Triple t = new Triple(p, (Pulse) pair.getPoint2(), null);
            triples.add(t);
            map.put(p, t);
          }
          // Complete the reference set of points
          for (final Coordinate c : fn) {
            final Pulse p = (Pulse) c;
            final Triple t = new Triple(p, null, null);
            triples.add(t);
            map.put(p, t);
          }

          // Add the unmatched points
          for (final Coordinate c : fp) {
            triples.add(new Triple(null, (Pulse) c, null));
          }
          for (final Coordinate c : fp2) {
            triples.add(new Triple(null, null, (Pulse) c));
          }

          // Add the results from the second match
          for (final PointPair pair : pairs2) {
            final Pulse p = (Pulse) pair.getPoint1();
            final Pulse pp = (Pulse) pair.getPoint2();
            final Triple triple = map.get(p);
            if (triple != null) {
              triple.p3 = pp;
            } else {
              triples.add(new Triple(null, null, pp));
            }
          }

          final List<? extends Triple> sortedTriples = sortTriples(triples);

          for (final Triple t : sortedTriples) {
            addTripleResult(wap.textWindow, t);
          }
        }
      }
    } else {
      if (writeHeader.compareAndSet(true, false)) {
        IJ.log(createResultsHeader());
      }
      resultsOutput = IJ::log;
    }

    final StringBuilder sb = new StringBuilder();
    addResult(resultsOutput, sb, settings.inputOption1, settings.inputOption2, distanceThreshold,
        result);
    if (result2 != null) {
      addResult(resultsOutput, sb, settings.inputOption1, settings.inputOption3, distanceThreshold,
          result2);
    }
  }

  private static WindowAndPainter refresh(AtomicReference<WindowAndPainter> ref, boolean pairs,
      TextWindow resultsWindow, MemoryPeakResults results1) {
    // Produce a pairs output
    final WindowAndPainter wap = ConcurrencyUtils.refresh(ref,
        // Test the window is showing
        w -> ImageJUtils.isShowing(w.textWindow),
        // Create
        () -> {
          final String title = TITLE + ((pairs) ? " Pairs" : " Triples");
          final String header = pairs ? createPairsHeader() : createTriplesHeader();
          final TextWindow window = new TextWindow(title, header, "", 900, 300);

          // Position relative to results window
          final Point p = resultsWindow.getLocation();
          p.y += resultsWindow.getHeight();
          window.setLocation(p);

          final CoordinateProvider coordinateProvider = line -> {
            // Extract the startT and x,y coordinates from the first pulse in the line
            final int[] index = {0, 5, 12};
            final String[] fields = line.split("\t");
            for (final int i : index) {
              if (i < fields.length) {
                if ("-".equals(fields[i])) {
                  continue;
                }
                final int startT = Integer.parseInt(fields[i]);
                final double x = Double.parseDouble(fields[i + 2]);
                final double y = Double.parseDouble(fields[i + 3]);
                return new double[] {startT, x, y};
              }
            }
            return null;
          };

          final ImageRoiPainter painter = new ImageRoiPainter(window.getTextPanel(),
              results1.getSource().getOriginal().getName(), coordinateProvider);
          final WindowAndPainter result = new WindowAndPainter(window, painter);

          // Free memory on close
          window.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosed(WindowEvent event) {
              ref.compareAndSet(result, null);
              super.windowClosed(event);
            }
          });

          return result;
        });

    wap.textWindow.getTextPanel().clear();
    wap.painter.setTitle(results1.getSource().getOriginal().getName());
    return wap;
  }

  @Nullable
  private static Pulse[] extractPulses(MemoryPeakResults results) {
    if (results == null) {
      return null;
    }
    final Pulse[] pulses = new Pulse[results.size()];
    final Counter i = new Counter();
    results.forEach(DistanceUnit.PIXEL,
        (XyrResultProcedure) (x, y, result) -> pulses[i.getAndIncrement()] =
            new Pulse(x, y, result.getFrame(), result.getEndFrame()));
    return pulses;
  }

  private static String toHeader(String[] names) {
    final StringBuilder sb = new StringBuilder(256);
    Arrays.stream(names).forEach(x -> sb.append(x).append('\t'));
    sb.setLength(sb.length() - 1);
    return sb.toString();
  }

  private static String createResultsHeader() {
    final String[] names = {
        // @formatter:off
        "Image 1",
        "Image 2",
        "Distance (px)",
        "N",
        "TP",
        "FP",
        "FN",
        "Jaccard",
        "Score",
        "Precision",
        "Recall",
        "F0.5",
        "F1",
        "F2",
        "F-beta",
        // @formatter:on
    };
    return toHeader(names);
  }

  private void addResult(Consumer<String> output, StringBuilder sb, String i1, String i2,
      double distanceThrehsold, MatchResult result) {
    sb.setLength(0);
    // @formatter:off
    sb.append(i1).append('\t')
      .append(i2).append('\t')
      .append(IJ.d2s(distanceThrehsold, 2)).append('\t')
      .append(result.getNumberPredicted()).append('\t')
      .append(result.getTruePositives()).append('\t')
      .append(result.getFalsePositives()).append('\t')
      .append(result.getFalseNegatives()).append('\t')
      .append(IJ.d2s(result.getJaccard(), 4)).append('\t')
      .append(IJ.d2s(result.getRmsd(), 4)).append('\t')
      .append(IJ.d2s(result.getPrecision(), 4)).append('\t')
      .append(IJ.d2s(result.getRecall(), 4)).append('\t')
      .append(IJ.d2s(result.getFScore(0.5), 4)).append('\t')
      .append(IJ.d2s(result.getFScore(1.0), 4)).append('\t')
      .append(IJ.d2s(result.getFScore(2.0), 4)).append('\t')
      .append(IJ.d2s(result.getFScore(settings.beta), 4));
    // @formatter:on

    output.accept(sb.toString());
  }

  private static String createPairsHeader() {
    final String[] names = {
        // @formatter:off
        "Start1",
        "End1",
        "X1",
        "Y1",
        "Z1",
        "Start2",
        "End2",
        "X2",
        "Y2",
        "Z2",
        "Distance",
        "Score",
        // @formatter:on
    };
    return toHeader(names);
  }

  private void addPairResult(TextWindow pairsWindow, PointPair pair) {
    final StringBuilder sb = new StringBuilder();
    final Pulse p1 = (Pulse) pair.getPoint1();
    final Pulse p2 = (Pulse) pair.getPoint2();
    addPoint(sb, p1);
    addPointPairResult(sb, p1, p2);
    pairsWindow.append(sb.toString());
  }

  private static void addPoint(StringBuilder sb, Pulse pulse) {
    if (pulse == null) {
      sb.append("-\t-\t-\t-\t-\t");
    } else {
      sb.append(pulse.getStart()).append('\t');
      sb.append(pulse.getEnd()).append('\t');
      sb.append(IJ.d2s(pulse.getX())).append('\t');
      sb.append(IJ.d2s(pulse.getY())).append('\t');
      sb.append(IJ.d2s(pulse.getZ())).append('\t');
    }
  }

  private static String createTriplesHeader() {
    final String[] names = {
        // @formatter:off
        "Start1",
        "End1",
        "X1",
        "Y1",
        "Z1",
        "Start2",
        "End2",
        "X2",
        "Y2",
        "Z2",
        "Distance",
        "Score",
        "Start3",
        "End3",
        "X3",
        "Y3",
        "Z3",
        "Distance2",
        "Score2",
        // @formatter:on
    };
    return toHeader(names);
  }

  private void addTripleResult(TextWindow triplesWindow, Triple triple) {
    final StringBuilder sb = new StringBuilder();
    final Pulse p1 = triple.p1;
    final Pulse p2 = triple.p2;
    final Pulse p3 = triple.p3;
    addPoint(sb, p1);
    addPointPairResult(sb, p1, p2);
    addPointPairResult(sb, p1, p3);
    triplesWindow.append(sb.toString());
  }

  private void addPointPairResult(StringBuilder sb, Pulse p1, Pulse p2) {
    addPoint(sb, p2);
    final PointPair pair = new PointPair(p1, p2);
    final double d = pair.getXyDistance();
    if (d >= 0) {
      sb.append(MathUtils.rounded(d, 4)).append('\t');
    } else {
      sb.append("-\t");
    }
    if (p1 != null && p2 != null) {
      sb.append(MathUtils.rounded(p1.score(p2, d * d, settings.distanceThreshold), 4)).append('\t');
    } else {
      sb.append("-\t");
    }
  }

  private List<? extends PointPair> sortPointPairs(List<PointPair> pairs) {
    if (settings.sortIndex == 1) {
      // Sort by time
      final ArrayList<TimeComparablePointPair> newPairs = new ArrayList<>(pairs.size());
      for (final PointPair pair : pairs) {
        newPairs.add(new TimeComparablePointPair(pair));
      }
      Collections.sort(newPairs, TimeComparablePointPair::compare);
      return newPairs;
    }
    // Already sorted by score
    return pairs;
  }

  private List<? extends Triple> sortTriples(List<Triple> triples) {
    if (settings.sortIndex == 1) {
      // Sort by time
      final List<TimeComparableTriple> sorted = new ArrayList<>(triples.size());
      for (final Triple t : triples) {
        sorted.add(new TimeComparableTriple(t));
      }
      Collections.sort(sorted, TimeComparableTriple::compare);
      return sorted;
    }
    final List<ScoreComparableTriple> sorted = new ArrayList<>(triples.size());
    for (final Triple t : triples) {
      sorted.add(new ScoreComparableTriple(t, settings.distanceThreshold));
    }
    Collections.sort(sorted, ScoreComparableTriple::compare);
    return sorted;
  }

  /**
   * A PointPair that can be compared using time.
   */
  private static class TimeComparablePointPair extends PointPair {
    int startT = Integer.MAX_VALUE;
    final Pulse p1;
    final Pulse p2;

    TimeComparablePointPair(PointPair pair) {
      super(pair.getPoint1(), pair.getPoint2());

      p1 = (Pulse) pair.getPoint1();
      p2 = (Pulse) pair.getPoint2();
      if (p2 != null) {
        startT = p2.getStart();
      }
      if (p1 != null) {
        startT = p1.getStart();
      }
    }

    static int compare(TimeComparablePointPair r1, TimeComparablePointPair r2) {
      // Significance of points are: p1, p2, p3
      // Sort by the earliest start time of the most significant point.
      if (r1.startT == r2.startT) {
        // Sort using the significant points first.
        int result = comparePulsesStartT(r1.p1, r2.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesStartT(r1.p2, r2.p2);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(r1.p1, r2.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(r1.p2, r2.p2);
        if (result != 0) {
          return result;
        }

        // Sort using coords
        result = comparePulsesCoords(r1.p1, r2.p1);
        if (result != 0) {
          return result;
        }
        return comparePulsesCoords(r1.p2, r2.p2);
      }
      return r1.startT - r2.startT;
    }
  }

  /**
   * Compare pulses using the start time.
   *
   * @param p1 pulse 1
   * @param p2 pulse 2
   * @return the result
   */
  static int comparePulsesStartT(Pulse p1, Pulse p2) {
    // Data for a point always beats no data
    if (p1 == null) {
      return (p2 == null) ? 0 : 1;
    }
    if (p2 == null) {
      return -1;
    }
    return p1.getStart() - p2.getStart();
  }

  /**
   * Compare pulses using the end time.
   *
   * @param p1 pulse 1
   * @param p2 pulse 2
   * @return the result
   */
  static int comparePulsesEndT(Pulse p1, Pulse p2) {
    // Data for a point always beats no data
    if (p1 == null) {
      return (p2 == null) ? 0 : 1;
    }
    if (p2 == null) {
      return -1;
    }
    return p1.getEnd() - p2.getEnd();
  }

  /**
   * Compare pulses using coordinates.
   *
   * @param p1 pulse 1
   * @param p2 pulse 2
   * @return the result
   */
  static int comparePulsesCoords(Pulse p1, Pulse p2) {
    // Data for a point always beats no data
    if (p1 == null) {
      return (p2 == null) ? 0 : 1;
    }
    if (p2 == null) {
      return -1;
    }
    if (p1.getX() < p2.getX()) {
      return -1;
    }
    if (p1.getX() > p2.getX()) {
      return 1;
    }
    if (p1.getY() < p2.getY()) {
      return -1;
    }
    if (p1.getY() > p2.getY()) {
      return 1;
    }
    return 0;
  }

  /**
   * A triple of Pulse results.
   */
  private static class Triple {
    final Pulse p1;
    final Pulse p2;
    Pulse p3;

    Triple(Pulse p1, Pulse p2, Pulse p3) {
      this.p1 = p1;
      this.p2 = p2;
      this.p3 = p3;
    }

    Triple(Triple triple) {
      this.p1 = triple.p1;
      this.p2 = triple.p2;
      this.p3 = triple.p3;
    }
  }

  /**
   * A triple that can be compared using time.
   */
  private static class TimeComparableTriple extends Triple {
    int startT = Integer.MAX_VALUE;

    TimeComparableTriple(Triple triple) {
      super(triple);

      if (p3 != null) {
        startT = p3.getStart();
      }
      if (p2 != null) {
        startT = p2.getStart();
      }
      if (p1 != null) {
        startT = p1.getStart();
      }
    }

    static int compare(TimeComparableTriple r1, TimeComparableTriple r2) {
      // Significance of points are: p1, p2, p3
      // Sort by the earliest start time of the most significant point.
      if (r1.startT == r2.startT) {
        // Make the order the same as for the pairs table
        int result = comparePulsesStartT(r1.p1, r2.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesStartT(r1.p2, r2.p2);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(r1.p1, r2.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(r1.p2, r2.p2);
        if (result != 0) {
          return result;
        }

        // Then sort using the third
        result = comparePulsesStartT(r1.p3, r2.p3);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(r1.p3, r2.p3);
        if (result != 0) {
          return result;
        }

        // Sort using coords
        result = comparePulsesCoords(r1.p1, r2.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesCoords(r1.p2, r2.p2);
        if (result != 0) {
          return result;
        }
        return comparePulsesCoords(r1.p3, r2.p3);
      }
      return r1.startT - r2.startT;
    }
  }

  /**
   * A triple that can be compared using score.
   */
  private static class ScoreComparableTriple extends Triple {
    double score;

    ScoreComparableTriple(Triple triple, double distanceThreshold) {
      super(triple);
      if (p1 != null) {
        if (p2 != null) {
          score = p1.score(p2, distanceThreshold);
        }
        if (p3 != null) {
          score = Math.max(score, p1.score(p2, distanceThreshold));
        }
      }
    }

    static int compare(ScoreComparableTriple r1, ScoreComparableTriple r2) {
      if (r1.score > r2.score) {
        return -1;
      }
      if (r1.score < r2.score) {
        return 1;
      }
      return 0;
    }
  }
}
