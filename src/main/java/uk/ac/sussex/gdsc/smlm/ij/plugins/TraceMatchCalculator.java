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

import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.MatchCalculator;
import uk.ac.sussex.gdsc.core.match.MatchResult;
import uk.ac.sussex.gdsc.core.match.PointPair;
import uk.ac.sussex.gdsc.core.match.Pulse;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageROIPainter;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.XYRResultProcedure;
import uk.ac.sussex.gdsc.smlm.utils.CoordinateProvider;

import ij.IJ;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import org.apache.commons.math3.util.FastMath;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

/**
 * Compares the coordinates in sets of traced results and computes the match statistics.
 */
public class TraceMatchCalculator implements PlugIn, CoordinateProvider {
  private static String TITLE = "Trace Match Calculator";

  private static String inputOption1 = "";
  private static String inputOption2 = "";
  private static String inputOption3 = "";
  private static double dThreshold = 1;
  private static double beta = 4;
  private static boolean showPairs;
  private static String[] SORT_OPTIONS = new String[] {"Score", "Time"};
  private static int sortIndex = 1;

  private static boolean writeHeader = true;
  private static TextWindow resultsWindow;
  private static TextWindow pairsWindow;
  private static TextWindow triplesWindow;
  private static ImageROIPainter pairPainter;
  private static ImageROIPainter triplePainter;

  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    // Load the results
    final MemoryPeakResults results1 =
        ResultsManager.loadInputResults(inputOption1, false, null, null);
    final MemoryPeakResults results2 =
        ResultsManager.loadInputResults(inputOption2, false, null, null);
    final MemoryPeakResults results3 =
        ResultsManager.loadInputResults(inputOption3, false, null, null);
    IJ.showStatus("");
    if (results1 == null || results1.size() == 0) {
      IJ.error(TITLE, "No results 1 could be loaded");
      return;
    }
    if (results2 == null || results2.size() == 0) {
      IJ.error(TITLE, "No results 2 could be loaded");
      return;
    }
    if (results1.getDistanceUnit() != results2.getDistanceUnit()) {
      IJ.error(TITLE, "Distance unit should be the same for the results 1 & 2");
      return;
    }
    if (results3 != null && results1.getDistanceUnit() != results3.getDistanceUnit()) {
      IJ.error(TITLE, "Distance unit should be the same for the results 1 & 3");
      return;
    }

    final long start = System.nanoTime();
    compareCoordinates(results1, results2, results3, dThreshold);
    final double seconds = (System.nanoTime() - start) / 1000000000.0;

    IJ.showStatus(String.format("%s = %ss", TITLE, MathUtils.rounded(seconds, 4)));
  }

  private static boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addMessage("Compare the points in two results sets\nand compute the match statistics");
    ResultsManager.addInput(gd, "Results1", inputOption1, InputSource.MEMORY_MULTI_FRAME);
    ResultsManager.addInput(gd, "Results2", inputOption2, InputSource.MEMORY_MULTI_FRAME);
    ResultsManager.addInput(gd, "Results3", inputOption3, InputSource.NONE,
        InputSource.MEMORY_MULTI_FRAME);
    gd.addNumericField("Distance", dThreshold, 2, 6, "px");

    gd.addNumericField("Beta", beta, 2);
    gd.addCheckbox("Show_pairs", showPairs);
    gd.addChoice("Sort_pairs", SORT_OPTIONS, SORT_OPTIONS[sortIndex]);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption1 = gd.getNextChoice();
    inputOption2 = gd.getNextChoice();
    inputOption3 = gd.getNextChoice();
    dThreshold = gd.getNextNumber();
    beta = gd.getNextNumber();
    showPairs = gd.getNextBoolean();
    sortIndex = gd.getNextChoiceIndex();

    // Check arguments
    try {
      Parameters.isAboveZero("Distance threshold", dThreshold);
      Parameters.isPositive("Beta", beta);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  @SuppressWarnings("null")
  private void compareCoordinates(MemoryPeakResults results1, MemoryPeakResults results2,
      MemoryPeakResults results3, double dThreshold) {
    final Pulse[] p1 = extractPulses(results1);
    final Pulse[] p2 = extractPulses(results2);
    final Pulse[] p3 = extractPulses(results3);

    final List<Pulse> TP = null;
    List<Pulse> FP = null;
    List<Pulse> FN = null;
    List<PointPair> pairs = null;

    final List<Pulse> TP2 = null;
    List<Pulse> FP2 = null;
    List<Pulse> FN2 = null;
    List<PointPair> pairs2 = null;

    if (showPairs) {
      pairs = new LinkedList<>();
      FP = new LinkedList<>();
      FN = new LinkedList<>();
      pairs2 = new LinkedList<>();
      FP2 = new LinkedList<>();
      FN2 = new LinkedList<>();
    }

    final MatchResult result =
        MatchCalculator.analyseResults2D(p1, p2, dThreshold, TP, FP, FN, pairs);
    final MatchResult result2 =
        MatchCalculator.analyseResults2D(p1, p3, dThreshold, TP2, FP2, FN2, pairs2);

    // Create output
    if (!java.awt.GraphicsEnvironment.isHeadless()) {
      if (resultsWindow == null || !resultsWindow.isShowing()) {
        resultsWindow = new TextWindow(TITLE + " Results", createResultsHeader(), "", 900, 300);
      }
      if (showPairs) {
        if (p3 == null) {
          // Produce a pairs output
          if (pairsWindow == null || !pairsWindow.isShowing()) {
            pairsWindow = new TextWindow(TITLE + " Pairs", createPairsHeader(), "", 900, 300);
            final Point p = resultsWindow.getLocation();
            p.y += resultsWindow.getHeight();
            pairsWindow.setLocation(p);
            pairPainter = new ImageROIPainter(pairsWindow.getTextPanel(),
                results1.getSource().getOriginal().getName(), this);
          }
          pairsWindow.getTextPanel().clear();
          pairPainter.setTitle(results1.getSource().getOriginal().getName());

          // Add the unmatched points
          WindowManager.getIDList();

          for (final Coordinate c : FN) {
            pairs.add(new PointPair(c, null));
          }
          for (final Coordinate c : FP) {
            pairs.add(new PointPair(null, c));
          }

          final List<? extends PointPair> sortedPairs = sort(pairs);

          for (final PointPair pair : sortedPairs) {
            addPairResult(pair);
          }
        } else {
          // Produce a triple output
          if (triplesWindow == null || !triplesWindow.isShowing()) {
            triplesWindow = new TextWindow(TITLE + " Triples", createTriplesHeader(), "", 900, 300);
            final Point p = resultsWindow.getLocation();
            p.y += resultsWindow.getHeight();
            triplesWindow.setLocation(p);
            triplePainter = new ImageROIPainter(triplesWindow.getTextPanel(),
                results1.getSource().getName(), this);
          }
          triplesWindow.getTextPanel().clear();
          triplePainter.setTitle(results1.getSource().getOriginal().getName());

          final HashMap<Pulse, Triple> map = new HashMap<>();
          final ArrayList<Triple> triples = new ArrayList<>(pairs.size());
          for (final PointPair pair : pairs) {
            final Pulse p = (Pulse) pair.getPoint1();
            final Triple t = new Triple(p, (Pulse) pair.getPoint2(), null);
            triples.add(t);
            map.put(p, t);
          }
          // Complete the reference set of points
          for (final Coordinate c : FN) {
            final Pulse p = (Pulse) c;
            final Triple t = new Triple(p, null, null);
            triples.add(t);
            map.put(p, t);
          }

          // Add the unmatched points
          for (final Coordinate c : FP) {
            triples.add(new Triple(null, (Pulse) c, null));
          }
          for (final Coordinate c : FP2) {
            triples.add(new Triple(null, null, (Pulse) c));
          }

          // Add the results from the second match
          for (final PointPair pair : pairs2) {
            final Pulse p = (Pulse) pair.getPoint1();
            final Pulse pp = (Pulse) pair.getPoint2();
            final Triple t = map.get(p);
            if (t != null) {
              t.p3 = pp;
            } else {
              triples.add(new Triple(null, null, pp));
            }
          }

          final List<? extends Triple> sortedTriples = sort(triples);

          for (final Triple t : sortedTriples) {
            addTripleResult(t);
          }
        }
      }
    } else if (writeHeader) {
      writeHeader = false;
      IJ.log(createResultsHeader());
    }

    addResult(inputOption1, inputOption2, dThreshold, result);
    if (p3 != null) {
      addResult(inputOption1, inputOption3, dThreshold, result2);
    }
  }

  private static Pulse[] extractPulses(MemoryPeakResults results) {
    if (results == null) {
      return null;
    }
    final Pulse[] pulses = new Pulse[results.size()];
    final Counter i = new Counter();
    results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure() {
      @Override
      public void executeXYR(float x, float y, PeakResult p) {
        pulses[i.getAndIncrement()] = new Pulse(x, y, p.getFrame(), p.getEndFrame());
      }
    });
    return pulses;
  }

  private static String createResultsHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Image 1\t");
    sb.append("Image 2\t");
    sb.append("Distance (px)\t");
    sb.append("N\t");
    sb.append("TP\t");
    sb.append("FP\t");
    sb.append("FN\t");
    sb.append("Jaccard\t");
    sb.append("Score\t");
    sb.append("Precision\t");
    sb.append("Recall\t");
    sb.append("F0.5\t");
    sb.append("F1\t");
    sb.append("F2\t");
    sb.append("F-beta");
    return sb.toString();
  }

  private static void addResult(String i1, String i2, double dThrehsold, MatchResult result) {
    final StringBuilder sb = new StringBuilder();
    sb.append(i1).append('\t');
    sb.append(i2).append('\t');
    sb.append(IJ.d2s(dThrehsold, 2)).append('\t');
    sb.append(result.getNumberPredicted()).append('\t');
    sb.append(result.getTruePositives()).append('\t');
    sb.append(result.getFalsePositives()).append('\t');
    sb.append(result.getFalseNegatives()).append('\t');
    sb.append(IJ.d2s(result.getJaccard(), 4)).append('\t');
    sb.append(IJ.d2s(result.getRmsd(), 4)).append('\t');
    sb.append(IJ.d2s(result.getPrecision(), 4)).append('\t');
    sb.append(IJ.d2s(result.getRecall(), 4)).append('\t');
    sb.append(IJ.d2s(result.getFScore(0.5), 4)).append('\t');
    sb.append(IJ.d2s(result.getFScore(1.0), 4)).append('\t');
    sb.append(IJ.d2s(result.getFScore(2.0), 4)).append('\t');
    sb.append(IJ.d2s(result.getFScore(beta), 4));

    if (java.awt.GraphicsEnvironment.isHeadless()) {
      IJ.log(sb.toString());
    } else {
      resultsWindow.append(sb.toString());
    }
  }

  private static String createPairsHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Start1\t");
    sb.append("End2\t");
    sb.append("X1\t");
    sb.append("Y1\t");
    sb.append("Z1\t");
    sb.append("Start2\t");
    sb.append("End2\t");
    sb.append("X2\t");
    sb.append("Y2\t");
    sb.append("Z2\t");
    sb.append("Distance\t");
    sb.append("Score\t");
    return sb.toString();
  }

  private static void addPairResult(PointPair pair) {
    final StringBuilder sb = new StringBuilder();
    final Pulse p1 = (Pulse) pair.getPoint1();
    final Pulse p2 = (Pulse) pair.getPoint2();
    addPoint(sb, p1);
    addPointPairResult(sb, p1, p2);
    pairsWindow.append(sb.toString());
  }

  private static void addPoint(StringBuilder sb, Pulse p) {
    if (p == null) {
      sb.append("-\t-\t-\t-\t-\t");
    } else {
      sb.append(p.getStart()).append('\t');
      sb.append(p.getEnd()).append('\t');
      sb.append(IJ.d2s(p.getX())).append('\t');
      sb.append(IJ.d2s(p.getY())).append('\t');
      sb.append(IJ.d2s(p.getZ())).append('\t');
    }
  }

  private static String createTriplesHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Start1\t");
    sb.append("End2\t");
    sb.append("X1\t");
    sb.append("Y1\t");
    sb.append("Z1\t");
    sb.append("Start2\t");
    sb.append("End2\t");
    sb.append("X2\t");
    sb.append("Y2\t");
    sb.append("Z2\t");
    sb.append("Distance\t");
    sb.append("Score\t");
    sb.append("Start3\t");
    sb.append("End3\t");
    sb.append("X3\t");
    sb.append("Y3\t");
    sb.append("Z3\t");
    sb.append("Distance2\t");
    sb.append("Score2\t");
    return sb.toString();
  }

  @Override
  public double[] getCoordinates(String line) {
    // Extract the startT and x,y coordinates from the first pulse in the line
    final int[] index = {0, 5, 12};
    final String[] fields = line.split("\t");
    for (final int i : index) {
      if (i < fields.length) {
        if (fields[i].equals("-")) {
          continue;
        }
        final int startT = Integer.valueOf(fields[i]);
        final double x = Double.valueOf(fields[i + 2]);
        final double y = Double.valueOf(fields[i + 3]);
        return new double[] {startT, x, y};
      }
    }
    return null;
  }

  private static void addTripleResult(Triple triple) {
    final StringBuilder sb = new StringBuilder();
    final Pulse p1 = triple.p1;
    final Pulse p2 = triple.p2;
    final Pulse p3 = triple.p3;
    addPoint(sb, p1);
    addPointPairResult(sb, p1, p2);
    addPointPairResult(sb, p1, p3);
    triplesWindow.append(sb.toString());
  }

  private static void addPointPairResult(StringBuilder sb, Pulse p1, Pulse p2) {
    addPoint(sb, p2);
    final PointPair pair = new PointPair(p1, p2);
    final double d = pair.getXyDistance();
    if (d >= 0) {
      sb.append(MathUtils.rounded(d, 4)).append('\t');
    } else {
      sb.append("-\t");
    }
    if (p1 != null && p2 != null) {
      sb.append(MathUtils.rounded(p1.score(p2, d * d, dThreshold), 4)).append('\t');
    } else {
      sb.append("-\t");
    }
  }

  private List<? extends PointPair> sort(List<PointPair> pairs) {
    switch (sortIndex) {
      case 1: // Sort by time
        final ArrayList<TimeComparablePointPair> newPairs = new ArrayList<>(pairs.size());
        for (final PointPair pair : pairs) {
          newPairs.add(new TimeComparablePointPair(pair));
        }
        Collections.sort(newPairs);
        return newPairs;

      default:
        // Already sorted by score
        return pairs;
    }
  }

  private List<? extends Triple> sort(ArrayList<Triple> triples) {
    if (sortIndex == 1) {
      final List<TimeComparableTriple> sorted = new ArrayList<>(triples.size());
      for (final Triple t : triples) {
        sorted.add(new TimeComparableTriple(t));
      }
      Collections.sort(sorted);
      return sorted;
    }
    final List<ScoreComparableTriple> sorted = new ArrayList<>(triples.size());
    for (final Triple t : triples) {
      sorted.add(new ScoreComparableTriple(t));
    }
    Collections.sort(sorted);
    return sorted;
  }

  private class TimeComparablePointPair extends PointPair
      implements Comparable<TimeComparablePointPair> {
    int startT = Integer.MAX_VALUE;
    public Pulse p1;
    public Pulse p2;

    public TimeComparablePointPair(PointPair pair) {
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

    @Override
    public int compareTo(TimeComparablePointPair o) {
      // Significance of points are: p1, p2, p3
      // Sort by the earliest start time of the most significant point.
      if (startT == o.startT) {
        // Sort using the significant points first.
        int result = comparePulsesStartT(p1, o.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesStartT(p2, o.p2);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(p1, o.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(p2, o.p2);
        if (result != 0) {
          return result;
        }

        // Sort using coords
        result = comparePulsesCoords(p1, o.p1);
        if (result != 0) {
          return result;
        }
        return comparePulsesCoords(p2, o.p2);
      }
      return startT - o.startT;
    }
  }

  private static int comparePulsesStartT(Pulse p1, Pulse p2) {
    // Data for a point always beats no data
    if (p1 == null) {
      return (p2 == null) ? 0 : 1;
    }
    if (p2 == null) {
      return -1;
    }
    return p1.getStart() - p2.getStart();
  }

  private static int comparePulsesEndT(Pulse p1, Pulse p2) {
    // Data for a point always beats no data
    if (p1 == null) {
      return (p2 == null) ? 0 : 1;
    }
    if (p2 == null) {
      return -1;
    }
    return p1.getEnd() - p2.getEnd();
  }

  private static int comparePulsesCoords(Pulse p1, Pulse p2) {
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

  private class Triple {
    public Pulse p1;
    public Pulse p2;
    public Pulse p3;

    public Triple(Pulse p1, Pulse p2, Pulse p3) {
      this.p1 = p1;
      this.p2 = p2;
      this.p3 = p3;
    }

    public Triple(Triple t) {
      this.p1 = t.p1;
      this.p2 = t.p2;
      this.p3 = t.p3;
    }
  }

  private class TimeComparableTriple extends Triple implements Comparable<TimeComparableTriple> {
    int startT = Integer.MAX_VALUE;

    public TimeComparableTriple(Triple t) {
      super(t);

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

    @Override
    public int compareTo(TimeComparableTriple o) {
      // Significance of points are: p1, p2, p3
      // Sort by the earliest start time of the most significant point.
      if (startT == o.startT) {
        // Sort using the significant points first.
        // int result = comparePulsesStartT(p1, o.p1);
        // if (result != 0)
        // return result;
        // result = comparePulsesStartT(p2, o.p2);
        // if (result != 0)
        // return result;
        // result = comparePulsesStartT(p3, o.p3);
        // if (result != 0)
        // return result;
        // result = comparePulsesEndT(p1, o.p1);
        // if (result != 0)
        // return result;
        // result = comparePulsesEndT(p2, o.p2);
        // if (result != 0)
        // return result;
        // return comparePulsesEndT(p3, o.p3);

        // Make the order the same as for the pairs table
        int result = comparePulsesStartT(p1, o.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesStartT(p2, o.p2);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(p1, o.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(p2, o.p2);
        if (result != 0) {
          return result;
        }

        // Then sort using the third
        result = comparePulsesStartT(p3, o.p3);
        if (result != 0) {
          return result;
        }
        result = comparePulsesEndT(p3, o.p3);
        if (result != 0) {
          return result;
        }

        // Sort using coords
        result = comparePulsesCoords(p1, o.p1);
        if (result != 0) {
          return result;
        }
        result = comparePulsesCoords(p2, o.p2);
        if (result != 0) {
          return result;
        }
        return comparePulsesCoords(p3, o.p3);
      }
      return startT - o.startT;
    }
  }

  private class ScoreComparableTriple extends Triple implements Comparable<ScoreComparableTriple> {
    double score;

    public ScoreComparableTriple(Triple t) {
      super(t);
      if (p1 != null) {
        if (p2 != null) {
          score = p1.score(p2, dThreshold);
        }
        if (p3 != null) {
          score = FastMath.max(score, p1.score(p2, dThreshold));
        }
      }
    }

    @Override
    public int compareTo(ScoreComparableTriple o) {
      if (score > o.score) {
        return -1;
      }
      if (score < o.score) {
        return 1;
      }
      return 0;
    }
  }
}
