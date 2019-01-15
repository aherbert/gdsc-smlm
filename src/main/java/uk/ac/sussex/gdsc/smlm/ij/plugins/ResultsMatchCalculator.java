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

import uk.ac.sussex.gdsc.core.data.utils.Rounder;
import uk.ac.sussex.gdsc.core.data.utils.RounderUtils;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.match.BasePoint;
import uk.ac.sussex.gdsc.core.match.Coordinate;
import uk.ac.sussex.gdsc.core.match.MatchCalculator;
import uk.ac.sussex.gdsc.core.match.MatchResult;
import uk.ac.sussex.gdsc.core.match.PointPair;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageROIPainter;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.set.hash.TIntHashSet;

import ij.IJ;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Compares the coordinates in two sets of results and computes the match statistics.
 */
public class ResultsMatchCalculator implements PlugIn {
  private static String TITLE = "Results Match Calculator";

  private static String inputOption1 = "";
  private static String inputOption2 = "";
  private static double distanceThreshold = 0.5;
  private static int increments = 5;
  private static double delta = 0.1;
  private static double beta = 4;
  private static boolean showTable = true;
  private static boolean showPairs;
  private static boolean saveClassifications;
  private static String classificationsFile = "";
  private static boolean idAnalysis;

  private static boolean writeHeader = true;
  private static TextWindow resultsWindow;
  private static TextWindow pairsWindow;
  private static ImageROIPainter pairPainter;

  private final Rounder rounder = RounderUtils.create(4);

  /**
   * A point that holds a reference to a PeakResult.
   */
  public static class PeakResultPoint extends BasePoint {
    /** The time. */
    int time;

    /** The peak result. */
    PeakResult peakResult;

    /**
     * Instantiates a new peak result point.
     *
     * @param time the time
     * @param x the x
     * @param y the y
     * @param z the z
     * @param peakResult the peak result
     */
    public PeakResultPoint(int time, float x, float y, float z, PeakResult peakResult) {
      super(x, y, z);
      this.time = time;
      this.peakResult = peakResult;
    }

    /**
     * Gets the time.
     *
     * @return the time
     */
    public int getTime() {
      return time;
    }

    /**
     * Gets the peak result.
     *
     * @return the peak result
     */
    public PeakResult getPeakResult() {
      return peakResult;
    }
  }

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
      IJ.error(TITLE, "Distance unit should be the same for the results");
      return;
    }

    final long start = System.nanoTime();
    runCompareCoordinates(results1, results2, distanceThreshold, increments, delta);
    final double seconds = (System.nanoTime() - start) / 1000000000.0;

    IJ.showStatus(String.format("%s = %ss", TITLE, MathUtils.rounded(seconds, 4)));
  }

  private static boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addMessage("Compare the points in two results sets\nand compute the match statistics");
    ResultsManager.addInput(gd, "Results1", inputOption1, InputSource.MEMORY);
    ResultsManager.addInput(gd, "Results2", inputOption2, InputSource.MEMORY);
    gd.addNumericField("Distance", distanceThreshold, 2);

    gd.addSlider("Increments", 0, 10, increments);
    gd.addNumericField("Delta", delta, 2);
    gd.addNumericField("Beta", beta, 2);
    gd.addCheckbox("Show_table", showTable);
    gd.addCheckbox("Show_pairs", showPairs);
    gd.addCheckbox("Save_classifications", saveClassifications);
    gd.addCheckbox("Id_analysis", idAnalysis);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption1 = gd.getNextChoice();
    inputOption2 = gd.getNextChoice();
    distanceThreshold = gd.getNextNumber();
    increments = (int) gd.getNextNumber();
    delta = gd.getNextNumber();
    beta = gd.getNextNumber();
    showTable = gd.getNextBoolean();
    showPairs = gd.getNextBoolean();
    saveClassifications = gd.getNextBoolean();
    idAnalysis = gd.getNextBoolean();

    if (!(showTable || showPairs || saveClassifications)) {
      IJ.error(TITLE, "No outputs specified");
      return false;
    }

    // Check arguments
    try {
      Parameters.isPositive("Distance threshold", distanceThreshold);
      Parameters.isPositive("Increments", increments);
      Parameters.isAboveZero("Delta", delta);
      Parameters.isPositive("Beta", beta);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  @SuppressWarnings("null")
  private void runCompareCoordinates(MemoryPeakResults results1, MemoryPeakResults results2,
      double maxDistanceThreshold, int increments, double delta) {
    final boolean requirePairs = showPairs || saveClassifications;

    final TextFilePeakResults fileResults = createFilePeakResults(results2);

    final List<PointPair> allMatches = new LinkedList<>();
    final List<PointPair> pairs = (requirePairs) ? new LinkedList<>() : null;

    final double maxDistance = maxDistanceThreshold + increments * delta;

    // Divide the results into time points
    final TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates = getCoordinates(results1);
    final TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates = getCoordinates(results2);

    int n1 = 0;
    int n2 = 0;

    // Process each time point
    for (final Integer t : getTimepoints(actualCoordinates, predictedCoordinates)) {
      final Coordinate[] actual = getCoordinates(actualCoordinates, t);
      final Coordinate[] predicted = getCoordinates(predictedCoordinates, t);

      final List<Coordinate> tp = null;
      List<Coordinate> fp = null;
      List<Coordinate> fn = null;
      final List<PointPair> matches = new LinkedList<>();
      if (requirePairs) {
        fp = new LinkedList<>();
        fn = new LinkedList<>();
      }

      MatchCalculator.analyseResults2D(actual, predicted, maxDistance, tp, fp, fn, matches);

      // Aggregate
      n1 += actual.length;
      n2 += predicted.length;

      allMatches.addAll(matches);
      if (showPairs) {
        pairs.addAll(matches);
        for (final Coordinate c : fn) {
          pairs.add(new PointPair(c, null));
        }
        for (final Coordinate c : fp) {
          pairs.add(new PointPair(null, c));
        }
      }
      if (fileResults != null) {
        // Matches are marked in the original value with 1 for true, 0 for false
        for (final PointPair pair : matches) {
          PeakResult result = ((PeakResultPoint) pair.getPoint2()).peakResult;
          result = result.copy();
          result.setOrigValue(1);
          fileResults.add(result);
        }
        for (final Coordinate c : fp) {
          PeakResult result = ((PeakResultPoint) c).peakResult;
          result = result.copy();
          result.setOrigValue(0);
          fileResults.add(result);
        }
      }
    }

    if (fileResults != null) {
      fileResults.end();
    }

    // XXX : DEBUGGING : Output for signal correlation and fitting analysis
    /*
     * try { OutputStreamWriter o = new OutputStreamWriter(new
     * FileOutputStream("/tmp/ResultsMatchCalculator.txt")); FilePeakResults r1 = new
     * FilePeakResults("/tmp/" + results1.getName() + ".txt", false); FilePeakResults r2 = new
     * FilePeakResults("/tmp/" + results2.getName() + ".txt", false); r1.begin(); r2.begin();
     * //OutputStreamWriter o2 = new OutputStreamWriter(new
     * FileOutputStream("/tmp/"+results1.getName()+".txt")); //OutputStreamWriter o3 = new
     * OutputStreamWriter(new FileOutputStream("/tmp/"+results2.getName()+".txt")); for (PointPair
     * pair : allMatches) { PeakResult p1 = ((PeakResultPoint) pair.getPoint1()).peakResult;
     * PeakResult p2 = ((PeakResultPoint) pair.getPoint2()).peakResult; r1.add(p1); r2.add(p2);
     * o.write(Float.toString(p1.getSignal())); o.write('\t');
     * o.write(Float.toString(p2.getSignal())); o.write('\n'); } o.close(); r1.end(); r2.end(); }
     * catch (Exception ex) { e.printStackTrace(); }
     */

    final boolean doIdAnalysis1 = (idAnalysis) ? haveIds(results1) : false;
    final boolean doIdAnalysis2 = (idAnalysis) ? haveIds(results2) : false;
    final boolean doIdAnalysis = doIdAnalysis1 || doIdAnalysis2;

    // Create output
    if (!java.awt.GraphicsEnvironment.isHeadless()) {
      final String header = createResultsHeader(doIdAnalysis);
      ImageJUtils.refreshHeadings(resultsWindow, header, true);

      if (showTable && (resultsWindow == null || !resultsWindow.isShowing())) {
        resultsWindow = new TextWindow(TITLE + " Results", header, "", 900, 300);
      }
      if (showPairs) {
        if (pairsWindow == null || !pairsWindow.isShowing()) {
          pairsWindow = new TextWindow(TITLE + " Pairs", createPairsHeader(), "", 900, 300);
          if (resultsWindow != null) {
            final Point p = resultsWindow.getLocation();
            p.y += resultsWindow.getHeight();
            pairsWindow.setLocation(p);
          }
          pairPainter = new ImageROIPainter(pairsWindow.getTextPanel(), "", line -> {
            // Extract the startT and x,y coordinates from the first pulse in the line
            final int[] index = {1, 4};
            final String[] fields = line.split("\t");
            final int startT = Integer.parseInt(fields[0]);
            for (final int i : index) {
              if (i < fields.length) {
                if (fields[i].equals("-")) {
                  continue;
                }
                final double x = Double.parseDouble(fields[i]);
                final double y = Double.parseDouble(fields[i + 1]);
                return new double[] {startT, x, y};
              }
            }
            return null;
          });
        }
        pairsWindow.getTextPanel().clear();
        String title = "Results 1";
        if (results1.getSource() != null
            && results1.getSource().getOriginal().getName().length() > 0) {
          title = results1.getSource().getOriginal().getName();
        }
        pairPainter.setTitle(title);
        IJ.showStatus("Writing pairs table");
        IJ.showProgress(0);
        int count = 0;
        final int total = pairs.size();
        final int step = ImageJUtils.getProgressInterval(total);
        final ArrayList<String> list = new ArrayList<>(total);
        boolean flush = true;
        for (final PointPair pair : pairs) {

          if (++count % step == 0) {
            IJ.showProgress(count, total);
          }
          list.add(addPairResult(pair));
          if (flush && count == 9) {
            pairsWindow.getTextPanel().append(list);
            list.clear();
            flush = false;
          }
        }
        pairsWindow.getTextPanel().append(list);
        IJ.showProgress(1);
      }
    } else if (writeHeader && showTable) {
      writeHeader = false;
      IJ.log(createResultsHeader(idAnalysis));
    }

    if (!showTable) {
      return;
    }

    // We have the results for the largest distance.
    // Now reduce the distance threshold and recalculate the results
    final double[] distanceThresholds = getDistances(maxDistanceThreshold, increments, delta);
    final double[] pairDistances = getPairDistances(allMatches);
    // Re-use storage for the ID analysis
    TIntHashSet id1 = null;
    TIntHashSet id2 = null;
    TIntHashSet matchId1 = null;
    TIntHashSet matchId2 = null;
    if (doIdAnalysis) {
      if (doIdAnalysis1) {
        id1 = getIds(results1);
        matchId1 = new TIntHashSet(id1.size());
      }
      if (doIdAnalysis2) {
        id2 = getIds(results2);
        matchId2 = new TIntHashSet(id2.size());
      }
    }
    for (final double distanceThreshold : distanceThresholds) {
      double rms = 0;
      int tp2 = 0;
      final double d2 = distanceThreshold * distanceThreshold;
      for (final double d : pairDistances) {
        if (d <= d2) {
          rms += d;
          tp2++;
        }
      }
      // All non-true positives must be added to the false totals.
      final int fp2 = n2 - tp2;
      final int fn2 = n1 - tp2;

      final MatchResult result =
          new MatchResult(tp2, fp2, fn2, (tp2 > 0) ? Math.sqrt(rms / tp2) : 0);

      MatchResult idResult1 = null;
      MatchResult idResult2 = null;
      if (doIdAnalysis) {
        if (doIdAnalysis1) {
          matchId1.clear();
        }
        if (doIdAnalysis2) {
          matchId2.clear();
        }
        int index = 0;
        for (final PointPair pair : allMatches) {
          if (pairDistances[index++] <= d2) {
            if (doIdAnalysis1) {
              matchId1.add(((PeakResultPoint) pair.getPoint1()).peakResult.getId());
            }
            if (doIdAnalysis2) {
              matchId2.add(((PeakResultPoint) pair.getPoint2()).peakResult.getId());
            }
          }
        }
        // Only the actual points are checked for Ids. For example these could be from the
        // Create Data plugin with actual fluorophore Ids.
        // => Only the recall will be valid: tp / (tp + fn)
        if (doIdAnalysis1) {
          idResult1 = new MatchResult(matchId1.size(), 0, id1.size() - matchId1.size(), 0);
        }
        if (doIdAnalysis2) {
          idResult2 = new MatchResult(matchId2.size(), 0, id2.size() - matchId2.size(), 0);
        }
      }

      addResult(inputOption1, inputOption2, distanceThreshold, result, idResult1, idResult2);
    }
  }

  @SuppressWarnings("unused")
  private static boolean haveIds(MemoryPeakResults results1, MemoryPeakResults results2) {
    return haveIds(results1) && haveIds(results2);
  }

  private static boolean haveIds(MemoryPeakResults results) {
    return results.hasId();
  }

  private static TextFilePeakResults createFilePeakResults(MemoryPeakResults results2) {
    if (!saveClassifications) {
      return null;
    }
    final String[] path = ImageJUtils.decodePath(classificationsFile);
    final OpenDialog chooser = new OpenDialog("Classifications_File", path[0], path[1]);
    if (chooser.getFileName() != null) {
      classificationsFile = chooser.getDirectory() + chooser.getFileName();
      final TextFilePeakResults r = new TextFilePeakResults(classificationsFile, false, false);
      r.copySettings(results2);
      r.begin();
      return r;
    }
    return null;
  }

  /**
   * Build a map between the peak id (time point) and a list of coordinates.
   *
   * @param results the results
   * @return the coordinates
   */
  public static TIntObjectHashMap<ArrayList<Coordinate>> getCoordinates(MemoryPeakResults results) {
    return getCoordinates(results, false);
  }

  /**
   * Build a map between the peak id (time point) and a list of coordinates.
   *
   * @param results the results
   * @param integerCoordinates True if the values should be rounded down to integers
   * @return the coordinates
   */
  public static TIntObjectHashMap<ArrayList<Coordinate>> getCoordinates(MemoryPeakResults results,
      final boolean integerCoordinates) {
    final TIntObjectHashMap<ArrayList<Coordinate>> coords = new TIntObjectHashMap<>();
    if (results.size() > 0) {
      // Do not use HashMap directly to build the coords object since there
      // will be many calls to getEntry(). Instead sort the results and use
      // a new list for each time point
      results.sort();
      final int minT = results.getFirstFrame();
      final int maxT = results.getLastFrame();

      // Create lists
      final ArrayList<ArrayList<Coordinate>> tmpCoords = new ArrayList<>(maxT - minT + 1);
      for (int t = minT; t <= maxT; t++) {
        tmpCoords.add(new ArrayList<Coordinate>());
      }

      // Add the results to the lists
      results.forEach((PeakResultProcedure) result -> {
        final float x;
        final float y;
        final float z;
        if (integerCoordinates) {
          x = (int) result.getXPosition();
          y = (int) result.getYPosition();
          z = (int) result.getZPosition();
        } else {
          x = result.getXPosition();
          y = result.getYPosition();
          z = result.getZPosition();
        }
        for (int t = result.getFrame() - minT, i = result.getEndFrame() - result.getFrame() + 1;
            i-- > 0; t++) {
          tmpCoords.get(t).add(new PeakResultPoint(t + minT, x, y, z, result));
        }
      });

      // Put in the map
      for (int t = minT, i = 0; t <= maxT; t++, i++) {
        coords.put(t, tmpCoords.get(i));
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
  public static Coordinate[] getCoordinates(TIntObjectHashMap<ArrayList<Coordinate>> coords,
      int time) {
    final ArrayList<Coordinate> tmp = coords.get(time);
    if (tmp != null) {
      return tmp.toArray(new Coordinate[tmp.size()]);
    }
    return new Coordinate[0];
  }

  /**
   * Merge the time points from each map into a single sorted list of unique time points.
   *
   * @param actualCoordinates the actual coordinates
   * @param predictedCoordinates the predicted coordinates
   * @return a list of time points
   */
  private static int[] getTimepoints(TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates,
      TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates) {

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

  /**
   * Merge the time points from the two sets into a single sorted list of unique time points.
   *
   * @param actualPoints the actual points
   * @param predictedPoints the predicted points
   * @return the timepoints
   */
  @SuppressWarnings("unused")
  private static int[] getTimepoints(List<PeakResult> actualPoints,
      List<PeakResult> predictedPoints) {
    final TIntHashSet set = new TIntHashSet();
    for (final PeakResult r : actualPoints) {
      set.add(r.getFrame());
    }
    for (final PeakResult r : predictedPoints) {
      set.add(r.getFrame());
    }
    final int[] t = set.toArray();
    Arrays.sort(t);
    return t;
  }

  private static String createResultsHeader(boolean idAnalysis) {
    final StringBuilder sb = new StringBuilder();
    sb.append("Image 1\t");
    sb.append("Image 2\t");
    sb.append("Distance (px)\t");
    sb.append("N\t");
    sb.append("TP\t");
    sb.append("FP\t");
    sb.append("FN\t");
    sb.append("Jaccard\t");
    sb.append("RMSD\t");
    sb.append("Precision\t");
    sb.append("Recall\t");
    sb.append("F0.5\t");
    sb.append("F1\t");
    sb.append("F2\t");
    sb.append("F-beta");
    if (idAnalysis) {
      sb.append("\tId1-N");
      sb.append("\tId1-TP");
      sb.append("\tId1-Recall");
      sb.append("\tId2-N");
      sb.append("\tId2-TP");
      sb.append("\tId2-Recall");
    }
    return sb.toString();
  }

  private void addResult(String i1, String i2, double distanceThrehsold, MatchResult result,
      MatchResult idResult1, MatchResult idResult2) {
    final StringBuilder sb = new StringBuilder();
    sb.append(i1).append('\t');
    sb.append(i2).append('\t');
    sb.append(rounder.round(distanceThrehsold)).append('\t');
    sb.append(result.getNumberPredicted()).append('\t');
    sb.append(result.getTruePositives()).append('\t');
    sb.append(result.getFalsePositives()).append('\t');
    sb.append(result.getFalseNegatives()).append('\t');
    sb.append(rounder.round(result.getJaccard())).append('\t');
    sb.append(rounder.round(result.getRmsd())).append('\t');
    sb.append(rounder.round(result.getPrecision())).append('\t');
    sb.append(rounder.round(result.getRecall())).append('\t');
    sb.append(rounder.round(result.getFScore(0.5))).append('\t');
    sb.append(rounder.round(result.getFScore(1.0))).append('\t');
    sb.append(rounder.round(result.getFScore(2.0))).append('\t');
    sb.append(rounder.round(result.getFScore(beta)));
    if (idResult1 != null) {
      sb.append('\t').append(idResult1.getNumberPredicted());
      sb.append('\t').append(idResult1.getTruePositives());
      sb.append('\t').append(rounder.round(idResult1.getRecall()));
    } else if (idResult2 != null) {
      sb.append("\t-\t-\t-");
    }
    if (idResult2 != null) {
      sb.append('\t').append(idResult2.getNumberPredicted());
      sb.append('\t').append(idResult2.getTruePositives());
      sb.append('\t').append(rounder.round(idResult2.getRecall()));
    } else if (idResult1 != null) {
      sb.append("\t-\t-\t-");
    }

    if (java.awt.GraphicsEnvironment.isHeadless()) {
      IJ.log(sb.toString());
    } else {
      resultsWindow.append(sb.toString());
    }
  }

  private static String createPairsHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("T\t");
    sb.append("X1\t");
    sb.append("Y1\t");
    sb.append("Z1\t");
    sb.append("X2\t");
    sb.append("Y2\t");
    sb.append("Z2\t");
    sb.append("Distance\t");
    return sb.toString();
  }

  private String addPairResult(PointPair pair) {
    final StringBuilder sb = new StringBuilder();
    final PeakResultPoint p1 = (PeakResultPoint) pair.getPoint1();
    final PeakResultPoint p2 = (PeakResultPoint) pair.getPoint2();
    final int t = (p1 != null) ? p1.getTime() : p2.getTime();
    sb.append(t).append('\t');
    addPoint(sb, p1);
    addPoint(sb, p2);
    final double d = pair.getXyDistance();
    if (d >= 0) {
      sb.append(rounder.round(d)).append('\t');
    } else {
      sb.append("-\t");
    }
    return sb.toString();
  }

  private void addPoint(StringBuilder sb, PeakResultPoint result) {
    if (result == null) {
      sb.append("-\t-\t-\t");
    } else {
      sb.append(rounder.round(result.getX())).append('\t');
      sb.append(rounder.round(result.getY())).append('\t');
      sb.append(rounder.round(result.getZ())).append('\t');
    }
  }

  private static TIntHashSet getIds(MemoryPeakResults results) {
    final TIntHashSet ids = new TIntHashSet(results.size());
    results.forEach((PeakResultProcedure) result -> ids.add(result.getId()));
    return ids;
  }

  private static double[] getDistances(double distanceThreshold, int increments, double delta) {
    final double[] d = new double[increments + 1];
    for (int i = 0; i <= increments; i++) {
      d[i] = distanceThreshold + i * delta;
    }
    return d;
  }

  private static double[] getPairDistances(List<PointPair> pairs) {
    final double[] d = new double[pairs.size()];
    int index = 0;
    for (final PointPair pair : pairs) {
      d[index++] = pair.getXyDistanceSquared();
    }
    return d;
  }

  /**
   * Compare the coordinates in two results sets.
   *
   * @param results1 the results 1
   * @param results2 the results 2
   * @param distance the distance
   * @return the match result
   */
  public static MatchResult compareCoordinates(MemoryPeakResults results1,
      MemoryPeakResults results2, double distance) {
    // Divide the results into time points
    final TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates = getCoordinates(results1);
    final TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates = getCoordinates(results2);

    return compareCoordinates(actualCoordinates, predictedCoordinates, distance);
  }

  /**
   * Compare the coordinates on a frame-by-frame basis.
   *
   * @param actualCoordinates the actual coordinates
   * @param predictedCoordinates the predicted coordinates
   * @param distance the distance
   * @return the match result
   */
  public static MatchResult compareCoordinates(
      TIntObjectHashMap<ArrayList<Coordinate>> actualCoordinates,
      TIntObjectHashMap<ArrayList<Coordinate>> predictedCoordinates, double distance) {
    int tp = 0;
    int fp = 0;
    int fn = 0;

    // Process each time point
    for (final Integer t : getTimepoints(actualCoordinates, predictedCoordinates)) {
      final Coordinate[] actual = getCoordinates(actualCoordinates, t);
      final Coordinate[] predicted = getCoordinates(predictedCoordinates, t);

      final MatchResult r = MatchCalculator.analyseResults2D(actual, predicted, distance);

      // Aggregate
      tp += r.getTruePositives();
      fp += r.getFalsePositives();
      fn += r.getFalseNegatives();
    }

    return new MatchResult(tp, fp, fn, 0);
  }
}
