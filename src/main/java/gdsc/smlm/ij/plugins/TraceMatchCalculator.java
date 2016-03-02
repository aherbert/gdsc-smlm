package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC Plugins for ImageJ
 * 
 * Copyright (C) 2011 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.utils.CoordinateProvider;
import gdsc.smlm.ij.utils.ImageROIPainter;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.match.Coordinate;
import gdsc.smlm.results.match.MatchCalculator;
import gdsc.smlm.results.match.MatchResult;
import gdsc.smlm.results.match.PointPair;
import gdsc.smlm.results.match.Pulse;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

/**
 * Compares the coordinates in sets of traced results and computes the match statistics.
 */
public class TraceMatchCalculator implements PlugIn, CoordinateProvider
{
	private static String TITLE = "Trace Match Calculator";

	private static String inputOption1 = "";
	private static String inputOption2 = "";
	private static String inputOption3 = "";
	private static double dThreshold = 1;
	private static double beta = 4;
	private static boolean showPairs = false;
	private static String[] SORT_OPTIONS = new String[] { "Score", "Time" };
	private static int sortIndex = 1;

	private static boolean writeHeader = true;
	private static TextWindow resultsWindow = null;
	private static TextWindow pairsWindow = null;
	private static TextWindow triplesWindow = null;
	private static ImageROIPainter pairPainter = null, triplePainter = null;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
		
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		if (!showDialog())
			return;

		// Load the results
		MemoryPeakResults results1 = ResultsManager.loadInputResults(inputOption1, false);
		MemoryPeakResults results2 = ResultsManager.loadInputResults(inputOption2, false);
		MemoryPeakResults results3 = ResultsManager.loadInputResults(inputOption3, false);
		if (results1 == null || results1.size() == 0 || results2 == null || results2.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		final long start = System.nanoTime();
		compareCoordinates(results1, results2, results3, dThreshold);
		double seconds = (System.nanoTime() - start) / 1000000000.0;

		IJ.showStatus(String.format("%s = %ss", TITLE, Utils.rounded(seconds, 4)));
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		gd.addMessage("Compare the points in two results sets\nand compute the match statistics");
		ResultsManager.addInput(gd, "Results1", inputOption1, InputSource.MEMORY_MULTI_FRAME);
		ResultsManager.addInput(gd, "Results2", inputOption2, InputSource.MEMORY_MULTI_FRAME);
		ResultsManager.addInput(gd, "Results3", inputOption3, InputSource.NONE, InputSource.MEMORY_MULTI_FRAME);
		gd.addNumericField("Distance", dThreshold, 2);

		gd.addNumericField("Beta", beta, 2);
		gd.addCheckbox("Show_pairs", showPairs);
		gd.addChoice("Sort_pairs", SORT_OPTIONS, SORT_OPTIONS[sortIndex]);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption1 = gd.getNextChoice();
		inputOption2 = gd.getNextChoice();
		inputOption3 = gd.getNextChoice();
		dThreshold = gd.getNextNumber();
		beta = gd.getNextNumber();
		showPairs = gd.getNextBoolean();
		sortIndex = gd.getNextChoiceIndex();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Distance threshold", dThreshold);
			Parameters.isPositive("Beta", beta);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}
		
		return true;
	}

	private void compareCoordinates(MemoryPeakResults results1, MemoryPeakResults results2, MemoryPeakResults results3,
			double dThreshold)
	{
		Pulse[] p1 = extractPulses(results1);
		Pulse[] p2 = extractPulses(results2);
		Pulse[] p3 = extractPulses(results3);

		List<Pulse> TP = null;
		List<Pulse> FP = null;
		List<Pulse> FN = null;
		List<PointPair> pairs = null;

		List<Pulse> TP2 = null;
		List<Pulse> FP2 = null;
		List<Pulse> FN2 = null;
		List<PointPair> pairs2 = null;

		if (showPairs)
		{
			pairs = new LinkedList<PointPair>();
			FP = new LinkedList<Pulse>();
			FN = new LinkedList<Pulse>();
			pairs2 = new LinkedList<PointPair>();
			FP2 = new LinkedList<Pulse>();
			FN2 = new LinkedList<Pulse>();
		}

		MatchResult result = MatchCalculator.analyseResults2D(p1, p2, dThreshold, TP, FP, FN, pairs);
		MatchResult result2 = MatchCalculator.analyseResults2D(p1, p3, dThreshold, TP2, FP2, FN2, pairs2);

		// Create output
		if (!java.awt.GraphicsEnvironment.isHeadless())
		{
			if (resultsWindow == null || !resultsWindow.isShowing())
			{
				resultsWindow = new TextWindow(TITLE + " Results", createResultsHeader(), "", 900, 300);
			}
			if (showPairs)
			{
				if (p3 == null)
				{
					// Produce a pairs output
					if (pairsWindow == null || !pairsWindow.isShowing())
					{
						pairsWindow = new TextWindow(TITLE + " Pairs", createPairsHeader(), "", 900, 300);
						Point p = resultsWindow.getLocation();
						p.y += resultsWindow.getHeight();
						pairsWindow.setLocation(p);
						pairPainter = new ImageROIPainter(pairsWindow.getTextPanel(), results1.getSource().getOriginal().getName(),
								this);
					}
					pairsWindow.getTextPanel().clear();
					pairPainter.setTitle(results1.getSource().getOriginal().getName());

					// Add the unmatched points
					for (Coordinate c : FN)
						pairs.add(new PointPair(c, null));
					for (Coordinate c : FP)
						pairs.add(new PointPair(null, c));

					List<? extends PointPair> sortedPairs = sort(pairs);

					for (PointPair pair : sortedPairs)
						addPairResult(pair);
				}
				else
				{
					// Produce a triple output
					if (triplesWindow == null || !triplesWindow.isShowing())
					{
						triplesWindow = new TextWindow(TITLE + " Triples", createTriplesHeader(), "", 900, 300);
						Point p = resultsWindow.getLocation();
						p.y += resultsWindow.getHeight();
						triplesWindow.setLocation(p);
						triplePainter = new ImageROIPainter(triplesWindow.getTextPanel(), results1.getSource()
								.getName(), this);
					}
					triplesWindow.getTextPanel().clear();
					triplePainter.setTitle(results1.getSource().getOriginal().getName());

					HashMap<Pulse, Triple> map = new HashMap<Pulse, Triple>();
					ArrayList<Triple> triples = new ArrayList<Triple>(pairs.size());
					for (PointPair pair : pairs)
					{
						Pulse p = (Pulse) pair.getPoint1();
						Triple t = new Triple(p, (Pulse) pair.getPoint2(), null);
						triples.add(t);
						map.put(p, t);
					}
					// Complete the reference set of points
					for (Coordinate c : FN)
					{
						Pulse p = (Pulse) c;
						Triple t = new Triple(p, null, null);
						triples.add(t);
						map.put(p, t);
					}

					// Add the unmatched points
					for (Coordinate c : FP)
						triples.add(new Triple(null, (Pulse) c, null));
					for (Coordinate c : FP2)
						triples.add(new Triple(null, null, (Pulse) c));

					// Add the results from the second match
					for (PointPair pair : pairs2)
					{
						Pulse p = (Pulse) pair.getPoint1();
						Pulse pp = (Pulse) pair.getPoint2();
						if (map.containsKey(p))
						{
							Triple t = map.get(p);
							t.p3 = pp;
						}
						else
						{
							triples.add(new Triple(null, null, pp));
						}
					}

					List<? extends Triple> sortedTriples = sort(triples);

					for (Triple t : sortedTriples)
						addTripleResult(t);
				}
			}
		}
		else
		{
			if (writeHeader)
			{
				writeHeader = false;
				IJ.log(createResultsHeader());
			}
		}

		addResult(inputOption1, inputOption2, dThreshold, result);
		if (p3 != null)
			addResult(inputOption1, inputOption3, dThreshold, result2);
	}

	private Pulse[] extractPulses(MemoryPeakResults results)
	{
		if (results == null)
			return null;
		results.getResults();
		Pulse[] pulses = new Pulse[results.size()];
		int i = 0;
		for (PeakResult p : results.getResults())
		{
			pulses[i++] = new Pulse(p.getXPosition(), p.getYPosition(), p.peak, p.getEndFrame());
		}
		return pulses;
	}

	private String createResultsHeader()
	{
		StringBuilder sb = new StringBuilder();
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

	private void addResult(String i1, String i2, double dThrehsold, MatchResult result)
	{
		StringBuilder sb = new StringBuilder();
		sb.append(i1).append("\t");
		sb.append(i2).append("\t");
		sb.append(IJ.d2s(dThrehsold, 2)).append("\t");
		sb.append(result.getNumberPredicted()).append("\t");
		sb.append(result.getTruePositives()).append("\t");
		sb.append(result.getFalsePositives()).append("\t");
		sb.append(result.getFalseNegatives()).append("\t");
		sb.append(IJ.d2s(result.getJaccard(), 4)).append("\t");
		sb.append(IJ.d2s(result.getRMSD(), 4)).append("\t");
		sb.append(IJ.d2s(result.getPrecision(), 4)).append("\t");
		sb.append(IJ.d2s(result.getRecall(), 4)).append("\t");
		sb.append(IJ.d2s(result.getFScore(0.5), 4)).append("\t");
		sb.append(IJ.d2s(result.getFScore(1.0), 4)).append("\t");
		sb.append(IJ.d2s(result.getFScore(2.0), 4)).append("\t");
		sb.append(IJ.d2s(result.getFScore(beta), 4));

		if (java.awt.GraphicsEnvironment.isHeadless())
		{
			IJ.log(sb.toString());
		}
		else
		{
			resultsWindow.append(sb.toString());
		}
	}

	private String createPairsHeader()
	{
		StringBuilder sb = new StringBuilder();
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

	private void addPairResult(PointPair pair)
	{
		StringBuilder sb = new StringBuilder();
		Pulse p1 = (Pulse) pair.getPoint1();
		Pulse p2 = (Pulse) pair.getPoint2();
		addPoint(sb, p1);
		addResult(sb, p1, p2);
		pairsWindow.append(sb.toString());
	}

	private void addPoint(StringBuilder sb, Pulse p)
	{
		if (p == null)
		{
			sb.append("-\t-\t-\t-\t-\t");
		}
		else
		{
			sb.append(p.getStart()).append("\t");
			sb.append(p.getEnd()).append("\t");
			sb.append(IJ.d2s(p.getX())).append("\t");
			sb.append(IJ.d2s(p.getY())).append("\t");
			sb.append(IJ.d2s(p.getZ())).append("\t");
		}
	}

	private String createTriplesHeader()
	{
		StringBuilder sb = new StringBuilder();
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.CoordinateProvider#getCoordinates(java.lang.String)
	 */
	public double[] getCoordinates(String line)
	{
		// Extract the startT and x,y coordinates from the first pulse in the line
		final int[] index = { 0, 5, 12 };
		String[] fields = line.split("\t");
		for (int i : index)
		{
			if (i < fields.length)
			{
				if (fields[i].equals("-"))
					continue;
				int startT = Integer.valueOf(fields[i]);
				double x = Double.valueOf(fields[i + 2]);
				double y = Double.valueOf(fields[i + 3]);
				return new double[] { startT, x, y };
			}
		}
		return null;
	}

	private void addTripleResult(Triple triple)
	{
		StringBuilder sb = new StringBuilder();
		Pulse p1 = triple.p1;
		Pulse p2 = triple.p2;
		Pulse p3 = triple.p3;
		addPoint(sb, p1);
		addResult(sb, p1, p2);
		addResult(sb, p1, p3);
		triplesWindow.append(sb.toString());
	}

	private void addResult(StringBuilder sb, Pulse p1, Pulse p2)
	{
		addPoint(sb, p2);
		PointPair pair = new PointPair(p1, p2);
		double d = pair.getXYDistance();
		if (d >= 0)
			sb.append(Utils.rounded(d, 4)).append("\t");
		else
			sb.append("-\t");
		if (p1 != null && p2 != null)
			sb.append(Utils.rounded(p1.score(p2, d * d, dThreshold), 4)).append("\t");
		else
			sb.append("-\t");
	}

	private List<? extends PointPair> sort(List<PointPair> pairs)
	{
		switch (sortIndex)
		{
			case 1: // Sort by time
				ArrayList<TimeComparablePointPair> newPairs = new ArrayList<TimeComparablePointPair>(pairs.size());
				for (PointPair pair : pairs)
				{
					newPairs.add(new TimeComparablePointPair(pair));
				}
				Collections.sort(newPairs);
				return newPairs;

			default:
				// Already sorted by score
				return pairs;
		}
	}

	private List<? extends Triple> sort(ArrayList<Triple> triples)
	{
		if (sortIndex == 1)
		{
			List<TimeComparableTriple> sorted = new ArrayList<TimeComparableTriple>(triples.size());
			for (Triple t : triples)
				sorted.add(new TimeComparableTriple(t));
			Collections.sort(sorted);
			return sorted;
		}
		else
		{
			List<ScoreComparableTriple> sorted = new ArrayList<ScoreComparableTriple>(triples.size());
			for (Triple t : triples)
				sorted.add(new ScoreComparableTriple(t));
			Collections.sort(sorted);
			return sorted;
		}
	}

	private class TimeComparablePointPair extends PointPair implements Comparable<TimeComparablePointPair>
	{
		int startT = Integer.MAX_VALUE;
		public Pulse p1, p2;

		public TimeComparablePointPair(PointPair pair)
		{
			super(pair.getPoint1(), pair.getPoint2());

			p1 = (Pulse) pair.getPoint1();
			p2 = (Pulse) pair.getPoint2();
			if (p2 != null)
				startT = p2.getStart();
			if (p1 != null)
				startT = p1.getStart();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		public int compareTo(TimeComparablePointPair o)
		{
			// Significance of points are: p1, p2, p3
			// Sort by the earliest start time of the most significant point.
			if (startT == o.startT)
			{
				// Sort using the significant points first.
				int result = comparePulsesStartT(p1, o.p1);
				if (result != 0)
					return result;
				result = comparePulsesStartT(p2, o.p2);
				if (result != 0)
					return result;
				result = comparePulsesEndT(p1, o.p1);
				if (result != 0)
					return result;
				result = comparePulsesEndT(p2, o.p2);
				if (result != 0)
					return result;

				// Sort using coords
				result = comparePulsesCoords(p1, o.p1);
				if (result != 0)
					return result;
				return comparePulsesCoords(p2, o.p2);
			}
			return startT - o.startT;
		}
	}

	private static int comparePulsesStartT(Pulse p1, Pulse p2)
	{
		// Data for a point always beats no data
		if (p1 == null)
			return (p2 == null) ? 0 : 1;
		if (p2 == null)
			return -1;
		return p1.getStart() - p2.getStart();
	}

	private static int comparePulsesEndT(Pulse p1, Pulse p2)
	{
		// Data for a point always beats no data
		if (p1 == null)
			return (p2 == null) ? 0 : 1;
		if (p2 == null)
			return -1;
		return p1.getEnd() - p2.getEnd();
	}

	private static int comparePulsesCoords(Pulse p1, Pulse p2)
	{
		// Data for a point always beats no data
		if (p1 == null)
			return (p2 == null) ? 0 : 1;
		if (p2 == null)
			return -1;
		if (p1.getX() < p2.getX())
			return -1;
		if (p1.getX() > p2.getX())
			return 1;
		if (p1.getY() < p2.getY())
			return -1;
		if (p1.getY() > p2.getY())
			return 1;
		return 0;
	}

	private class Triple
	{
		public Pulse p1, p2, p3;

		public Triple(Pulse p1, Pulse p2, Pulse p3)
		{
			this.p1 = p1;
			this.p2 = p2;
			this.p3 = p3;
		}

		public Triple(Triple t)
		{
			this.p1 = t.p1;
			this.p2 = t.p2;
			this.p3 = t.p3;
		}
	}

	private class TimeComparableTriple extends Triple implements Comparable<TimeComparableTriple>
	{
		int startT = Integer.MAX_VALUE;

		public TimeComparableTriple(Triple t)
		{
			super(t);

			if (p3 != null)
				startT = p3.getStart();
			if (p2 != null)
				startT = p2.getStart();
			if (p1 != null)
				startT = p1.getStart();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		public int compareTo(TimeComparableTriple o)
		{
			// Significance of points are: p1, p2, p3
			// Sort by the earliest start time of the most significant point.
			if (startT == o.startT)
			{
				// Sort using the significant points first.
				//int result = comparePulsesStartT(p1, o.p1);
				//if (result != 0)
				//	return result;
				//result = comparePulsesStartT(p2, o.p2);
				//if (result != 0)
				//	return result;
				//result = comparePulsesStartT(p3, o.p3);
				//if (result != 0)
				//	return result;
				//result = comparePulsesEndT(p1, o.p1);
				//if (result != 0)
				//	return result;
				//result = comparePulsesEndT(p2, o.p2);
				//if (result != 0)
				//	return result;
				//return comparePulsesEndT(p3, o.p3);

				// Make the order the same as for the pairs table
				int result = comparePulsesStartT(p1, o.p1);
				if (result != 0)
					return result;
				result = comparePulsesStartT(p2, o.p2);
				if (result != 0)
					return result;
				result = comparePulsesEndT(p1, o.p1);
				if (result != 0)
					return result;
				result = comparePulsesEndT(p2, o.p2);
				if (result != 0)
					return result;

				// Then sort using the third
				result = comparePulsesStartT(p3, o.p3);
				if (result != 0)
					return result;
				result = comparePulsesEndT(p3, o.p3);
				if (result != 0)
					return result;

				// Sort using coords
				result = comparePulsesCoords(p1, o.p1);
				if (result != 0)
					return result;
				result = comparePulsesCoords(p2, o.p2);
				if (result != 0)
					return result;
				return comparePulsesCoords(p3, o.p3);
			}
			return startT - o.startT;
		}
	}

	private class ScoreComparableTriple extends Triple implements Comparable<ScoreComparableTriple>
	{
		double score = 0;

		public ScoreComparableTriple(Triple t)
		{
			super(t);
			if (p1 != null)
			{
				if (p2 != null)
					score = p1.score(p2, dThreshold);
				if (p3 != null)
					score = FastMath.max(score, p1.score(p2, dThreshold));
			}
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		public int compareTo(ScoreComparableTriple o)
		{
			if (score > o.score)
				return -1;
			if (score < o.score)
				return 1;
			return 0;
		}
	}
}
