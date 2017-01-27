package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.TreeSet;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import gdsc.core.clustering.optics.ClusteringResult;
import gdsc.core.clustering.optics.DBSCANResult;
import gdsc.core.clustering.optics.OPTICSCluster;
import gdsc.core.clustering.optics.OPTICSManager;
import gdsc.core.clustering.optics.OPTICSManager.Option;
import gdsc.core.clustering.optics.OPTICSResult;
import gdsc.core.clustering.optics.SampleMode;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.match.RandIndex;
import gdsc.core.utils.ConvexHull;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NotImplementedException;
import gdsc.core.utils.Settings;
import gdsc.core.utils.Sort;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.OPTICSSettings;
import gdsc.smlm.ij.settings.OPTICSSettings.OPTICSMode;
import gdsc.smlm.ij.settings.OPTICSSettings.OutlineMode;
import gdsc.smlm.ij.settings.OPTICSSettings.ClusteringMode;
import gdsc.smlm.ij.settings.OPTICSSettings.ImageMode;
import gdsc.smlm.ij.settings.OPTICSSettings.SpanningTreeMode;
import gdsc.smlm.ij.settings.OPTICSSettings.PlotMode;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.process.LUTHelper;
import ij.process.LUTHelper.LUTMapper;
import ij.process.LUTHelper.LutColour;

/**
 * Run the OPTICS algorithm on the peak results.
 * <p>
 * This is an implementation of the OPTICS method. Mihael Ankerst, Markus M Breunig, Hans-Peter Kriegel, and Jorg
 * Sander. Optics: ordering points to identify the clustering structure. In ACM Sigmod Record, volume 28, pages
 * 49–60. ACM, 1999.
 */
public class OPTICS implements PlugIn
{
	private static final String TITLE_OPTICS = "OPTICS";
	private static final String TITLE_DBSCAN = "DBSCAN";

	private static LUT clusterLut = LUTHelper.getColorModel();
	private static LUT valueLut = LUTHelper.createLUT(LutColour.FIRE);
	// Need to be able to see all colours against white (plot) or black (image) background 
	//private static LUT clusterDepthLut = LUTHelper.createLUT(LutColour.FIRE_LIGHT, true);
	private static LUT clusterDepthLut = LUTHelper.createLUT(LutColour.FIRE_GLOW, true);
	private static LUT clusterOrderLut = LUTHelper.createLUT(LutColour.FIRE_GLOW, true);

	private String TITLE;
	/**
	 * Delay (in milliseconds) used when entering new values in the dialog before the preview is processed
	 */
	private long DELAY = 500;
	private boolean isPreview = false;

	private GlobalSettings globalSettings;
	private OPTICSSettings inputSettings;

	private boolean extraOptions, debug;

	private static class Work
	{
		long time = 0;
		OPTICSSettings inputSettings;
		Settings settings;

		public Work(OPTICSSettings inputSettings, Settings settings)
		{
			this.inputSettings = inputSettings;
			this.settings = settings;
		}

		public Work(OPTICSSettings inputSettings, Object... settings)
		{
			this.inputSettings = inputSettings;
			this.settings = new Settings(settings);
		}
	}

	/**
	 * Allow work to be added to a FIFO stack in a synchronised manner
	 * 
	 * @author Alex Herbert
	 */
	private class WorkStack
	{
		// We only support a stack size of 1
		private Work work = null;

		synchronized void setWork(Work work)
		{
			this.work = work;
		}

		synchronized void addWork(Work work)
		{
			this.work = work;
			this.notify();
		}

		@SuppressWarnings("unused")
		synchronized void close()
		{
			this.work = null;
			this.notify();
		}

		synchronized Work getWork()
		{
			Work work = this.work;
			this.work = null;
			return work;
		}

		boolean isEmpty()
		{
			return work == null;
		}
	}

	// Stack to which the work is first added
	private WorkStack inputStack = new WorkStack();

	private abstract class Worker implements Runnable
	{
		private boolean running = true;
		private Work lastWork = null;
		private Work result = null;
		private WorkStack inbox;
		private WorkStack[] outbox;

		public void run()
		{
			while (running)
			{
				try
				{
					Work work = null;
					synchronized (inbox)
					{
						if (inbox.isEmpty())
						{
							debug("Inbox empty, waiting ...");
							inbox.wait();
						}
						work = inbox.getWork();
						if (work != null)
							debug(" Found work");
					}
					if (work == null)
					{
						debug(" No work, stopping");
						break;
					}

					// Commented out. To stop the worker from doing any more work just close the inbox
					// (so passing in null work).
					//if (!running)
					//{
					//	debug(" Not running, stopping");
					//	break;
					//}

					// Delay processing the work. Allows the work to be updated before we process it.
					if (work.time != 0)
					{
						long time = work.time;
						while (System.currentTimeMillis() < time)
						{
							debug(" Delaying");
							Thread.sleep(50);
							// Assume new work can be added to the inbox. Here we are peaking at the inbox
							// so we do not take ownership with synchronized
							if (inbox.work != null)
								time = inbox.work.time;
						}
						// If we intend to modify the inbox then we should take ownership
						synchronized (inbox)
						{
							if (!inbox.isEmpty())
							{
								work = inbox.getWork();
								debug(" Found updated work");
							}
						}
						// No further delay checking
						work.time = 0;
					}

					if (!equals(work, lastWork))
					{
						// Create a new result
						debug(" Creating new result");
						result = createResult(work);
					}
					else
					{
						// Pass through the new settings with the existing results
						debug(" Updating existing result");
						result = new Work(work.inputSettings, result.settings);
					}
					lastWork = work;
					// Add the result to the output
					if (outbox != null)
					{
						debug(" Posting result");
						for (int i = outbox.length; i-- > 0;)
							outbox[i].addWork(result);
					}
				}
				catch (InterruptedException e)
				{
					debug(" Interrupted, stopping");
					break;
				}
			}
		}

		private void debug(String msg)
		{
			if (debug)
				System.out.println(this.getClass().getSimpleName() + msg);
		}

		boolean equals(Work work, Work lastWork)
		{
			if (lastWork == null)
				return false;

			// We must selectively compare this as not all settings changes matter.
			if (!equals(work.inputSettings, lastWork.inputSettings))
				return false;

			// We can compare these here using object references. 
			// Any new results passed in will trigger equals to fail.
			boolean result = work.settings.equals(lastWork.settings);
			if (!result)
				newResults();
			return result;
		}

		/**
		 * Compare the settings and return false if any settings that the work depends on have changed
		 * 
		 * @param current
		 * @param previous
		 * @return
		 */
		abstract boolean equals(OPTICSSettings current, OPTICSSettings previous);

		abstract Work createResult(Work work);

		/**
		 * Called when there are new results in the current work.
		 */
		void newResults()
		{
		}
	}

	private class InputWorker extends Worker
	{
		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			// Nothing in the settings effects if we have to create a new OPTICS manager
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			// The first item should be the memory peak results 
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			// Convert results to coordinates
			float[] x, y;
			int size = results.size();
			x = new float[size];
			y = new float[size];
			ArrayList<PeakResult> list = (ArrayList<PeakResult>) results.getResults();
			for (int i = 0; i < size; i++)
			{
				PeakResult p = list.get(i);
				x[i] = p.getXPosition();
				y[i] = p.getYPosition();
			}
			Rectangle bounds = results.getBounds(true);
			OPTICSManager opticsManager = new OPTICSManager(x, y, bounds);
			opticsManager.setTracker(new IJTrackProgress());
			opticsManager.setOptions(Option.CACHE);
			return new Work(work.inputSettings, results, opticsManager);
		}
	}

	private class OpticsWorker extends Worker
	{
		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.minPoints != previous.minPoints)
				return false;
			if (current.getOPTICSMode() != previous.getOPTICSMode())
				return false;
			if (current.getOPTICSMode() == OPTICSMode.OPTICS)
			{
				if (current.generatingDistance != previous.generatingDistance)
					return false;
			}
			else
			{
				if (current.numberOfSplitSets != previous.numberOfSplitSets)
					return false;
				if (extraOptions)
				{
					if (current.useRandomVectors != previous.useRandomVectors)
						return false;
					if (current.saveApproximateSets != previous.saveApproximateSets)
						return false;
					if (current.getSampleMode() != previous.getSampleMode())
						return false;
				}
			}
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			// The first item should be the memory peak results 
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			// The second item should be the OPTICS manager
			OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);

			int minPts = work.inputSettings.minPoints;

			OPTICSResult opticsResult;
			if (work.inputSettings.getOPTICSMode() == OPTICSMode.FAST_OPTICS)
			{
				int n = work.inputSettings.numberOfSplitSets;
				// Q. Should these be options
				boolean useRandomVectors = false;
				boolean saveApproximateSets = false;
				SampleMode sampleMode = SampleMode.RANDOM;
				if (extraOptions)
				{
					useRandomVectors = work.inputSettings.useRandomVectors;
					saveApproximateSets = work.inputSettings.saveApproximateSets;
					sampleMode = work.inputSettings.getSampleMode();
				}
				synchronized (opticsManager)
				{
					opticsManager.setNumberOfThreads(Prefs.getThreads());
					opticsResult = opticsManager.fastOptics(minPts, n, n, useRandomVectors, saveApproximateSets,
							sampleMode);
				}
			}
			else
			{
				double distance = work.inputSettings.generatingDistance;
				if (distance > 0)
				{
					// Convert generating distance to pixels
					double nmPerPixel = getNmPerPixel(results);
					if (nmPerPixel != 1)
					{
						double newDistance = distance / nmPerPixel;
						Utils.log(TITLE + ": Converting generating distance %s nm to %s pixels",
								Utils.rounded(distance), Utils.rounded(newDistance));
						distance = newDistance;
					}
				}
				else
				{
					double nmPerPixel = getNmPerPixel(results);
					if (nmPerPixel != 1)
					{
						Utils.log(TITLE + ": Default generating distance %s nm",
								Utils.rounded(opticsManager.computeGeneratingDistance(minPts) * nmPerPixel));
					}
				}

				synchronized (opticsManager)
				{
					opticsResult = opticsManager.optics((float) distance, minPts);
				}
			}
			// It may be null if cancelled. However return null Work will close down the next thread
			return new Work(work.inputSettings, results, opticsManager, opticsResult);
		}
	}

	private class OpticsClusterWorker extends Worker
	{
		int clusterCount = 0;

		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.getClusteringMode() != previous.getClusteringMode())
				return false;
			if (current.getClusteringMode() == ClusteringMode.XI)
			{
				if (current.xi != previous.xi)
					return false;
				if (current.topLevel != previous.topLevel)
					return false;
			}
			else
			{
				if (current.clusteringDistance != previous.clusteringDistance)
					return false;
				if (current.core != previous.core)
					return false;
			}
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);
			OPTICSResult opticsResult = (OPTICSResult) work.settings.get(2);
			// It may be null if cancelled.
			if (opticsResult != null)
			{
				synchronized (opticsResult)
				{
					if (work.inputSettings.getClusteringMode() == ClusteringMode.XI)
					{
						int options = (work.inputSettings.topLevel) ? OPTICSResult.XI_OPTION_TOP_LEVEL : 0;
						opticsResult.extractClusters(work.inputSettings.xi, options);
					}
					else
					{
						double nmPerPixel = getNmPerPixel(results);

						double distance;
						if (work.inputSettings.getOPTICSMode() == OPTICSMode.FAST_OPTICS)
						{
							if (work.inputSettings.clusteringDistance > 0)
								distance = work.inputSettings.clusteringDistance;
							else
							{
								distance = opticsManager.computeGeneratingDistance(work.inputSettings.minPoints) *
										nmPerPixel;
								if (nmPerPixel != 1)
								{
									Utils.log(TITLE + ": Default clustering distance %s nm", Utils.rounded(distance));
								}
							}
						}
						else
						{
							// Ensure that the distance is valid
							distance = opticsResult.generatingDistance * nmPerPixel;
							if (work.inputSettings.clusteringDistance > 0)
								distance = Math.min(work.inputSettings.clusteringDistance, distance);
						}

						if (nmPerPixel != 1)
						{
							double newDistance = distance / nmPerPixel;
							Utils.log(TITLE + ": Converting clustering distance %s nm to %s pixels",
									Utils.rounded(distance), Utils.rounded(newDistance));
							distance = newDistance;
						}

						opticsResult.extractDBSCANClustering((float) distance, work.inputSettings.core);
					}
				}
				// We created a new clustering
				clusterCount++;
			}
			return new Work(work.inputSettings, results, opticsManager, opticsResult, clusterCount);
		}
	}

	private class ResultsWorker extends Worker
	{
		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			// Only depends on if the clustering results are new. This is triggered 
			// in the default comparison of the Settings object.
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			// The result is in position 2.
			// It may be null if cancelled.
			if (work.settings.get(2) == null)
			{
				// Only log here so it happens once
				IJ.log(TITLE + ": No results to display");
			}
			return work;
		}
	}

	private class ClusterResult
	{
		final int n1, n2;
		final int[] c1, c2;

		ClusterResult(int[] clusters, int[] topClusters)
		{
			// The original set of clusters does not need to be compacted
			n1 = Maths.max(clusters) + 1;
			this.c1 = clusters;
			if (topClusters != null)
			{
				// The top clusters from OPTICS may contain non-sequential integers
				n2 = RandIndex.compact(topClusters);
				this.c2 = topClusters;
			}
			else
			{
				n2 = 0;
				c2 = null;
			}
		}
	}

	private class RandIndexWorker extends Worker
	{
		Queue<ClusterResult> queue = new LinkedList<ClusterResult>();

		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (!current.inputOption.equals(previous.inputOption))
			{
				// We only cache results for the same set of raw results, i.e. the input coordinates.
				queue.clear();
				return false;
			}

			return true;
		}

		@Override
		Work createResult(Work work)
		{
			ClusteringResult clusteringResult = (ClusteringResult) work.settings.get(2);
			if (clusteringResult == null)
				return work;

			int[] clusters, topClusters = null;
			synchronized (clusteringResult)
			{
				clusters = clusteringResult.getClusters();
				if (clusteringResult instanceof OPTICSResult)
				{
					topClusters = ((OPTICSResult) clusteringResult).getTopLevelClusters(false);
				}
			}

			ClusterResult current = new ClusterResult(clusters, topClusters);

			// Compare to previous results
			if (!queue.isEmpty())
			{
				int i = -queue.size();
				StringBuilder sb = new StringBuilder();
				sb.append("Cluster comparison: RandIndex (AdjustedRandIndex)\n");
				for (Iterator<ClusterResult> it = queue.iterator(); it.hasNext();)
				{
					ClusterResult previous = it.next();
					sb.append("[").append(i++).append("] ");
					compare(sb, "Clusters", current.c1, current.n1, previous.c1, previous.n1);
					if (current.c2 != null)
					{
						sb.append(" : ");
						compare(sb, "Top-level clusters", current.c2, current.n2, previous.c2, previous.n2);
					}
					sb.append('\n');
				}
				IJ.log(sb.toString());
			}

			queue.add(current);

			// Limit size
			if (queue.size() > 2)
				queue.poll();

			return work;
		}

		private void compare(StringBuilder sb, String title, int[] set1, int n1, int[] set2, int n2)
		{
			RandIndex ri = new RandIndex();
			ri.compute(set1, n1, set2, n2);
			double r = ri.getRandIndex();
			double ari = ri.getAdjustedRandIndex();
			sb.append(title);
			sb.append(" ");
			sb.append(Utils.rounded(r));
			sb.append(" (");
			sb.append(Utils.rounded(ari));
			sb.append(")");

		}
	}

	private class MemoryResultsWorker extends Worker
	{
		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			// Only depends on if the clustering results are new. This is triggered 
			// in the default comparison of the Settings object.
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			//OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);
			ClusteringResult clusteringResult = (ClusteringResult) work.settings.get(2);
			// It may be null if cancelled.
			if (clusteringResult != null)
			{
				int[] clusters;
				synchronized (clusteringResult)
				{
					clusters = clusteringResult.getClusters();
				}
				int max = Maths.max(clusters);

				// Save the clusters to memory
				Trace[] traces = new Trace[max + 1];
				for (int i = 0; i <= max; i++)
				{
					traces[i] = new Trace();
					traces[i].setId(i);
				}
				ArrayList<PeakResult> list = (ArrayList<PeakResult>) results.getResults();
				for (int i = 0, size = results.size(); i < size; i++)
				{
					PeakResult r = list.get(i);
					traces[clusters[i++]].add(r);
				}
				TraceMolecules.saveResults(results, traces, TITLE);
			}

			// We have not created anything new so return the current object
			return work;
		}
	}

	private class ReachabilityResultsWorker extends Worker
	{
		// TODO - options to colour the reachability plot by the cluster depth

		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.getPlotMode() != previous.getPlotMode())
				return false;
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			OPTICSResult opticsResult = (OPTICSResult) work.settings.get(2);
			// It may be null if cancelled.
			if (opticsResult == null)
			{
				return work;
			}

			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			double nmPerPixel = getNmPerPixel(results);

			// Draw the reachability profile
			PlotMode mode = work.inputSettings.getPlotMode();
			if (mode != PlotMode.OFF)
			{
				double[] profile;
				synchronized (opticsResult)
				{
					profile = opticsResult.getReachabilityDistanceProfile(true);
				}
				String units = " (px)";
				if (nmPerPixel != 1)
				{
					units = " (nm)";
					for (int i = 0; i < profile.length; i++)
						profile[i] *= nmPerPixel;
				}

				double[] order = Utils.newArray(profile.length, 1.0, 1.0);
				String title = TITLE + " Reachability Distance";
				Plot plot = new Plot(title, "Order", "Reachability" + units);
				double[] limits = Maths.limits(profile);
				// plot to zero
				limits[0] = 0;

				ArrayList<OPTICSCluster> clusters = null;
				LUT lut = clusterLut;
				int maxClusterId = 0;
				int maxLevel = 0;

				if (mode.requiresClusters())
				{
					synchronized (opticsResult)
					{
						clusters = opticsResult.getAllClusters();
					}

					for (OPTICSCluster cluster : clusters)
					{
						if (maxLevel < cluster.getLevel())
							maxLevel = cluster.getLevel();
						if (maxClusterId < cluster.clusterId)
							maxClusterId = cluster.clusterId;
					}
				}

				if (work.inputSettings.getOPTICSMode() == OPTICSMode.FAST_OPTICS)
				{
					// The profile may be very high. Compute the outliers and remove.
					Percentile p = new Percentile();
					p.setData(profile);
					double max;
					boolean useIQR = true;
					if (useIQR)
					{
						double lq = p.evaluate(25);
						double uq = p.evaluate(75);
						max = (uq - lq) * 2 + uq;
					}
					else
					{
						// Remove top 2%
						max = p.evaluate(98);
					}
					if (limits[1] > max)
						limits[1] = max;
				}
				else
				{
					// Show extra at the top
					limits[1] *= 1.05;
				}

				double distance = -1;
				if (inputSettings.getClusteringMode() == ClusteringMode.DBSCAN)
				{
					if (work.inputSettings.getOPTICSMode() == OPTICSMode.FAST_OPTICS)
					{
						if (work.inputSettings.clusteringDistance > 0)
							distance = work.inputSettings.clusteringDistance;
						else
						{
							OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);
							distance = opticsManager.computeGeneratingDistance(work.inputSettings.minPoints) *
									nmPerPixel;
						}
					}
					else
					{
						// Ensure that the distance is valid
						distance = opticsResult.generatingDistance * nmPerPixel;
						if (work.inputSettings.clusteringDistance > 0)
							distance = Math.min(work.inputSettings.clusteringDistance, distance);
					}

					if (distance > limits[1])
						limits[1] = distance * 1.05;
				}

				// Draw the clusters using lines underneath
				if (mode.isDrawClusters() && maxClusterId > 0)
				{
					// Get a good distance to start the lines, and the separation
					// Make the clusters fill 1/3 of the plot.
					double range = 0.5 * limits[1];
					double separation = range / (maxLevel + 2);
					double start = -separation;

					LUTMapper mapper = new LUTHelper.NonZeroLUTMapper(1, maxClusterId);
					for (OPTICSCluster cluster : clusters)
					{
						int level = cluster.getLevel();
						double y = start - (maxLevel - level) * separation;
						plot.setColor(mapper.getColour(lut, cluster.clusterId));
						plot.drawLine(cluster.start, y, cluster.end, y);
					}

					// Update the limits if we are plotting lines underneath for the clusters
					limits[0] = -range; //start - (maxLevel + 1) * separation;
				}

				plot.setLimits(1, order.length, limits[0], limits[1]);

				// Create the colour for each point on the line
				int[] profileColour = new int[profile.length];

				//plot.setColor(Color.black);
				//plot.addPoints(order, profile, Plot.LINE);

				// Create a colour to match the LUT of the image
				LUTMapper mapper = null;

				// Colour the reachability plot line if it is in a cluster. Use a default colour
				if (mode.isColourProfile())
				{
					if (mode.isColourProfileByOrder())
					{
						lut = clusterOrderLut;
						mapper = new LUTHelper.NonZeroLUTMapper(0, profileColour.length - 1);
						for (int i = 0; i < profileColour.length; i++)
							profileColour[i] = mapper.map(i);
						// Ensure we correctly get colours for each value 
						mapper = new LUTHelper.DefaultLUTMapper(0, 255);
					}
					else
					{
						// Do all clusters so rank by level
						Collections.sort(clusters, new Comparator<OPTICSCluster>()
						{
							public int compare(OPTICSCluster o1, OPTICSCluster o2)
							{
								return o1.getLevel() - o2.getLevel();
							}
						});

						final boolean useLevel = mode.isColourProfileByDepth();
						if (useLevel)
							lut = clusterDepthLut;

						for (OPTICSCluster cluster : clusters)
						{
							Arrays.fill(profileColour, cluster.start, cluster.end + 1,
									(useLevel) ? cluster.getLevel() + 1 : cluster.clusterId);
						}
					}
				}
				else if (mode.isHighlightProfile())
				{
					for (OPTICSCluster cluster : clusters)
					{
						// Only do top level clusters
						if (cluster.getLevel() != 0)
							continue;
						Arrays.fill(profileColour, cluster.start, cluster.end + 1, 1);
					}
				}

				// Now draw the line
				int maxColour = Maths.max(profileColour);
				if (mapper == null)
					// Make zero black
					mapper = new LUTHelper.NonZeroLUTMapper(1, maxColour);

				// Cache all the colours
				Color[] colors = new Color[maxColour + 1];
				if (mode.isColourProfile())
				{
					for (int c = mapper.getMin(); c <= maxColour; c++)
						colors[c] = mapper.getColour(lut, c);
				}
				else if (maxColour == 1)
				{
					colors[1] = Color.BLUE;
				}
				if (colors[0] == null)
					colors[0] = Color.BLACK;

				int from = 0;
				for (int i = 1; i < profileColour.length; i++)
				{
					if (profileColour[from] != profileColour[i])
					{
						// Draw the line on the plot
						int to = i + 1;
						double[] order1 = Arrays.copyOfRange(order, from, to);
						double[] profile1 = Arrays.copyOfRange(profile, from, to);
						plot.setColor(colors[profileColour[from]]);
						plot.addPoints(order1, profile1, Plot.LINE);

						from = i;
					}
				}

				// Draw the final line
				if (from != profileColour.length - 1)
				{
					int to = profileColour.length;
					double[] order1 = Arrays.copyOfRange(order, from, to);
					double[] profile1 = Arrays.copyOfRange(profile, from, to);
					plot.setColor(colors[profileColour[from]]);
					plot.addPoints(order1, profile1, Plot.LINE);
				}

				// Add the DBSCAN clustering distance
				if (distance > -1)
				{
					plot.setColor(Color.red);
					plot.drawLine(1, distance, order.length, distance);
				}

				Utils.display(title, plot);
			}
			else
			{
				// We could close an existing plot here. 
				// However we leave it as the user may wish to keep it for something.
			}

			// We have not created anything new so return the current object
			return work;
		}
	}

	private class OrderProvider
	{
		int getOrder(int i)
		{
			return 0;
		}
	}

	private class RealOrderProvider extends OrderProvider
	{
		final int[] order;

		RealOrderProvider(int[] order)
		{
			this.order = order;
		}

		@Override
		int getOrder(int i)
		{
			return order[i];
		}
	}

	private class ImageResultsWorker extends Worker
	{
		IJImagePeakResults image = null;
		OutlineMode lastOutlineMode = null;
		Overlay outline = null;
		SpanningTreeMode lastSpanningTreeMode = null;
		Overlay spanningTree = null;

		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.imageScale != previous.imageScale)
			{
				// Clear all the cached results
				clearCache(true);
				return false;
			}
			boolean result = true;
			if (current.getOutlineMode() != previous.getOutlineMode())
			{
				if (current.getOutlineMode().isOutline() && current.getOutlineMode() != lastOutlineMode)
					outline = null;
				result = false;
			}
			if (current.getSpanningTreeMode() != previous.getSpanningTreeMode())
			{
				if (current.getSpanningTreeMode().isSpanningTree() &&
						current.getSpanningTreeMode() != lastSpanningTreeMode)
					spanningTree = null;
				result = false;
			}
			if (current.getImageMode() != previous.getImageMode() ||
					getDisplayFlags(current) != getDisplayFlags(previous))
			{
				// We can only cache the image if the display mode is the same
				image = null;
				result = false;
			}
			return result;
		}

		@Override
		void newResults()
		{
			// We can keep the image but should clear the overlays
			clearCache(false);
		}

		private void clearCache(boolean clearImage)
		{
			// Clear cache
			if (clearImage)
				image = null;
			lastOutlineMode = null;
			outline = null;
			lastSpanningTreeMode = null;
			spanningTree = null;
		}

		@Override
		Work createResult(Work work)
		{
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);
			ClusteringResult clusteringResult = (ClusteringResult) work.settings.get(2);
			int clusterCount = (Integer) work.settings.get(3);
			// It may be null if cancelled.
			if (clusteringResult == null)
			{
				clearCache(true);
				return new Work(work.inputSettings, results, opticsManager, clusteringResult, clusterCount, image);
			}

			int[] clusters = null, order = null;
			int max = 0; // max cluster value
			float[] x = null, y = null;
			ImageMode mode = work.inputSettings.getImageMode();

			if (work.inputSettings.imageScale > 0)
			{
				if (image == null)
				{
					// Display the results ...

					Rectangle bounds = results.getBounds();
					image = new IJImagePeakResults(results.getName() + " " + TITLE, bounds,
							(float) work.inputSettings.imageScale);
					// Options to control rendering
					image.copySettings(results);
					image.setDisplayFlags(getDisplayFlags(work.inputSettings));
					image.setLiveImage(false);
					image.begin();
					ImagePlus imp = image.getImagePlus();
					imp.setOverlay(null);
					if (mode != ImageMode.NONE)
					{
						int[] map = null; // Used to map clusters to a display value
						int max2 = 0; // max mapped cluster value

						synchronized (clusteringResult)
						{
							clusters = clusteringResult.getClusters();
							max = Maths.max(clusters);
							if (clusteringResult instanceof OPTICSResult)
							{
								OPTICSResult opticsResult = (OPTICSResult) clusteringResult;
								if (mode == ImageMode.CLUSTER_ORDER)
								{
									order = opticsResult.getOrder();
								}
								else if (mode == ImageMode.CLUSTER_DEPTH)
								{
									ArrayList<OPTICSCluster> allClusters = opticsResult.getAllClusters();
									map = new int[max + 1];
									for (OPTICSCluster c : allClusters)
										map[c.clusterId] = c.getLevel() + 1;
									max2 = Maths.max(map);
								}
							}
							if (map == null)
							{
								map = Utils.newArray(max + 1, 0, 1);
								max2 = max;
							}
						}

						// Draw each cluster in a new colour
						LUT lut;
						LUTMapper mapper = null;
						if (mode.isMapped())
						{
							switch (mode)
							{
								//@formatter:off
								case CLUSTER_ORDER: 
									lut = clusterOrderLut;
									mapper = new LUTHelper.NonZeroLUTMapper(1, results.size());																		
									break;								
								case CLUSTER_ID: lut = clusterLut; break;
								case CLUSTER_DEPTH: lut = clusterDepthLut; break;
								default:
									throw new NotImplementedException();
								//@formatter:on
							}
							if (mapper == null)
								mapper = new LUTHelper.NonZeroLUTMapper(1, max2);
						}
						else
						{
							lut = valueLut;
							mapper = new LUTHelper.NullLUTMapper();
						}
						image.getImagePlus().getProcessor().setColorModel(lut);

						// Add in a single batch
						float[] v;
						x = new float[results.size()];
						y = new float[x.length];
						v = new float[x.length];
						ArrayList<PeakResult> list = (ArrayList<PeakResult>) results.getResults();
						OrderProvider op = (order == null) ? new OrderProvider() : new RealOrderProvider(order);
						for (int i = 0, size = results.size(); i < size; i++)
						{
							PeakResult r = list.get(i);
							x[i] = r.getXPosition();
							y[i] = r.getYPosition();
							v[i] = mapper.mapf(mode.getValue(r.getSignal(), map[clusters[i]], op.getOrder(i)));
						}
						image.add(x, y, v);
					}
					image.end();
					if (mode.isMapped())
					{
						// Convert already mapped image to 8-bit (so the values are fixed)
						imp.setProcessor(imp.getProcessor().convertToByteProcessor(false));
					}
					imp.getWindow().toFront();
				}
			}
			else
			{
				// We could close an image here. 
				// However we leave it as the user may wish to keep it for something.
			}

			// Note: If the image scale is set to zero then the image cache will be cleared and the image will be null.
			// This will prevent computing an overlay even if the 'outline' setting is enabled.

			if (image != null)
			{
				ImagePlus imp = image.getImagePlus();
				Overlay overlay = null;

				int[] map = null; // Used to map clusters to a display value
				int max2 = 0; // max mapped cluster value

				if (work.inputSettings.getOutlineMode().isOutline())
				{
					if (outline == null)
					{
						lastOutlineMode = work.inputSettings.getOutlineMode();
						if (clusters == null)
						{
							synchronized (clusteringResult)
							{
								clusters = clusteringResult.getClusters();
								max = Maths.max(clusters);
							}
						}

						max2 = max;
						map = Utils.newArray(max + 1, 0, 1);

						LUT lut = clusterLut;

						if (clusteringResult instanceof OPTICSResult)
						{
							if (work.inputSettings.getOutlineMode().isColourByDepth())
							{
								lut = clusterDepthLut;
								synchronized (clusteringResult)
								{
									OPTICSResult opticsResult = (OPTICSResult) clusteringResult;
									ArrayList<OPTICSCluster> allClusters = opticsResult.getAllClusters();
									Arrays.fill(map, 0);
									for (OPTICSCluster c : allClusters)
										map[c.clusterId] = c.getLevel() + 1;
									max2 = Maths.max(map);
								}
							}
						}

						outline = new Overlay();
						ConvexHull[] hulls = new ConvexHull[max + 1];
						synchronized (clusteringResult)
						{
							clusteringResult.computeConvexHulls();
							for (int c = 1; c <= max; c++)
							{
								hulls[c] = clusteringResult.getConvexHull(c);
							}
						}

						// Create a colour to match the LUT of the image
						LUTMapper mapper = new LUTHelper.NonZeroLUTMapper(1, max2);

						// Cache all the colours
						Color[] colors = new Color[max2 + 1];
						for (int c = 1; c <= max2; c++)
							colors[c] = mapper.getColour(lut, c);

						// Extract the ConvexHull of each cluster
						for (int c = 1; c <= max; c++)
						{
							ConvexHull hull = hulls[c];
							if (hull != null)
							{
								// Convert the Hull to the correct image scale.
								float[] x2 = hull.x.clone();
								float[] y2 = hull.y.clone();
								for (int i = 0; i < x2.length; i++)
								{
									x2[i] = image.mapX(x2[i]);
									y2[i] = image.mapY(y2[i]);
								}
								PolygonRoi roi = new PolygonRoi(x2, y2, Roi.POLYGON);
								roi.setStrokeColor(colors[map[c]]);
								// TODO: Options to set a fill colour?
								outline.add(roi);
							}
						}
					}
					overlay = outline;
				}

				if (work.inputSettings.getSpanningTreeMode().isSpanningTree() &&
						clusteringResult instanceof OPTICSResult)
				{
					if (spanningTree == null)
					{
						OPTICSResult opticsResult = (OPTICSResult) clusteringResult;
						lastSpanningTreeMode = work.inputSettings.getSpanningTreeMode();

						int[] predecessor, topLevelClusters;
						synchronized (opticsResult)
						{
							predecessor = opticsResult.getPredecessor();
							if (order == null)
								order = opticsResult.getOrder();
							topLevelClusters = opticsResult.getTopLevelClusters(false);

							if (clusters == null)
							{
								clusters = opticsResult.getClusters();
								max = Maths.max(clusters);
							}
						}

						if (map == null)
						{
							max2 = max;
							map = Utils.newArray(max + 1, 0, 1);
						}

						LUT lut = clusterLut;

						if (lastSpanningTreeMode.isColourByOrder())
						{
							lut = clusterOrderLut;
						}
						else if (lastSpanningTreeMode.isColourByDepth() && max == max2)
						{
							lut = clusterDepthLut;
							synchronized (clusteringResult)
							{
								ArrayList<OPTICSCluster> allClusters = opticsResult.getAllClusters();
								Arrays.fill(map, 0);
								for (OPTICSCluster c : allClusters)
									map[c.clusterId] = c.getLevel() + 1;
								max2 = Maths.max(map);
							}
						}

						// Get the coordinates
						if (x == null)
						{
							int size = results.size();
							x = new float[size];
							y = new float[x.length];
							ArrayList<PeakResult> list = (ArrayList<PeakResult>) results.getResults();
							for (int i = 0; i < size; i++)
							{
								PeakResult r = list.get(i);
								x[i] = r.getXPosition();
								y[i] = r.getYPosition();
							}
						}

						spanningTree = new Overlay();

						// Cache all the colours
						Color[] colors;
						// Create a colour to match the LUT of the image
						LUTMapper mapper;

						boolean useMap;
						if (lastSpanningTreeMode.isColourByOrder())
						{
							// We will use the order for the colour
							useMap = false;
							mapper = new LUTHelper.DefaultLUTMapper(0, 255);
							colors = new Color[256];
							for (int c = 1; c < colors.length; c++)
								colors[c] = mapper.getColour(lut, c);
							mapper = new LUTHelper.NonZeroLUTMapper(1, clusters.length);
						}
						else
						{
							// Alternative is to colour by cluster Id/Depth using a map
							useMap = true;
							mapper = new LUTHelper.NonZeroLUTMapper(1, max2);
							colors = new Color[max2 + 1];
							for (int c = 1; c <= max2; c++)
								colors[c] = mapper.getColour(lut, c);
						}

						for (int i = 1; i < predecessor.length; i++)
						{
							if (clusters[i] == 0 || predecessor[i] < 0)
								continue;

							int j = predecessor[i];

							// The spanning tree can jump across hierachical clusters.
							// Prevent jumps across top-level clusters
							if (topLevelClusters[i] != topLevelClusters[i])
								continue;

							float xi = image.mapX(x[i]);
							float yi = image.mapY(y[i]);
							float xj = image.mapX(x[j]);
							float yj = image.mapY(y[j]);

							Line roi = new Line(xi, yi, xj, yj);
							if (useMap)
								roi.setStrokeColor(colors[map[clusters[i]]]);
							else
								roi.setStrokeColor(colors[mapper.map(order[i])]);
							spanningTree.add(roi);
						}
					}
					if (overlay == null)
					{
						overlay = spanningTree;
					}
					else
					{
						// Merge the two
						overlay = new Overlay();
						for (int i = outline.size(); i-- > 0;)
							overlay.add(outline.get(i));
						for (int i = spanningTree.size(); i-- > 0;)
							overlay.add(spanningTree.get(i));
					}
				}

				imp.setOverlay(overlay);
			}

			return new Work(work.inputSettings, results, opticsManager, clusteringResult, clusterCount, image);
		}

		private int getDisplayFlags(OPTICSSettings inputSettings)
		{
			int displayFlags = 0;
			if (inputSettings.getImageMode().canBeWeighted())
			{
				if (inputSettings.weighted)
					displayFlags |= IJImagePeakResults.DISPLAY_WEIGHTED;
				if (inputSettings.equalised)
					displayFlags |= IJImagePeakResults.DISPLAY_EQUALIZED;
			}
			if (inputSettings.getImageMode() == ImageMode.CLUSTER_ID)
				displayFlags = IJImagePeakResults.DISPLAY_REPLACE;
			else if (inputSettings.getImageMode() == ImageMode.CLUSTER_DEPTH ||
					inputSettings.getImageMode() == ImageMode.CLUSTER_ORDER)
				displayFlags = IJImagePeakResults.DISPLAY_MAX;
			return displayFlags;
		}
	}

	private class KNNWorker extends Worker
	{
		double[] profile = null;

		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.minPoints != previous.minPoints)
			{
				newResults();
				return false;
			}
			if (current.samples != previous.samples || current.sampleFraction != previous.sampleFraction)
			{
				newResults();
				return false;
			}
			if (current.fractionNoise != previous.fractionNoise)
				return false;
			if (clusteringDistanceChange(current.clusteringDistance, previous.clusteringDistance))
				return false;

			return true;
		}

		@Override
		void newResults()
		{
			// Clear cache
			profile = null;
		}

		@Override
		Work createResult(Work work)
		{
			// The first item should be the memory peak results 
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			// The second item should be the OPTICS manager
			OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);

			int minPts = work.inputSettings.minPoints;
			int k = minPts - 1; // Since min points includes the actual point
			double fractionNoise = work.inputSettings.fractionNoise;

			double nmPerPixel = getNmPerPixel(results);

			// Create a profile of the K-Nearest Neighbour distances
			if (profile == null)
			{
				synchronized (opticsManager)
				{
					int samples = work.inputSettings.samples;
					if (samples > 1 || work.inputSettings.sampleFraction > 0)
					{
						// Ensure we take a reasonable amount of samples (min=100)
						samples = Maths.max(100, samples,
								(int) Math.ceil(opticsManager.getSize() * work.inputSettings.sampleFraction));
					}
					float[] d = opticsManager.nearestNeighbourDistance(k, samples, true);
					profile = new double[d.length];
					for (int i = d.length; i-- > 0;)
						profile[i] = d[i];
				}
				Arrays.sort(profile);
				Sort.reverse(profile);
				if (nmPerPixel != 1)
				{
					for (int i = 0; i < profile.length; i++)
						profile[i] *= nmPerPixel;
				}
			}

			String units = (nmPerPixel != 1) ? " (nm)" : " (px)";

			double[] order = Utils.newArray(profile.length, 1.0, 1.0);
			String title = TITLE + " KNN Distance";
			Plot plot = new Plot(title, "Sample", k + "-NN Distance" + units);
			double[] limits = new double[] { profile[profile.length - 1], profile[0] };

			plot.setLimits(1, order.length, limits[0], limits[1] * 1.05);

			plot.setColor(Color.black);
			plot.addPoints(order, profile, Plot.LINE);

			// Add the DBSCAN clustering distance
			double distance = work.inputSettings.clusteringDistance;
			plot.setColor(Color.red);
			plot.drawLine(1, distance, order.length, distance);

			// Find the clustering distance using a % noise in the KNN distance samples
			distance = findClusteringDistance(profile, fractionNoise);
			plot.setColor(Color.blue);
			plot.drawDottedLine(1, distance, order.length, distance, 2);

			Utils.display(title, plot);

			if (work.inputSettings.clusteringDistance <= 0)
			{
				// Set this distance into the settings if there is no clustering distance
				// Use a negative value to show it is an auto-distance
				work.inputSettings.clusteringDistance = -distance;
			}

			// We have not created anything new so return the current object
			return work;
		}
	}

	private boolean clusteringDistanceChange(double newD, double oldD)
	{
		if (newD <= 0 && oldD <= 0)
			// Auto-distance
			return false;

		return newD != oldD;
	}

	/**
	 * Find the clustering distance using a sorted profile of the KNN distance.
	 *
	 * @param profile
	 *            the profile (sorted high to low)
	 * @param fractionNoise
	 *            the fraction noise
	 * @return the clustering distance
	 */
	public static double findClusteringDistance(double[] profile, double fractionNoise)
	{
		// Return the next distance after the fraction has been achieved
		int n = Maths.clip(0, profile.length - 1, (int) Math.ceil(profile.length * fractionNoise));
		return profile[n];
	}

	private class DBSCANWorker extends Worker
	{
		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.minPoints != previous.minPoints)
				return false;
			if (clusteringDistanceChange(current.clusteringDistance, previous.clusteringDistance))
				return false;
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			// The first item should be the memory peak results 
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			// The second item should be the OPTICS manager
			OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);

			double clusteringDistance = Math.abs(work.inputSettings.clusteringDistance);
			int minPts = work.inputSettings.minPoints;
			if (clusteringDistance > 0)
			{
				// Convert clustering distance to pixels
				double nmPerPixel = getNmPerPixel(results);
				if (nmPerPixel != 1)
				{
					double newGeneratingDistance = clusteringDistance / nmPerPixel;
					Utils.log(TITLE + ": Converting clustering distance %s nm to %s pixels",
							Utils.rounded(clusteringDistance), Utils.rounded(newGeneratingDistance));
					clusteringDistance = newGeneratingDistance;
				}
			}
			else
			{
				// Note: This should not happen since the clustering distance is set using the KNN distance samples

				double nmPerPixel = getNmPerPixel(results);
				if (nmPerPixel != 1)
				{
					Utils.log(TITLE + ": Default clustering distance %s nm",
							Utils.rounded(opticsManager.computeGeneratingDistance(minPts) * nmPerPixel));
				}
			}

			DBSCANResult dbscanResult;
			synchronized (opticsManager)
			{
				dbscanResult = opticsManager.dbscan((float) clusteringDistance, minPts);
			}
			// It may be null if cancelled. However return null Work will close down the next thread
			return new Work(work.inputSettings, results, opticsManager, dbscanResult);
		}
	}

	private class DBSCANClusterWorker extends Worker
	{
		int clusterCount = 0;

		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.core != previous.core)
				return false;
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			MemoryPeakResults results = (MemoryPeakResults) work.settings.get(0);
			OPTICSManager opticsManager = (OPTICSManager) work.settings.get(1);
			DBSCANResult dbscanResult = (DBSCANResult) work.settings.get(2);
			// It may be null if cancelled.
			if (dbscanResult != null)
			{
				synchronized (dbscanResult)
				{
					dbscanResult.extractClusters(work.inputSettings.core);
				}
				// We created a new clustering
				clusterCount++;
			}
			return new Work(work.inputSettings, results, opticsManager, dbscanResult, clusterCount);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		extraOptions = Utils.isExtraOptions();

		globalSettings = SettingsManager.loadSettings();
		inputSettings = globalSettings.getOPTICSSettings();

		IJ.showStatus("");

		if ("dbscan".equals(arg))
		{
			runDBSCAN();
		}
		else
		{
			runOPTICS();
		}

		IJ.showStatus(TITLE + " finished");

		// Update the settings
		SettingsManager.saveSettings(globalSettings);
	}

	private void runOPTICS()
	{
		TITLE = TITLE_OPTICS;

		// Create the working threads, connected in a chain
		ArrayList<Worker> workers = new ArrayList<Worker>();
		add(workers, new InputWorker());
		add(workers, new OpticsWorker());
		add(workers, new OpticsClusterWorker());
		// The following could operate in parallel
		//		add(workers, new ResultsWorker());
		//		add(workers, new MemoryResultsWorker());
		//		add(workers, new ReachabilityResultsWorker());
		//		add(workers, new ImageResultsWorker());
		Worker previous = workers.get(workers.size() - 1);
		add(workers, new ResultsWorker(), previous);
		add(workers, new RandIndexWorker(), previous);
		add(workers, new MemoryResultsWorker(), previous);
		add(workers, new ReachabilityResultsWorker(), previous);
		add(workers, new ImageResultsWorker(), previous);

		ArrayList<Thread> threads = startWorkers(workers);

		boolean cancelled = !showDialog(false);

		finishWorkers(workers, threads, cancelled);
	}

	private ArrayList<Thread> startWorkers(ArrayList<Worker> workers)
	{
		ArrayList<Thread> threads = new ArrayList<Thread>();
		for (Worker w : workers)
		{
			Thread t = new Thread(w);
			t.setDaemon(true);
			t.start();
			threads.add(t);
		}
		return threads;
	}

	private void finishWorkers(ArrayList<Worker> workers, ArrayList<Thread> threads, boolean cancelled)
	{
		// Finish work
		for (int i = 0; i < threads.size(); i++)
		{
			Thread t = threads.get(i);
			Worker w = workers.get(i);

			if (cancelled)
			{
				// Stop immediately any running worker
				try
				{
					t.interrupt();
				}
				catch (SecurityException e)
				{
					// We should have permission to interrupt this thread.
					e.printStackTrace();
				}
			}
			else
			{
				// Stop after the current work in the inbox
				w.running = false;

				// Notify a workers waiting on the inbox.
				// Q. How to check if the worker is sleeping?
				synchronized (w.inbox)
				{
					w.inbox.notify();
				}

				// Leave to finish their current work
				try
				{
					t.join(0);
				}
				catch (InterruptedException e)
				{
				}
			}
		}
	}

	/**
	 * Gets the nm per pixel.
	 *
	 * @param results
	 *            the results
	 * @return the nm per pixel
	 */
	public double getNmPerPixel(MemoryPeakResults results)
	{
		if (results.getCalibration() != null && results.getCalibration().getNmPerPixel() > 0)
			return results.getCalibration().getNmPerPixel();
		return 1;
	}

	/**
	 * Adds the worker. Connect the inbox to the previous worker outbox, or the primary input.
	 *
	 * @param workers
	 *            the workers
	 * @param worker
	 *            the worker
	 */
	private void add(ArrayList<Worker> workers, Worker worker)
	{
		if (workers.isEmpty())
		{
			add(workers, worker, null);
		}
		else
		{
			// Chain together
			Worker previous = workers.get(workers.size() - 1);
			add(workers, worker, previous);
		}
	}

	/**
	 * Adds the worker. Connect the inbox to the previous worker outbox, or the primary input.
	 *
	 * @param workers
	 *            the workers
	 * @param worker
	 *            the worker
	 * @param previous
	 *            the previous worker from which to take work
	 */
	private void add(ArrayList<Worker> workers, Worker worker, Worker previous)
	{
		if (previous == null)
		{
			// Take the primary input
			worker.inbox = inputStack;
		}
		else
		{
			// Chain together
			int size = (previous.outbox == null) ? 0 : previous.outbox.length;
			if (size == 0)
				previous.outbox = new WorkStack[1];
			else
				previous.outbox = Arrays.copyOf(previous.outbox, size + 1);
			previous.outbox[size] = new WorkStack();
			worker.inbox = previous.outbox[size];
		}
		workers.add(worker);
	}

	private Object[] imageModeArray;
	private Object[] outlineModeArray;

	private boolean showDialog(boolean isDBSCAN)
	{
		logReferences(isDBSCAN);

		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputSettings.inputOption, InputSource.MEMORY);

		//globalSettings = SettingsManager.loadSettings();
		//settings = globalSettings.getClusteringSettings();

		gd.addMessage("--- " + TITLE + " ---");
		gd.addNumericField("Min_points", inputSettings.minPoints, 0);
		if (isDBSCAN)
		{
			// Add fields to auto-compute the clustering distance from the K-nearest neighbour distance profile
			gd.addSlider("Noise (%)", 0, 50, inputSettings.fractionNoise * 100);
			gd.addNumericField("Samples", inputSettings.samples, 0);
			gd.addSlider("Sample_fraction (%)", 0, 15, inputSettings.sampleFraction * 100);
			gd.addNumericField("Clustering_distance", inputSettings.clusteringDistance, 2, 6, "nm");
		}
		else
		{
			String[] opticsModes = SettingsManager.getNames((Object[]) OPTICSMode.values());
			gd.addChoice("OPTICS_mode", opticsModes, opticsModes[inputSettings.getOPTICSModeOridinal()]);
			gd.addNumericField("Number_of_splits", inputSettings.numberOfSplitSets, 0);
			if (extraOptions)
			{
				gd.addCheckbox("Random_vectors", inputSettings.useRandomVectors);
				gd.addCheckbox("Approx_sets", inputSettings.saveApproximateSets);
				String[] sampleModes = SettingsManager.getNames((Object[]) SampleMode.values());
				gd.addChoice("Sample_mode", sampleModes, sampleModes[inputSettings.getSampleModeOridinal()]);
			}
			gd.addNumericField("Generating_distance", inputSettings.generatingDistance, 2, 6, "nm");
		}
		gd.addMessage("--- Clustering ---");
		if (isDBSCAN)
		{
			gd.addCheckbox("Core_points", inputSettings.core);
		}
		else
		{
			String[] clusteringModes = SettingsManager.getNames((Object[]) ClusteringMode.values());
			gd.addChoice("Clustering_mode", clusteringModes,
					clusteringModes[inputSettings.getClusteringModeOridinal()]);
			gd.addMessage(ClusteringMode.XI.toString() + " options:\n" + ClusteringMode.XI.toString() +
					" controls the change in reachability (profile steepness) to define a cluster");
			gd.addNumericField("Xi", inputSettings.xi, 4);
			gd.addCheckbox("Top_clusters", inputSettings.topLevel);
			gd.addMessage(ClusteringMode.DBSCAN.toString() + " options:");
			gd.addNumericField("Clustering_distance", inputSettings.clusteringDistance, 4);
			gd.addCheckbox("Core_points", inputSettings.core);
		}
		gd.addMessage("--- Image ---");
		gd.addSlider("Image_Scale", 0, 15, inputSettings.imageScale);
		TreeSet<ImageMode> imageModeSet = new TreeSet<ImageMode>();
		imageModeSet.addAll(Arrays.asList(ImageMode.values()));
		if (isDBSCAN)
		{
			imageModeSet.remove(ImageMode.CLUSTER_DEPTH);
			imageModeSet.remove(ImageMode.CLUSTER_ORDER);
		}
		imageModeArray = imageModeSet.toArray();
		String[] imageModes = SettingsManager.getNames(imageModeArray);
		gd.addChoice("Image_mode", imageModes, imageModes[inputSettings.getImageModeOridinal()]);

		gd.addCheckboxGroup(1, 2, new String[] { "Weighted", "Equalised" },
				new boolean[] { inputSettings.weighted, inputSettings.equalised }, new String[] { "Image" });
		//gd.addCheckbox("Weighted", inputSettings.weighted);
		//gd.addCheckbox("Equalised", inputSettings.equalised);

		TreeSet<OutlineMode> outlineModeSet = new TreeSet<OutlineMode>();
		outlineModeSet.addAll(Arrays.asList(OutlineMode.values()));
		if (isDBSCAN)
		{
			outlineModeSet.remove(OutlineMode.COLOURED_BY_DEPTH);
		}
		outlineModeArray = outlineModeSet.toArray();
		String[] outlineModes = SettingsManager.getNames(outlineModeArray);
		gd.addChoice("Outline", outlineModes, outlineModes[inputSettings.getOutlineModeOridinal()]);

		if (!isDBSCAN)
		{
			String[] spanningTreeModes = SettingsManager.getNames((Object[]) SpanningTreeMode.values());
			gd.addChoice("Spanning_tree", spanningTreeModes,
					spanningTreeModes[inputSettings.getSpanningTreeModeOridinal()]);

			gd.addMessage("--- Reachability Plot ---");
			String[] plotModes = SettingsManager.getNames((Object[]) PlotMode.values());
			gd.addChoice("Plot_mode", plotModes, plotModes[inputSettings.getPlotModeOridinal()]);
		}

		// Start disabled so the user can choose settings to update
		gd.addCheckbox("Preview", false);

		if (extraOptions)
			gd.addCheckbox("Debug", false);

		// Everything is done within the dialog listener
		if (isDBSCAN)
			gd.addDialogListener(new DBSCANDialogListener());
		else
			gd.addDialogListener(new OPTICSDialogListener());

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		if (!isPreview)
		{
			// The dialog was OK'd so run if work was stashed in the input stack.
			inputStack.addWork(inputStack.work);
		}

		return true;
	}

	private static byte logged = 0;
	private static final byte LOG_DBSCAN = 0x01;
	private static final byte LOG_OPTICS = 0x02;

	private static void logReferences(boolean isDBSCAN)
	{
		int width = 80;
		StringBuilder sb = new StringBuilder();
		if (isDBSCAN && (logged & LOG_DBSCAN) != LOG_DBSCAN)
		{
			logged |= LOG_DBSCAN;
			sb.append("DBSCAN: ");
			sb.append(TextUtils.wrap(
					"Ester, et al (1996). 'A density-based algorithm for discovering clusters in large spatial databases with noise'. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD-96). AAAI Press. pp. 226–231.",
					width)).append('\n');
			;
		}
		else if ((logged & LOG_OPTICS) != LOG_OPTICS)
		{
			logged |= LOG_OPTICS;
			sb.append("OPTICS: ");
			sb.append(TextUtils.wrap(
					"Kriegel, et al (2011). 'Density-based clustering'. Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery. 1 (3): 231–240.",
					width)).append('\n');
			;
			sb.append("FastOPTICS: ");
			sb.append(TextUtils.wrap(
					"Schneider, et al (2013). 'Fast parameterless density-based clustering via random projections'. 22nd ACM International Conference on Information and Knowledge Management(CIKM). ACM.",
					width)).append('\n');
			;
		}
		if (sb.length() > 0)
			IJ.log(sb.toString());
	}

	private class OPTICSDialogListener implements DialogListener
	{
		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			//		if (e == null)
			//		{
			//			// This happens when the dialog is first shown and can be ignored.
			//			// It also happens when called from within a macro. In this case we should run.
			//			if (!Utils.isMacro())
			//				return true;
			//		}

			if (debug)
				System.out.println("dialogItemChanged: " + e);

			// A previous run may have been cancelled so we have to handle this.
			if (Utils.isInterrupted())
			{
				if (Utils.isMacro())
					return true;

				// Q. Should we ask if the user wants to restart?
				IJ.resetEscape();
			}

			inputSettings.inputOption = ResultsManager.getInputSource(gd);

			// Load the results
			MemoryPeakResults results = ResultsManager.loadInputResults(inputSettings.inputOption, true);
			if (results == null || results.size() == 0)
			{
				IJ.error(TITLE, "No results could be loaded");
				return false;
			}

			inputSettings.minPoints = (int) gd.getNextNumber();
			inputSettings.setOPTICSMode(gd.getNextChoiceIndex());
			inputSettings.numberOfSplitSets = (int) gd.getNextNumber();
			if (extraOptions)
			{
				inputSettings.useRandomVectors = gd.getNextBoolean();
				inputSettings.saveApproximateSets = gd.getNextBoolean();
				inputSettings.setSampleMode(gd.getNextChoiceIndex());
			}
			inputSettings.generatingDistance = gd.getNextNumber();
			inputSettings.setClusteringMode(gd.getNextChoiceIndex());
			inputSettings.xi = gd.getNextNumber();
			inputSettings.topLevel = gd.getNextBoolean();
			inputSettings.clusteringDistance = gd.getNextNumber();
			inputSettings.core = gd.getNextBoolean();
			inputSettings.imageScale = gd.getNextNumber();
			inputSettings.setImageMode((ImageMode) imageModeArray[gd.getNextChoiceIndex()]);
			inputSettings.weighted = gd.getNextBoolean();
			inputSettings.equalised = gd.getNextBoolean();
			inputSettings.setOutlineMode((OutlineMode) outlineModeArray[gd.getNextChoiceIndex()]);
			inputSettings.setSpanningTreeMode(gd.getNextChoiceIndex());
			inputSettings.setPlotMode(gd.getNextChoiceIndex());
			boolean preview = gd.getNextBoolean();
			if (extraOptions)
				debug = gd.getNextBoolean();

			if (gd.invalidNumber())
				return false;

			// Check arguments
			try
			{
				Parameters.isAboveZero("Xi", inputSettings.xi);
				Parameters.isBelow("Xi", inputSettings.xi, 1);
			}
			catch (IllegalArgumentException ex)
			{
				Utils.log(TITLE + ": " + ex.getMessage());
				return false;
			}

			// Clone so that the work queue has it's own unique reference
			Work work = new Work(inputSettings.clone(), results);
			if (preview)
			{
				// Queue the settings
				if (debug)
					System.out.println("Adding work");
				if (isPreview)
					// Use a delay next time. This prevents delay when the preview is first switched on. 
					work.time = System.currentTimeMillis() + DELAY;
				else
					isPreview = true;
				inputStack.addWork(work);
			}
			else
			{
				// Preview is off
				isPreview = false;
				// Stash the work (this does not notify the input worker)
				inputStack.setWork(work);
			}

			return true;
		}
	}

	private class DBSCANDialogListener implements DialogListener
	{
		public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
		{
			//		if (e == null)
			//		{
			//			// This happens when the dialog is first shown and can be ignored.
			//			// It also happens when called from within a macro. In this case we should run.
			//			if (!Utils.isMacro())
			//				return true;
			//		}

			if (debug)
				System.out.println("dialogItemChanged: " + e);

			// A previous run may have been cancelled so we have to handle this.
			if (Utils.isInterrupted())
			{
				if (Utils.isMacro())
					return true;

				// Q. Should we ask if the user wants to restart?
				IJ.resetEscape();
			}

			inputSettings.inputOption = ResultsManager.getInputSource(gd);

			// Load the results
			MemoryPeakResults results = ResultsManager.loadInputResults(inputSettings.inputOption, true);
			if (results == null || results.size() == 0)
			{
				IJ.error(TITLE, "No results could be loaded");
				return false;
			}

			inputSettings.minPoints = (int) gd.getNextNumber();
			inputSettings.fractionNoise = gd.getNextNumber() / 100;
			inputSettings.samples = (int) gd.getNextNumber();
			inputSettings.sampleFraction = gd.getNextNumber() / 100;
			inputSettings.clusteringDistance = gd.getNextNumber();
			inputSettings.core = gd.getNextBoolean();
			inputSettings.imageScale = gd.getNextNumber();
			inputSettings.setImageMode((ImageMode) imageModeArray[gd.getNextChoiceIndex()]);
			inputSettings.weighted = gd.getNextBoolean();
			inputSettings.equalised = gd.getNextBoolean();
			inputSettings.setOutlineMode((OutlineMode) outlineModeArray[gd.getNextChoiceIndex()]);
			boolean preview = gd.getNextBoolean();
			if (extraOptions)
				debug = gd.getNextBoolean();

			if (gd.invalidNumber())
				return false;

			// Clone so that the work queue has it's own unique reference
			Work work = new Work(inputSettings.clone(), results);
			if (preview)
			{
				// Queue the settings
				if (debug)
					System.out.println("Adding work");
				if (isPreview)
					// Use a delay next time. This prevents delay when the preview is first switched on. 
					work.time = System.currentTimeMillis() + DELAY;
				else
					isPreview = true;
				inputStack.addWork(work);
			}
			else
			{
				// Preview is off
				isPreview = false;
				// Stash the work (this does not notify the input worker)
				inputStack.setWork(work);
			}

			return true;
		}
	}

	private void runDBSCAN()
	{
		TITLE = TITLE_DBSCAN;

		// Create the working threads, connected in a chain
		ArrayList<Worker> workers = new ArrayList<Worker>();
		add(workers, new InputWorker());
		add(workers, new KNNWorker());
		add(workers, new DBSCANWorker());
		add(workers, new DBSCANClusterWorker());
		// The following could operate in parallel
		Worker previous = workers.get(workers.size() - 1);
		add(workers, new ResultsWorker(), previous);
		add(workers, new RandIndexWorker(), previous);
		add(workers, new MemoryResultsWorker(), previous);
		add(workers, new ImageResultsWorker(), previous);

		ArrayList<Thread> threads = startWorkers(workers);

		boolean cancelled = !showDialog(true);

		finishWorkers(workers, threads, cancelled);
	}
}
