package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;

import gdsc.core.clustering.optics.ClusteringResult;
import gdsc.core.clustering.optics.DBSCANResult;
import gdsc.core.clustering.optics.OPTICSCluster;
import gdsc.core.clustering.optics.OPTICSManager;
import gdsc.core.clustering.optics.OPTICSResult;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.utils.ConvexHull;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Settings;
import gdsc.core.utils.Sort;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.OPTICSSettings;
import gdsc.smlm.ij.settings.OPTICSSettings.ClusteringMode;
import gdsc.smlm.ij.settings.OPTICSSettings.ImageMode;
import gdsc.smlm.ij.settings.OPTICSSettings.PlotMode;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.process.LUTHelper;
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

	private String TITLE;
	/**
	 * Delay (in milliseconds) used when entering new values in the dialog before the preview is processed
	 */
	private long DELAY = 500;
	private boolean isPreview = false;

	private GlobalSettings globalSettings;
	private OPTICSSettings inputSettings;

	private boolean extraOptions;

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
		private WorkStack inbox, outbox;

		public Worker()
		{
			this(null, null);
		}

		public Worker(WorkStack inbox, WorkStack outbox)
		{
			this.inbox = inbox;
			this.outbox = outbox;
		}

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
						outbox.addWork(result);
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
			if (extraOptions)
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
			return new Work(work.inputSettings, results, opticsManager);
		}
	}

	private class OpticsWorker extends Worker
	{
		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.generatingDistance != previous.generatingDistance)
				return false;
			if (current.minPoints != previous.minPoints)
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

			double distance = work.inputSettings.generatingDistance;
			int minPts = work.inputSettings.minPoints;
			if (distance > 0)
			{
				// Convert generating distance to pixels
				double nmPerPixel = getNmPerPixel(results);
				if (nmPerPixel != 1)
				{
					double newDistance = distance / nmPerPixel;
					Utils.log(TITLE + ": Converting generating distance %s nm to %s pixels", Utils.rounded(distance),
							Utils.rounded(newDistance));
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

			OPTICSResult opticsResult = opticsManager.optics((float) distance, minPts);
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
			if (current.clusteringMode != previous.clusteringMode)
				return false;
			if (current.clusteringMode == ClusteringMode.XI)
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
					if (work.inputSettings.clusteringMode == ClusteringMode.XI)
					{
						int options = (work.inputSettings.topLevel) ? OPTICSResult.XI_OPTION_TOP_LEVEL : 0;
						opticsResult.extractClusters(work.inputSettings.xi, options);
					}
					else
					{
						double nmPerPixel = getNmPerPixel(results);

						// Ensure that the distance is valid
						double distance = opticsResult.generatingDistance * nmPerPixel;
						if (work.inputSettings.clusteringDistance > 0)
							distance = Math.min(work.inputSettings.clusteringDistance, distance);

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

				ArrayList<OPTICSCluster> clusters = null;
				LUT lut = Utils.getColorModel();
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

				// Draw the clusters using lines underneath
				if (mode.isDrawClusters())
				{
					// TODO - Get a good distance to start the lines, and the separation
					double start = -limits[0];
					double separation = limits[0];

					for (OPTICSCluster cluster : clusters)
					{
						int level = cluster.getLevel();
						double y = start - (maxLevel - level) * separation;
						Color color = LUTHelper.getColour(lut, cluster.clusterId, 0f, maxClusterId);
						plot.setColor(color);
						plot.drawLine(cluster.start, y, cluster.end, y);
					}

					// Update the limits if we are plotting lines underneath for the clusters
					limits[0] = start - (maxLevel + 1) * separation;
				}
				else
				{
					// Just plot to zero
					limits[0] = 0;
				}

				plot.setLimits(1, order.length, limits[0], limits[1] * 1.05);

				plot.setColor(Color.black);
				plot.addPoints(order, profile, Plot.LINE);

				// Colour the reachability plot line if it is in a cluster. Use a default colour 
				plot.setColor(Color.blue);
				if (mode.isHighlightProfile())
				{
					for (OPTICSCluster cluster : clusters)
					{
						// Only do top level clusters
						if (cluster.getLevel() != 0)
							continue;

						if (mode.isColourProfile())
							plot.setColor(LUTHelper.getColour(lut, cluster.clusterId, 0f, maxClusterId));
						int from = cluster.start;
						int to = cluster.end + 1;
						double[] order1 = Arrays.copyOfRange(order, from, to);
						double[] profile1 = Arrays.copyOfRange(profile, from, to);
						plot.addPoints(order1, profile1, Plot.LINE);
					}
				}

				// Add the DBSCAN clustering distance
				if (inputSettings.clusteringMode == ClusteringMode.DBSCAN)
				{
					// Ensure that the distance is valid
					double distance = opticsResult.generatingDistance * nmPerPixel;
					if (work.inputSettings.clusteringDistance > 0)
						distance = Math.min(work.inputSettings.clusteringDistance, distance);
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

	private class ImageResultsWorker extends Worker
	{
		IJImagePeakResults image = null;
		Overlay o = null;

		@Override
		boolean equals(OPTICSSettings current, OPTICSSettings previous)
		{
			if (current.imageScale != previous.imageScale)
			{
				// Clear all the cached results
				newResults();
				return false;
			}
			if (current.outline != previous.outline)
				return false;
			if (current.getImageMode() != previous.getImageMode() ||
					getDisplayFlags(current) != getDisplayFlags(previous))
			{
				// We can only cache the image if the display mode is the same
				image = null;
				return false;
			}
			return true;
		}

		@Override
		void newResults()
		{
			// Clear cache
			image = null;
			o = null;
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
				image = null;
				o = null;
				return new Work(work.inputSettings, results, opticsManager, clusteringResult, clusterCount, image);
			}

			int[] clusters = null;

			if (work.inputSettings.imageScale > 0)
			{
				if (image == null)
				{
					synchronized (clusteringResult)
					{
						clusters = clusteringResult.getClusters();
					}

					// Display the results ...

					// TODO: Options to not draw the points

					Rectangle bounds = results.getBounds();
					image = new IJImagePeakResults(results.getName() + " " + TITLE, bounds,
							(float) work.inputSettings.imageScale);
					// TODO - options to control rendering
					ImageMode mode = work.inputSettings.getImageMode();
					image.setDisplayFlags(getDisplayFlags(work.inputSettings));
					image.setLiveImage(false);
					image.begin();
					ImagePlus imp = image.getImagePlus();
					imp.setOverlay(null);
					if (mode != ImageMode.NONE)
					{
						// Draw each cluster in a new colour. Set get an appropriate LUT.
						LUT lut = (mode == ImageMode.CLUSTER_ID) ? Utils.getColorModel()
								: LUTHelper.createLUT(LutColour.FIRE);
						image.getImagePlus().getProcessor().setColorModel(lut);

						// Add in a single batch
						float[] x, y, v;
						x = new float[results.size()];
						y = new float[x.length];
						v = new float[x.length];
						ArrayList<PeakResult> list = (ArrayList<PeakResult>) results.getResults();
						for (int i = 0, size = results.size(); i < size; i++)
						{
							PeakResult r = list.get(i);
							x[i] = r.getXPosition();
							y[i] = r.getYPosition();
							v[i] = mode.getValue(r.getSignal(), clusters[i]);
						}
						image.add(x, y, v);
					}
					image.end();
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
				if (work.inputSettings.outline)
				{
					if (o == null)
					{
						if (clusters == null)
						{
							synchronized (clusteringResult)
							{
								clusters = clusteringResult.getClusters();
							}
						}
						int max = Maths.max(clusters);

						// The current overlay may be fine. 
						// If the result has cached convex hulls then assume we have constructed an overlay
						// for these results.
						if (!clusteringResult.hasConvexHulls() || o == null)
						{
							// We need to recompute
							o = new Overlay();
							ConvexHull[] hulls = new ConvexHull[max + 1];
							synchronized (clusteringResult)
							{
								clusteringResult.computeConvexHulls();
								for (int c = 1; c <= max; c++)
								{
									hulls[c] = clusteringResult.getConvexHull(c);
								}
							}

							// Create a colour to match the LUT
							LUT lut = Utils.getColorModel();
							// Extract the ConvexHull of each cluster
							for (int c = 1; c <= max; c++)
							{
								ConvexHull hull = hulls[c];
								if (hull != null)
								{
									// Convert the Hull to the correct image scale.
									float[] x = hull.x.clone();
									float[] y = hull.y.clone();
									for (int i = 0; i < x.length; i++)
									{
										x[i] = image.mapX(x[i]);
										y[i] = image.mapY(y[i]);
									}
									PolygonRoi roi = new PolygonRoi(x, y, Roi.POLYGON);
									Color color = LUTHelper.getColour(lut, c, 0f, max);
									roi.setStrokeColor(color);
									// TODO: Options to set a fill colour?
									o.add(roi);
								}
							}
						}
					}
					imp.setOverlay(o);
				}
				else
				{
					imp.setOverlay(null);
				}
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
			if (current.samples != previous.samples)
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
					if (samples != -1)
						samples = Math.max(100, samples);
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

			DBSCANResult dbscanResult = opticsManager.dbscan((float) clusteringDistance, minPts);
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
		add(workers, new ResultsWorker());
		add(workers, new MemoryResultsWorker());
		add(workers, new ReachabilityResultsWorker());
		add(workers, new ImageResultsWorker());

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
			// Take the primary input
			worker.inbox = inputStack;
		}
		else
		{
			// Chain together
			Worker previous = workers.get(workers.size() - 1);
			previous.outbox = new WorkStack();
			worker.inbox = previous.outbox;
		}
		workers.add(worker);
	}

	private boolean showDialog(boolean isDBSCAN)
	{
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
			gd.addNumericField("Clustering_distance", inputSettings.clusteringDistance, 2, 6, "nm");
		}
		else
		{
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
			gd.addMessage(ClusteringMode.XI.toString() + " options:");
			gd.addMessage("Xi controls the change in reachability (profile steepness) to define a cluster");
			gd.addNumericField("Xi", inputSettings.xi, 4);
			gd.addCheckbox("Top_clusters", inputSettings.topLevel);
			gd.addMessage(ClusteringMode.DBSCAN.toString() + " options:");
			gd.addNumericField("Clustering_distance", inputSettings.clusteringDistance, 4);
			gd.addCheckbox("Core_points", inputSettings.core);
		}
		gd.addMessage("--- Image ---");
		gd.addSlider("Image_Scale", 0, 15, inputSettings.imageScale);
		String[] imageModes = SettingsManager.getNames((Object[]) ImageMode.values());
		gd.addChoice("Image_mode", imageModes, imageModes[inputSettings.getImageModeOridinal()]);
		gd.addCheckbox("Weighted", inputSettings.weighted);
		gd.addCheckbox("Equalised", inputSettings.equalised);
		gd.addCheckbox("Outline", inputSettings.outline);
		if (!isDBSCAN)
		{
			gd.addMessage("--- Reachability Plot ---");
			String[] plotModes = SettingsManager.getNames((Object[]) PlotMode.values());
			gd.addChoice("Plot_mode", plotModes, plotModes[inputSettings.getPlotModeOridinal()]);
		}

		// Start disabled so the user can choose settings to update
		gd.addCheckbox("Preview", false);

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

			if (extraOptions)
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
			inputSettings.generatingDistance = gd.getNextNumber();
			inputSettings.setClusteringMode(gd.getNextChoiceIndex());
			inputSettings.xi = gd.getNextNumber();
			inputSettings.topLevel = gd.getNextBoolean();
			inputSettings.clusteringDistance = gd.getNextNumber();
			inputSettings.core = gd.getNextBoolean();
			inputSettings.imageScale = gd.getNextNumber();
			inputSettings.setImageMode(gd.getNextChoiceIndex());
			inputSettings.weighted = gd.getNextBoolean();
			inputSettings.equalised = gd.getNextBoolean();
			inputSettings.outline = gd.getNextBoolean();
			inputSettings.setPlotMode(gd.getNextChoiceIndex());
			boolean preview = gd.getNextBoolean();

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
				IJ.error(TITLE, ex.getMessage());
				return false;
			}

			// Clone so that the work queue has it's own unique reference
			Work work = new Work(inputSettings.clone(), results);
			if (preview)
			{
				// Queue the settings
				if (extraOptions)
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

			if (extraOptions)
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
			inputSettings.clusteringDistance = gd.getNextNumber();
			inputSettings.core = gd.getNextBoolean();
			inputSettings.imageScale = gd.getNextNumber();
			inputSettings.setImageMode(gd.getNextChoiceIndex());
			inputSettings.weighted = gd.getNextBoolean();
			inputSettings.equalised = gd.getNextBoolean();
			inputSettings.outline = gd.getNextBoolean();
			boolean preview = gd.getNextBoolean();

			if (gd.invalidNumber())
				return false;

			// Clone so that the work queue has it's own unique reference
			Work work = new Work(inputSettings.clone(), results);
			if (preview)
			{
				// Queue the settings
				if (extraOptions)
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
		add(workers, new ResultsWorker());
		add(workers, new MemoryResultsWorker());
		add(workers, new ImageResultsWorker());

		ArrayList<Thread> threads = startWorkers(workers);

		boolean cancelled = !showDialog(true);

		finishWorkers(workers, threads, cancelled);
	}
}
