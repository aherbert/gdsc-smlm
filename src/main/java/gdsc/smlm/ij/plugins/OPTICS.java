package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;

import gdsc.core.clustering.optics.OPTICSCluster;
import gdsc.core.clustering.optics.OPTICSManager;
import gdsc.core.clustering.optics.OPTICSResult;
import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.utils.ConvexHull;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Settings;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.results.IJImagePeakResults;
import gdsc.smlm.ij.settings.ClusteringSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
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
public class OPTICS implements PlugIn, DialogListener
{
	private String TITLE = "OPTICS";

	private enum ImageMode
	{
		//@formatter:off
		CLUSTER_ID {
			@Override
			public String getName() { return "Cluster Id"; };
			@Override
			float getValue(float value, int clusterId) { return clusterId; }
		},
		VALUE {
			@Override
			public String getName() { return "Value"; };
			@Override
			boolean canBeWeighted() { return true; }
			@Override
			float getValue(float value, int clusterId) { return value; }
		},
		COUNT {
			@Override
			public String getName() { return "Count"; };
			@Override
			boolean canBeWeighted() { return true; }
			@Override
			float getValue(float value, int clusterId) { return 1f; }
		},
		NONE {
			@Override
			public String getName() { return "None"; };
			@Override
			float getValue(float value, int clusterId) { return 0; }
		};
		//@formatter:on

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();

		/**
		 * Return the value to draw
		 * 
		 * @param value
		 *            The value of the cluster point
		 * @param clusterId
		 *            The cluster Id of the cluster point
		 * @return The value
		 */
		abstract float getValue(float value, int clusterId);

		/**
		 * Can be weighted.
		 *
		 * @return true, if successful
		 */
		boolean canBeWeighted()
		{
			return false;
		}

		@Override
		public String toString()
		{
			return getName();
		}
	}

	private enum PlotMode
	{
		//@formatter:off
		ON {
			@Override
			public String getName() { return "On"; };
		},
		WITH_CLUSTERS {
			@Override
			public String getName() { return "With clusters"; };
			@Override
			boolean isDrawClusters() { return true; }
		},
		HIGHLIGHTED {
			@Override
			public String getName() { return "Highlighted"; };
			@Override
			boolean isHighlightProfile() { return true; }
		},
		HIGHLIGHTED_WITH_CLUSTERS {
			@Override
			public String getName() { return "Highlighted with clusters"; };
			@Override
			boolean isHighlightProfile() { return true; }
			@Override
			boolean isDrawClusters() { return true; }
		},
		COLOURED_WITH_CLUSTERS {
			@Override
			public String getName() { return "Coloured with clusters"; };
			@Override
			boolean isHighlightProfile() { return true; }
			@Override
			boolean isColourProfile() { return true; }
			@Override
			boolean isDrawClusters() { return true; }
		},
		OFF {
			@Override
			public String getName() { return "Off"; };
		};
		//@formatter:on

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();

		/**
		 * @return True if the profile should be highlighted for top-cluster regions
		 */
		boolean isHighlightProfile()
		{
			return false;
		}
		
		/**
		 * @return True if the profile should be coloured using the cluster colour for top-cluster regions
		 */
		boolean isColourProfile()
		{
			return false;
		}

		/**
		 * @return If clusters should be drawn on the plot
		 */
		boolean isDrawClusters()
		{
			return false;
		}
		
		/**
		 * @return True if the clusters are needed
		 */
		boolean requiresClusters()
		{
			return isDrawClusters() || isHighlightProfile() || isColourProfile();
		}

		@Override
		public String toString()
		{
			return getName();
		}
	}

	private static class InputSettings
	{
		// Affect creating the OPTICS manager
		String inputOption = "";

		// Affect running OPTICS
		double generatingDistance = 0;
		int minPoints = 5;

		// Affect running OPTICS Xi
		double xi = 0.03;
		boolean topLevel = false;

		// Affect display of results
		double imageScale = 2;
		ImageMode imageMode = ImageMode.CLUSTER_ID;
		boolean weighted = false;
		boolean equalised = false;

		PlotMode plotMode = PlotMode.COLOURED_WITH_CLUSTERS;

		boolean outline = true;

		public int getImageMode()
		{
			if (imageMode == null)
				return 0;
			return imageMode.ordinal();
		}

		public void setImageMode(int mode)
		{
			ImageMode[] values = ImageMode.values();
			if (mode < 0 || mode >= values.length)
				mode = 0;
			this.imageMode = values[mode];
		}

		public int getPlotMode()
		{
			if (plotMode == null)
				return 0;
			return plotMode.ordinal();
		}

		public void setPlotMode(int mode)
		{
			PlotMode[] values = PlotMode.values();
			if (mode < 0 || mode >= values.length)
				mode = 0;
			this.plotMode = values[mode];
		}
	}

	private static InputSettings inputSettings = new InputSettings();

	// TODO - store the settings between sessions
	private GlobalSettings globalSettings;
	private ClusteringSettings settings;

	private boolean extraOptions;

	private static class Work
	{
		long time = 0;
		InputSettings inputSettings;
		Settings settings;

		public Work(InputSettings inputSettings, Settings settings)
		{
			this.inputSettings = inputSettings;
			this.settings = settings;
		}

		public Work(InputSettings inputSettings, Object... settings)
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
					if (!running)
					{
						debug(" Not running, stopping");
						break;
					}

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
		abstract boolean equals(InputSettings current, InputSettings previous);

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
		boolean equals(InputSettings current, InputSettings previous)
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
			synchronized (results)
			{
				int size = results.size();
				x = new float[size];
				y = new float[size];
				int i = 0;
				for (PeakResult p : results)
				{
					x[i] = p.getXPosition();
					y[i] = p.getYPosition();
					i++;
				}
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
		boolean equals(InputSettings current, InputSettings previous)
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
			OPTICSResult opticsResult = opticsManager.optics((float) work.inputSettings.generatingDistance,
					work.inputSettings.minPoints);
			// It may be null if cancelled. However return null Work will close down the next thread
			return new Work(work.inputSettings, results, opticsManager, opticsResult);
		}
	}

	private class OpticsClusterWorker extends Worker
	{
		int clusterCount = 0;

		@Override
		boolean equals(InputSettings current, InputSettings previous)
		{
			if (current.xi != previous.xi)
				return false;
			if (current.topLevel != previous.topLevel)
				return false;
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
				int options = (work.inputSettings.topLevel) ? OPTICSResult.XI_OPTION_TOP_LEVEL : 0;
				synchronized (opticsResult)
				{
					opticsResult.extractClusters(work.inputSettings.xi, options);
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
		boolean equals(InputSettings current, InputSettings previous)
		{
			// Only depends on if the clustering results are new. This is triggered 
			// in the default comparison of the Settings object.
			return true;
		}

		@Override
		Work createResult(Work work)
		{
			OPTICSResult opticsResult = (OPTICSResult) work.settings.get(2);
			// It may be null if cancelled.
			if (opticsResult == null)
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
		boolean equals(InputSettings current, InputSettings previous)
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
			OPTICSResult opticsResult = (OPTICSResult) work.settings.get(2);
			// It may be null if cancelled.
			if (opticsResult != null)
			{
				int[] clusters;
				synchronized (opticsResult)
				{
					clusters = opticsResult.getClusters();
				}
				int max = Maths.max(clusters);

				// Save the clusters to memory
				Trace[] traces = new Trace[max + 1];
				for (int i = 0; i <= max; i++)
				{
					traces[i] = new Trace();
					traces[i].setId(i);
				}
				int i = 0;
				synchronized (results)
				{
					for (PeakResult r : results.getResults())
					{
						traces[clusters[i++]].add(r);
					}
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
		boolean equals(InputSettings current, InputSettings previous)
		{
			if (current.plotMode != previous.plotMode)
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

			// Draw the reachability profile
			PlotMode mode = work.inputSettings.plotMode;
			if (mode != PlotMode.OFF)
			{
				double[] profile;
				synchronized (opticsResult)
				{
					profile = opticsResult.getReachabilityDistanceProfile(true);
				}
				double[] order = Utils.newArray(profile.length, 1.0, 1.0);
				String title = TITLE + " Reachability Distance";
				Plot plot = new Plot(title, "Order", "Reachability");
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
		boolean equals(InputSettings current, InputSettings previous)
		{
			if (current.imageScale != previous.imageScale)
			{
				// Clear all the cached results
				newResults();
				return false;
			}
			if (current.outline != previous.outline)
				return false;
			if (getDisplayFlags(current) != getDisplayFlags(previous))
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
			OPTICSResult opticsResult = (OPTICSResult) work.settings.get(2);
			int clusterCount = (Integer) work.settings.get(3);
			// It may be null if cancelled.
			if (opticsResult == null)
			{
				image = null;
				o = null;
				return new Work(work.inputSettings, results, opticsManager, opticsResult, clusterCount, image);
			}

			int[] clusters = null;

			if (work.inputSettings.imageScale > 0)
			{
				if (image == null)
				{
					synchronized (opticsResult)
					{
						clusters = opticsResult.getClusters();
					}

					// Display the results ...

					// TODO: Options to not draw the points

					Rectangle bounds = results.getBounds();
					image = new IJImagePeakResults(results.getName() + " " + TITLE, bounds,
							(float) work.inputSettings.imageScale);
					// TODO - options to control rendering
					ImageMode mode = work.inputSettings.imageMode;
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
						synchronized (results)
						{
							int i = 0;
							x = new float[results.size()];
							y = new float[x.length];
							v = new float[x.length];
							for (PeakResult r : results.getResults())
							{
								x[i] = r.getXPosition();
								y[i] = r.getYPosition();
								v[i] = mode.getValue(r.getSignal(), clusters[i]);
								i++;
							}
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
							synchronized (opticsResult)
							{
								clusters = opticsResult.getClusters();
							}
						}
						int max = Maths.max(clusters);

						// The current overlay may be fine. 
						// If the result has cached convex hulls then assume we have constructed an overlay
						// for these results.
						if (!opticsResult.hasConvexHulls() || o == null)
						{
							// We need to recompute
							o = new Overlay();
							ConvexHull[] hulls = new ConvexHull[max + 1];
							synchronized (opticsResult)
							{
								opticsResult.computeConvexHulls();
								clusters = opticsResult.getClusters();
								for (int c = 1; c <= max; c++)
								{
									hulls[c] = opticsResult.getConvexHull(c);
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
										y[i] = image.mapX(y[i]);
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

			return new Work(work.inputSettings, results, opticsManager, opticsResult, clusterCount, image);
		}

		private int getDisplayFlags(InputSettings inputSettings)
		{
			int displayFlags = 0;
			if (inputSettings.imageMode.canBeWeighted())
			{
				if (inputSettings.weighted)
					displayFlags |= IJImagePeakResults.DISPLAY_WEIGHTED;
				if (inputSettings.equalised)
					displayFlags |= IJImagePeakResults.DISPLAY_EQUALIZED;
			}
			if (inputSettings.imageMode == ImageMode.CLUSTER_ID)
				displayFlags = IJImagePeakResults.DISPLAY_REPLACE;
			return displayFlags;
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

		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		extraOptions = Utils.isExtraOptions();

		// Create the working threads, connected in a chain
		ArrayList<Worker> workers = new ArrayList<Worker>();
		add(workers, new InputWorker());
		add(workers, new OpticsWorker());
		add(workers, new OpticsClusterWorker());
		add(workers, new ResultsWorker());
		add(workers, new MemoryResultsWorker());
		add(workers, new ReachabilityResultsWorker());
		add(workers, new ImageResultsWorker());

		ArrayList<Thread> threads = new ArrayList<Thread>();

		for (Worker w : workers)
		{
			Thread t = new Thread(w);
			t.setDaemon(true);
			t.start();
			threads.add(t);
		}

		boolean cancelled = !showDialog();

		// Finish work
		for (int i = 0; i < threads.size(); i++)
		{
			Thread t = threads.get(i);
			Worker w = workers.get(i);

			// Wait for threads to end
			if (cancelled)
			{
				w.running = false;
				w.inbox.close();

				// Stop immediately
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
				// Leave to finish
				try
				{
					t.join(0);
				}
				catch (InterruptedException e)
				{
				}

				w.running = false;
				w.inbox.close();
			}
		}

		IJ.showStatus(TITLE + " finished");
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

	private boolean showDialog()
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		ResultsManager.addInput(gd, inputSettings.inputOption, InputSource.MEMORY);

		//globalSettings = SettingsManager.loadSettings();
		//settings = globalSettings.getClusteringSettings();

		gd.addNumericField("Generating_distance", inputSettings.generatingDistance, 2);
		gd.addNumericField("Min_points", inputSettings.minPoints, 0);
		gd.addNumericField("Xi", inputSettings.xi, 4);
		gd.addCheckbox("Top_clusters", inputSettings.topLevel);
		gd.addSlider("Image_Scale", 0, 15, inputSettings.imageScale);
		String[] imageModes = SettingsManager.getNames((Object[]) ImageMode.values());
		gd.addChoice("Image_mode", imageModes, imageModes[inputSettings.getImageMode()]);
		gd.addCheckbox("Weighted", inputSettings.weighted);
		gd.addCheckbox("Equalised", inputSettings.equalised);
		gd.addCheckbox("Outline", inputSettings.outline);
		String[] plotModes = SettingsManager.getNames((Object[]) PlotMode.values());
		gd.addChoice("Plot_mode", plotModes, plotModes[inputSettings.getPlotMode()]);

		// Start disabled so the user can choose settings to update
		gd.addCheckbox("Preview", false);

		// Everything is done within the dialog listener
		gd.addDialogListener(this);

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

	private boolean isPreview = false;

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

		InputSettings settings = new InputSettings();

		settings.inputOption = ResultsManager.getInputSource(gd);

		// Load the results
		MemoryPeakResults results = ResultsManager.loadInputResults(settings.inputOption, true);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return false;
		}

		settings.generatingDistance = gd.getNextNumber();
		settings.minPoints = (int) gd.getNextNumber();
		settings.xi = gd.getNextNumber();
		settings.topLevel = gd.getNextBoolean();
		settings.imageScale = gd.getNextNumber();
		settings.setImageMode(gd.getNextChoiceIndex());
		settings.weighted = gd.getNextBoolean();
		settings.equalised = gd.getNextBoolean();
		settings.outline = gd.getNextBoolean();
		settings.setPlotMode(gd.getNextChoiceIndex());
		boolean preview = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;

		// Check arguments
		try
		{
			Parameters.isAboveZero("Xi", settings.xi);
			Parameters.isBelow("Xi", settings.xi, 1);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		// Store for next time
		inputSettings = settings;

		Work work = new Work(settings, results);
		if (preview)
		{
			// Queue the settings
			if (extraOptions)
				System.out.println("Adding work");
			if (isPreview)
				// Use a delay next time. This prevents delay when the preview is first switched on. 
				work.time = System.currentTimeMillis() + 100;
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
