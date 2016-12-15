package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.ArrayList;

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
		boolean outline = true;
	}

	private static InputSettings inputSettings = new InputSettings();

	// TODO - store the settings between sessions
	private GlobalSettings globalSettings;
	private ClusteringSettings settings;

	private boolean extraOptions;

	private static class Work
	{
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
			this.inbox = null;
			this.outbox = null;
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
						debug(" Found work");
						work = inbox.getWork();
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
			return work.settings.equals(lastWork.settings);
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
			int size = results.size();
			float[] x = new float[size];
			float[] y = new float[size];
			int i = 0;
			for (PeakResult p : results)
			{
				x[i] = p.getXPosition();
				y[i] = p.getYPosition();
				i++;
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
				opticsResult.extractClusters(work.inputSettings.xi, options);
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
				int[] clusters = opticsResult.getClusters();
				int max = Maths.max(clusters);

				// Save the clusters to memory
				Trace[] traces = new Trace[max + 1];
				for (int i = 0; i <= max; i++)
				{
					traces[i] = new Trace();
					traces[i].setId(i);
				}
				int i = 0;
				for (PeakResult r : results.getResults())
				{
					traces[clusters[i++]].add(r);
				}
				TraceMolecules.saveResults(results, traces, "");
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
			// Only depends on if the OPTICS results are new. This is triggered 
			// in the default comparison of the Settings object.

			// TODO - add options to control the plot output

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
			double[] profile = opticsResult.getReachabilityDistanceProfile(true);
			double[] order = Utils.newArray(profile.length, 1.0, 1.0);
			String title = TITLE + " Reachability Distance";
			Plot plot = new Plot(title, "Order", "Reachability");
			double[] limits = Maths.limits(profile);
			// We should update the limits if we are plotting lines underneath for the clusters
			limits[0] = 0;
			plot.setLimits(1, order.length, limits[0], limits[1] * 1.05);
			plot.addPoints(order, profile, Plot.LINE);
			// TODO - Draw the clusters using lines underneath

			Utils.display(title, plot);

			// We have not created anything new so return the current object
			return work;
		}
	}

	private class ImageResultsWorker extends Worker
	{
		IJImagePeakResults image = null;

		@Override
		boolean equals(InputSettings current, InputSettings previous)
		{
			if (current.imageScale != previous.imageScale)
				return false;
			// TODO - options for other image output ...
			return true;
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
				return new Work(work.inputSettings, results, opticsManager, opticsResult, clusterCount, image);
			}

			int[] clusters = opticsResult.getClusters();

			if (work.inputSettings.imageScale > 0)
			{
				// Display the results ...

				// Working coordinates
				float[] x, y;
				{
					int i = 0;
					x = new float[results.size()];
					y = new float[results.size()];
					for (PeakResult r : results.getResults())
					{
						x[i] = r.getXPosition();
						y[i++] = r.getYPosition();
					}
				}

				Rectangle bounds = results.getBounds();
				image = new IJImagePeakResults(results.getName() + " " + TITLE, bounds,
						(float) work.inputSettings.imageScale);
				int displayFlags = IJImagePeakResults.DISPLAY_REPLACE; // | IJImagePeakResults.DISPLAY_WEIGHTED;
				image.setDisplayFlags(displayFlags);
				image.begin();
				// TODO - Draw each cluster in a new colour. Set get an appropriate LUT.
				for (int i = 0; i < x.length; i++)
				{
					// TODO: Options to not draw the points
					image.add(0, x[i], y[i], clusters[i]);
				}
				image.end();
			}

			return new Work(work.inputSettings, results, opticsManager, opticsResult, clusterCount, image);
		}
	}

	private class OverlayResultsWorker extends Worker
	{
		Overlay o = null;

		@Override
		boolean equals(InputSettings current, InputSettings previous)
		{
			if (current.imageScale != previous.imageScale)
			{
				// We can cache the overlay only if the image is the same size
				o = null;
				return false;
			}
			if (current.outline != previous.outline)
				return false;
			// TODO - options for other results output ...
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
			IJImagePeakResults image = (IJImagePeakResults) work.settings.get(4);
			if (image == null)
			{
				return work;
			}
			ImagePlus imp = image.getImagePlus();

			if (work.inputSettings.outline)
			{
				int[] clusters = opticsResult.getClusters();
				int max = Maths.max(clusters);
				
				// The current overlay may be fine. 
				// If the result has cached convex hulls then assume we have constructed an overlay
				// for these results.
				if (!opticsResult.hasConvexHulls() || o == null)
				{
					// We need to recompute
					o = new Overlay();
					opticsResult.computeConvexHulls();

					// Extract the ConvexHull of each cluster
					for (int c = 1; c <= max; c++)
					{
						ConvexHull hull = opticsResult.getConvexHull(c);
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
							// TODO: Create a colour to match the LUT
							//roi.setStrokeColor(color);
							// TODO: Options to set a fill colour?
							o.add(roi);
						}
					}
				}
				imp.setOverlay(o);
			}
			else
			{
				imp.setOverlay(null);
			}

			return work;
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
		add(workers, new OverlayResultsWorker());

		ArrayList<Thread> threads = new ArrayList<Thread>();

		for (Worker w : workers)
		{
			Thread t = new Thread(w);
			t.setDaemon(true);
			t.start();
			threads.add(t);
		}

		boolean cancelled = showDialog();

		// Finish work
		for (int i = 0; i < threads.size(); i++)
		{
			Thread t = threads.get(i);
			Worker w = workers.get(i);

			w.running = false;
			w.inbox.close();

			// Wait for threads to end
			if (cancelled)
			{
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
		gd.addCheckbox("Outline", inputSettings.outline);

		// Start disabled so the user can choose settings to update
		gd.addCheckbox("Preview", false);

		// Everything is done within the dialog listener
		gd.addDialogListener(this);

		gd.showDialog();

		return gd.wasCanceled();
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		if (e == null)
		{
			// This happens when the dialog is first shown and can be ignored.
			// It also happens when called from within a macro. In this case we should run.
			if (!Utils.isMacro())
				return true;
		}

		if (extraOptions)
			System.out.println("dialogItemChanged: " + e);

		// A previous run may have been cancelled so we have to handle this.
		if (Utils.isInterrupted())
		{
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
		settings.outline = gd.getNextBoolean();
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

		if (preview)
		{
			// Queue the settings
			if (extraOptions)
				System.out.println("Adding work:" + e);
			inputStack.addWork(new Work(settings, results));
		}

		return true;
	}
}
