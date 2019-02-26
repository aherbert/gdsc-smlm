package gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.ij.ImageJUtils; import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils; import uk.ac.sussex.gdsc.core.utils.TextUtils; import uk.ac.sussex.gdsc.core.utils.MathUtils;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;

import java.awt.Choice;
import java.awt.Label;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Arrays;

/**
 * Produces a summary table of the results that are stored in memory.
 */
public class OverlayResults implements PlugIn, ItemListener, ImageListener
{
	private static final String TITLE = "Overlay Results";
	private static String name = "";
	private String[] names;
	private int[] ids;
	private Choice choice;
	private Label label;

	private int currentIndex = 0;
	private int currentSlice = -1;

	private class Job
	{
		final int index;

		Job(int index)
		{
			this.index = index;
		}
	}

	private class InBox
	{
		private Job job = null;

		synchronized void add(int index)
		{
			this.job = new Job(index);
			this.notify();
		}

		synchronized void close()
		{
			this.job = null;
			this.notify();
		}

		synchronized Job next()
		{
			Job job = this.job;
			this.job = null;
			return job;
		}

		boolean isEmpty()
		{
			return job == null;
		}
	}

	private InBox inbox = new InBox();

	private class Worker implements Runnable
	{
		private boolean running = true;
		private boolean[] error = new boolean[ids.length];

		public void run()
		{
			while (running)
			{
				try
				{
					Job job = null;
					synchronized (inbox)
					{
						if (inbox.isEmpty())
							inbox.wait();
						job = inbox.next();
					}
					if (job == null || !running)
						break;
					if (job.index == 0)
					{
						// This may be selection of no image
						clearOldOverlay();
						continue;
					}

					// Check name of the image
					if (currentIndex != job.index)
						clearOldOverlay();

					currentIndex = job.index;
					drawOverlay();
				}
				catch (InterruptedException e)
				{
					break;
				}
			}
			clearOldOverlay();
		}

		private void clearOldOverlay()
		{
			if (currentIndex != 0)
			{
				ImagePlus oldImp = WindowManager.getImage(ids[currentIndex]);
				if (oldImp != null)
					oldImp.setOverlay(null);
			}
			currentSlice = -1;
			currentIndex = 0;
		}

		/**
		 * Draw the overlay.
		 * <p>
		 * This is only called when index > 0.
		 */
		private void drawOverlay()
		{
			ImagePlus imp = WindowManager.getImage(ids[currentIndex]);
			String name = names[currentIndex];

			if (imp == null)
			{
				// Image has been closed.
				logError("Image not available", name);
				return;
			}

			// Check slice
			int newSlice = imp.getCurrentSlice();
			if (currentSlice == newSlice)
				// No change from last time
				return;
			currentSlice = newSlice;

			MemoryPeakResults results = MemoryPeakResults.getResults(name);
			if (results == null)
			{
				// Results have been cleared from memory (or renamed).
				logError("Results not available", name);
				return;
			}
			clearError();

			float[] ox = new float[100];
			float[] oy = new float[100];
			int points = 0;
			for (PeakResult r : results.getResults())
			{
				if (r.peak != currentSlice)
					continue;
				if (points == ox.length)
				{
					ox = Arrays.copyOf(ox, (int) (points * 1.5));
					oy = Arrays.copyOf(oy, ox.length);
				}
				ox[points] = r.getXPosition();
				oy[points] = r.getYPosition();
				points++;
			}
			PointRoi roi = new PointRoi(ox, oy, points);
			imp.getWindow().toFront();
			imp.setOverlay(new Overlay(roi));
		}

		private void logError(String msg, String name)
		{
			if (!error[currentIndex])
			{
				ImageJUtils.log("%s Error: %s for results '%s'", TITLE, msg, name);
				label.setText("Error: " + msg + ". Restart this plugin to refresh.");
			}
			error[currentIndex] = true;
		}

		private void clearError()
		{
			error[currentIndex] = false;
			if (!TextUtils.isNullOrEmpty(label.getText()))
				label.setText("");
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
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		names = new String[MemoryPeakResults.getResultNames().size() + 1];
		ids = new int[names.length];
		int c = 0;
		names[c++] = "(None)";
		for (MemoryPeakResults results : MemoryPeakResults.getAllResults())
		{
			if (results.getSource().getOriginal() instanceof IJImageSource)
			{
				IJImageSource source = (IJImageSource) (results.getSource().getOriginal());
				ImagePlus imp = WindowManager.getImage(source.getName());
				if (imp != null)
				{
					ids[c] = imp.getID();
					names[c++] = results.getName();
				}
			}
		}
		if (c == 0)
		{
			IJ.error(TITLE, "There are no result images available");
			return;
		}
		Arrays.copyOf(names, c);

		Thread t = null;
		Worker w = null;
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addMessage("Overlay results on current image frame");
		gd.addChoice("Results", names, name);
		gd.addMessage("");
		gd.addHelp(About.HELP_URL);
		gd.hideCancelButton();
		gd.setOKLabel("Close");
		if (!(IJ.isMacro() || java.awt.GraphicsEnvironment.isHeadless()))
		{
			choice = (Choice) gd.getChoices().get(0);
			choice.addItemListener(this);
			label = (Label) gd.getMessage();

			// Listen for changes to an image
			ImagePlus.addImageListener(this);

			show();

			t = new Thread(w = new Worker());
			t.setDaemon(true);
			t.start();
		}
		gd.showDialog();
		if (!(IJ.isMacro() || java.awt.GraphicsEnvironment.isHeadless()))
			ImagePlus.removeImageListener(this);
		if (!gd.wasCanceled())
		{
			name = gd.getNextChoice();
		}
		if (t != null)
		{
			w.running = false;
			inbox.close();
			try
			{
				t.join(0);
			}
			catch (InterruptedException e)
			{
			}
			t = null;
		}
	}

	public void itemStateChanged(ItemEvent e)
	{
		show();
	}

	public void imageClosed(ImagePlus arg0)
	{
	}

	public void imageOpened(ImagePlus arg0)
	{
	}

	public void imageUpdated(ImagePlus imp)
	{
		if (imp == null)
			return;
		if (ids[currentIndex] == imp.getID())
		{
			show();
		}
	}

	private void show()
	{
		inbox.add(choice.getSelectedIndex());
	}
}
