/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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
package gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
import java.awt.Choice;
import java.awt.Label;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.Arrays;

import gdsc.core.data.utils.TypeConverter;
import gdsc.core.ij.Utils;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.results.IJTablePeakResults;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultView;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageWindow;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.text.TextPanel;
import ij.text.TextWindow;

/**
 * Produces a summary table of the results that are stored in memory.
 */
public class OverlayResults implements PlugIn, ItemListener, ImageListener
{
	private static final String TITLE = "Overlay Results";
	private static String name = "";
	private static boolean showTable = false;

	private String[] names;
	private int[] ids;
	private Choice choice;
	private Checkbox checkbox;
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
		// The results text window (so we can close it)
		private TextWindow tw = null;
		private Rectangle windowBounds = null;

		TFloatArrayList ox = new TFloatArrayList(100);
		TFloatArrayList oy = new TFloatArrayList(100);
		PeakResultView view = null;
		TypeConverter<DistanceUnit> converter;

		@Override
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
			closeTextWindow();
		}

		private void clearOldOverlay()
		{
			if (currentIndex != 0)
			{
				ImagePlus oldImp = WindowManager.getImage(ids[currentIndex]);
				if (oldImp != null)
					oldImp.setOverlay(null);
			}
			view = null;
			currentSlice = -1;
			currentIndex = 0;
		}

		private void closeTextWindow()
		{
			if (tw != null)
			{
				windowBounds = tw.getBounds();
				tw.close();
				tw = null;
			}
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
			{
				boolean isShowing = tw != null;
				if (showTable == isShowing)
					// No change from last time
					return;
			}
			currentSlice = newSlice;

			MemoryPeakResults results = MemoryPeakResults.getResults(name);
			if (results == null)
			{
				// Results have been cleared from memory (or renamed).
				logError("Results not available", name);
				return;
			}
			clearError();

			final IJTablePeakResults table;
			TIntHashSet selectedId = null;
			if (showTable)
			{
				boolean hasId = results.hasId();

				// Old selected item
				boolean is3D = false;
				if (hasId && tw != null)
				{
					TextPanel tp = tw.getTextPanel();
					int idColumn = Utils.getColumn(tp, "Id");
					int start = tp.getSelectionStart();
					if (start != -1 && idColumn != -1)
					{
						selectedId = new TIntHashSet();
						int end = tp.getSelectionEnd();
						for (int index = start; index <= end; index++)
						{
							String text = tp.getLine(index).split("\t")[idColumn];
							selectedId.add(Integer.parseInt(text));
						}
					}
					// Keep the z column to avoid table redraw
					is3D = tp.getColumnHeadings().contains("\tZ");
				}

				// New table
				is3D = is3D || results.is3D();
				table = new IJTablePeakResults(false);
				table.setTableTitle(TITLE);
				table.copySettings(results);
				table.setClearAtStart(true);
				table.setAddCounter(true);
				table.setHideSourceText(true);
				table.setShowZ(is3D);
				table.setShowId(hasId);
				//table.setShowFittingData(true);
				//table.setShowNoiseData(true);
				table.begin();

				tw = table.getResultsWindow();
				if (windowBounds != null)
				{
					tw.setBounds(windowBounds);
				}
				else
				{
					// Position under the window
					ImageWindow win = imp.getWindow();
					Point p = win.getLocation();
					p.y += win.getHeight();
					tw.setLocation(p);
				}
			}
			else
			{
				table = null;
				closeTextWindow();
			}

			ox.resetQuick();
			oy.resetQuick();
			if (view == null)
			{
				view = results.getSnapshotPeakResultView();
				converter = results.getDistanceConverter(DistanceUnit.PIXEL);
			}
			int select = -1;
			PeakResult[] frameResults = view.getResultsByFrame(currentSlice);
			for (int i = 0; i < frameResults.length; i++)
			{
				PeakResult r = frameResults[i];
				ox.add(converter.convert(r.getXPosition()));
				oy.add(converter.convert(r.getYPosition()));
				if (table != null)
				{
					table.add(r);
					if (selectedId != null)
					{
						if (selectedId.contains(r.getId()))
						{
							// For now just preserve the first selected ID.
							// This at least allows tracking a single localisation.
							select = i;
							selectedId = null;
						}
					}
				}
			}
			// Old method without the cached view
			//			results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure()
			//			{
			//				public void executeXYR(float x, float y, PeakResult r)
			//				{
			//					if (r.getFrame() == currentSlice)
			//					{
			//						ox.add(x);
			//						oy.add(y);
			//						if (table != null)
			//							table.add(r);
			//					}
			//				}
			//			});
			PointRoi roi = new PointRoi(ox.toArray(), oy.toArray());
			roi.setPointType(3);
			// Leave to ImageJ default. Then the user can change it using the options.
			//Color c = Color.GREEN;
			//roi.setStrokeColor(c);
			//roi.setFillColor(c);
			imp.getWindow().toFront();
			imp.setOverlay(new Overlay(roi));

			if (table != null)
			{
				table.end();
				TextWindow tw = table.getResultsWindow();
				TextPanel tp = tw.getTextPanel();
				tp.scrollToTop();

				// Reselect the same Id
				if (select != -1)
					table.select(select);
			}
		}

		private void logError(String msg, String name)
		{
			if (!error[currentIndex])
			{
				Utils.log("%s Error: %s for results '%s'", TITLE, msg, name);
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
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
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
			if (results.getSource() != null && results.getSource().getOriginal() instanceof IJImageSource)
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
		if (c == 1)
		{
			IJ.error(TITLE, "There are no result images available");
			return;
		}
		names = Arrays.copyOf(names, c);

		Thread t = null;
		Worker w = null;
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(TITLE);
		gd.addMessage("Overlay results on current image frame");
		gd.addChoice("Results", names, (name == null) ? "" : name);
		gd.addCheckbox("Show_table", showTable);
		gd.addMessage("");
		gd.addHelp(About.HELP_URL);
		gd.hideCancelButton();
		gd.setOKLabel("Close");
		if (!(IJ.isMacro() || java.awt.GraphicsEnvironment.isHeadless()))
		{
			choice = (Choice) gd.getChoices().get(0);
			choice.addItemListener(this);
			checkbox = (Checkbox) gd.getCheckboxes().get(0);
			checkbox.addItemListener(this);
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
			showTable = gd.getNextBoolean();
		}
		if (t != null && w != null)
		{
			w.running = false;
			inbox.close();
			try
			{
				t.join(0);
			}
			catch (InterruptedException e)
			{
				// Ignore
			}
			t = null;
		}
	}

	@Override
	public void itemStateChanged(ItemEvent e)
	{
		show();
	}

	@Override
	public void imageClosed(ImagePlus arg0)
	{
		// Ignore
	}

	@Override
	public void imageOpened(ImagePlus arg0)
	{
		// Ignore
	}

	@Override
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
		showTable = checkbox.getState();
		inbox.add(choice.getSelectedIndex());
	}
}
