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
/*
 * 
 */
package gdsc.smlm.ij.plugins;

import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.Macro;
import ij.WindowManager;
import ij.gui.GUI;
import ij.macro.Interpreter;
import ij.plugin.frame.Recorder;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Component;
import java.awt.Dialog;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.List;
import java.awt.Panel;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Locale;

/**
 * Shows a list of all the results sets held in memory, allowing multiple results to be selected
 */
public class MultiDialog extends Dialog
		implements ActionListener, KeyListener, WindowListener, MouseListener, ItemListener
{
	private static final long serialVersionUID = -881270633231897572L;

	private java.util.List<String> selected;
	private boolean selectAll = false;

	private Button cancel, okay, all, none;
	private boolean wasCanceled;
	private List list;
	private String macroOptions;
	private boolean macro;

	/**
	 * Interface to allow a list of any type to be shown in the MultiDialog
	 * 
	 * @author Alex Herbert
	 */
	public interface Items
	{
		/**
		 * Get the number of items to display.
		 *
		 * @return the size
		 */
		int size();

		/**
		 * Gets the formatted name of the result for display in the dialog.
		 *
		 * @param i
		 *            the result i
		 * @return the formatted name
		 */
		String getFormattedName(int i);

		/**
		 * Removes the formatting from the name. The plain name will be in the list returned by
		 * {@link MultiDialog#getSelectedResults()}.
		 *
		 * @param formattedName
		 *            the formatted name
		 * @return the plain name string
		 */
		String removeFormatting(String formattedName);
	}

	/**
	 * Base class for default implementation of the Items interface
	 * 
	 * @author Alex Herbert
	 */
	public static abstract class BaseItems implements Items
	{
		/**
		 * Returns the same formatted name.
		 * <p>
		 * {@inheritDoc}
		 * 
		 * @see gdsc.smlm.ij.plugins.MultiDialog.Items#removeFormatting(java.lang.String)
		 */
		@Override
		public String removeFormatting(String formattedName)
		{
			return formattedName;
		}
	}

	/**
	 * Interface to allow resulst to populate the items in the multi dialog
	 */
	public interface MemoryResultsFilter
	{
		/**
		 * Accept the results.
		 *
		 * @param results
		 *            the results
		 * @return true, if successful
		 */
		public boolean accept(MemoryPeakResults results);
	}

	private static class NullMemoryResultsFilter implements MemoryResultsFilter
	{
		@Override
		public boolean accept(MemoryPeakResults results)
		{
			return true;
		}
	}

	/**
	 * Class that allows the current results held in memory to be shown in the dialog
	 * 
	 * @author Alex Herbert
	 */
	public static class MemoryResultsItems implements Items
	{
		private String[] names;
		private int size;

		public MemoryResultsItems()
		{
			this(new NullMemoryResultsFilter());
		}

		public MemoryResultsItems(MemoryResultsFilter filter)
		{
			Collection<MemoryPeakResults> allResults = MemoryPeakResults.getAllResults();
			names = new String[allResults.size()];
			size = 0;
			for (MemoryPeakResults results : allResults)
				if (filter.accept(results))
					names[size++] = ResultsManager.getName(results);
		}

		@Override
		public int size()
		{
			return size;
		}

		@Override
		public String getFormattedName(int i)
		{
			return names[i];
		}

		@Override
		public String removeFormatting(String formattedName)
		{
			return ResultsManager.removeFormatting(formattedName);
		}
	}

	private Items items;

	public MultiDialog(String title, Items items)
	{
		super(WindowManager.getCurrentImage() != null ? (Frame) WindowManager.getCurrentImage().getWindow()
				: IJ.getInstance() != null ? IJ.getInstance() : new Frame(), title, true);
		addKeyListener(this);
		addWindowListener(this);
		macroOptions = Macro.getOptions();
		macro = macroOptions != null;
		this.items = items;
	}

	public void addSelected(java.util.List<String> selected)
	{
		this.selected = selected;
	}

	public boolean isSelectAll()
	{
		return selectAll;
	}

	public void setSelectAll(boolean selectAll)
	{
		this.selectAll = selectAll;
	}

	public void showDialog()
	{
		// Detect if running in a macro and just collect the input options
		if (macro)
		{
			dispose();
		}
		else
		{
			add(buildPanels());
			this.addKeyListener(this);
			if (IJ.isMacintosh())
				setResizable(false);
			pack();
			GUI.center(this);
			setVisible(true);
			IJ.wait(50); // work around for Sun/WinNT bug
		}
	}

	protected Panel buildPanels()
	{
		Panel p = new Panel();
		BorderLayout layout = new BorderLayout();
		layout.setVgap(3);
		p.setLayout(layout);
		p.add(buildResultsList(), BorderLayout.NORTH, 0);
		p.add(buildButtonPanel(), BorderLayout.CENTER, 1);
		return p;
	}

	protected Component buildResultsList()
	{
		final int MAX_SIZE = 30;
		int size = items.size();
		int rows = (size < MAX_SIZE) ? size : MAX_SIZE;
		list = new List(rows, true);
		for (int n = 0; n < size; n++)
		{
			String formattedName = items.getFormattedName(n);
			list.add(formattedName);
			// Initial selection
			if (selectAll || (selected != null && selected.contains(items.removeFormatting(formattedName))))
			{
				list.select(n);
			}
		}

		list.addMouseListener(this);
		list.addItemListener(this);
		list.addKeyListener(this);

		return list;
	}

	protected Panel buildButtonPanel()
	{
		Panel buttons = new Panel();
		buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 0));
		all = new Button("All");
		all.addActionListener(this);
		all.addKeyListener(this);
		buttons.add(all);
		none = new Button("None");
		none.addActionListener(this);
		none.addKeyListener(this);
		buttons.add(none);
		okay = new Button("OK");
		okay.addActionListener(this);
		okay.addKeyListener(this);
		buttons.add(okay);
		cancel = new Button("Cancel");
		cancel.addActionListener(this);
		cancel.addKeyListener(this);
		buttons.add(cancel);
		return buttons;
	}

	public boolean wasCanceled()
	{
		return wasCanceled;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	@Override
	public void actionPerformed(ActionEvent e)
	{
		Object source = e.getSource();
		if (source == okay || source == cancel)
		{
			wasCanceled = source == cancel;
			dispose();
		}
		else if (source == all)
		{
			for (int i = 0; i < list.getItemCount(); i++)
				list.select(i);
		}
		else if (source == none)
		{
			for (int i = 0; i < list.getItemCount(); i++)
				list.deselect(i);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.KeyListener#keyTyped(java.awt.event.KeyEvent)
	 */
	@Override
	public void keyTyped(KeyEvent paramKeyEvent)
	{
	}

	@Override
	public void keyPressed(KeyEvent e)
	{
		int keyCode = e.getKeyCode();
		IJ.setKeyDown(keyCode);
		if (keyCode == KeyEvent.VK_ENTER)
		{
			Object source = e.getSource();
			if (source == okay || source == cancel || source == list)
			{
				wasCanceled = source == cancel;
				dispose();
			}
			else if (source == all)
			{
				for (int i = 0; i < list.getItemCount(); i++)
					list.select(i);
			}
			else if (source == none)
			{
				for (int i = 0; i < list.getItemCount(); i++)
					list.deselect(i);
			}
		}
		else if (keyCode == KeyEvent.VK_ESCAPE)
		{
			wasCanceled = true;
			dispose();
			IJ.resetEscape();
		}
		else if (keyCode == KeyEvent.VK_W &&
				(e.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0)
		{
			wasCanceled = true;
			dispose();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.KeyListener#keyReleased(java.awt.event.KeyEvent)
	 */
	@Override
	public void keyReleased(KeyEvent paramKeyEvent)
	{
	}

	public ArrayList<String> getSelectedResults()
	{
		ArrayList<String> selected;

		// Get the selected names
		if (macro)
		{
			selected = new ArrayList<String>();
			String name = getValue("input");
			while (name != null)
			{
				selected.add(name);
				name = getValue("input" + selected.size());
			}
		}
		else
		{
			final int[] listIndexes = list.getSelectedIndexes();
			selected = new ArrayList<String>(listIndexes.length);
			if (listIndexes.length > 0)
			{
				for (int index : listIndexes)
				{
					selected.add(items.removeFormatting(list.getItem(index)));
				}
			}
		}

		// Record as if we use the multiple_inputs option
		if ((macro && Recorder.record && Recorder.recordInMacros) || Recorder.record)
		{
			if (!selected.isEmpty())
			{
				Recorder.recordOption("Input", selected.get(0));
				if (selected.size() > 1)
				{
					Recorder.recordOption("Multiple_inputs");
					for (int n = 1; n < selected.size(); ++n)
					{
						Recorder.recordOption("Input" + n, selected.get(n));
					}
				}
			}
		}

		return selected;
	}

	/**
	 * Get a value from the macro options. Adapted from ij.gui.GenericDialog.
	 * 
	 * @param label
	 * @return The value (or null)
	 */
	private String getValue(String label)
	{
		String theText = Macro.getValue(macroOptions, label, null);
		if (theText != null && (theText.startsWith("&") || label.toLowerCase(Locale.US).startsWith(theText)))
		{
			// Is the value a macro variable?
			if (theText.startsWith("&"))
				theText = theText.substring(1);
			Interpreter interp = Interpreter.getInstance();
			String s = interp != null ? interp.getVariableAsString(theText) : null;
			if (s != null)
				theText = s;
		}
		return theText;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.WindowListener#windowClosing(java.awt.event.WindowEvent)
	 */
	@Override
	public void windowClosing(WindowEvent e)
	{
		wasCanceled = true;
		dispose();
	}

	//@formatter:off
	@Override
	public void windowActivated(WindowEvent e)
	{
	}

	@Override
	public void windowOpened(WindowEvent e)
	{
	}

	@Override
	public void windowClosed(WindowEvent e)
	{
	}

	@Override
	public void windowIconified(WindowEvent e)
	{
	}

	@Override
	public void windowDeiconified(WindowEvent e)
	{
	}

	@Override
	public void windowDeactivated(WindowEvent e)
	{
	}

	@Override
	public void mousePressed(MouseEvent paramMouseEvent)
	{
	}

	@Override
	public void mouseReleased(MouseEvent paramMouseEvent)
	{
	}

	@Override
	public void mouseEntered(MouseEvent paramMouseEvent)
	{
	}

	@Override
	public void mouseExited(MouseEvent paramMouseEvent)
	{
	}

	//@formatter:on

	int lastIndex;
	int modifiers;
	int lastEvent = -1;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseClicked(java.awt.event.MouseEvent)
	 */
	@Override
	public void mouseClicked(MouseEvent paramMouseEvent)
	{
		modifiers = paramMouseEvent.getModifiers();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ItemListener#itemStateChanged(java.awt.event.ItemEvent)
	 */
	@Override
	public void itemStateChanged(ItemEvent paramItemEvent)
	{
		int index = (Integer) paramItemEvent.getItem();
		int event = paramItemEvent.getStateChange();

		// If we have the shift key down, support multiple select/deselect
		if (event == lastEvent && (modifiers & InputEvent.SHIFT_MASK) != 0 &&
				(event == ItemEvent.SELECTED || event == ItemEvent.DESELECTED))
		{
			if (lastIndex != index)
			{
				int top = Math.max(index, lastIndex);
				int bottom = Math.min(index, lastIndex);
				for (int i = bottom + 1; i < top; i++)
				{
					if (event == ItemEvent.SELECTED)
						list.select(i);
					else
						list.deselect(i);
				}
			}
		}

		lastEvent = event;
		lastIndex = index;
	}
}
