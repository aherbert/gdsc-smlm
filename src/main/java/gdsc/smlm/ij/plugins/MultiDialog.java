package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

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
public class MultiDialog extends Dialog implements ActionListener, KeyListener, WindowListener, MouseListener,
		ItemListener
{
	private static final long serialVersionUID = -881270633231897572L;

	private ArrayList<String> selected;

	private Button cancel, okay, all, none;
	private boolean wasCanceled;
	private List list;
	private String macroOptions;
	private boolean macro;

	public MultiDialog(String title)
	{
		super(WindowManager.getCurrentImage() != null ? (Frame) WindowManager.getCurrentImage().getWindow() : IJ
				.getInstance() != null ? IJ.getInstance() : new Frame(), title, true);
		addKeyListener(this);
		addWindowListener(this);
		macroOptions = Macro.getOptions();
		macro = macroOptions != null;
	}

	public void addSelected(ArrayList<String> selected)
	{
		this.selected = selected;
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
		Collection<MemoryPeakResults> alResults = MemoryPeakResults.getAllResults();
		final int MAX_SIZE = 30;
		int size;
		if (alResults.size() < MAX_SIZE)
		{
			size = alResults.size();
		}
		else
		{
			size = MAX_SIZE;
		}
		list = new List(size, true);
		int n = 0;
		for (MemoryPeakResults results : alResults)
		{
			String formattedName = ResultsManager.getName(results);
			list.add(formattedName);
			// Select the same as last time
			if (selected != null && selected.contains(results.getName()))
			{
				list.select(n);
			}
			n++;
		}

		list.addMouseListener(this);
		list.addItemListener(this);
		list.addKeyListener(this);

		return (Component) list;
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
					selected.add(ResultsManager.removeFormatting(list.getItem(index)));
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

	@Override
	public void windowClosing(WindowEvent e)
	{
		wasCanceled = true;
		dispose();
	}

	//@formatter:off
    public void windowActivated(WindowEvent e) {}
    public void windowOpened(WindowEvent e) {}
    public void windowClosed(WindowEvent e) {}
    public void windowIconified(WindowEvent e) {}
    public void windowDeiconified(WindowEvent e) {}
    public void windowDeactivated(WindowEvent e) {}
	public void mousePressed(MouseEvent paramMouseEvent) {}
	public void mouseReleased(MouseEvent paramMouseEvent) {}
	public void mouseEntered(MouseEvent paramMouseEvent) {}
	public void mouseExited(MouseEvent paramMouseEvent) {}
	//@formatter:on

	int lastIndex;
	int modifiers;
	int lastEvent = -1;

	@Override
	public void mouseClicked(MouseEvent paramMouseEvent)
	{
		modifiers = paramMouseEvent.getModifiers();
	}

	@Override
	public void itemStateChanged(ItemEvent paramItemEvent)
	{
		int index = (int) paramItemEvent.getItem();
		int event = paramItemEvent.getStateChange();

		// If we have the shift key down, support multiple select/deselect
		if (event == lastEvent && (modifiers & MouseEvent.SHIFT_MASK) != 0 &&
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