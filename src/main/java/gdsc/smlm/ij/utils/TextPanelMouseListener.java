package gdsc.smlm.ij.utils;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import ij.text.TextPanel;

/**
 * Attaches to a text panel and listens for mouse events. Responds to double click on a single line or a selection of
 * multiple lines.
 */
public abstract class TextPanelMouseListener implements MouseListener
{
	protected TextPanel textPanel;

	/**
	 * Instantiates a new text panel mouse listener.
	 */
	public TextPanelMouseListener()
	{
	}

	/**
	 * Instantiates a new text panel mouse listener.
	 *
	 * @param textPanel
	 *            The text panel to listen to for mouse events
	 */
	public TextPanelMouseListener(TextPanel textPanel)
	{
		setTextPanel(textPanel);
	}

	/**
	 * Sets the text panel.
	 *
	 * @param textPanel
	 *            the new text panel
	 */
	public void setTextPanel(TextPanel textPanel)
	{
		if (this.textPanel != null)
			this.textPanel.removeMouseListener(this);
		this.textPanel = textPanel;
		if (this.textPanel != null)
			this.textPanel.addMouseListener(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseClicked(java.awt.event.MouseEvent)
	 */
	public void mouseClicked(MouseEvent e)
	{
		// Show the result that was double clicked in the result table
		if (e.getClickCount() > 1)
		{
			selected(textPanel.getSelectionStart());
		}
	}

	/**
	 * Trigger that a single line from the panel has been selected.
	 *
	 * @param selectedIndex
	 *            the selected index
	 */
	protected abstract void selected(int selectedIndex);

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mousePressed(java.awt.event.MouseEvent)
	 */
	public void mousePressed(MouseEvent e)
	{
		int index = textPanel.getSelectionStart();
		int index2 = textPanel.getSelectionEnd();
		if (index == index2)
			return;
		selected(textPanel.getSelectionStart(), textPanel.getSelectionEnd());
	}

	/**
	 * Trigger that multiple lines from the panel have been selected.
	 *
	 * @param selectionStart
	 *            the selection start
	 * @param selectionEnd
	 *            the selection end
	 */
	protected abstract void selected(int selectionStart, int selectionEnd);

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseReleased(java.awt.event.MouseEvent)
	 */
	public void mouseReleased(MouseEvent e)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseEntered(java.awt.event.MouseEvent)
	 */
	public void mouseEntered(MouseEvent e)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.MouseListener#mouseExited(java.awt.event.MouseEvent)
	 */
	public void mouseExited(MouseEvent e)
	{
		// Ignore
	}
}
