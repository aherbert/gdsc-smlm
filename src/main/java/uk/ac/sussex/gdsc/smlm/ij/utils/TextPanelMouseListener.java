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
package uk.ac.sussex.gdsc.smlm.ij.utils;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import ij.text.TextPanel;

/**
 * Attaches to a text panel and listens for mouse events. Responds to double click on a single line or a selection of
 * multiple lines.
 */
public abstract class TextPanelMouseListener implements MouseListener
{
	/** The text panel. */
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
	@Override
	public void mouseClicked(MouseEvent e)
	{
		// Show the result that was double clicked in the result table
		if (e.getClickCount() > 1)
			selected(textPanel.getSelectionStart());
	}

	/**
	 * Triggered when a single line from the panel has been selected.
	 *
	 * @param selectedIndex
	 *            the selected index
	 */
	public abstract void selected(int selectedIndex);

	/*
	 * (non-Javadoc)
	 *
	 * @see java.awt.event.MouseListener#mousePressed(java.awt.event.MouseEvent)
	 */
	@Override
	public void mousePressed(MouseEvent e)
	{
		final int index = textPanel.getSelectionStart();
		final int index2 = textPanel.getSelectionEnd();
		if (index == index2)
			return;
		selected(textPanel.getSelectionStart(), textPanel.getSelectionEnd());
	}

	/**
	 * Triggered when multiple lines from the panel have been selected.
	 *
	 * @param selectionStart
	 *            the selection start
	 * @param selectionEnd
	 *            the selection end
	 */
	public abstract void selected(int selectionStart, int selectionEnd);

	/*
	 * (non-Javadoc)
	 *
	 * @see java.awt.event.MouseListener#mouseReleased(java.awt.event.MouseEvent)
	 */
	@Override
	public void mouseReleased(MouseEvent e)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.awt.event.MouseListener#mouseEntered(java.awt.event.MouseEvent)
	 */
	@Override
	public void mouseEntered(MouseEvent e)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.awt.event.MouseListener#mouseExited(java.awt.event.MouseEvent)
	 */
	@Override
	public void mouseExited(MouseEvent e)
	{
		// Ignore
	}
}
