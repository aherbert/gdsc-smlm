/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import ij.text.TextPanel;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

/**
 * Attaches to a text panel and listens for mouse events. Responds to double click on a single line
 * or a selection of multiple lines.
 */
public abstract class TextPanelMouseListener extends MouseAdapter {
  /** The text panel. */
  protected TextPanel textPanel;

  /**
   * Instantiates a new text panel mouse listener.
   *
   * @param textPanel The text panel to listen to for mouse events
   */
  public TextPanelMouseListener(TextPanel textPanel) {
    // Note: The TextPanel.addMouseListener passes through the listener to
    // the backing TextCanvas. There is no corresponding TextPanel.removeMouseListener
    // method. So this listener is bound for the entire the lifetime of the TextPanel.
    textPanel.addMouseListener(this);
    this.textPanel = textPanel;
  }

  @Override
  public void mouseClicked(MouseEvent event) {
    // Show the result that was double clicked in the result table
    if (event.getClickCount() > 1) {
      selected(textPanel.getSelectionStart());
    }
  }

  /**
   * Triggered when a single line from the panel has been selected.
   *
   * <p>The default implementation calls {@link #selected(int, int)} with the provided selected
   * index for the start and end.
   *
   * @param selectedIndex the selected index
   */
  public void selected(int selectedIndex) {
    selected(selectedIndex, selectedIndex);
  }

  /**
   * Triggered when multiple lines from the panel have been selected.
   *
   * @param selectionStart the selection start
   * @param selectionEnd the selection end
   */
  public abstract void selected(int selectionStart, int selectionEnd);

  @Override
  public void mousePressed(MouseEvent event) {
    final int index = textPanel.getSelectionStart();
    final int index2 = textPanel.getSelectionEnd();
    if (index == index2) {
      return;
    }
    selected(textPanel.getSelectionStart(), textPanel.getSelectionEnd());
  }
}
