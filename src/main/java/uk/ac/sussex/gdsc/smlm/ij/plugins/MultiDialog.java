/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

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
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Locale;
import java.util.function.Predicate;

/**
 * Shows a list of all the results sets held in memory, allowing multiple results to be selected.
 */
public class MultiDialog extends Dialog {
  private static final long serialVersionUID = 21012019L;

  /**
   * Maximum number of items to show in the displayed list.
   */
  private static final int MAX_SIZE = 30;

  private java.util.List<String> selected;
  private boolean selectAll;

  private Button cancel;
  private Button okay;
  private Button all;
  private Button none;
  private boolean wasCanceled;
  private List list;
  private final String macroOptions;
  private final boolean macro;

  private final Items items;

  /**
   * The last index from {@link ItemEvent#getItem()} captured in
   * {@link ItemListener#itemStateChanged(ItemEvent)}.
   */
  protected int lastIndex;

  /**
   * The modifiers captured in from {@link MouseListener#mouseClicked(MouseEvent)}.
   */
  protected int modifiers;

  /**
   * The last event from {@link ItemEvent#getStateChange()} captured in
   * {@link ItemListener#itemStateChanged(ItemEvent)}.
   */
  protected int lastEvent = -1;

  private LocalActionListener actionListener = new LocalActionListener();
  private LocalKeyListener keyListener = new LocalKeyListener();
  private LocalWindowAdpater windowAdpater = new LocalWindowAdpater();
  private LocalMouseAdpater mouseAdpater = new LocalMouseAdpater();
  private LocalItemListener itemListener = new LocalItemListener();

  private class LocalActionListener implements ActionListener {
    @Override
    public void actionPerformed(ActionEvent event) {
      final Object source = event.getSource();
      if (source == okay || source == cancel) {
        wasCanceled = source == cancel;
        dispose();
      } else if (source == all) {
        for (int i = 0; i < list.getItemCount(); i++) {
          list.select(i);
        }
      } else if (source == none) {
        for (int i = 0; i < list.getItemCount(); i++) {
          list.deselect(i);
        }
      }
    }
  }

  private class LocalKeyListener extends KeyAdapter {
    @Override
    public void keyPressed(KeyEvent event) {
      final int keyCode = event.getKeyCode();
      IJ.setKeyDown(keyCode);
      if (keyCode == KeyEvent.VK_ENTER) {
        final Object source = event.getSource();
        if (source == okay || source == cancel || source == list) {
          wasCanceled = source == cancel;
          dispose();
        } else if (source == all) {
          for (int i = 0; i < list.getItemCount(); i++) {
            list.select(i);
          }
        } else if (source == none) {
          for (int i = 0; i < list.getItemCount(); i++) {
            list.deselect(i);
          }
        }
      } else if (keyCode == KeyEvent.VK_ESCAPE) {
        wasCanceled = true;
        dispose();
        IJ.resetEscape();
      } else if (keyCode == KeyEvent.VK_W
          && (event.getModifiers() & Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()) != 0) {
        wasCanceled = true;
        dispose();
      }
    }
  }

  private class LocalWindowAdpater extends WindowAdapter {
    @Override
    public void windowClosing(WindowEvent event) {
      wasCanceled = true;
      dispose();
    }
  }

  private class LocalMouseAdpater extends MouseAdapter {
    @Override
    public void mouseClicked(MouseEvent event) {
      modifiers = event.getModifiers();
    }
  }

  private class LocalItemListener implements ItemListener {
    @Override
    public void itemStateChanged(ItemEvent paramItemEvent) {
      final int index = (Integer) paramItemEvent.getItem();
      final int event = paramItemEvent.getStateChange();

      // If we have the shift key down, support multiple select/deselect
      if (event == lastEvent && (modifiers & InputEvent.SHIFT_MASK) != 0
          && (event == ItemEvent.SELECTED || event == ItemEvent.DESELECTED) && lastIndex != index) {
        final int top = Math.max(index, lastIndex);
        final int bottom = Math.min(index, lastIndex);
        for (int i = bottom + 1; i < top; i++) {
          if (event == ItemEvent.SELECTED) {
            list.select(i);
          } else {
            list.deselect(i);
          }
        }
      }

      lastEvent = event;
      lastIndex = index;
    }
  }

  /**
   * Interface to allow a list of any type to be shown in the MultiDialog.
   */
  public interface Items {
    /**
     * Get the number of items to display.
     *
     * @return the size
     */
    int size();

    /**
     * Gets the formatted name of the result for display in the dialog.
     *
     * @param index the result index
     * @return the formatted name
     */
    String getFormattedName(int index);

    /**
     * Removes the formatting from the name. The plain name will be in the list returned by
     * {@link MultiDialog#getSelectedResults()}.
     *
     * @param formattedName the formatted name
     * @return the plain name string
     */
    String removeFormatting(String formattedName);
  }

  /**
   * Base class for default implementation of the Items interface.
   */
  public abstract static class BaseItems implements Items {
    /**
     * Returns the same formatted name.
     */
    @Override
    public String removeFormatting(String formattedName) {
      return formattedName;
    }
  }

  private static class NullMemoryResultsFilter implements Predicate<MemoryPeakResults> {
    @Override
    public boolean test(MemoryPeakResults results) {
      return true;
    }
  }

  /**
   * Class that allows the current results held in memory to be shown in the dialog.
   */
  public static class MemoryResultsItems implements Items {
    private final String[] names;
    private int size;

    /**
     * Instantiates a new memory results items.
     */
    public MemoryResultsItems() {
      this(new NullMemoryResultsFilter());
    }

    /**
     * Instantiates a new memory results items.
     *
     * @param filter the filter
     */
    public MemoryResultsItems(Predicate<MemoryPeakResults> filter) {
      final Collection<MemoryPeakResults> allResults = MemoryPeakResults.getAllResults();
      names = new String[allResults.size()];
      size = 0;
      for (final MemoryPeakResults results : allResults) {
        if (filter.test(results)) {
          names[size++] = ResultsManager.getName(results);
        }
      }
    }

    @Override
    public int size() {
      return size;
    }

    @Override
    public String getFormattedName(int index) {
      return names[index];
    }

    @Override
    public String removeFormatting(String formattedName) {
      return ResultsManager.removeFormatting(formattedName);
    }
  }

  /**
   * Instantiates a new multi dialog.
   *
   * @param title the title
   * @param items the items
   */
  public MultiDialog(String title, Items items) {
    super(getDialogOwner(), title, true);
    addKeyListener(keyListener);
    addWindowListener(windowAdpater);
    macroOptions = Macro.getOptions();
    macro = macroOptions != null;
    this.items = items;
  }

  private static Frame getDialogOwner() {
    if (WindowManager.getCurrentImage() != null) {
      return WindowManager.getCurrentImage().getWindow();
    }
    return IJ.getInstance();
  }

  /**
   * Adds the list of selected items.
   *
   * @param selected the selected
   */
  public void addSelected(java.util.List<String> selected) {
    this.selected = selected;
  }

  /**
   * Checks if select all items in the list.
   *
   * @return true, if is select all
   */
  public boolean isSelectAll() {
    return selectAll;
  }

  /**
   * Sets the select all flag.
   *
   * @param selectAll Set to true to select all items in the list
   */
  public void setSelectAll(boolean selectAll) {
    this.selectAll = selectAll;
  }

  /**
   * Show the dialog.
   */
  public void showDialog() {
    // Detect if running in a macro and just collect the input options
    if (macro) {
      dispose();
    } else {
      add(buildPanel());
      this.addKeyListener(keyListener);
      if (IJ.isMacintosh()) {
        setResizable(false);
      }
      pack();
      GUI.center(this);
      setVisible(true);
      IJ.wait(50); // work around for Sun/WinNT bug
    }
  }

  /**
   * Builds the main panel for the dialog.
   *
   * @return the panel
   */
  protected Panel buildPanel() {
    final Panel p = new Panel();
    final BorderLayout layout = new BorderLayout();
    layout.setVgap(3);
    p.setLayout(layout);
    p.add(buildResultsList(), BorderLayout.NORTH, 0);
    p.add(buildButtonPanel(), BorderLayout.CENTER, 1);
    return p;
  }

  /**
   * Builds the results list component for the dialog.
   *
   * @return the component
   */
  protected Component buildResultsList() {
    final int size = items.size();
    final int rows = (size < MAX_SIZE) ? size : MAX_SIZE;
    list = new List(rows, true);
    for (int n = 0; n < size; n++) {
      final String formattedName = items.getFormattedName(n);
      list.add(formattedName);
      // Initial selection
      if (selectAll
          || (selected != null && selected.contains(items.removeFormatting(formattedName)))) {
        list.select(n);
      }
    }

    list.addMouseListener(mouseAdpater);
    list.addItemListener(itemListener);
    list.addKeyListener(keyListener);

    return list;
  }

  /**
   * Builds the button panel for the dialog.
   *
   * @return the panel
   */
  protected Panel buildButtonPanel() {
    final Panel buttons = new Panel();
    buttons.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 0));
    all = new Button("All");
    all.addActionListener(actionListener);
    all.addKeyListener(keyListener);
    buttons.add(all);
    none = new Button("None");
    none.addActionListener(actionListener);
    none.addKeyListener(keyListener);
    buttons.add(none);
    okay = new Button("OK");
    okay.addActionListener(actionListener);
    okay.addKeyListener(keyListener);
    buttons.add(okay);
    cancel = new Button("Cancel");
    cancel.addActionListener(actionListener);
    cancel.addKeyListener(keyListener);
    buttons.add(cancel);
    return buttons;
  }

  /**
   * Check if the dialog was cancelled.
   *
   * @return true, if cancelled
   */
  public boolean wasCancelled() {
    return wasCanceled;
  }

  /**
   * Gets the selected results from the dialog.
   *
   * @return the selected results
   */
  public java.util.List<String> getSelectedResults() {
    ArrayList<String> selectedResults;

    // Get the selected names
    if (macro) {
      selectedResults = new ArrayList<>();
      String name = getValue("input");
      while (name != null) {
        selectedResults.add(name);
        name = getValue("input" + selectedResults.size());
      }
    } else {
      final int[] listIndexes = list.getSelectedIndexes();
      selectedResults = new ArrayList<>(listIndexes.length);
      if (listIndexes.length > 0) {
        for (final int index : listIndexes) {
          selectedResults.add(items.removeFormatting(list.getItem(index)));
        }
      }
    }

    // Record as if we use the multiple_inputs option
    if (((macro && Recorder.record && Recorder.recordInMacros) || Recorder.record)
        && !selectedResults.isEmpty()) {
      Recorder.recordOption("Input", selectedResults.get(0));
      if (selectedResults.size() > 1) {
        Recorder.recordOption("Multiple_inputs");
        for (int n = 1; n < selectedResults.size(); ++n) {
          Recorder.recordOption("Input" + n, selectedResults.get(n));
        }
      }
    }

    return selectedResults;
  }

  /**
   * Get a value from the macro options. Adapted from ij.gui.GenericDialog.
   *
   * @param label the label
   * @return The value (or null)
   */
  private String getValue(String label) {
    String theText = Macro.getValue(macroOptions, label, null);
    if (theText != null
        && (theText.startsWith("&") || label.toLowerCase(Locale.US).startsWith(theText))) {
      // Is the value a macro variable?
      if (theText.startsWith("&")) {
        theText = theText.substring(1);
      }
      final Interpreter interp = Interpreter.getInstance();
      final String s = interp != null ? interp.getVariableAsString(theText) : null;
      if (s != null) {
        theText = s;
      }
    }
    return theText;
  }
}
