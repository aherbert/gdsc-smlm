/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.gui;

import ij.WindowManager;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Arrays;
import java.util.Collections;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.RowSorter;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ScreenDimensionHelper;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SummariseResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.ArrayPeakResultStore;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * A frame that shows a TraceDatasTableModel.
 */
public class TraceDataTableModelFrame extends JFrame {
  private static final long serialVersionUID = 20230906L;

  private final TraceDataTableModelJTable table;
  private JMenuItem resultsSave;
  private JMenuItem resultsShowInfo;
  private JCheckBoxMenuItem editReadOnly;
  private JMenuItem editDelete;
  private JMenuItem editDeleteAll;
  private JMenuItem editSelectAll;
  private JMenuItem editSelectNone;
  private JMenuItem editUnsort;
  private JMenuItem editTableSettings;
  private JMenuItem sourceAttachImage;
  private JMenuItem sourceShowInfo;
  private JMenuItem sourceShowImage;
  private JMenuItem sourceOverlay;
  private String saveName;

  /**
   * Create an instance.
   *
   * @param model the model
   */
  public TraceDataTableModelFrame(TraceDataTableModel model) {
    this(model, null, null);
  }

  /**
   * Create an instance.
   *
   * @param model the model
   * @param columnModel the column model
   * @param selectionModel the selection model
   */
  private TraceDataTableModelFrame(final TraceDataTableModel model, TableColumnModel columnModel,
      ListSelectionModel selectionModel) {
    setJMenuBar(createMenuBar());

    // This is required to get the column sizes for the model data.
    model.setLive(true);

    final int[] indices =
        selectionModel == null ? null : ListSelectionModelHelper.getSelectedIndices(selectionModel);

    table = new TraceDataTableModelJTable(model, columnModel, selectionModel);

    if (indices != null) {
      ListSelectionModelHelper.setSelectedIndices(selectionModel, indices);
    }

    final JScrollPane scroll = new JScrollPane(table);

    final ScreenDimensionHelper helper = new ScreenDimensionHelper();
    helper.setMinHeight(300);
    helper.setup(scroll);

    add(scroll);
    pack();

    // If the window is never set visible then we cannot ensure this happens in
    // the closing event so do it here and again when opened.
    model.setLive(false);

    addWindowListener(new WindowAdapter() {
      @Override
      public void windowOpened(WindowEvent event) {
        model.setLive(true);
        WindowManager.addWindow(TraceDataTableModelFrame.this);
      }

      @Override
      public void windowClosing(WindowEvent event) {
        model.setLive(false);
        WindowManager.removeWindow(TraceDataTableModelFrame.this);
      }
    });
  }

  private JMenuBar createMenuBar() {
    final JMenuBar menubar = new JMenuBar();
    menubar.add(createResultsMenu());
    menubar.add(createEditMenu());
    menubar.add(createSourceMenu());
    return menubar;
  }

  private JMenu createResultsMenu() {
    final JMenu menu = new JMenu("Results");
    menu.setMnemonic(KeyEvent.VK_R);
    menu.add(resultsSave = add("Save ...", KeyEvent.VK_S, "ctrl pressed S"));
    menu.add(resultsShowInfo = add("Show Info", KeyEvent.VK_I, null));
    return menu;
  }

  private JMenu createEditMenu() {
    final JMenu menu = new JMenu("Edit");
    menu.setMnemonic(KeyEvent.VK_E);
    menu.add(editReadOnly = addToggle("Read-only", KeyEvent.VK_R, null, false));
    menu.add(editDelete = add("Delete", KeyEvent.VK_D, "DELETE"));
    menu.add(editDeleteAll = add("Delete All", KeyEvent.VK_A, null));
    menu.addSeparator();
    menu.add(editSelectNone = add("Select None", KeyEvent.VK_N, "ctrl shift pressed A"));
    menu.add(editSelectAll = add("Select All", KeyEvent.VK_S, null));
    menu.add(editUnsort = add("Unsort", KeyEvent.VK_U, null));
    menu.addSeparator();
    menu.add(editTableSettings = add("Table Settings ...", KeyEvent.VK_T, "ctrl pressed T"));

    // Read-only by default
    editReadOnly.setSelected(true);
    editDelete.setEnabled(false);
    editDeleteAll.setEnabled(false);
    editReadOnly.addActionListener(event -> {
      final boolean allowDelete = !isReadOnly();
      editDelete.setEnabled(allowDelete);
      editDeleteAll.setEnabled(allowDelete);
    });

    return menu;
  }

  private JMenu createSourceMenu() {
    final JMenu menu = new JMenu("Source");
    menu.setMnemonic(KeyEvent.VK_S);
    menu.add(sourceAttachImage = add("Attach Image ...", KeyEvent.VK_A, null));
    menu.add(sourceShowInfo = add("Show Info", KeyEvent.VK_I, "ctrl pressed I"));
    menu.add(sourceShowImage = add("Show Image", KeyEvent.VK_S, "ctrl shift pressed I"));
    menu.add(sourceOverlay = add("Overlay", KeyEvent.VK_O, "ctrl pressed Y"));
    return menu;
  }

  private JMenuItem add(String text, int mnemonic, String keyStroke) {
    final JMenuItem item = new JMenuItem(text, mnemonic);
    if (keyStroke != null) {
      item.setAccelerator(KeyStroke.getKeyStroke(keyStroke));
    }
    item.addActionListener(this::actionPerformed);
    return item;
  }

  private JCheckBoxMenuItem addToggle(String text, int mnemonic, String keyStroke,
      boolean selected) {
    final JCheckBoxMenuItem item = new JCheckBoxMenuItem(text, selected);
    item.setMnemonic(mnemonic);
    if (keyStroke != null) {
      item.setAccelerator(KeyStroke.getKeyStroke(keyStroke));
    }
    item.addActionListener(this::actionPerformed);
    return item;
  }

  private void actionPerformed(ActionEvent event) {
    final Object src = event.getSource();
    if (src == resultsSave) {
      doResultsSave();
    } else if (src == resultsShowInfo) {
      doResultsShowInfo();
    } else if (src == editDelete) {
      doEditDelete();
    } else if (src == editDeleteAll) {
      doEditDeleteAll();
    } else if (src == editSelectNone) {
      doEditSelectNone();
    } else if (src == editSelectAll) {
      doEditSelectAll();
    } else if (src == editUnsort) {
      doEditUnsort();
    } else if (src == editTableSettings) {
      doEditTableSettings();
    } else if (src == sourceAttachImage) {
      doSourceAttachImage();
    } else if (src == sourceShowInfo) {
      doSourceShowInfo();
    } else if (src == sourceShowImage) {
      doSourceShowImage();
    } else if (src == sourceOverlay) {
      doSourceOverlay();
    }
  }

  private void doResultsSave() {
    final TraceDataTableModel model = getModel();
    if (model == null || model.getRowCount() == 0) {
      return;
    }
    saveName = TableHelper.saveResults(getTitle(), saveName, model::toMemoryPeakResults);
  }

  private void doResultsShowInfo() {
    // Delegate to Summarise Results
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    final MemoryPeakResults results = model.toMemoryPeakResults();
    results.setName("Table data: " + getTitle());
    SummariseResults.showSummary(Collections.singletonList(results));
  }

  /**
   * Checks if is read only.
   *
   * @return true, if is read only
   */
  private boolean isReadOnly() {
    return editReadOnly.isSelected();
  }

  private void doEditDelete() {
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    final int[] indices = table.getSelectedRows();
    table.convertRowIndexToModel(indices);
    model.remove(this, indices);
  }

  private void doEditDeleteAll() {
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    model.clear(this);
  }

  private void doEditSelectNone() {
    table.clearSelection();
  }

  private void doEditSelectAll() {
    table.selectAll();
  }

  private void doEditUnsort() {
    final RowSorter<?> rs = table.getRowSorter();
    if (rs != null) {
      rs.setSortKeys(null);
    }
  }

  private void doEditTableSettings() {
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    final ResultsTableSettings.Builder tableSettings = model.getTableSettings().toBuilder();
    // Copied from ResultsManager.addTableResultsOptions
    final ExtendedGenericDialog egd = new ExtendedGenericDialog("Table Settings", this);
    egd.addChoice("Table_distance_unit", SettingsManager.getDistanceUnitNames(),
        tableSettings.getDistanceUnitValue());
    egd.addCheckbox("Table_show_precision", tableSettings.getShowTracePrecision());
    egd.addCheckbox("Table_show_D", tableSettings.getShowTraceDiffusionCoefficient());
    egd.addSlider("Table_precision", 0, 10, tableSettings.getRoundingPrecision());
    egd.addCheckbox("Table_show_counter", tableSettings.getShowRowCounter());
    egd.showDialog();
    if (egd.wasCanceled()) {
      return;
    }
    tableSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
    tableSettings.setShowTracePrecision(egd.getNextBoolean());
    tableSettings.setShowTraceDiffusionCoefficient(egd.getNextBoolean());
    tableSettings.setRoundingPrecision((int) egd.getNextNumber());
    tableSettings.setShowRowCounter(egd.getNextBoolean());
    model.setTableSettings(tableSettings.build());
  }

  private void doSourceAttachImage() {
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    TableHelper.updateSource(model::setSource);
  }

  private void doSourceShowInfo() {
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    TableHelper.showInfo(getTitle(), model.getSource());
  }

  private void doSourceShowImage() {
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    TableHelper.showImage(model.getSource());
  }

  private void doSourceOverlay() {
    final TraceDataTableModel model = getModel();
    if (model == null) {
      return;
    }
    final TraceData[] list = table.getSelectedData();
    if (list.length == 0) {
      return;
    }
    final ImageSource source = model.getSource();
    if (source == null) {
      return;
    }
    final ArrayPeakResultStore store = new ArrayPeakResultStore(100);
    Arrays.stream(list).map(t -> t.getTrace().getPoints()).forEach(store::addStore);
    TableHelper.showOverlay(model.getSource(), model.getCalibration(), store.toArray());
  }

  private TraceDataTableModel getModel() {
    final TableModel model = table.getModel();
    return (model instanceof TraceDataTableModel) ? (TraceDataTableModel) model : null;
  }
}
