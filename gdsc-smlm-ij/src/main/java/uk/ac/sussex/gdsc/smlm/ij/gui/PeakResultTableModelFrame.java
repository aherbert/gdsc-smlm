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
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
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
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * A frame that shows a PeakResultsTableModel.
 */
public class PeakResultTableModelFrame extends JFrame implements ActionListener {
  private static final long serialVersionUID = -3671174621388288975L;

  private final PeakResultTableModelJTable table;
  private JMenuItem resultsSave;
  private JMenuItem resultsSaveSelection;
  private JMenuItem resultsShowInfo;
  private JCheckBoxMenuItem editReadOnly;
  private JMenuItem editDelete;
  private JMenuItem editDeleteAll;
  private JMenuItem editSelectAll;
  private JMenuItem editSelectNone;
  private JMenuItem editSelectInverse;
  private JMenuItem editUnsort;
  private JMenuItem editTableSettings;
  private JMenuItem sourceAttachImage;
  private JMenuItem sourceShowInfo;
  private JMenuItem sourceShowImage;
  private JMenuItem sourceOverlay;
  private String saveName;

  /**
   * Instantiates a new peak result table model frame.
   *
   * @param model the model
   */
  public PeakResultTableModelFrame(PeakResultTableModel model) {
    this(model, null, null);
  }

  /**
   * Instantiates a new peak result table model frame.
   *
   * @param model the model
   * @param selectionModel the selection model
   */
  public PeakResultTableModelFrame(PeakResultTableModel model, ListSelectionModel selectionModel) {
    this(model, null, selectionModel);
  }

  /**
   * Instantiates a new peak result table model frame.
   *
   * @param model the model
   * @param columnModel the column model
   * @param selectionModel the selection model
   */
  public PeakResultTableModelFrame(final PeakResultTableModel model, TableColumnModel columnModel,
      ListSelectionModel selectionModel) {
    setJMenuBar(createMenuBar());

    // This is required to get the column sizes for the model data.
    model.setLive(true);

    final int[] indices =
        selectionModel == null ? null : ListSelectionModelHelper.getSelectedIndices(selectionModel);

    table = new PeakResultTableModelJTable(model, columnModel, selectionModel);

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
        WindowManager.addWindow(PeakResultTableModelFrame.this);
      }

      @Override
      public void windowClosing(WindowEvent event) {
        model.setLive(false);
        WindowManager.removeWindow(PeakResultTableModelFrame.this);
      }
    });
  }

  /**
   * Clean up this table. This should only be called when the table is no longer required as it
   * removes the JTable from the model listeners.
   */
  public void cleanUp() {
    // Since the models may be shared
    table.getModel().removeTableModelListener(table);
    table.getColumnModel().removeColumnModelListener(table);
    table.getSelectionModel().removeListSelectionListener(table);
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
    menu.add(
        resultsSaveSelection = add("Save Selection ...", KeyEvent.VK_E, "ctrl shift pressed S"));
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
    menu.add(editSelectInverse = add("Select Inverse", KeyEvent.VK_I, null));
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
    item.addActionListener(this);
    return item;
  }

  private JCheckBoxMenuItem addToggle(String text, int mnemonic, String keyStroke,
      boolean selected) {
    final JCheckBoxMenuItem item = new JCheckBoxMenuItem(text, selected);
    item.setMnemonic(mnemonic);
    if (keyStroke != null) {
      item.setAccelerator(KeyStroke.getKeyStroke(keyStroke));
    }
    item.addActionListener(this);
    return item;
  }

  @Override
  public void actionPerformed(ActionEvent event) {
    final Object src = event.getSource();
    if (src == resultsSave) {
      doResultsSave();
    } else if (src == resultsSaveSelection) {
      doResultsSaveSelection();
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
    } else if (src == editSelectInverse) {
      doEditSelectInverse();
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
    final PeakResultTableModel model = getModel();
    if (model == null || model.getRowCount() == 0) {
      return;
    }
    saveName = TableHelper.saveResults(getTitle(), saveName, model::toMemoryPeakResults);
  }

  private void doResultsSaveSelection() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final int[] selected = table.getSelectedRows();
    if (selected.length == 0) {
      return;
    }
    table.convertRowIndexToModel(selected);
    saveName =
        TableHelper.saveResults(getTitle(), saveName, () -> model.toMemoryPeakResults(selected));
  }

  private void doResultsShowInfo() {
    // Delegate to Summarise Results
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final MemoryPeakResults results = model.toMemoryPeakResults();
    results.setName("Table data: " + getTitle());
    SummariseResults.showSummary(Collections.singletonList(results));
  }

  /**
   * Sets the read only.
   *
   * @param readOnly the new read only
   */
  public void setReadOnly(boolean readOnly) {
    editReadOnly.setSelected(readOnly);
  }

  /**
   * Checks if is read only.
   *
   * @return true, if is read only
   */
  public boolean isReadOnly() {
    return editReadOnly.isSelected();
  }

  private void doEditDelete() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final int[] indices = table.getSelectedRows();
    table.convertRowIndexToModel(indices);
    model.remove(this, indices);
  }

  private void doEditDeleteAll() {
    final PeakResultTableModel model = getModel();
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

  private void doEditSelectInverse() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    ListSelectionModelHelper.invertSelection(model.getRowCount(), table.getSelectionModel());
  }

  private void doEditUnsort() {
    final RowSorter<?> rs = table.getRowSorter();
    if (rs != null) {
      rs.setSortKeys(null);
    }
  }

  private void doEditTableSettings() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final ResultsTableSettings.Builder tableSettings = model.getTableSettings().toBuilder();
    // Copied from ResultsManager.addTableResultsOptions
    final ExtendedGenericDialog egd = new ExtendedGenericDialog("Table Settings", this);
    egd.addChoice("Table_distance_unit", SettingsManager.getDistanceUnitNames(),
        tableSettings.getDistanceUnitValue());
    egd.addChoice("Table_intensity_unit", SettingsManager.getIntensityUnitNames(),
        tableSettings.getIntensityUnitValue());
    egd.addChoice("Table_angle_unit", SettingsManager.getAngleUnitNames(),
        tableSettings.getAngleUnitValue());
    egd.addCheckbox("Table_show_fitting_data", tableSettings.getShowFittingData());
    egd.addCheckbox("Table_show_noise_data", tableSettings.getShowNoiseData());
    egd.addCheckbox("Table_show_precision", tableSettings.getShowPrecision());
    egd.addCheckbox("Table_show_deviations", model.isShowDeviations());
    egd.addCheckbox("Table_show_end_frame", model.isShowEndFrame());
    egd.addSlider("Table_precision", 0, 10, tableSettings.getRoundingPrecision());
    egd.addCheckbox("Table_show_counter", tableSettings.getShowRowCounter());
    egd.showDialog();
    if (egd.wasCanceled()) {
      return;
    }
    tableSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
    tableSettings.setIntensityUnitValue(egd.getNextChoiceIndex());
    tableSettings.setAngleUnitValue(egd.getNextChoiceIndex());
    tableSettings.setShowFittingData(egd.getNextBoolean());
    tableSettings.setShowNoiseData(egd.getNextBoolean());
    tableSettings.setShowPrecision(egd.getNextBoolean());
    model.setShowDeviations(egd.getNextBoolean());
    model.setShowEndFrame(egd.getNextBoolean());
    tableSettings.setRoundingPrecision((int) egd.getNextNumber());
    tableSettings.setShowRowCounter(egd.getNextBoolean());
    model.setTableSettings(tableSettings.build());
  }

  private void doSourceAttachImage() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    TableHelper.updateSource(model::setSource);
  }

  private void doSourceShowInfo() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    TableHelper.showInfo(getTitle(), model.getSource());
  }

  private void doSourceShowImage() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    TableHelper.showImage(model.getSource());
  }

  private void doSourceOverlay() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final ImageSource source = model.getSource();
    if (source == null) {
      return;
    }
    TableHelper.showOverlay(source, model.getCalibration(), table.getSelectedData());
  }

  private PeakResultTableModel getModel() {
    final TableModel model = table.getModel();
    return (model instanceof PeakResultTableModel) ? (PeakResultTableModel) model : null;
  }

  /**
   * Maps the index of the row in terms of the view to the underlying {@code TableModel}. If the
   * contents of the model are not sorted the model and view indices are the same.
   *
   * @param viewRowIndex the index of the row in the view
   * @return the index of the corresponding row in the model
   * @throws IndexOutOfBoundsException if sorting is enabled and passed an index outside the range
   *         of the {@code JTable} as determined by the method {@code getRowCount}
   * @see javax.swing.table.TableRowSorter
   * @see javax.swing.JTable#getRowCount()
   */
  public int convertRowIndexToModel(int viewRowIndex) {
    return table.convertRowIndexToModel(viewRowIndex);
  }

  /**
   * Returns the location of {@code index} in terms of the underlying model. That is, for the row
   * {@code index} in the coordinates of the view this returns the row index in terms of the
   * underlying model.
   *
   * @param indices the indices (updated in-place)
   * @throws IndexOutOfBoundsException if {@code index} is outside the range of the view
   */
  public void convertRowIndexToModel(int[] indices) {
    table.convertRowIndexToModel(indices);
  }

  /**
   * Maps the index of the row in terms of the {@code TableModel} to the view. If the contents of
   * the model are not sorted the model and view indices are the same.
   *
   * @param modelRowIndex the index of the row in terms of the model
   * @return the index of the corresponding row in the view, or -1 if the row isn't visible
   * @throws IndexOutOfBoundsException if sorting is enabled and passed an index outside the number
   *         of rows of the {@code TableModel}
   * @see javax.swing.table.TableRowSorter
   */
  public int convertRowIndexToView(int modelRowIndex) {
    return table.convertRowIndexToView(modelRowIndex);
  }

  /**
   * Returns the location of {@code index} in terms of the view. That is, for the row {@code index}
   * in the coordinates of the underlying model this returns the row index in terms of the view.
   *
   * @param indices the indices (updated in-place)
   * @throws IndexOutOfBoundsException if {@code index} is outside the range of the model
   */
  public void convertRowIndexToView(int[] indices) {
    table.convertRowIndexToView(indices);
  }
}
