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

package uk.ac.sussex.gdsc.smlm.ij.gui;

import gnu.trove.list.array.TFloatArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Arrays;
import javax.swing.DefaultListSelectionModel;
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
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ScreenDimensionHelper;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.XmlUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SplitMix;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.SeriesImageSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.TiffSeriesViewer.TiffSeriesVirtualStack;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.ArrayPeakResultStore;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.ImageSource.ReadHint;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList;
import uk.ac.sussex.gdsc.smlm.results.sort.FramePeakResultComparator;

/**
 * A frame that shows a PeakResultsTableModel.
 */
public class PeakResultTableModelFrame extends JFrame implements ActionListener {
  private static final long serialVersionUID = -3671174621388288975L;

  private final PeakResultTableModelJTable table;
  private JMenuItem fileSave;
  private JCheckBoxMenuItem editReadOnly;
  private JMenuItem editDelete;
  private JMenuItem editDeleteAll;
  private JMenuItem editSelectAll;
  private JMenuItem editSelectNone;
  private JMenuItem editUnsort;
  private JMenuItem editTableSettings;
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
        (selectionModel != null) ? ListSelectionModelHelper.getSelectedIndices(selectionModel)
            : null;

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
    menubar.add(createFileMenu());
    menubar.add(createEditMenu());
    menubar.add(createSourceMenu());
    return menubar;
  }

  private JMenu createFileMenu() {
    final JMenu menu = new JMenu("File");
    menu.setMnemonic(KeyEvent.VK_F);
    menu.add(fileSave = add("Save ...", KeyEvent.VK_S, "ctrl pressed S"));
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
    if (src == fileSave) {
      doFileSave();
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
    } else if (src == sourceShowInfo) {
      doSourceShowInfo();
    } else if (src == sourceShowImage) {
      doSourceShowImage();
    } else if (src == sourceOverlay) {
      doSourceOverlay();
    }
  }

  private void doFileSave() {
    final PeakResultTableModel model = getModel();
    if (model == null || model.getRowCount() == 0) {
      return;
    }
    final ExtendedGenericDialog gd = new ExtendedGenericDialog("Save Results", this);
    if (TextUtils.isNullOrEmpty(saveName)) {
      saveName = getTitle();
    }
    gd.addStringField("Results_set_name", saveName, 30);
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    saveName = gd.getNextString();
    if (TextUtils.isNullOrEmpty(saveName)) {
      IJ.error("No results set name");
      return;
    }
    final MemoryPeakResults results = model.toMemoryPeakResults();
    results.setName(saveName);
    MemoryPeakResults.addResults(results);
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
    tableSettings.setRoundingPrecision((int) egd.getNextNumber());
    tableSettings.setShowRowCounter(egd.getNextBoolean());
    model.setTableSettings(tableSettings.build());
  }

  private void doSourceShowInfo() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final ImageSource source = model.getSource();
    String text = getTitle() + " source";
    if (source == null) {
      text += " = NA";
    } else {
      text += "\n" + XmlUtils.prettyPrintXml(source.toXml());
    }
    // Note:
    // If a raw path is printed to the ImageJ log double-clicking it will open the image.
    // We could separate these onto multiple lines:
    // <path>/path/to/image.tif</path>
    // <string>/path/to/image.tif</string>
    IJ.log(text);
  }

  private void doSourceShowImage() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final ImageSource source = model.getSource();
    if (source == null) {
      return;
    }
    // Check if already open
    final ImagePlus imp = WindowManager.getImage(source.getName());
    if (imp != null) {
      imp.getWindow().toFront();
      return;
    }

    // Check if an ImageJ image source
    if (source instanceof IJImageSource) {
      final IJImageSource imageSource = (IJImageSource) source;
      final String path = imageSource.getPath();
      if (path != null && new File(path).exists()) {
        IJ.showStatus("Opening image ...");
        IJ.open(path);
        IJ.showStatus("");
      } else {
        IJ.log("Cannot find the image source: " + path);
      }
      return;
    }
    // Open a SeriesImageSource.
    if (source instanceof SeriesImageSource) {
      final SeriesImageSource imageSource = (SeriesImageSource) source;
      imageSource.setBufferLimit(0); // No memory buffer
      imageSource.setReadHint(ReadHint.NONSEQUENTIAL);
      if (!source.open()) {
        IJ.log("Cannot open the series image source");
        return;
      }
      new TiffSeriesVirtualStack(imageSource).show();
    }
  }

  private void doSourceOverlay() {
    final PeakResultTableModel model = getModel();
    if (model == null) {
      return;
    }
    final PeakResult[] list = table.getSelectedData();
    if (list.length == 0) {
      return;
    }
    final ImageSource source = model.getSource();
    if (source == null) {
      return;
    }

    final String title = source.getOriginal().getName();
    final ImagePlus imp = WindowManager.getImage(title);
    if (imp == null) {
      return;
    }
    // Assumes 3D stack (no channel/time)
    if (imp.getNDimensions() > 3) {
      return;
    }
    try {
      final TypeConverter<DistanceUnit> converter =
          CalibrationHelper.getDistanceConverter(model.getCalibration(), DistanceUnit.PIXEL);
      final Overlay o = new Overlay();

      if (list.length == 1) {
        final PeakResult p = list[0];
        final PointRoi roi =
            new PointRoi(converter.convert(p.getXPosition()), converter.convert(p.getYPosition()));
        roi.setPointType(3);
        roi.setPosition(p.getFrame());
        o.add(roi);
      } else {
        Arrays.sort(list, FramePeakResultComparator.INSTANCE);
        final TFloatArrayList ox = new TFloatArrayList(list.length);
        final TFloatArrayList oy = new TFloatArrayList(list.length);
        int frame = list[0].getFrame() - 1;
        for (int i = 0; i < list.length; i++) {
          if (frame != list[i].getFrame()) {
            if (ox.size() > 0) {
              final PointRoi roi = new PointRoi(ox.toArray(), oy.toArray());
              roi.setPointType(3);
              roi.setPosition(frame);
              ox.resetQuick();
              oy.resetQuick();
              o.add(roi);
            }
            frame = list[i].getFrame();
          }
          ox.add(converter.convert(list[i].getXPosition()));
          oy.add(converter.convert(list[i].getYPosition()));
        }
        if (ox.size() > 0) {
          final PointRoi roi = new PointRoi(ox.toArray(), oy.toArray());
          roi.setPointType(3);
          roi.setPosition(frame);
          o.add(roi);
        }
      }
      imp.setOverlay(o);
      final PeakResult p = list[0];
      imp.setSlice(p.getFrame());
      ImageJUtils.adjustSourceRect(imp, 0, (int) converter.convert(p.getXPosition()),
          (int) converter.convert(p.getYPosition()));
      imp.getWindow().toFront();
    } catch (final ConversionException ex) {
      // Ignore
    }
  }

  private PeakResultTableModel getModel() {
    final TableModel model = table.getModel();
    return (model instanceof PeakResultTableModel) ? (PeakResultTableModel) model : null;
  }

  /**
   * Maps the index of the row in terms of the view to the underlying <code>TableModel</code>. If
   * the contents of the model are not sorted the model and view indices are the same.
   *
   * @param viewRowIndex the index of the row in the view
   * @return the index of the corresponding row in the model
   * @throws IndexOutOfBoundsException if sorting is enabled and passed an index outside the range
   *         of the <code>JTable</code> as determined by the method <code>getRowCount</code>
   * @see javax.swing.table.TableRowSorter
   * @see javax.swing.JTable#getRowCount()
   * @since 1.6
   */
  public int convertRowIndexToModel(int viewRowIndex) {
    return table.convertRowIndexToModel(viewRowIndex);
  }

  /**
   * Returns the location of <code>index</code> in terms of the underlying model. That is, for the
   * row <code>index</code> in the coordinates of the view this returns the row index in terms of
   * the underlying model.
   *
   * @param indices the indices (updated in-place)
   * @throws IndexOutOfBoundsException if <code>index</code> is outside the range of the view
   */
  public void convertRowIndexToModel(int[] indices) {
    table.convertRowIndexToModel(indices);
  }

  /**
   * Maps the index of the row in terms of the <code>TableModel</code> to the view. If the contents
   * of the model are not sorted the model and view indices are the same.
   *
   * @param modelRowIndex the index of the row in terms of the model
   * @return the index of the corresponding row in the view, or -1 if the row isn't visible
   * @throws IndexOutOfBoundsException if sorting is enabled and passed an index outside the number
   *         of rows of the <code>TableModel</code>
   * @see javax.swing.table.TableRowSorter
   * @since 1.6
   */
  public int convertRowIndexToView(int modelRowIndex) {
    return table.convertRowIndexToView(modelRowIndex);
  }

  /**
   * Returns the location of <code>index</code> in terms of the view. That is, for the row
   * <code>index</code> in the coordinates of the underlying model this returns the row index in
   * terms of the view.
   *
   * @param indices the indices (updated in-place)
   * @throws IndexOutOfBoundsException if <code>index</code> is outside the range of the model
   */
  public void convertRowIndexToView(int[] indices) {
    table.convertRowIndexToView(indices);
  }

  /**
   * Launch the application.
   *
   * @param args the arguments
   */
  public static void main(String[] args) {
    final SplitMix r = SplitMix.new64(System.currentTimeMillis());
    final int n = 20;

    final ListSelectionModel selectionModel = new DefaultListSelectionModel();

    EventQueue.invokeLater((Runnable) () -> {
      try {
        final PeakResultStoreList store = new ArrayPeakResultStore(10);
        for (int i = n; i-- > 0;) {
          store.add(new PeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(),
              r.nextDouble(), r.nextFloat(), r.nextFloat(), PeakResult.createParams(r.nextFloat(),
                  r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat()),
              null));
        }

        final CalibrationWriter cw = new CalibrationWriter();
        cw.setNmPerPixel(100);
        cw.setCountPerPhoton(10);
        cw.setDistanceUnit(DistanceUnit.PIXEL);
        cw.setIntensityUnit(IntensityUnit.COUNT);

        final ResultsTableSettings.Builder tableSettings = ResultsTableSettings.newBuilder();
        tableSettings.setDistanceUnit(DistanceUnit.NM);
        tableSettings.setIntensityUnit(IntensityUnit.PHOTON);
        tableSettings.setShowFittingData(true);
        tableSettings.setShowNoiseData(true);
        tableSettings.setShowPrecision(true);
        tableSettings.setRoundingPrecision(4);

        final PeakResultTableModel model =
            new PeakResultTableModel(store, cw.getCalibration(), null, tableSettings.build());

        final PeakResultTableModelFrame d =
            new PeakResultTableModelFrame(model, null, selectionModel);
        d.setTitle("D");
        d.setDefaultCloseOperation(EXIT_ON_CLOSE);
        d.setVisible(true);

        // Selecting in one list activates the other list

        final PeakResultTableModelFrame d2 =
            new PeakResultTableModelFrame(model, null, selectionModel);
        d2.setTitle("D2");
        // Since we have the same selection model we need the same row sorter,
        // otherwise the selection is scrambled by sorting.
        // The alternative would be to get the source for the selection event (the table)
        // and get the row sorter to do the mapping.
        // However this breaks deletion of data as the row sorter double processes the deletion.
        // Basically only one table can use the same selection model when sorting is desired.
        // d2.table.setRowSorter(d.table.getRowSorter())
        d2.setDefaultCloseOperation(EXIT_ON_CLOSE);
        d2.setVisible(true);
      } catch (final Exception ex) {
        ex.printStackTrace();
      }
    });
  }
}
