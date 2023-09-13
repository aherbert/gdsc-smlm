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

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import uk.ac.sussex.gdsc.core.data.utils.Rounder;
import uk.ac.sussex.gdsc.core.data.utils.RounderUtils;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.results.ArrayPeakResultStore;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

/**
 * Stores traced peak results and allows event propagation to listeners of the model.
 */
public class TraceDataTableModel extends AbstractTableModel {
  private static final long serialVersionUID = 20230906L;

  /** Count used to count how many PeakResultTableModelFrame are displaying this model. */
  private final AtomicInteger liveCount = new AtomicInteger();

  /** The flag used to indicate the columns require updating. */
  private final AtomicBoolean columnsComputed = new AtomicBoolean(false);

  private final List<TraceData> data;
  private final Calibration calibration;
  private final PSF psf;
  private ResultsTableSettings tableSettings;

  // These depend on the source results
  private ImageSource source;
  private String configuration = "";
  private final boolean showCategory;
  private final boolean showZ;

  // Used for the columns
  /** The rounder for the numeric cell entries. */
  Rounder rounder;
  private Object[] values;
  private Class<?>[] valuesTypes;
  private String[] names;
  private boolean rowCounter;

  private Consumer<ResultsTableSettings> settingsUpdatedAction;


  /**
   * Create an instance from the resultsSource.
   *
   * @param resultsSource the results source
   * @param tableSettings the table settings
   */
  public TraceDataTableModel(MemoryPeakResults resultsSource, ResultsTableSettings tableSettings) {
    if (tableSettings == null) {
      tableSettings = ResultsProtosHelper.DefaultResultsSettings.INSTANCE.getResultsTableSettings();
    }
    this.tableSettings = tableSettings;

    this.calibration = resultsSource.getCalibration();
    this.psf = resultsSource.getPsf();

    this.showZ = resultsSource.is3D();
    this.showCategory = resultsSource.hasCategory();
    this.source = resultsSource.getSource();
    this.configuration = resultsSource.getConfiguration();

    data = extractTraceData(resultsSource, showZ);
  }

  private static List<TraceData> extractTraceData(MemoryPeakResults results, boolean showZ) {
    final Trace[] traces = TraceManager.convert(results);
    return Arrays.stream(traces).map(t -> new TraceData(t, showZ)).collect(Collectors.toList());
  }

  // *************************************************************************/
  // Table column appearance
  // *************************************************************************/

  /**
   * Gets the table settings.
   *
   * @return the table settings
   */
  ResultsTableSettings getTableSettings() {
    return tableSettings;
  }

  /**
   * Sets the table settings.
   *
   * @param tableSettings the new table settings
   */
  void setTableSettings(ResultsTableSettings tableSettings) {
    if (tableSettings == null) {
      tableSettings = ResultsProtosHelper.DefaultResultsSettings.INSTANCE.getResultsTableSettings();
    }
    if (equals(this.tableSettings, tableSettings)) {
      return;
    }
    this.tableSettings = tableSettings;
    createTableStructure(true);
    if (settingsUpdatedAction != null) {
      settingsUpdatedAction.accept(tableSettings);
    }
  }

  /**
   * Add an action to perform if the table settings have changed. This is combined with any previous
   * actions.
   *
   * @param action the action
   */
  public void addSettingsUpdatedAction(Consumer<ResultsTableSettings> action) {
    if (settingsUpdatedAction == null) {
      this.settingsUpdatedAction = action;
    } else {
      settingsUpdatedAction = settingsUpdatedAction.andThen(action);
    }
  }

  /**
   * Check if the settings are equal.
   *
   * @param t1 the t1
   * @param t2 the t2
   * @return true, if successful
   */
  private static boolean equals(ResultsTableSettings t1, ResultsTableSettings t2) {
    // Adapted from ResultsTableSettings.equals() to only use the settings of interest
    boolean result = t1.getDistanceUnitValue() == t2.getDistanceUnitValue();
    result = result && (t1.getShowTracePrecision() == t2.getShowTracePrecision());
    result =
        result && (t1.getShowTraceDiffusionCoefficient() == t2.getShowTraceDiffusionCoefficient());
    result = result && (t1.getRoundingPrecision() == t2.getRoundingPrecision());
    result = result && (t1.getShowRowCounter() == t2.getShowRowCounter());
    return result;
  }

  /**
   * Called when the structure of the table should be created. If the structure has not changed or
   * no live table are attached then this does nothing.
   *
   * @param changed Set to true if a property controlling the structure has changed
   */
  private void createTableStructure(boolean changed) {
    if (changed) {
      columnsComputed.set(false);
    }

    // If no live attached tables or nothing has changed then return
    // Note: We use a live count instead of the listener count (see listenerList.getListenerCount())
    // so it can be turned off by PeakResultTableModelFrame.
    if (liveCount.get() == 0 || columnsComputed.get()) {
      return;
    }

    columnsComputed.set(true);

    rounder = RounderUtils.create(tableSettings.getRoundingPrecision());

    // Create the converters
    // Handle null calibration
    final CalibrationReader cr =
        new CalibrationReader(calibration == null ? Calibration.getDefaultInstance() : calibration);
    final TypeConverter<DistanceUnit> distanceConverter =
        cr.getDistanceConverterSafe(tableSettings.getDistanceUnit());
    // Use the same values for precision and D
    final TypeConverter<DistanceUnit> nmConverter = cr.getDistanceConverterSafe(DistanceUnit.NM);
    final TypeConverter<DistanceUnit> umConverter = cr.getDistanceConverterSafe(DistanceUnit.UM);
    final TypeConverter<TimeUnit> secConverter = cr.getTimeConverterSafe(TimeUnit.SECOND);

    // Organise the data columns.
    final LocalList<Function<TraceData, Object>> valuesList = new LocalList<>();
    final LocalList<Class<?>> typesList = new LocalList<>();
    final LocalList<String> namesList = new LocalList<>();

    rowCounter = tableSettings.getShowRowCounter();
    if (rowCounter) {
      addColumn(valuesList, typesList, namesList, null, Integer.class, "#");
    }

    addColumn(valuesList, typesList, namesList, t -> t.getTrace().getId(), Integer.class, "id");
    addColumn(valuesList, typesList, namesList, t -> t.getTrace().size(), Integer.class, "size");
    if (showCategory) {
      addColumn(valuesList, typesList, namesList, t -> t.getCategoryString(), String.class,
          "category");
    }
    addColumn(valuesList, typesList, namesList, t -> distanceConverter.convert(t.getCx()),
        Double.class, "cx", UnitHelper.getShortName(distanceConverter.to()));
    addColumn(valuesList, typesList, namesList, t -> distanceConverter.convert(t.getCy()),
        Double.class, "cy", UnitHelper.getShortName(distanceConverter.to()));
    if (showZ) {
      addColumn(valuesList, typesList, namesList, t -> distanceConverter.convert(t.getCz()),
          Double.class, "cz", UnitHelper.getShortName(distanceConverter.to()));
    }
    addColumn(valuesList, typesList, namesList, t -> t.getTrace().getHead().getFrame(),
        Integer.class, "start");
    addColumn(valuesList, typesList, namesList, t -> t.getTrace().getTail().getEndFrame(),
        Integer.class, "end");

    if (tableSettings.getShowTracePrecision()) {
      addColumn(valuesList, typesList, namesList, t -> nmConverter.convert(t.getSx()), Double.class,
          "sx", UnitHelper.getShortName(nmConverter.to()));
      addColumn(valuesList, typesList, namesList, t -> nmConverter.convert(t.getSy()), Double.class,
          "sy", UnitHelper.getShortName(nmConverter.to()));
      if (showZ) {
        addColumn(valuesList, typesList, namesList, t -> nmConverter.convert(t.getSz()),
            Double.class, "sz", UnitHelper.getShortName(nmConverter.to()));
      }
    }

    if (tableSettings.getShowTraceDiffusionCoefficient()) {
      addColumn(valuesList, typesList, namesList,
          t -> secConverter.convert(umConverter.convert(umConverter.convert(t.getMsd() / 4))),
          Double.class, "D", UnitHelper.getShortName(umConverter.to()) + "^2/"
              + UnitHelper.getShortName(secConverter.to()));
    }

    // unchecked cast
    values = valuesList.toArray();
    valuesTypes = typesList.toArray(new Class<?>[0]);
    names = namesList.toArray(new String[0]);

    fireTableStructureChanged();
  }

  private static void addColumn(LocalList<Function<TraceData, Object>> valuesList,
      LocalList<Class<?>> typesList, LocalList<String> namesList, Function<TraceData, Object> fun,
      Class<?> cls, String name) {
    valuesList.add(fun);
    typesList.add(cls);
    namesList.add(name);
  }

  private static void addColumn(LocalList<Function<TraceData, Object>> valuesList,
      LocalList<Class<?>> typesList, LocalList<String> namesList, Function<TraceData, Object> fun,
      Class<?> cls, String name, String unitName) {
    valuesList.add(fun);
    typesList.add(cls);
    addName(namesList, name, unitName);
  }

  private static void addName(LocalList<String> namesList, String name, String unitName) {
    if (!TextUtils.isNullOrEmpty(unitName)) {
      name += " (" + unitName + ")";
    }
    namesList.add(name);
  }

  // *************************************************************************/
  // Table model methods
  // *************************************************************************/

  @Override
  public int getRowCount() {
    return data.size();
  }

  @Override
  public String getColumnName(int columnIndex) {
    return names[columnIndex];
  }

  @Override
  public Class<?> getColumnClass(int columnIndex) {
    return valuesTypes[columnIndex];
  }

  @Override
  public int getColumnCount() {
    return values.length;
  }

  @SuppressWarnings("unchecked")
  @Override
  public Object getValueAt(int rowIndex, int columnIndex) {
    if (rowCounter && columnIndex == 0) {
      return Integer.valueOf(rowIndex + 1);
    }
    final TraceData r = get(rowIndex);
    // unchecked cast is allowed as we created the array
    return ((Function<TraceData, Object>) values[columnIndex]).apply(r);
  }

  // *************************************************************************/
  // Data management
  // *************************************************************************/

  /**
   * To memory peak results.
   *
   * @return the memory peak results
   */
  MemoryPeakResults toMemoryPeakResults() {
    final ArrayPeakResultStore store = new ArrayPeakResultStore(100);
    data.stream().map(t -> t.getTrace().getPoints()).forEach(store::addStore);
    return toMemoryPeakResults(store);
  }

  /**
   * To memory peak results.
   *
   * @param indices the indices
   * @return the memory peak results
   */
  MemoryPeakResults toMemoryPeakResults(int[] indices) {
    final ArrayPeakResultStore store = new ArrayPeakResultStore(100);
    for (final int i : indices) {
      store.addStore(data.get(i).getTrace().getPoints());
    }
    return toMemoryPeakResults(store);
  }

  private MemoryPeakResults toMemoryPeakResults(final ArrayPeakResultStore store) {
    final MemoryPeakResults results = new MemoryPeakResults(store);
    results.setPsf(psf);
    results.setCalibration(calibration);
    results.setSource(source);
    results.setConfiguration(configuration);
    return results;
  }

  /**
   * Gets the results for the given index.
   *
   * @param index the index
   * @return the peak result
   */
  TraceData get(int index) {
    return data.get(index);
  }

  /**
   * Removes the result.
   *
   * @param source the source
   * @param index the index
   */
  void remove(Object source, int index) {
    if (index < 0 || index >= data.size()) {
      return;
    }
    data.remove(index);
    fireTableRowsDeleted(index, index);
  }

  /**
   * Removes the result.
   *
   * @param source the source
   * @param indices the indices
   */
  void remove(Object source, int[] indices) {
    if (indices == null || indices.length == 0) {
      return;
    }

    if (indices.length == 1) {
      remove(source, indices[0]);
      return;
    }

    int size = 0;
    for (final int index : indices) {
      if (index < 0 || index >= data.size()) {
        continue;
      }
      indices[size++] = index;
    }

    if (size == 0) {
      return;
    }
    if (size < indices.length) {
      indices = Arrays.copyOf(indices, size);
    }

    final int[] pairs = SimpleArrayUtils.getRanges(indices);

    size = pairs.length;
    final int firstRow = pairs[0];
    final int lastRow = pairs[size - 1];

    // Remove ranges starting at the end (to preserve the list order)
    for (int i = size - 1; i > 0; i -= 2) {
      data.subList(pairs[i - 1], pairs[i] + 1).clear();
    }

    fireTableRowsDeleted(firstRow, lastRow);
  }

  /**
   * Clear the results.
   *
   * @param source the source
   */
  void clear(Object source) {
    final int index1 = data.size() - 1;
    if (index1 >= 0) {
      data.clear();
      fireTableRowsDeleted(0, index1);
    }
  }

  /**
   * Deletes the components at the specified range of indexes. The removal is inclusive, so
   * specifying a range of (1,5) removes the component at index 1 and the component at index 5, as
   * well as all components in between.
   *
   * <p>Throws an {@code ArrayIndexOutOfBoundsException} if the index was invalid.
   *
   * @param source the source
   * @param fromIndex the index of the lower end of the range
   * @param toIndex the index of the upper end of the range
   */
  void removeRange(Object source, int fromIndex, int toIndex) {
    data.subList(fromIndex, toIndex + 1).clear();
    fireTableRowsDeleted(fromIndex, toIndex);
  }

  /**
   * Performs the given action for each element of the data.
   *
   * @param action the action
   */
  void forEach(Consumer<? super TraceData> action) {
    data.forEach(action);
  }

  // *************************************************************************/
  // Properties
  // *************************************************************************/

  /**
   * Sets the model to the live state. This creates all the table layout information and causes it
   * to update when properties are changed.
   *
   * @param isLive the new live
   */
  void setLive(boolean isLive) {
    if (isLive) {
      liveCount.getAndIncrement();
    } else {
      liveCount.getAndDecrement();
    }
    createTableStructure(false);
  }

  /**
   * Checks if is live.
   *
   * @return true, if is live
   */
  boolean isLive() {
    return liveCount.get() != 0;
  }

  /**
   * Checks if showing the results category in the table.
   *
   * @return If true show the results category in the table.
   */
  boolean isShowCategory() {
    return showCategory;
  }

  /**
   * Checks if showing the Z column.
   *
   * @return true, if is show Z
   */
  boolean isShowZ() {
    return showZ;
  }

  /**
   * Sets the source.
   *
   * @param source the new source
   */
  void setSource(ImageSource source) {
    // This does not require a change to the table structure
    this.source = source;
  }

  /**
   * Gets the source.
   *
   * @return the source
   */
  ImageSource getSource() {
    return source;
  }

  /**
   * Gets the calibration.
   *
   * @return the calibration
   */
  Calibration getCalibration() {
    return calibration;
  }

  /**
   * Gets the psf.
   *
   * @return the psf
   */
  PSF getPsf() {
    return psf;
  }

  // *************************************************************************/
  // Table cell appearance (Rendering)
  // *************************************************************************/

  /**
   * Gets the double renderer.
   *
   * @return the double renderer
   */
  TableCellRenderer getDoubleRenderer() {
    return PeakResultTableModel.createTableCellRenderer(new DefaultTableCellRenderer() {
      private static final long serialVersionUID = 20230906L;

      @Override
      protected void setValue(Object value) {
        setText((value == null) ? "" : rounder.toString((Double) value));
      }
    });
  }

  /**
   * Gets the integer renderer.
   *
   * @return the integer renderer
   */
  TableCellRenderer getIntegerRenderer() {
    return PeakResultTableModel.createTableCellRenderer(new DefaultTableCellRenderer() {
      private static final long serialVersionUID = 20230906L;
    });
  }
}
