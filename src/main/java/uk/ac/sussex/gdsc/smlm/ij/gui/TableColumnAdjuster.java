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

package uk.ac.sussex.gdsc.smlm.ij.gui;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.HashMap;
import java.util.Map;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

// This class has been taken from https://tips4java.wordpress.com/2008/11/10/table-column-adjuster/

/**
 * Class to manage the widths of columns in a table.
 *
 * Various properties control how the width of the column is calculated. Another property controls
 * whether column width calculation should be dynamic. Finally, various Actions will be added to the
 * table to allow the user to customize the functionality.
 *
 * This class was designed to be used with tables that use an auto resize mode of AUTO_RESIZE_OFF.
 * With all other modes you are constrained as the width of the columns must fit inside the table.
 * So if you increase one column, one or more of the other columns must decrease. Because of this
 * the resize mode of RESIZE_ALL_COLUMNS will work the best.
 *
 * @see <a href=
 *      "https://tips4java.wordpress.com/2008/11/10/table-column-adjuster/">https://tips4java.wordpress.com/2008/11/10/table-column-adjuster/</a>
 */
public class TableColumnAdjuster implements PropertyChangeListener, TableModelListener {
  private final JTable table;
  private final int spacing;
  private boolean isColumnHeaderIncluded;
  private boolean isColumnDataIncluded;
  private boolean isOnlyAdjustLarger;
  private boolean isDynamicAdjustment;
  private final Map<TableColumn, Integer> columnSizes = new HashMap<>();

  private int maxRows;

  /**
   * Instantiates a new table column adjuster and use default spacing.
   *
   * @param table the table
   */
  public TableColumnAdjuster(JTable table) {
    this(table, 6);
  }

  /**
   * Instantiates a new table column adjuster.
   *
   * @param table the table
   * @param spacing the spacing
   */
  public TableColumnAdjuster(JTable table, int spacing) {
    this(table, 6, true);
  }

  /**
   * Instantiates a new table column adjuster.
   *
   * @param table the table
   * @param spacing the spacing
   * @param installActions the install actions flag
   */
  public TableColumnAdjuster(JTable table, int spacing, boolean installActions) {
    this.table = table;
    this.spacing = spacing;
    setColumnHeaderIncluded(true);
    setColumnDataIncluded(true);
    setOnlyAdjustLarger(false);
    setDynamicAdjustment(false);
    if (installActions) {
      installActions();
    }
  }

  /**
   * Adjust the widths of all the columns in the table.
   */
  public void adjustColumns() {
    final TableColumnModel tcm = table.getColumnModel();

    for (int i = 0; i < tcm.getColumnCount(); i++) {
      adjustColumn(i);
    }
  }

  /**
   * Adjust the width of the specified column in the table.
   *
   * @param column the column
   */
  public void adjustColumn(final int column) {
    final TableColumn tableColumn = table.getColumnModel().getColumn(column);

    if (!tableColumn.getResizable()) {
      return;
    }

    final int columnHeaderWidth = getColumnHeaderWidth(column);
    final int columnDataWidth = getColumnDataWidth(column);
    final int preferredWidth = Math.max(columnHeaderWidth, columnDataWidth);

    updateTableColumn(column, preferredWidth);
  }

  /*
   * Calculated the width based on the column name
   */
  private int getColumnHeaderWidth(int column) {
    if (!isColumnHeaderIncluded) {
      return 0;
    }

    final TableColumn tableColumn = table.getColumnModel().getColumn(column);
    final Object value = tableColumn.getHeaderValue();
    TableCellRenderer renderer = tableColumn.getHeaderRenderer();

    if (renderer == null) {
      renderer = table.getTableHeader().getDefaultRenderer();
    }

    final Component c =
        renderer.getTableCellRendererComponent(table, value, false, false, -1, column);
    return c.getPreferredSize().width;
  }

  /*
   * Calculate the width based on the widest cell renderer for the given column.
   */
  private int getColumnDataWidth(int column) {
    if (!isColumnDataIncluded) {
      return 0;
    }

    int preferredWidth = 0;
    final int maxWidth = table.getColumnModel().getColumn(column).getMaxWidth();

    // First N
    final int firstN = getRowCount();
    for (int row = firstN; row-- > 0;) {
      preferredWidth = Math.max(preferredWidth, getCellDataWidth(row, column));

      // We've exceeded the maximum width, no need to check other rows

      if (preferredWidth >= maxWidth) {
        return maxWidth;
      }
    }

    // Last N
    if (maxRows > 0) {
      for (int row = table.getRowCount(), i = maxRows; row-- > firstN && i-- > 0;) {
        preferredWidth = Math.max(preferredWidth, getCellDataWidth(row, column));

        // We've exceeded the maximum width, no need to check other rows

        if (preferredWidth >= maxWidth) {
          return maxWidth;
        }
      }
    }

    return preferredWidth;
  }

  private int getRowCount() {
    if (maxRows > 0) {
      return Math.min(maxRows, table.getRowCount());
    }
    return table.getRowCount();
  }

  /**
   * Sets the max rows.
   *
   * @param maxRows the new max rows
   */
  public void setMaxRows(int maxRows) {
    this.maxRows = maxRows;
  }

  /**
   * Gets the max rows.
   *
   * @return the max rows
   */
  public int getMaxRows() {
    return maxRows;
  }

  /*
   * Get the preferred width for the specified cell
   */
  private int getCellDataWidth(int row, int column) {
    // Inovke the renderer for the cell to calculate the preferred width

    final TableCellRenderer cellRenderer = table.getCellRenderer(row, column);
    final Component c = table.prepareRenderer(cellRenderer, row, column);
    final int width = c.getPreferredSize().width + table.getIntercellSpacing().width;

    return width;
  }

  /*
   * Update the TableColumn with the newly calculated width
   */
  private void updateTableColumn(int column, int width) {
    final TableColumn tableColumn = table.getColumnModel().getColumn(column);

    if (!tableColumn.getResizable()) {
      return;
    }

    width += spacing;

    // Don't shrink the column width

    if (isOnlyAdjustLarger) {
      width = Math.max(width, tableColumn.getPreferredWidth());
    }

    columnSizes.put(tableColumn, new Integer(tableColumn.getWidth()));

    table.getTableHeader().setResizingColumn(tableColumn);
    tableColumn.setWidth(width);
  }

  /**
   * Restore the widths of the columns in the table to its previous width.
   */
  public void restoreColumns() {
    final TableColumnModel tcm = table.getColumnModel();

    for (int i = 0; i < tcm.getColumnCount(); i++) {
      restoreColumn(i);
    }
  }

  /**
   * Restore the width of the specified column to its previous width.
   *
   * @param column the column
   */
  private void restoreColumn(int column) {
    final TableColumn tableColumn = table.getColumnModel().getColumn(column);
    final Integer width = columnSizes.get(tableColumn);

    if (width != null) {
      table.getTableHeader().setResizingColumn(tableColumn);
      tableColumn.setWidth(width.intValue());
    }
  }

  /**
   * Indicates whether to include the header in the width calculation.
   *
   * @param isColumnHeaderIncluded the new column header included
   */
  public void setColumnHeaderIncluded(boolean isColumnHeaderIncluded) {
    this.isColumnHeaderIncluded = isColumnHeaderIncluded;
  }

  /**
   * Indicates whether to include the model data in the width calculation.
   *
   * @param isColumnDataIncluded the new column data included
   */
  public void setColumnDataIncluded(boolean isColumnDataIncluded) {
    this.isColumnDataIncluded = isColumnDataIncluded;
  }

  /**
   * Indicates whether columns can only be increased in size.
   *
   * @param isOnlyAdjustLarger the new only adjust larger
   */
  public void setOnlyAdjustLarger(boolean isOnlyAdjustLarger) {
    this.isOnlyAdjustLarger = isOnlyAdjustLarger;
  }

  /**
   * Indicate whether changes to the model should cause the width to be dynamically recalculated.
   *
   * @param isDynamicAdjustment the new dynamic adjustment
   */
  public void setDynamicAdjustment(boolean isDynamicAdjustment) {
    // May need to add or remove the TableModelListener when changed

    if (this.isDynamicAdjustment != isDynamicAdjustment) {
      if (isDynamicAdjustment) {
        table.addPropertyChangeListener(this);
        table.getModel().addTableModelListener(this);
      } else {
        table.removePropertyChangeListener(this);
        table.getModel().removeTableModelListener(this);
      }
    }

    this.isDynamicAdjustment = isDynamicAdjustment;
  }

  //
  // Implement the PropertyChangeListener
  //
  @Override
  public void propertyChange(PropertyChangeEvent event) {
    // When the TableModel changes we need to update the listeners
    // and column widths

    if ("model".equals(event.getPropertyName())) {
      TableModel model = (TableModel) event.getOldValue();
      model.removeTableModelListener(this);

      model = (TableModel) event.getNewValue();
      model.addTableModelListener(this);
      adjustColumns();
    }
  }

  //
  // Implement the TableModelListener
  //
  @Override
  public void tableChanged(final TableModelEvent event) {
    if (!isColumnDataIncluded) {
      return;
    }

    // Needed when table is sorted.

    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        // A cell has been updated

        final int column = table.convertColumnIndexToView(event.getColumn());

        if (event.getType() == TableModelEvent.UPDATE && column != -1) {
          // Only need to worry about an increase in width for this cell

          if (isOnlyAdjustLarger) {
            final int row = event.getFirstRow();
            final TableColumn tableColumn = table.getColumnModel().getColumn(column);

            if (tableColumn.getResizable()) {
              final int width = getCellDataWidth(row, column);
              updateTableColumn(column, width);
            }
          } else {
            adjustColumn(column);
          }
        } else {
          adjustColumns();
        }
      }
    });
  }

  /*
   * Install Actions to give user control of certain functionality.
   */
  private void installActions() {
    installColumnAction(true, true, "adjustColumn", "control ADD");
    installColumnAction(false, true, "adjustColumns", "control shift ADD");
    installColumnAction(true, false, "restoreColumn", "control SUBTRACT");
    installColumnAction(false, false, "restoreColumns", "control shift SUBTRACT");

    installToggleAction(true, false, "toggleDynamic", "control MULTIPLY");
    installToggleAction(false, true, "toggleLarger", "control DIVIDE");
  }

  /*
   * Update the input and action maps with a new ColumnAction
   */
  private void installColumnAction(boolean isSelectedColumn, boolean isAdjust, String key,
      String keyStroke) {
    final Action action = new ColumnAction(isSelectedColumn, isAdjust);
    final KeyStroke ks = KeyStroke.getKeyStroke(keyStroke);
    table.getInputMap().put(ks, key);
    table.getActionMap().put(key, action);
  }

  /*
   * Update the input and action maps with new ToggleAction
   */
  private void installToggleAction(boolean isToggleDynamic, boolean isToggleLarger, String key,
      String keyStroke) {
    final Action action = new ToggleAction(isToggleDynamic, isToggleLarger);
    final KeyStroke ks = KeyStroke.getKeyStroke(keyStroke);
    table.getInputMap().put(ks, key);
    table.getActionMap().put(key, action);
  }

  /**
   * Action to adjust or restore the width of a single column or all columns.
   */
  private class ColumnAction extends AbstractAction {
    private static final long serialVersionUID = -6176797789918378461L;
    private final boolean isSelectedColumn;
    private final boolean isAdjust;

    public ColumnAction(boolean isSelectedColumn, boolean isAdjust) {
      this.isSelectedColumn = isSelectedColumn;
      this.isAdjust = isAdjust;
    }

    @Override
    public void actionPerformed(ActionEvent event) {
      // Handle selected column(s) width change actions

      if (isSelectedColumn) {
        final int[] columns = table.getSelectedColumns();

        for (int i = 0; i < columns.length; i++) {
          if (isAdjust) {
            adjustColumn(columns[i]);
          } else {
            restoreColumn(columns[i]);
          }
        }
      } else if (isAdjust) {
        adjustColumns();
      } else {
        restoreColumns();
      }
    }
  }

  /**
   * Toggle properties of the TableColumnAdjuster so the user can customize the functionality to
   * their preferences
   */
  private class ToggleAction extends AbstractAction {
    private static final long serialVersionUID = 8147471317715640019L;
    private final boolean isToggleDynamic;
    private final boolean isToggleLarger;

    public ToggleAction(boolean isToggleDynamic, boolean isToggleLarger) {
      this.isToggleDynamic = isToggleDynamic;
      this.isToggleLarger = isToggleLarger;
    }

    @Override
    public void actionPerformed(ActionEvent event) {
      if (isToggleDynamic) {
        setDynamicAdjustment(!isDynamicAdjustment);
        return;
      }

      if (isToggleLarger) {
        setOnlyAdjustLarger(!isOnlyAdjustLarger);
        return;
      }
    }
  }
}
