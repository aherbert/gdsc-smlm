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
package gdsc.smlm.ij.gui;

import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.RowSorter;
import javax.swing.SwingUtilities;
import javax.swing.event.TableModelEvent;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

import gdsc.smlm.results.PeakResult;

/**
 * Stores peak results and allows event propagation to listeners of the model.
 */
public class PeakResultTableModelJTable extends JTable
{
	private static final long serialVersionUID = 7144289957208169053L;

	private TableColumnAdjuster tca;

	public PeakResultTableModelJTable(PeakResultTableModel model, TableColumnModel cm,
			ListSelectionModel selectionModel)
	{
		super(model, cm, selectionModel);
		updateRenderer();

		// Make all the columns show the full data.
		setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

		// Note that this is not dynamic and so must manually be called when columns change
		tca = new TableColumnAdjuster(this, 6, false);
		//  Only process 10 rows max.
		tca.setMaxRows(5);
		tca.adjustColumns();

		setAutoCreateRowSorter(true);
	}

	@Override
	public void tableChanged(final TableModelEvent e)
	{
		if (e.getType() == PeakResultTableModel.RENDERER)
		{
			// Special event when the rendering has changed,
			// e.g. the rounding precision has changed
			updateRenderer();
			return;
		}

		super.tableChanged(e);

		if (e.getFirstRow() == TableModelEvent.HEADER_ROW)
		{
			// The whole thing changed so resize the columns
			SwingUtilities.invokeLater(new Runnable()
			{
				@Override
				public void run()
				{
					// This is null when the table is first created
					if (tca != null)
						tca.adjustColumns();
				}
			});
		}
	}

	private void updateRenderer()
	{
		// For rounding
		TableModel m = getModel();
		if (m instanceof PeakResultTableModel)
		{
			PeakResultTableModel model = (PeakResultTableModel) m;
			setDefaultRenderer(Float.class, model.getFloatRenderer());
			setDefaultRenderer(Double.class, model.getDoubleRenderer());
			setDefaultRenderer(Integer.class, model.getIntegerRenderer());
		}
		else
		{
			// Reset to null?
		}
	}

	/**
	 * Returns the data of all selected rows. This maps the indices from the view to the data model.
	 *
	 * @return an array containing the data of all selected rows,
	 *         or an empty array if no row is selected
	 * @see #getSelectedRow
	 */
	public PeakResult[] getSelectedData()
	{
		if (dataModel instanceof PeakResultTableModel)
		{
			PeakResultTableModel model = (PeakResultTableModel) dataModel;
			int iMin = selectionModel.getMinSelectionIndex();
			int iMax = selectionModel.getMaxSelectionIndex();

			if ((iMin == -1) || (iMax == -1))
			{
				return new PeakResult[0];
			}

			PeakResult[] rvTmp = new PeakResult[1 + (iMax - iMin)];
			int n = 0;

			RowSorter<?> sorter = getRowSorter();
			if (sorter != null)
			{
				for (int i = iMin; i <= iMax; i++)
				{
					if (selectionModel.isSelectedIndex(i))
					{
						rvTmp[n++] = model.get(sorter.convertRowIndexToModel(i));
					}
				}
			}
			else
			{
				for (int i = iMin; i <= iMax; i++)
				{
					if (selectionModel.isSelectedIndex(i))
					{
						rvTmp[n++] = model.get(i);
					}
				}
			}
			if (n == rvTmp.length)
				return rvTmp;
			PeakResult[] rv = new PeakResult[n];
			System.arraycopy(rvTmp, 0, rv, 0, n);
			return rv;
		}
		return new PeakResult[0];
	}

	/**
	 * Returns the location of <code>index</code> in terms of the
	 * underlying model. That is, for the row <code>index</code> in
	 * the coordinates of the view this returns the row index in terms
	 * of the underlying model.
	 *
	 * @param indices
	 *            the indices (updated in-place)
	 * @throws IndexOutOfBoundsException
	 *             if <code>index</code> is outside the
	 *             range of the view
	 */
	public void convertRowIndexToModel(int[] indices)
	{
		RowSorter<?> sorter = getRowSorter();
		if (sorter != null)
		{
			for (int i = 0; i < indices.length; i++)
				indices[i] = sorter.convertRowIndexToModel(indices[i]);
		}
	}

	/**
	 * Returns the location of <code>index</code> in terms of the
	 * view. That is, for the row <code>index</code> in the
	 * coordinates of the underlying model this returns the row index
	 * in terms of the view.
	 *
	 * @param indices
	 *            the indices (updated in-place)
	 * @throws IndexOutOfBoundsException
	 *             if <code>index</code> is outside
	 *             the range of the model
	 */
	public void convertRowIndexToView(int[] indices)
	{
		RowSorter<?> sorter = getRowSorter();
		if (sorter != null)
		{
			for (int i = 0; i < indices.length; i++)
				indices[i] = sorter.convertRowIndexToView(indices[i]);
		}
	}
}
