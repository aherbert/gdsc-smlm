package gdsc.smlm.gui;

import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.event.TableModelEvent;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

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
		tca.setMaxRows(10);
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

		if (e == null || e.getFirstRow() == TableModelEvent.HEADER_ROW)
		{
			// The whole thing changed so resize the columns
			SwingUtilities.invokeLater(new Runnable()
			{
				public void run()
				{
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
}
