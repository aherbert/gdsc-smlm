package gdsc.smlm.gui;

import java.util.Arrays;

import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.TableModelEvent;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.Converter;
import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtosHelper;

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

import gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.CalibrationProtosHelper;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.ResultsProtosHelper;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.results.ArrayPeakResultStore;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultConversionHelper;
import gdsc.smlm.results.PeakResultData;
import gdsc.smlm.results.PeakResultStoreList;
import gdsc.smlm.results.data.*;

/**
 * Stores peak results and allows event propagation to listeners of the model.
 */
public class PeakResultTableModelJTable extends JTable
{
	private static final long serialVersionUID = 7144289957208169053L;

	public PeakResultTableModelJTable(PeakResultTableModel model, TableColumnModel cm,
			ListSelectionModel selectionModel)
	{
		super(model, cm, selectionModel);
		updateRenderer();
	}

	@Override
	public void tableChanged(TableModelEvent e)
	{
		if (e.getType() == PeakResultTableModel.RENDERER)
		{
			updateRenderer();
			return;
		}
		if (e.getType() == TableModelEvent.ALL_COLUMNS)
		{
			updateRenderer();
		}
		super.tableChanged(e);
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
			// Reset?
		}
	}
}
