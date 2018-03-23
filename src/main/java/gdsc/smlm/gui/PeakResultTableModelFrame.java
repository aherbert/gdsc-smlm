package gdsc.smlm.gui;

import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Arrays;

import javax.swing.DefaultListSelectionModel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.utils.TextUtils;
import gdsc.core.utils.XmlUtils;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.ij.settings.SettingsManager;

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

import gdsc.smlm.results.ArrayPeakResultStore;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultStoreList;
import ij.IJ;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ScreenDimensionHelper;

/**
 * A frame that shows a PeakResultsTableModel
 * 
 * @author Alex Herbert
 */
public class PeakResultTableModelFrame extends JFrame implements ActionListener
{
	private static final long serialVersionUID = -3671174621388288975L;

	private JTable table;
	private JMenuItem fileSave;
	private JMenuItem fileShowSource;
	private JMenuItem editDelete;
	private JMenuItem editClear;
	private JMenuItem editTableSettings;
	private String saveName;

	/**
	 * Instantiates a new peak result table model frame.
	 *
	 * @param model
	 *            the model
	 */
	public PeakResultTableModelFrame(PeakResultTableModel model)
	{
		this(model, null, null);
	}

	/**
	 * Instantiates a new peak result table model frame.
	 *
	 * @param model
	 *            the model
	 * @param selectionModel
	 *            the selection model
	 */
	public PeakResultTableModelFrame(PeakResultTableModel model, ListSelectionModel selectionModel)
	{
		this(model, null, selectionModel);
	}

	/**
	 * Instantiates a new peak result table model frame.
	 *
	 * @param model
	 *            the model
	 * @param columnModel
	 *            the column model
	 * @param selectionModel
	 *            the selection model
	 */
	public PeakResultTableModelFrame(final PeakResultTableModel model, TableColumnModel columnModel,
			ListSelectionModel selectionModel)
	{
		setJMenuBar(createMenuBar());

		// This is required to get the column sizes for the model data.
		model.setLive(true);

		table = new PeakResultTableModelJTable(model, columnModel, selectionModel);

		final JScrollPane scroll = new JScrollPane(table);

		ScreenDimensionHelper helper = new ScreenDimensionHelper();
		helper.setMinHeight(300);
		helper.setup(scroll);

		add(scroll);
		pack();

		// If the window is never set visible then we cannot ensure this happens in 
		// the closing event so do it here and again when opened.
		model.setLive(false);

		addWindowListener(new WindowAdapter()
		{
			public void windowOpened(WindowEvent e)
			{
				model.setLive(true);
				WindowManager.addWindow(PeakResultTableModelFrame.this);
			}

			public void windowClosed(WindowEvent e)
			{
				model.setLive(false);
				WindowManager.removeWindow(PeakResultTableModelFrame.this);
			}
		});
	}

	/**
	 * Clean up this table. This should only be called when the table is no longer required as it removes the JTable
	 * from the model listeners.
	 */
	public void cleanUp()
	{
		// Since the models may be shared
		table.getModel().removeTableModelListener(table);
		table.getColumnModel().removeColumnModelListener(table);
		table.getSelectionModel().removeListSelectionListener(table);
	}

	private JMenuBar createMenuBar()
	{
		JMenuBar menubar = new JMenuBar();
		menubar.add(createFileMenu());
		menubar.add(createEditMenu());
		return menubar;
	}

	private JMenu createFileMenu()
	{
		final JMenu menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		menu.add(fileSave = add(menu, "Save ...", KeyEvent.VK_S, "ctrl pressed S"));
		menu.add(fileShowSource = add(menu, "Show source", KeyEvent.VK_W, null));
		return menu;
	}

	private JMenu createEditMenu()
	{
		final JMenu menu = new JMenu("Edit");
		menu.setMnemonic(KeyEvent.VK_E);
		menu.add(editDelete = add(menu, "Delete", KeyEvent.VK_D, null));
		menu.add(editClear = add(menu, "Clear", KeyEvent.VK_C, null));
		menu.add(editTableSettings = add(menu, "Table Settings ...", KeyEvent.VK_T, "ctrl pressed T"));
		return menu;
	}

	private JMenuItem add(JMenu menu, String text, int mnemonic, String keyStroke)
	{
		JMenuItem item = new JMenuItem(text, mnemonic);
		if (keyStroke != null)
			item.setAccelerator(KeyStroke.getKeyStroke(keyStroke));
		item.addActionListener(this);
		return item;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e)
	{
		final Object src = e.getSource();
		if (src == fileSave)
			doSave();
		else if (src == fileShowSource)
			doShowSource();
		else if (src == editDelete)
			doDelete();
		else if (src == editClear)
			doClear();
		else if (src == editTableSettings)
			doEditTableSettings();
	}

	private void doSave()
	{
		PeakResultTableModel model = getModel();
		if (model == null || model.getRowCount() == 0)
			return;
		ExtendedGenericDialog gd = new ExtendedGenericDialog("Save Results", this);
		if (TextUtils.isNullOrEmpty(saveName))
			saveName = getTitle();
		gd.addStringField("Results_set_name", saveName, 30);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		saveName = gd.getNextString();
		if (TextUtils.isNullOrEmpty(saveName))
		{
			IJ.error("No results set name");
			return;
		}
		MemoryPeakResults results = model.toMemoryPeakResults();
		results.setName(saveName);
		MemoryPeakResults.addResults(results);
	}

	private void doShowSource()
	{
		PeakResultTableModel model = getModel();
		if (model == null)
			return;
		ImageSource source = model.getSource();
		String text = getTitle() + " source";
		if (source == null)
		{
			text += " = NA";
		}
		else
		{
			text += "\n" + XmlUtils.prettyPrintXml(source.toXML());
		}
		IJ.log(text);
	}

	private void doDelete()
	{
		PeakResultTableModel model = getModel();
		if (model == null)
			return;
		int[] indices = table.getSelectedRows();
		model.remove(this, indices);
	}

	private void doClear()
	{
		PeakResultTableModel model = getModel();
		if (model == null)
			return;
		model.clear();
	}

	private void doEditTableSettings()
	{
		PeakResultTableModel model = getModel();
		if (model == null)
			return;
		ResultsTableSettings.Builder tableSettings = model.getTableSettings().toBuilder();
		// Copied from ResultsManager.addTableResultsOptions
		ExtendedGenericDialog egd = new ExtendedGenericDialog("Table Settings", this);
		egd.addChoice("Table_distance_unit", SettingsManager.getDistanceUnitNames(),
				tableSettings.getDistanceUnit().getNumber());
		egd.addChoice("Table_intensity_unit", SettingsManager.getIntensityUnitNames(),
				tableSettings.getIntensityUnit().getNumber());
		egd.addChoice("Table_angle_unit", SettingsManager.getAngleUnitNames(),
				tableSettings.getAngleUnit().getNumber());
		egd.addCheckbox("Table_show_fitting_data", tableSettings.getShowFittingData());
		egd.addCheckbox("Table_show_noise_data", tableSettings.getShowNoiseData());
		egd.addCheckbox("Table_show_precision", tableSettings.getShowPrecision());
		egd.addSlider("Table_precision", 0, 10, tableSettings.getRoundingPrecision());
		egd.showDialog();
		if (egd.wasCanceled())
			return;
		tableSettings.setDistanceUnitValue(egd.getNextChoiceIndex());
		tableSettings.setIntensityUnitValue(egd.getNextChoiceIndex());
		tableSettings.setAngleUnitValue(egd.getNextChoiceIndex());
		tableSettings.setShowFittingData(egd.getNextBoolean());
		tableSettings.setShowNoiseData(egd.getNextBoolean());
		tableSettings.setShowPrecision(egd.getNextBoolean());
		tableSettings.setRoundingPrecision((int) egd.getNextNumber());
		model.setTableSettings(tableSettings.build());
	}

	private PeakResultTableModel getModel()
	{
		TableModel model = table.getModel();
		return (model instanceof PeakResultTableModel) ? (PeakResultTableModel) model : null;
	}

	/**
	 * Launch the application.
	 * 
	 * @throws InterruptedException
	 */
	public static void main(String[] args) throws InterruptedException
	{
		final RandomGenerator r = new Well19937c();
		final int n = 20;

		final ListSelectionModel selectionModel = new DefaultListSelectionModel();

		// The way to interact with a model from many parts of the same GUI is through
		// the same selection model. Each JList updates using the same selection.
		selectionModel.addListSelectionListener(new ListSelectionListener()
		{
			public void valueChanged(ListSelectionEvent e)
			{
				if (e.getValueIsAdjusting())
					return;
				System.out.printf("Model Selected %d-%d [%b] : %s\n", e.getFirstIndex(), e.getLastIndex(),
						e.getValueIsAdjusting(),
						Arrays.toString(ListSelectionModelHelper.getSelectedIndices(selectionModel)));
			}
		});

		EventQueue.invokeLater(new Runnable()
		{
			public void run()
			{
				try
				{
					final PeakResultStoreList store = new ArrayPeakResultStore(10);
					for (int i = n; i-- > 0;)
					{
						store.add(new PeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(), r.nextDouble(),
								r.nextFloat(), PeakResult.createParams(r.nextFloat(), r.nextFloat(), r.nextFloat(),
										r.nextFloat(), r.nextFloat()),
								null));
					}

					CalibrationWriter cw = new CalibrationWriter();
					cw.setNmPerPixel(100);
					cw.setCountPerPhoton(10);
					cw.setDistanceUnit(DistanceUnit.PIXEL);
					cw.setIntensityUnit(IntensityUnit.COUNT);

					ResultsTableSettings.Builder tableSettings = ResultsTableSettings.newBuilder();
					tableSettings.setDistanceUnit(DistanceUnit.NM);
					tableSettings.setIntensityUnit(IntensityUnit.PHOTON);
					tableSettings.setShowFittingData(true);
					tableSettings.setShowNoiseData(true);
					tableSettings.setShowPrecision(true);
					tableSettings.setRoundingPrecision(4);

					final PeakResultTableModel model = new PeakResultTableModel(store, cw.getCalibration(), null,
							tableSettings.build());

					final PeakResultTableModelFrame d = new PeakResultTableModelFrame(model, null, selectionModel);
					d.setTitle("D");
					//					d.addListSelectionListener(new ListSelectionListener()
					//					{
					//						public void valueChanged(ListSelectionEvent e)
					//						{
					//							// Only process the event if the value is not adjusting.
					//							// Then to determine what has changed only process the 
					//							// indices between the first and last index. 
					//
					//							if (e.getValueIsAdjusting())
					//								return;
					//							System.out.printf("D Selected %d-%d [%b] : %s\n", e.getFirstIndex(), e.getLastIndex(),
					//									e.getValueIsAdjusting(), Arrays.toString(d.list.getSelectedIndices()));
					//						}
					//					});
					d.setDefaultCloseOperation(EXIT_ON_CLOSE);
					d.setVisible(true);

					// Selecting in one list activates the other list

					final PeakResultTableModelFrame d2 = new PeakResultTableModelFrame(model, null, selectionModel);
					d2.setTitle("D2");
					//					d2.addListSelectionListener(new ListSelectionListener()
					//					{
					//						public void valueChanged(ListSelectionEvent e)
					//						{
					//							if (e.getValueIsAdjusting())
					//								return;
					//							int[] indices = d2.list.getSelectedIndices();
					//							System.out.printf("D2 Selected %d-%d [%b] : %s\n", e.getFirstIndex(), e.getLastIndex(),
					//									e.getValueIsAdjusting(), Arrays.toString(indices));
					//							//d.list.setSelectedIndices(indices);
					//						}
					//					});
					d2.setDefaultCloseOperation(EXIT_ON_CLOSE);
					d2.setVisible(true);
				}
				catch (Exception e)
				{
					e.printStackTrace();
				}
			}
		});
	}
}
