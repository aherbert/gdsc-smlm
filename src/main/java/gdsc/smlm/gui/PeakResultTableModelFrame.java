package gdsc.smlm.gui;

import java.awt.EventQueue;
import java.util.Arrays;

import javax.swing.DefaultListSelectionModel;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableColumnModel;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;

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
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultStoreList;

/**
 * A frame that shows a PeakResultsTableModel
 * 
 * @author Alex Herbert
 */
public class PeakResultTableModelFrame extends JFrame
{
	private static final long serialVersionUID = -3671174621388288975L;

	private JTable table;

	public PeakResultTableModelFrame(PeakResultTableModel model)
	{
		this(model, null, null);
	}

	public PeakResultTableModelFrame(PeakResultTableModel model, TableColumnModel cm, ListSelectionModel selectionModel)
	{
		table = new PeakResultTableModelJTable(model, cm, selectionModel);
		
		// TODO - sort this
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

		final JScrollPane scroll = new JScrollPane(table);
		add(scroll);
		pack();
	}
	
	/**
	 * Launch the application.
	 * 
	 * @throws InterruptedException
	 */
	public static void main(String[] args) throws InterruptedException
	{
		final RandomGenerator r = new Well19937c();
		final int n = 10;

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
						e.getValueIsAdjusting(), Arrays.toString(getSelectedIndices(selectionModel)));
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
					
					final PeakResultTableModel model = new PeakResultTableModel(store, cw.getCalibration(), null, tableSettings.build());

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

		//		// This doesn't work as the frames are frozen
		//		Thread.sleep(1000);
		//		
		//		// Random selections ...
		//		EventQueue.invokeLater(new Runnable()
		//		{
		//			public void run()
		//			{
		//				try
		//				{
		//					while (true)
		//					{
		//						Thread.sleep(5000);
		//						int k = r.nextInt(n);
		//						int[] indices = Random.sample(k, n, r);
		//						selectionModel.clearSelection();
		//						for (int index : indices)
		//							selectionModel.addSelectionInterval(index, index);
		//					}
		//				}
		//				catch (InterruptedException e)
		//				{
		//					e.printStackTrace();
		//				}
		//			}
		//		});
	}

	/**
	 * Gets the selected indices from the selection model.
	 * <p>
	 * Copied from javax.swing.JList
	 *
	 * @param sm
	 *            the sm
	 * @return the selected indices
	 */
	public static int[] getSelectedIndices(ListSelectionModel sm)
	{
		int iMin = sm.getMinSelectionIndex();
		int iMax = sm.getMaxSelectionIndex();

		if ((iMin < 0) || (iMax < 0))
		{
			return new int[0];
		}

		int[] rvTmp = new int[1 + (iMax - iMin)];
		int n = 0;
		for (int i = iMin; i <= iMax; i++)
		{
			if (sm.isSelectedIndex(i))
			{
				rvTmp[n++] = i;
			}
		}
		int[] rv = new int[n];
		System.arraycopy(rvTmp, 0, rv, 0, n);
		return rv;
	}
}
