package gdsc.smlm.gui;

import java.awt.Component;
import java.awt.EventQueue;
import java.util.Arrays;

import javax.swing.DefaultListSelectionModel;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JScrollPane;
import javax.swing.ListCellRenderer;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
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
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultStoreList;

/**
 * A frame that shows a PeakResultsListModel
 * 
 * @author Alex Herbert
 */
public class PeakResultListModelFrame extends JFrame
{
	private static final long serialVersionUID = -1530205032042929260L;

	private class MyCellRenderer extends JLabel implements ListCellRenderer<PeakResult>
	{
		// This is the only method defined by ListCellRenderer.
		// We just reconfigure the JLabel each time we're called.
		private static final long serialVersionUID = 1998620838894273028L;

		public Component getListCellRendererComponent(JList<? extends PeakResult> list, // the list
				PeakResult value, // value to display
				int index, // cell index
				boolean isSelected, // is the cell selected
				boolean cellHasFocus) // does the cell have focus
		{
			// TODO - Make this a better representation of the Peak Result.
			// Build a configurable layout using the TableResults settings.
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < PeakResult.STANDARD_PARAMETERS; i++)
			{
				if (sb.length() != 0)
					sb.append(' ');
				sb.append(PeakResult.getParameterName(i)).append('=').append(value.getParameter(i));
			}

			String s = sb.toString();
			setText(s);
			if (isSelected)
			{
				setBackground(list.getSelectionBackground());
				setForeground(list.getSelectionForeground());
			}
			else
			{
				setBackground(list.getBackground());
				setForeground(list.getForeground());
			}
			setEnabled(list.isEnabled());
			setFont(list.getFont());
			setOpaque(true);
			return this;
		}
	}

	private JList<PeakResult> list;

	public PeakResultListModelFrame(PeakResultListModel model)
	{
		this(model, null);
	}

	public PeakResultListModelFrame(PeakResultListModel model, ListSelectionModel selectionModel)
	{
		list = new JList<PeakResult>(model);
		list.setPrototypeCellValue(new ExtendedPeakResult(1, 1, 1, 1));
		list.setCellRenderer(new MyCellRenderer());
		if (selectionModel != null)
			list.setSelectionModel(selectionModel);
		final JScrollPane scroll = new JScrollPane(list);
		add(scroll);
		pack();

	}

	public void addListSelectionListener(ListSelectionListener listener)
	{
		list.addListSelectionListener(listener);
	}

	public void removeListSelectionListener(ListSelectionListener listener)
	{
		list.removeListSelectionListener(listener);
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
					final PeakResultListModel model = new PeakResultListModel(store);

					final PeakResultListModelFrame d = new PeakResultListModelFrame(model, selectionModel);
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

					final PeakResultListModelFrame d2 = new PeakResultListModelFrame(model, selectionModel);
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
