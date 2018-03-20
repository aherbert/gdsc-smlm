package gdsc.smlm.gui;

import java.awt.Component;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JScrollPane;
import javax.swing.ListCellRenderer;
import javax.swing.event.ListSelectionListener;

import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.PeakResult;

public class PeakResultListDialog extends JDialog
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
				sb.append(PeakResult.getParameterName(i)).append('=').append(value.getParameter(i));

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

	public PeakResultListDialog(PeakResultModel model)
	{
		list = new JList<PeakResult>(model);
		list.setPrototypeCellValue(new ExtendedPeakResult(1, 1, 1, 1));
		list.setCellRenderer(new MyCellRenderer());
		final JScrollPane scroll = new JScrollPane(list);
		add(scroll);

	}

	public void addListSelectionListener(ListSelectionListener listener)
	{
		list.addListSelectionListener(listener);
	}

	public void removeListSelectionListener(ListSelectionListener listener)
	{
		list.removeListSelectionListener(listener);
	}
}
