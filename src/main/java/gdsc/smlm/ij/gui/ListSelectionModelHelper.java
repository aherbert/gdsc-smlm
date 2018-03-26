package gdsc.smlm.ij.gui;

import javax.swing.ListSelectionModel;

/**
 * A helper class for the ListSelectionModel
 */
public class ListSelectionModelHelper
{
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

	public static void setSelectedIndices(ListSelectionModel sm, int[] indices)
	{
		if (indices == null || indices.length == 0)
			return;

		sm.setValueIsAdjusting(true);
		sm.setSelectionInterval(indices[0], indices[0]);
		for (int i = 1; i < indices.length; i++)
			sm.addSelectionInterval(indices[i], indices[i]);
		sm.setValueIsAdjusting(false);
	}
}
