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

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingConstants;
import javax.swing.event.TableModelEvent;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.Converter;
import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import gdsc.smlm.data.config.ResultsProtosHelper;
import gdsc.smlm.results.ArrayPeakResultStore;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.ImageSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultConversionHelper;
import gdsc.smlm.results.PeakResultData;
import gdsc.smlm.results.PeakResultStoreList;
import gdsc.smlm.results.data.PeakResultDataEndFrame;
import gdsc.smlm.results.data.PeakResultDataError;
import gdsc.smlm.results.data.PeakResultDataFloat;
import gdsc.smlm.results.data.PeakResultDataFrame;
import gdsc.smlm.results.data.PeakResultDataId;
import gdsc.smlm.results.data.PeakResultDataOrigValue;
import gdsc.smlm.results.data.PeakResultDataOrigX;
import gdsc.smlm.results.data.PeakResultDataOrigY;
import gdsc.smlm.results.data.PeakResultDataParameterConverter;
import gdsc.smlm.results.data.PeakResultDataParameterDeviationConverter;
import gdsc.smlm.results.data.PeakResultDataPrecision;
import gdsc.smlm.results.data.PeakResultDataSNR;
import gnu.trove.list.array.TIntArrayList;

/**
 * Stores peak results and allows event propagation to listeners of the model.
 */
public class PeakResultTableModel extends AbstractTableModel
{
	private static final long serialVersionUID = 3600737169000976938L;

	/** Identifies the cell renderer has changed. */
	static final int RENDERER = -99;

	/** Count used to count how many PeakResultTableModelFrame are displaying this model. */
	private AtomicInteger liveCount = new AtomicInteger();

	/** The flag used to indicate the columns require updating. */
	private AtomicBoolean columnsComputed = new AtomicBoolean(false);

	private final PeakResultStoreList data;
	private final Calibration calibration;
	private final PSF psf;
	private ResultsTableSettings tableSettings;
	private boolean checkForDuplicates = false;

	// These depend on the source results
	private ImageSource source;
	private String configuration = "";
	private boolean showDeviations = false;
	private boolean showEndFrame = false;
	private boolean showId = false;
	private boolean showZ = false;

	// Used for the columns 
	private Rounder rounder;
	private PeakResultData<?>[] values;
	private String[] names;
	private boolean rowCounter;

	/**
	 * Instantiates a new peak result model using settings from the resultsSource.
	 *
	 * @param resultsSource
	 *            the results source
	 * @param copyData
	 *            Set to true to copy the data from the results source
	 * @param tableSettings
	 *            the table settings
	 */
	public PeakResultTableModel(MemoryPeakResults resultsSource, boolean copyData, ResultsTableSettings tableSettings)
	{
		if (tableSettings == null)
			tableSettings = ResultsProtosHelper.defaultResultsSettings.getResultsTableSettings();
		this.tableSettings = tableSettings;

		data = (copyData) ? new ArrayPeakResultStore(resultsSource.toArray()) : new ArrayPeakResultStore(10);
		this.calibration = resultsSource.getCalibration();
		this.psf = resultsSource.getPSF();

		setShowDeviations(resultsSource.hasDeviations());
		setShowZ(resultsSource.is3D());
		setShowId(resultsSource.hasId());
		setShowEndFrame(resultsSource.hasEndFrame());
		setSource(resultsSource.getSource());
		setConfiguration(resultsSource.getConfiguration());
	}

	/**
	 * Instantiates a new peak result model using the store.
	 *
	 * @param results
	 *            the results
	 * @param calibration
	 *            the calibration
	 * @param psf
	 *            the psf
	 * @param tableSettings
	 *            the table settings
	 */
	public PeakResultTableModel(PeakResultStoreList results, Calibration calibration, PSF psf,
			ResultsTableSettings tableSettings)
	{
		if (results == null)
			results = new ArrayPeakResultStore(10);
		if (calibration == null)
			calibration = Calibration.getDefaultInstance();
		if (psf == null)
			psf = PSF.getDefaultInstance();
		if (tableSettings == null)
			tableSettings = ResultsProtosHelper.defaultResultsSettings.getResultsTableSettings();
		this.data = results;
		this.calibration = calibration;
		this.psf = psf;
		this.tableSettings = tableSettings;
	}

	//*************************************************************************/
	// Table column appearance
	//*************************************************************************/

	/**
	 * Gets the table settings.
	 *
	 * @return the table settings
	 */
	public ResultsTableSettings getTableSettings()
	{
		return tableSettings;
	}

	/**
	 * Sets the table settings.
	 *
	 * @param tableSettings
	 *            the new table settings
	 */
	public void setTableSettings(ResultsTableSettings tableSettings)
	{
		if (tableSettings == null)
			tableSettings = ResultsProtosHelper.defaultResultsSettings.getResultsTableSettings();
		if (equals(this.tableSettings, tableSettings))
			return;
		this.tableSettings = tableSettings;
		createTableStructure(true);
	}

	/**
	 * Check if the settings are equal
	 *
	 * @param t1
	 *            the t1
	 * @param t2
	 *            the t2
	 * @return true, if successful
	 */
	private boolean equals(ResultsTableSettings t1, ResultsTableSettings t2)
	{
		// Adapted from ResultsTableSettings.equals() to only use the settings of interest
		boolean result = true;
		result = result && t1.getDistanceUnitValue() == t2.getDistanceUnitValue();
		result = result && t1.getIntensityUnitValue() == t2.getIntensityUnitValue();
		result = result && t1.getAngleUnitValue() == t2.getAngleUnitValue();
		result = result && (t1.getShowPrecision() == t2.getShowPrecision());
		result = result && (t1.getShowFittingData() == t2.getShowFittingData());
		result = result && (t1.getShowNoiseData() == t2.getShowNoiseData());
		result = result && (t1.getRoundingPrecision() == t2.getRoundingPrecision());
		result = result && (t1.getShowRowCounter() == t2.getShowRowCounter());
		return result;
	}

	/**
	 * Called when the structure of the table should be created. If the structure has not changed or no live table are
	 * attached then this does nothing.
	 *
	 * @param changed
	 *            Set to true if a property controlling the structure has changed
	 */
	private void createTableStructure(boolean changed)
	{
		if (changed)
			columnsComputed.set(false);

		// If no live attached tables or nothing has changed then return
		// Note: We use a live count instead of the listener count (see listenerList.getListenerCount())
		// so it can be turned off by PeakResultTableModelFrame.
		if (liveCount.get() == 0 || columnsComputed.get())
			return;

		columnsComputed.set(true);

		rounder = RounderFactory.create(tableSettings.getRoundingPrecision());

		// Create the converters
		PeakResultConversionHelper helper = new PeakResultConversionHelper(calibration, psf);
		helper.setIntensityUnit(tableSettings.getIntensityUnit());
		helper.setDistanceUnit(tableSettings.getDistanceUnit());
		helper.setAngleUnit(tableSettings.getAngleUnit());
		final Converter[] converters = helper.getConverters();
		final Converter ic = converters[PeakResult.INTENSITY];

		String[] paramNames = helper.getNames();
		String[] unitNames = helper.getUnitNames();

		//Calibration tableCalibration = (helper.isCalibrationChanged()) ? helper.getCalibration() : calibration;

		// Organise the data columns.
		// This is done as per the IJTablePeakResults for consistency
		TurboList<PeakResultData<?>> valuesList = new TurboList<PeakResultData<?>>();
		TurboList<String> namesList = new TurboList<String>();

		if (rowCounter = tableSettings.getShowRowCounter())
		{
			valuesList.add(new PeakResultDataFrame());
			namesList.add("#");
		}

		valuesList.add(new PeakResultDataFrame());
		addName(valuesList, namesList);
		if (showEndFrame)
		{
			valuesList.add(new PeakResultDataEndFrame());
			addName(valuesList, namesList);
		}
		if (showId)
		{
			valuesList.add(new PeakResultDataId());
			addName(valuesList, namesList);
		}
		if (tableSettings.getShowFittingData())
		{
			valuesList.add(new PeakResultDataOrigX());
			addName(valuesList, namesList);
			valuesList.add(new PeakResultDataOrigY());
			addName(valuesList, namesList);
			valuesList.add(new PeakResultDataOrigValue());
			addName(valuesList, namesList);
			valuesList.add(new PeakResultDataError());
			addName(valuesList, namesList);
		}
		if (tableSettings.getShowNoiseData())
		{
			// Must be converted
			valuesList.add(new PeakResultDataFloat()
			{
				@Override
				public Float getValue(PeakResult result)
				{
					return ic.convert(result.getNoise());
				}
			});
			addName("Noise", namesList, unitNames[PeakResult.INTENSITY]);
			valuesList.add(new PeakResultDataFloat()
			{
				@Override
				public Float getValue(PeakResult result)
				{
					return ic.convert(result.getMeanIntensity());
				}
			});
			addName("Mean" + paramNames[PeakResult.INTENSITY], namesList, unitNames[PeakResult.INTENSITY]);
			valuesList.add(new PeakResultDataSNR());
			addName(valuesList, namesList);
		}

		int[] outIndices = SimpleArrayUtils.newArray(converters.length, 0, 1);
		if (!showZ)
		{
			TIntArrayList list = new TIntArrayList(outIndices);
			list.remove(PeakResult.Z);
			outIndices = list.toArray();
		}

		for (int i : outIndices)
		{
			// Must be converted
			valuesList.add(new PeakResultDataParameterConverter(converters[i], i));
			addName(paramNames[i], namesList, unitNames[i]);
			if (showDeviations)
			{
				valuesList.add(new PeakResultDataParameterDeviationConverter(converters[i], i));
				namesList.add("+/-");
			}
		}

		if (tableSettings.getShowPrecision())
		{
			PeakResultDataPrecision p = null;

			try
			{
				final Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(getPSF(),
						calibration, Gaussian2DPeakResultHelper.LSE_PRECISION);
				p = new PeakResultDataPrecision()
				{
					@Override
					public Double getValue(PeakResult result)
					{
						if (result.hasPrecision())
							return result.getPrecision();
						if (calculator != null)
							return calculator.getLSEPrecision(result.getParameters(), result.getNoise());
						return 0.0;
					}
				};
			}
			catch (ConfigurationException e)
			{
			}
			catch (ConversionException e)
			{
			}
			if (p == null)
				p = new PeakResultDataPrecision();
			valuesList.add(p);
			namesList.add("Precision (nm)");
		}

		values = valuesList.toArray(new PeakResultData<?>[valuesList.size()]);
		names = namesList.toArray(new String[namesList.size()]);

		fireTableStructureChanged();
	}

	private static void addName(TurboList<PeakResultData<?>> valuesList, TurboList<String> namesList)
	{
		namesList.add(valuesList.get(valuesList.size() - 1).getValueName());
	}

	private static void addName(String name, TurboList<String> namesList, String unitName)
	{
		if (!TextUtils.isNullOrEmpty(unitName))
			name += " (" + unitName + ")";
		namesList.add(name);
	}

	//*************************************************************************/
	// Table model methods
	//*************************************************************************/

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.table.TableModel#getRowCount()
	 */
	@Override
	public int getRowCount()
	{
		return data.size();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.table.AbstractTableModel#getColumnName(int)
	 */
	@Override
	public String getColumnName(int column)
	{
		return names[column]; // values[column].getValueName();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.table.AbstractTableModel#getColumnClass(int)
	 */
	@Override
	public Class<?> getColumnClass(int columnIndex)
	{
		return values[columnIndex].getValueClass();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.table.TableModel#getColumnCount()
	 */
	@Override
	public int getColumnCount()
	{
		return values.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.table.TableModel#getValueAt(int, int)
	 */
	@Override
	public Object getValueAt(int rowIndex, int columnIndex)
	{
		if (rowCounter)
		{
			if (columnIndex == 0)
				return new Integer(rowIndex + 1);
		}
		PeakResult r = get(rowIndex);
		return values[columnIndex].getValue(r);
	}

	//*************************************************************************/
	// Data management
	//*************************************************************************/

	/**
	 * To memory peak results.
	 *
	 * @return the memory peak results
	 */
	public MemoryPeakResults toMemoryPeakResults()
	{
		ArrayPeakResultStore store = new ArrayPeakResultStore(data.size());
		store.addArray(data.toArray());
		MemoryPeakResults results = new MemoryPeakResults(store);
		results.setPSF(psf);
		results.setCalibration(calibration);
		results.setSource(source);
		results.setConfiguration(configuration);
		return results;
	}

	/**
	 * Convert the model to an array.
	 *
	 * @return the peak result array
	 */
	public PeakResult[] toArray()
	{
		return data.toArray();
	}

	/**
	 * Gets the results for the given index.
	 *
	 * @param index
	 *            the index
	 * @return the peak result
	 */
	public PeakResult get(int index)
	{
		return data.get(index);
	}

	/**
	 * Returns the index of the first occurrence of the specified result
	 * in this store, or -1 if this list does not contain the element.
	 * More formally, returns the lowest index <tt>i</tt> such that
	 * <tt>(result==null&nbsp;?&nbsp;get(i)==null&nbsp;:&nbsp;result.equals(get(i)))</tt>,
	 * or -1 if there is no such index.
	 *
	 * @param result
	 *            the result
	 * @return the index (or -1)
	 */
	public int indexOf(PeakResult result)
	{
		return data.indexOf(result);
	}

	/**
	 * Adds the results. Duplicates can be avoided using the check for duplicates property.
	 * 
	 * @param source
	 *            the source
	 * @param peakResults
	 *            the peak results
	 * @see {@link #isCheckDuplicates()}
	 */
	public void add(Object source, PeakResult... peakResults)
	{
		if (peakResults.length == 0)
			return;
		int index0 = data.size();
		if (checkForDuplicates)
		{
			int size = 0;
			for (int i = 0; i < peakResults.length; i++)
			{
				if (!data.contains(peakResults[i]))
					peakResults[size++] = peakResults[i];
			}
			if (size == 0)
				return;
			if (size != peakResults.length)
				peakResults = Arrays.copyOf(peakResults, size);
		}
		data.addArray(peakResults);
		int index1 = data.size() - 1;

		fireTableRowsInserted(index0, index1);
	}

	/**
	 * Removes the result
	 *
	 * @param source
	 *            the source
	 * @param peakResult
	 *            the peak result
	 */
	public void remove(Object source, PeakResult peakResult)
	{
		remove(source, data.indexOf(peakResult));
	}

	/**
	 * Removes the result.
	 *
	 * @param source
	 *            the source
	 * @param index
	 *            the index
	 */
	public void remove(Object source, int index)
	{
		if (index < 0 || index >= data.size())
			return;
		data.remove(index);
		fireTableRowsDeleted(index, index);
	}

	/**
	 * Removes the result.
	 *
	 * @param source
	 *            the source
	 * @param indices
	 *            the indices
	 */
	public void remove(Object source, int[] indices)
	{
		if (indices == null || indices.length == 0)
			return;

		if (indices.length == 1)
		{
			remove(source, indices[0]);
			return;
		}

		int size = 0;
		for (int i = 0; i < indices.length; i++)
		{
			int index = indices[i];
			if (index < 0 || index >= data.size())
				continue;
			indices[size++] = index;
		}

		if (size == 0)
			return;
		if (size < indices.length)
			indices = Arrays.copyOf(indices, size);

		int[] pairs = SimpleArrayUtils.getRanges(indices);

		size = pairs.length;
		int firstRow = pairs[0];
		int lastRow = pairs[size - 1];

		// Remove ranges starting at the end (to preserve the list order)
		for (int i = size - 1; i > 0; i -= 2)
		{
			data.remove(pairs[i - 1], pairs[i]);
		}

		fireTableRowsDeleted(firstRow, lastRow);
	}

	/**
	 * Removes the results.
	 *
	 * @param source
	 *            the source
	 * @param peakResults
	 *            the peak results
	 */
	public void remove(Object source, PeakResult... peakResults)
	{
		if (peakResults.length == 0)
			return;
		if (peakResults.length == 1)
		{
			remove(source, peakResults[0]);
			return;
		}
		int[] indices = new int[peakResults.length];
		int size = 0;
		for (int i = 0; i < peakResults.length; i++)
		{
			int j = data.indexOf(peakResults[i]);
			if (j >= 0)
				indices[size++] = j;
		}
		if (size == 0)
			return;
		if (size < peakResults.length)
			indices = Arrays.copyOf(indices, size);

		int[] pairs = SimpleArrayUtils.getRanges(indices);

		size = pairs.length;
		int firstRow = pairs[0];
		int lastRow = pairs[size - 1];

		// Remove ranges starting at the end (to preserve the list order)
		for (int i = size - 1; i > 0; i -= 2)
		{
			data.remove(pairs[i - 1], pairs[i]);
		}

		fireTableRowsDeleted(firstRow, lastRow);
	}

	/**
	 * Clear the results.
	 *
	 * @param source
	 *            the source
	 */
	public void clear(Object source)
	{
		int index1 = data.size() - 1;
		if (index1 >= 0)
		{
			data.clear();
			fireTableRowsDeleted(0, index1);
		}
	}

	/**
	 * Deletes the components at the specified range of indexes.
	 * The removal is inclusive, so specifying a range of (1,5)
	 * removes the component at index 1 and the component at index 5,
	 * as well as all components in between.
	 * <p>
	 * Throws an <code>ArrayIndexOutOfBoundsException</code>
	 * if the index was invalid.
	 *
	 * @param source
	 *            the source
	 * @param fromIndex
	 *            the index of the lower end of the range
	 * @param toIndex
	 *            the index of the upper end of the range
	 */
	public void removeRange(Object source, int fromIndex, int toIndex)
	{
		data.remove(fromIndex, toIndex);
		fireTableRowsDeleted(fromIndex, toIndex);
	}

	//*************************************************************************/
	// Properties
	//*************************************************************************/

	/**
	 * Sets the model to the live state. This creates all the table layout information and causes it to update when
	 * properties are changed.
	 */
	void setLive(boolean isLive)
	{
		if (isLive)
			liveCount.getAndDecrement();
		else
			liveCount.getAndDecrement();
		createTableStructure(false);
	}

	/**
	 * Checks if is live.
	 *
	 * @return true, if is live
	 */
	boolean isLive()
	{
		return liveCount.get() != 0;
	}

	/**
	 * Sets the source.
	 *
	 * @param source
	 *            the new source
	 */
	public void setSource(ImageSource source)
	{
		boolean changed = (this.source != source);
		this.source = source;
		createTableStructure(changed);
	}

	/**
	 * Gets the source.
	 *
	 * @return the source
	 */
	public ImageSource getSource()
	{
		return source;
	}

	/**
	 * Sets the configuration.
	 *
	 * @param configuration
	 *            the new configuration
	 */
	public void setConfiguration(String configuration)
	{
		this.configuration = configuration;
	}

	/**
	 * Gets the configuration.
	 *
	 * @return the configuration
	 */
	public String getConfiguration()
	{
		return configuration;
	}

	/**
	 * @return If true show the results deviations in the table
	 */
	public boolean isShowDeviations()
	{
		return showDeviations;
	}

	/**
	 * @param showDeviations
	 *            If true show the results deviations in the table
	 */
	public void setShowDeviations(boolean showDeviations)
	{
		boolean changed = this.showDeviations != showDeviations;
		this.showDeviations = showDeviations;
		createTableStructure(changed);
	}

	/**
	 * @return If true show the results end frame in the table
	 */
	public boolean isShowEndFrame()
	{
		return showEndFrame;
	}

	/**
	 * @param showEndFrame
	 *            If true show the results end frame in the table
	 */
	public void setShowEndFrame(boolean showEndFrame)
	{
		boolean changed = this.showEndFrame != showEndFrame;
		this.showEndFrame = showEndFrame;
		createTableStructure(changed);
	}

	/**
	 * @return If true show the results Id in the table
	 */
	public boolean isShowId()
	{
		return showId;
	}

	/**
	 * @param showId
	 *            If true show the results Id in the table
	 */
	public void setShowId(boolean showId)
	{
		boolean changed = this.showId != showId;
		this.showId = showId;
		createTableStructure(changed);
	}

	/**
	 * Checks if showing the Z column.
	 *
	 * @return true, if is show Z
	 */
	public boolean isShowZ()
	{
		return showZ;
	}

	/**
	 * Set to true to show the Z column.
	 *
	 * @param showZ
	 *            the new show Z
	 */
	public void setShowZ(boolean showZ)
	{
		boolean changed = this.showZ != showZ;
		this.showZ = showZ;
		createTableStructure(changed);
	}

	/**
	 * If true then all results will be checked against the current contents before
	 * allowing an addition.
	 *
	 * @return true, if checking duplicates
	 */
	public boolean isCheckDuplicates()
	{
		return checkForDuplicates;
	}

	/**
	 * Sets the check duplicates flag. If true then all results will be checked against the current contents before
	 * allowing an addition.
	 *
	 * @param checkForDuplicates
	 *            the new check duplicates flag
	 */
	public void setCheckDuplicates(boolean checkForDuplicates)
	{
		this.checkForDuplicates = checkForDuplicates;
	}

	/**
	 * Sets the rounding precision.
	 *
	 * @param roundingPrecision
	 *            the new rounding precision
	 */
	public void setRoundingPrecision(int roundingPrecision)
	{
		rounder = RounderFactory.create(roundingPrecision);
		fireTableChanged(new TableModelEvent(this, 0, 0, 0, RENDERER));
	}

	/**
	 * Gets the calibration.
	 *
	 * @return the calibration
	 */
	public Calibration getCalibration()
	{
		return calibration;
	}

	/**
	 * Gets the psf.
	 *
	 * @return the psf
	 */
	public PSF getPSF()
	{
		return psf;
	}
	//*************************************************************************/
	// Table cell appearance (Rendering)
	//*************************************************************************/

	/**
	 * Gets the float renderer.
	 *
	 * @return the float renderer
	 */
	TableCellRenderer getFloatRenderer()
	{
		return createTableCellRenderer(new DefaultTableCellRenderer()
		{
			private static final long serialVersionUID = 6315482677057619115L;

			@Override
			protected void setValue(Object value)
			{
				setText((value == null) ? "" : rounder.toString((Float) value));
			}
		});
	}

	/**
	 * Gets the double renderer.
	 *
	 * @return the double renderer
	 */
	TableCellRenderer getDoubleRenderer()
	{
		return createTableCellRenderer(new DefaultTableCellRenderer()
		{
			private static final long serialVersionUID = 476334804032977020L;

			@Override
			protected void setValue(Object value)
			{
				setText((value == null) ? "" : rounder.toString((Double) value));
			}
		});
	}

	/**
	 * Gets the double renderer.
	 *
	 * @return the double renderer
	 */
	TableCellRenderer getIntegerRenderer()
	{
		return createTableCellRenderer(new DefaultTableCellRenderer()
		{
			private static final long serialVersionUID = 9099181306853421196L;
		});
	}

	private TableCellRenderer createTableCellRenderer(DefaultTableCellRenderer renderer)
	{
		renderer.setHorizontalAlignment(SwingConstants.TRAILING);
		return renderer;
	}
}
