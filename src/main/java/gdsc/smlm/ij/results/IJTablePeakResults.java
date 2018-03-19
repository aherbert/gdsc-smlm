package gdsc.smlm.ij.results;

import java.awt.Frame;
import java.awt.event.MouseListener;
import java.util.HashMap;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.Converter;
import gdsc.core.data.utils.IdentityTypeConverter;
import gdsc.core.data.utils.Rounder;
import gdsc.core.data.utils.RounderFactory;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.data.config.UnitProtos.AngleUnit;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.ij.utils.CoordinateProvider;
import gdsc.smlm.ij.utils.ImageROIPainter;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultsHelper;
import gnu.trove.list.array.TIntArrayList;
import ij.WindowManager;
import ij.text.TextPanel;
import ij.text.TextWindow;

/**
 * Saves the fit results to an ImageJ results table.
 * <p>
 * The table supports mouse click events to draw the selected coordinates on the original source image using the
 * ImageROIPainter.
 */
public class IJTablePeakResults extends IJAbstractPeakResults implements CoordinateProvider
{
	/**
	 * Converter to change the distances to pixels. It is created in {@link #begin()} and may be an identity converter.
	 */
	protected TypeConverter<DistanceUnit> toPixelConverter;

	/** The calculator used to compute precision. */
	private Gaussian2DPeakResultCalculator calculator;

	private boolean canComputePrecision = false;

	private PeakResultsHelper helper;
	private Converter[] converters;
	private Converter ic;
	private int[] outIndices;

	private DistanceUnit distanceUnit = null;
	private IntensityUnit intensityUnit = null;
	private AngleUnit angleUnit = null;
	private boolean showPrecision = false;
	private boolean computePrecision = false;
	private int roundingPrecision = 0;

	// Store the ROI painters that have been attached to TextPanels so they can be updated
	// with a new image source
	private static HashMap<TextPanel, ImageROIPainter> map = new HashMap<TextPanel, ImageROIPainter>();

	private boolean showDeviations = false;
	private boolean showEndFrame = false;
	private boolean showId = false;
	private boolean showFittingData = false;
	private boolean showNoiseData = false;
	private boolean showZ = false;
	private boolean clearAtStart = false;
	private boolean hideSourceText = false;
	private String frameColumnName = "T";
	private String source = null;
	private String sourceText = null;
	private String tableTitle = "Fit Results";
	private boolean newWindow;
	private TextWindow resultsWindow;
	private TextPanel tp;
	private ImageROIPainter roiPainter;
	private boolean addCounter = false;
	protected boolean tableActive = false;
	private int nextRepaintSize = 0;
	private double repaintInterval = 0.1;
	private Rounder rounder;

	private int indexT = -1, indexX = -1, indexY = -1;

	private int size = 0;

	public IJTablePeakResults(boolean showDeviations)
	{
		this.showDeviations = showDeviations;
	}

	public IJTablePeakResults(boolean showDeviations, String source)
	{
		this.showDeviations = showDeviations;
		this.source = source;
	}

	public IJTablePeakResults(boolean showDeviations, String source, boolean clearAtStart)
	{
		this.showDeviations = showDeviations;
		this.source = source;
		this.clearAtStart = clearAtStart;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#begin()
	 */
	public void begin()
	{
		tableActive = false;

		// Set-up unit processing that requires the calibration 
		toPixelConverter = new IdentityTypeConverter<DistanceUnit>(null);
		calculator = null;
		canComputePrecision = false;
		rounder = RounderFactory.create(roundingPrecision);

		// We must correctly convert all the PSF parameter types
		helper = new PeakResultsHelper(getCalibration(), getPSF());
		helper.setIntensityUnit(intensityUnit);
		helper.setDistanceUnit(distanceUnit);
		helper.setAngleUnit(angleUnit);
		converters = helper.getConverters();

		if (hasCalibration())
		{
			if (showPrecision)
			{
				if (computePrecision)
				{
					try
					{
						calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibrationReader(),
								Gaussian2DPeakResultHelper.LSE_PRECISION);
						canComputePrecision = true;
					}
					catch (ConfigurationException e)
					{
					}
					catch (ConversionException e)
					{
					}
				}
			}

			try
			{
				if (helper.hasDistanceConverter())
				{
					toPixelConverter = UnitConverterFactory.createConverter(distanceUnit, DistanceUnit.PIXEL,
							getCalibrationReader().getNmPerPixel());
				}
			}
			catch (ConversionException e)
			{
				// Gracefully fail so ignore this
			}
		}

		ic = converters[PeakResult.INTENSITY];
		outIndices = SimpleArrayUtils.newArray(converters.length, 0, 1);
		if (!showZ)
		{
			TIntArrayList list = new TIntArrayList(outIndices);
			list.remove(PeakResult.Z);
			outIndices = list.toArray();
		}
		// Update the calibration if converters were created
		if (helper.isCalibrationChanged())
			setCalibration(helper.getCalibration());

		createSourceText();
		createResultsWindow();
		if (clearAtStart)
		{
			tp.clear();
		}
		size = 0;
		// Let some results appear before drawing.
		// ImageJ will auto-layout columns if it has less than 10 rows
		nextRepaintSize = 9;
		tableActive = true;
	}

	/**
	 * Clear the table contents.
	 */
	public void clear()
	{
		tp.clear();
		size = 0;
		// Let some results appear before drawing.
		// ImageJ will auto-layout columns if it has less than 10 rows
		nextRepaintSize = 9;
	}

	/**
	 * Create the result window (if it is not available)
	 */
	private void createResultsWindow()
	{
		String header = createResultsHeader();

		roiPainter = null;
		for (Frame f : WindowManager.getNonImageWindows())
		{
			if (f != null && tableTitle.equals(f.getTitle()) && f instanceof TextWindow)
			{
				resultsWindow = (TextWindow) f;

				// Check if the existing table matches the desired header
				String currentHeader = resultsWindow.getTextPanel().getColumnHeadings();
				if (!currentHeader.startsWith(header))
				{
					resultsWindow = null;
					continue;
				}

				roiPainter = map.get(resultsWindow.getTextPanel());
				break;
			}
		}

		newWindow = false;
		if (resultsWindow == null || !resultsWindow.isShowing())
		{
			newWindow = true;
			resultsWindow = new TextWindow(tableTitle, header, "", 800, 300);
			roiPainter = new ImageROIPainter(resultsWindow.getTextPanel(), "", this);

			// The ROI painter adds itself to the TextPanel as a mouse listener. However
			// the TextPanel addMouseListener() adds to the private TextCanvas object so it 
			// cannot be retrieved. Store the painter in a global lookup table.
			map.put(resultsWindow.getTextPanel(), roiPainter);
		}

		tp = resultsWindow.getTextPanel();

		if (roiPainter != null && getSource() != null)
		{
			roiPainter.setTitle(getSource().getOriginal().getName());

			// Update the coordinate provider (avoids memory leaks with old objects lying around)
			roiPainter.setCoordProvider(this);

			// Get the headings for extracting the coordinates 
			String[] headings = tp.getColumnHeadings().split("\t");
			indexT = indexX = indexY = -1;
			for (int i = 0; i < headings.length; i++)
			{
				if (headings[i].equals(frameColumnName))
				{
					indexT = i;
					continue;
				}
				// Allow for units
				if (headings[i].equals("X") || headings[i].startsWith("X ("))
				{
					indexX = i;
					continue;
				}
				if (headings[i].equals("Y") || headings[i].startsWith("Y ("))
				{
					indexY = i;
					continue;
				}
			}
		}
	}

	/**
	 * Checks if is new window.
	 *
	 * @return true, if is new window
	 */
	public boolean isNewWindow()
	{
		return newWindow;
	}

	private String createResultsHeader()
	{
		String[] names = helper.getNames();
		String[] unitNames = helper.getUnitNames();

		StringBuilder sb = new StringBuilder();
		if (addCounter)
			sb.append("#\t");
		if (sourceText != null)
			sb.append("Source\t");
		sb.append(frameColumnName);
		if (showEndFrame)
			sb.append("\tEnd ").append(frameColumnName);
		if (showId)
			sb.append("\tId");
		if (showFittingData)
		{
			sb.append("\torigX");
			sb.append("\torigY");
			sb.append("\torigValue");
			sb.append("\tError");
		}
		if (showNoiseData)
		{
			sb.append("\tNoise");
			if (!TextUtils.isNullOrEmpty(unitNames[PeakResult.INTENSITY]))
				sb.append(" (").append(unitNames[PeakResult.INTENSITY]).append(')');
			sb.append("\tSNR");
		}

		for (int i = 0; i < outIndices.length; i++)
		{
			sb.append('\t').append(names[outIndices[i]]);
			if (!TextUtils.isNullOrEmpty(unitNames[outIndices[i]]))
				sb.append(" (").append(unitNames[outIndices[i]]).append(')');
			addDeviation(sb);
		}
		if (showPrecision)
		{
			sb.append("\tPrecision (nm)");
		}
		return sb.toString();
	}

	private void createSourceText()
	{
		if (hideSourceText)
		{
			sourceText = null;
			return;
		}
		StringBuilder sb = new StringBuilder();
		if (source != null)
			sb.append(source);
		else if (getSource() != null)
			sb.append(getSource().getName());
		if (getBounds() != null)
		{
			if (sb.length() > 0)
				sb.append(": ");
			sb.append(getBoundsString());
		}
		if (sb.length() > 0)
		{
			sb.append('\t');
			sourceText = sb.toString();
		}
	}

	private void addDeviation(StringBuilder sb)
	{
		if (showDeviations)
			sb.append("\t+/-");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.AbstractPeakResults#add(int, int, int, float, double, float, float[], float[])
	 */
	public void add(int frame, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsDev)
	{
		addPeak(frame, frame, 0, origX, origY, origValue, error, noise, params, paramsDev, -1);
	}

	private void addPeak(int frame, int endFrame, int id, int origX, int origY, float origValue, double error,
			float noise, float[] params, float[] paramsStdDev, double precision)
	{
		if (!tableActive)
			return;

		final float snr = (noise > 0) ? params[PeakResult.INTENSITY] / noise : 0;
		StringBuilder sb = addStandardData(frame, endFrame, id, origX, origY, origValue, error, noise, snr);
		if (isShowDeviations())
		{
			if (paramsStdDev != null)
			{
				for (int i = 0; i < outIndices.length; i++)
				{
					add(sb, converters[outIndices[i]].convert(params[outIndices[i]]));
					add(sb, converters[outIndices[i]].convert(paramsStdDev[outIndices[i]]));
				}
			}
			else
			{
				for (int i = 0; i < outIndices.length; i++)
				{
					add(sb, converters[outIndices[i]].convert(params[outIndices[i]]));
					sb.append("\t0");
				}
			}
		}
		else
		{
			for (int i = 0; i < outIndices.length; i++)
				add(sb, converters[outIndices[i]].convert(params[outIndices[i]]));
		}
		if (isShowPrecision())
		{
			// The default precision in a peak result is NaN so this compare will be false
			if (precision >= 0)
				addPrecision(sb, precision, false);
			else if (canComputePrecision)
				addPrecision(sb, calculator.getLSEPrecision(params, noise), true);
			else
				sb.append("\t0");
		}
		append(sb.toString());
	}

	private StringBuilder addStandardData(int frame, int endFrame, int id, int origX, int origY, float origValue,
			double error, float noise, float snr)
	{
		StringBuilder sb = new StringBuilder();
		if (addCounter)
			sb.append(size + 1).append('\t');
		if (sourceText != null)
			sb.append(sourceText);
		// Do not calibrate the original values		
		//if (showCalibratedValues)
		//	sb.append(frame).append(String.format("\t%g", origX)).append(String.format("\t%g", origY));
		//else
		sb.append(frame);
		if (showEndFrame)
			sb.append('\t').append(endFrame);
		if (showId)
			sb.append('\t').append(id);
		if (showFittingData)
		{
			sb.append('\t').append(origX).append('\t').append(origY);
			add(sb, origValue);
			add(sb, error);
		}
		if (showNoiseData)
		{
			add(sb, ic.convert(noise)); // This should be converted
			add(sb, snr);
		}
		return sb;
	}

	private void add(StringBuilder sb, float value)
	{
		sb.append('\t').append(rounder.toString(value));
	}

	private void add(StringBuilder sb, double value)
	{
		sb.append('\t').append(rounder.toString(value));
	}

	private void addPrecision(StringBuilder sb, double value, boolean computed)
	{
		add(sb, value);
		if (computed)
			sb.append('*');
	}

	private void append(String result)
	{
		// Support for periodic refresh
		synchronized (tp)
		{
			addResult(result);
		}
		updateTable();
	}

	private void addResult(String result)
	{
		size++;
		tp.appendWithoutUpdate(result);
	}

	private void updateTable()
	{
		if (size < nextRepaintSize)
			return;

		if (!resultsWindow.isShowing())
		{
			//System.out.println("Table has been closed");
			tableActive = false;
			return;
		}

		drawTable();
	}

	private void drawTable()
	{
		synchronized (tp)
		{
			//System.out.printf("Refreshing table: size = %d\n", size);
			nextRepaintSize = (int) (size + size * repaintInterval);
			tp.updateDisplay();
		}
	}

	public void add(PeakResult result)
	{
		addPeak(result.getFrame(), result.getEndFrame(), result.getId(), result.getOrigX(), result.getOrigY(),
				result.getOrigValue(), result.getError(), result.getNoise(), result.getParameters(),
				result.getParameterDeviations(), result.getPrecision());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResults#addAll(gdsc.smlm.results.PeakResult[])
	 */
	public void addAll(PeakResult[] results)
	{
		if (!tableActive)
			return;
		int n = 0;
		for (PeakResult result : results)
		{
			addPeak(result.getFrame(), result.getEndFrame(), result.getId(), result.getOrigX(), result.getOrigY(),
					result.getOrigValue(), result.getError(), result.getNoise(), result.getParameters(),
					result.getParameterDeviations(), result.getPrecision());
			if (n++ > 31)
			{
				if (!tableActive)
					return;
				n = 0;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#size()
	 */
	public int size()
	{
		return size;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#end()
	 */
	public void end()
	{
		tableActive = false;
		drawTable();
	}

	/**
	 * Forces the table to be updated with the current contents.
	 */
	public void flush()
	{
		drawTable();
	}

	/**
	 * @return the name of the frame column
	 */
	public String getFrameColumnName()
	{
		return frameColumnName;
	}

	/**
	 * @param frameColumnName
	 *            the name of the frame column
	 */
	public void setFrameColumnName(String frameColumnName)
	{
		this.frameColumnName = frameColumnName;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		return tableActive;
	}

	/**
	 * @return the table title
	 */
	public String getTableTitle()
	{
		return tableTitle;
	}

	/**
	 * Use to set the title of the table. If an existing table exists with the same title then it will be appended,
	 * otherwise a new table is created.
	 * 
	 * @param tableTitle
	 *            the table title
	 */
	public void setTableTitle(String tableTitle)
	{
		if (tableTitle != null && tableTitle.length() > 0)
			this.tableTitle = tableTitle;
	}

	/**
	 * @return True if the deviations of the parameters should be shown
	 */
	public boolean isShowDeviations()
	{
		return showDeviations;
	}

	/**
	 * @param showDeviations
	 *            True if the deviations of the parameters should be shown
	 */
	public void setShowDeviations(boolean showDeviations)
	{
		this.showDeviations = showDeviations;
	}

	/**
	 * @return True if the table should be cleared in {@link #begin()}
	 */
	public boolean isClearAtStart()
	{
		return clearAtStart;
	}

	/**
	 * @param clearAtStart
	 *            True if the table should be cleared in {@link #begin()}
	 */
	public void setClearAtStart(boolean clearAtStart)
	{
		this.clearAtStart = clearAtStart;
	}

	/**
	 * @return the addCounter
	 */
	public boolean isAddCounter()
	{
		return addCounter;
	}

	/**
	 * @param addCounter
	 *            the addCounter to set
	 */
	public void setAddCounter(boolean addCounter)
	{
		this.addCounter = addCounter;
	}

	/**
	 * Checks if the source text will be added to each entry.
	 *
	 * @return true, if hiding the source text
	 */
	public boolean isHideSourceText()
	{
		return hideSourceText;
	}

	/**
	 * Sets the hide source text flag.
	 *
	 * @param hideSourceText
	 *            the new hide source text flag
	 */
	public void setHideSourceText(boolean hideSourceText)
	{
		this.hideSourceText = hideSourceText;
	}

	/**
	 * @return the resultsWindow
	 */
	public TextWindow getResultsWindow()
	{
		return resultsWindow;
	}

	/**
	 * @param line
	 * @return
	 */
	public double[] getCoordinates(String line)
	{
		// Extract the startT and x,y coordinates from the PeakResult line
		String[] fields = line.split("\t");
		try
		{
			int startT = Integer.valueOf(fields[indexT]);
			double x = Double.valueOf(fields[indexX]);
			double y = Double.valueOf(fields[indexY]);
			return new double[] { startT, toPixelConverter.convert(x), toPixelConverter.convert(y) };
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			// Will happen if any index is still at the default of -1 or if there are not enough fields
		}
		catch (NumberFormatException e)
		{
			// In case any field is not a number
		}
		return null;
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
		this.showEndFrame = showEndFrame;
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
		this.showId = showId;
	}

	/**
	 * @return If true then show the fitting data (original pixel data and fit error) in the table
	 */
	public boolean isShowFittingData()
	{
		return showFittingData;
	}

	/**
	 * @param showFittingData
	 *            If true then show the fitting data (original pixel data and fit error) in the table
	 */
	public void setShowFittingData(boolean showFittingData)
	{
		this.showFittingData = showFittingData;
	}

	/**
	 * @return If true then show the noise and SNR in the table
	 */
	public boolean isShowNoiseData()
	{
		return showNoiseData;
	}

	/**
	 * @param showNoiseData
	 *            If true then show the noise and SNR in the table
	 */
	public void setShowNoiseData(boolean showNoiseData)
	{
		this.showNoiseData = showNoiseData;
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
		this.showZ = showZ;
	}

	/**
	 * Image will be repainted when a fraction of new results have been added.
	 * 
	 * @param repaintInterval
	 *            the repaintInterval to set (range 0.001-1)
	 */
	public void setRepaintInterval(double repaintInterval)
	{
		if (repaintInterval < 0.001)
			repaintInterval = 0.001;
		if (repaintInterval > 1)
			repaintInterval = 1;

		this.repaintInterval = repaintInterval;
	}

	/**
	 * @return the repaintInterval
	 */
	public double getRepaintInterval()
	{
		return repaintInterval;
	}

	/**
	 * Gets the distance unit.
	 *
	 * @return the distance unit
	 */
	public DistanceUnit getDistanceUnit()
	{
		return distanceUnit;
	}

	/**
	 * Sets the distance unit.
	 *
	 * @param distanceUnit
	 *            the new distance unit
	 */
	public void setDistanceUnit(DistanceUnit distanceUnit)
	{
		this.distanceUnit = distanceUnit;
	}

	/**
	 * Gets the intensity unit.
	 *
	 * @return the intensity unit
	 */
	public IntensityUnit getIntensityUnit()
	{
		return intensityUnit;
	}

	/**
	 * Sets the intensity unit.
	 *
	 * @param intensityUnit
	 *            the new intensity unit
	 */
	public void setIntensityUnit(IntensityUnit intensityUnit)
	{
		this.intensityUnit = intensityUnit;
	}

	/**
	 * Gets the angle unit.
	 *
	 * @return the angle unit
	 */
	public AngleUnit getAngleUnit()
	{
		return angleUnit;
	}

	/**
	 * Sets the angle unit.
	 *
	 * @param angleUnit
	 *            the new angle unit
	 */
	public void setAngleUnit(AngleUnit angleUnit)
	{
		this.angleUnit = angleUnit;
	}

	/**
	 * Checks if the precision will be computed.
	 *
	 * @return true, if the precision will be computed
	 */
	public boolean isShowPrecision()
	{
		return showPrecision;
	}

	/**
	 * Sets the show precision flag.
	 *
	 * @param showPrecision
	 *            set to true to show the precision and write to the output
	 */
	public void setShowPrecision(boolean showPrecision)
	{
		this.showPrecision = showPrecision;
	}

	/**
	 * Checks if the precision will be computed if needed. This is only relevant if show precision is true (see
	 * {@link #isShowPrecision()}).
	 *
	 * @return true, if the precision will be computed
	 */
	public boolean isComputePrecision()
	{
		return computePrecision;
	}

	/**
	 * Sets the compute precision flag. This is only relevant if show precision is true (see
	 * {@link #isShowPrecision()}).
	 *
	 * @param computePrecision
	 *            set to true to compute the precision if needed
	 */
	public void setComputePrecision(boolean computePrecision)
	{
		this.computePrecision = computePrecision;
	}

	/**
	 * Gets the rounding precision.
	 *
	 * @return the rounding precision
	 */
	public int getRoundingPrecision()
	{
		return roundingPrecision;
	}

	/**
	 * Sets the rounding precision.
	 *
	 * @param roundingPrecision
	 *            the new rounding precision
	 */
	public void setRoundingPrecision(int roundingPrecision)
	{
		this.roundingPrecision = roundingPrecision;
	}

	/**
	 * Select an index from the text panel.
	 *
	 * @param selectedIndex
	 *            the selected index
	 */
	public void select(int selectedIndex)
	{
		if (selectedIndex < 0 || selectedIndex >= tp.getLineCount())
			return;
		tp.setSelection(selectedIndex, selectedIndex);
		if (roiPainter != null)
			roiPainter.selected(selectedIndex);
	}

	/**
	 * Select a range of indices from the text panel.
	 *
	 * @param selectionStart
	 *            the selection start
	 * @param selectionEnd
	 *            the selection end
	 */
	public void select(int selectionStart, int selectionEnd)
	{
		if (selectionStart < 0 || selectionStart >= tp.getLineCount())
			return;
		if (selectionEnd < selectionStart || selectionEnd >= tp.getLineCount())
			return;
		tp.setSelection(selectionStart, selectionEnd);
		if (roiPainter != null)
			roiPainter.selected(selectionStart, selectionEnd);
	}

	/**
	 * Adds the mouse listener to the text panel.
	 *
	 * @param listener
	 *            the listener
	 */
	public void addTextPanelMouseListener(MouseListener listener)
	{
		tp.addMouseListener(listener);
	}

	/**
	 * Removes the mouse listener to the text panel.
	 *
	 * @param listener
	 *            the listener
	 */
	public void removeTextPanelMouseListener(MouseListener listener)
	{
		tp.removeMouseListener(listener);
	}
	
	/**
	 * Gets the text panel.
	 *
	 * @return the text panel
	 */
	public TextPanel getTextPanel()
	{
		return tp;
	}

	// TODO - Extend the selection functionality so that selectionListeners can be added
	// to receive notification whan items from the table are selected. 
	//	public void addSelectionListener()
	//	{
	//		
	//	}
	//	
	//	public void removeSelectionListener()
	//	{
	//		
	//	}
}
