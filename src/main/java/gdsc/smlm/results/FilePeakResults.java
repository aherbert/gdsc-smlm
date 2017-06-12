package gdsc.smlm.results;

import java.io.FileOutputStream;
import java.io.IOException;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

import gdsc.smlm.data.units.DistanceUnit;
import gdsc.smlm.data.units.IntensityUnit;
import gdsc.smlm.data.units.ConversionException;
import gdsc.smlm.data.units.TypeConverter;
import gdsc.smlm.data.units.UnitConverterFactory;

/**
 * Saves the fit results to file
 */
public abstract class FilePeakResults extends AbstractPeakResults
{
	/** Converter to change the distances to nm. It is created in {@link #begin()} but may be null. */
	protected TypeConverter<DistanceUnit> toNMConverter;
	/** The nm per pixel if calibrated. */
	protected double nmPerPixel;
	/** Converter to change the intensity to photons. It is created in {@link #begin()} but may be null. */
	protected TypeConverter<IntensityUnit> toPhotonConverter;

	// Only write to a single results file
	protected FileOutputStream fos = null;

	protected String filename;
	private boolean sortAfterEnd = false;

	protected int size = 0;

	public FilePeakResults(String filename)
	{
		this.filename = filename;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.PeakResults#begin()
	 */
	public void begin()
	{
		createPrecisionConverters();

		fos = null;
		size = 0;
		try
		{
			fos = new FileOutputStream(filename);
			openOutput();
			write(createResultsHeader());
		}
		catch (Exception e)
		{
			// TODO - Add better handling of errors
			e.printStackTrace();
			closeOutput();
		}
	}

	/**
	 * Open the required output from the open file output stream
	 */
	protected abstract void openOutput();

	/**
	 * Write the data to the output.
	 *
	 * @param data
	 *            the data
	 */
	protected abstract void write(String data);

	/**
	 * Write the result and increment the size by the count.
	 *
	 * @param count
	 *            the count
	 * @param result
	 *            the result
	 */
	protected synchronized void writeResult(int count, String result)
	{
		// In case another thread caused the output to close
		if (fos == null)
			return;
		size += count;
		write(result);
	}

	/**
	 * Creates the standard converters for computing distance in nm and intensity in photons for use in a precision
	 * computation.
	 */
	protected void createPrecisionConverters()
	{
		toNMConverter = null;
		nmPerPixel = 0;
		toPhotonConverter = null;

		if (calibration != null)
		{
			// Create converters 
			if (calibration.hasNmPerPixel())
			{
				nmPerPixel = calibration.getNmPerPixel();
				if (calibration.hasDistanceUnit())
				{
					try
					{
						toNMConverter = UnitConverterFactory.createConverter(calibration.getDistanceUnit(),
								DistanceUnit.NM, nmPerPixel);
					}
					catch (ConversionException e)
					{
						// Gracefully fail so ignore this
					}
				}
			}
			if (calibration.hasIntensityUnit())
			{
				if (calibration.hasGain())
				{
					try
					{
						toPhotonConverter = UnitConverterFactory.createConverter(calibration.getIntensityUnit(),
								IntensityUnit.PHOTON, calibration.getGain());
					}
					catch (ConversionException e)
					{
						// Gracefully fail so ignore this
					}
				}
			}
		}
	}

	protected String createResultsHeader()
	{
		StringBuilder sb = new StringBuilder();

		addComment(sb, getHeaderTitle());
		sb.append(String.format("#FileVersion %s\n", getVersion()));

		// Add the standard details
		if (name != null)
			sb.append(String.format("#Name %s\n", singleLine(name)));
		if (source != null)
			sb.append(String.format("#Source %s\n", singleLine(source.toXML())));
		if (bounds != null)
			sb.append(String.format("#Bounds x%d y%d w%d h%d\n", bounds.x, bounds.y, bounds.width, bounds.height));
		if (calibration != null)
			sb.append(String.format("#Calibration %s\n", singleLine(toXML(calibration))));
		if (configuration != null && configuration.length() > 0)
			sb.append(String.format("#Configuration %s\n", singleLine(configuration)));

		// Add any extra comments
		String[] comments = getHeaderComments();
		if (comments != null)
		{
			for (String comment : comments)
				addComment(sb, comment);
		}

		// Output the field names
		String[] fields = getFieldNames();
		if (fields != null)
		{
			sb.append("#");
			for (int i = 0; i < fields.length; i++)
			{
				if (i != 0)
					sb.append("\t");
				sb.append(fields[i]);
			}
			sb.append('\n');
		}

		addComment(sb, getHeaderEnd());

		return sb.toString();
	}

	private static XStream xs = null;

	/**
	 * Convert the calibration to XML.
	 * <p>
	 * This method was added in version 3 to support an update to the calibration
	 *
	 * @param calibration
	 *            the calibration
	 * @return the string
	 */
	private static String toXML(Calibration calibration)
	{
		if (xs == null)
		{
			xs = new XStream(new DomDriver());
			// Control what is serialised
			xs.omitField(Calibration.class, "fields");
			// This field is deprecated
			xs.omitField(Calibration.class, "emCCD");
			// Do not ignore missing fields. This is so XStream deserialisation
			// will read the invalid values that are set in the object.
		}
		if (xs != null)
		{
			try
			{
				return xs.toXML(calibration);
			}
			catch (XStreamException ex)
			{
				//ex.printStackTrace();
			}
		}
		return null;
	}

	private void addComment(StringBuilder sb, String comment)
	{
		if (comment != null)
			sb.append("#").append(comment).append('\n');
	}

	/**
	 * @return The first line added to the header
	 */
	protected String getHeaderTitle()
	{
		return "Localisation Results File";
	}

	/**
	 * @return The last line added to the header (e.g. a header end tag)
	 */
	protected String getHeaderEnd()
	{
		return null;
	}

	/**
	 * @return A line containing the file format version
	 */
	protected abstract String getVersion();

	/**
	 * @return Any comment lines to add to the header after the standard output of source, name, bounds, etc.
	 */
	protected String[] getHeaderComments()
	{
		return null;
	}

	/**
	 * @return The names of the fields in each record. Will be the last comment of the header
	 */
	protected abstract String[] getFieldNames();

	protected String singleLine(String text)
	{
		return text.replaceAll("\n *", "");
	}

	protected void closeOutput()
	{
		if (fos == null)
			return;

		try
		{
			fos.close();
		}
		catch (Exception e)
		{
			// Ignore exception
		}
		finally
		{
			fos = null;
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
		if (fos == null)
			return;

		// Close the file.
		try
		{
			closeOutput();

			if (!isSortAfterEnd())
				return;

			sort();
		}
		catch (IOException e)
		{
			// ignore
		}
		finally
		{
			fos = null;
		}
	}

	/**
	 * Sort the data file records. This is called once the file has been closed for input.
	 *
	 * @throws IOException
	 *             if an IO error occurs
	 */
	protected abstract void sort() throws IOException;

	/**
	 * @param sortAfterEnd
	 *            True if the results should be sorted after the {@link #end()} method
	 */
	public void setSortAfterEnd(boolean sortAfterEnd)
	{
		this.sortAfterEnd = sortAfterEnd;
	}

	/**
	 * @return True if the results should be sorted after the {@link #end()} method
	 */
	public boolean isSortAfterEnd()
	{
		return sortAfterEnd;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.utils.fitting.results.PeakResults#isActive()
	 */
	public boolean isActive()
	{
		return fos != null;
	}

	/**
	 * @return true if the records are stored as binary data
	 */
	public boolean isBinary()
	{
		return false;
	}
}
