package gdsc.smlm.ij.plugins;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Comparator;

import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.NamedObject;
import gdsc.smlm.data.config.UnitConverterFactory;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.TimeUnit;
import gdsc.smlm.ij.plugins.MultiDialog.MemoryResultsItems;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.PeakResultPredicate;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

/**
 * Plugin to export traced datasets.
 */
public class TraceExporter implements PlugIn
{
	private enum ExportFormat implements NamedObject
	{
		SPOT_ON("Spot-On");

		private String name;

		ExportFormat(String name)
		{
			this.name = name;
		}

		public String getName()
		{
			return name;
		}

		public String getShortName()
		{
			return name;
		}
	}

	private static final String TITLE = "Trace Exporter";
	private static ArrayList<String> selected;
	private static String directory = "";
	private static int minLength = 2;

	private static Comparator<PeakResult> comp;
	private static String[] FORMAT_NAMES;
	private static int format = 0;

	private ExportFormat exportFormat;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		MemoryResultsItems items = new MemoryResultsItems(new MultiDialog.MemoryResultsFilter()
		{
			public boolean accept(MemoryPeakResults results)
			{
				return results.hasId();
			}
		});

		if (items.size() == 0)
		{
			IJ.error(TITLE, "No traced localisations in memory");
			return;
		}

		// Get input options
		if (!showDialog())
			return;

		ArrayList<MemoryPeakResults> allResults = new ArrayList<MemoryPeakResults>();

		// Pick multiple input datasets together using a list box.
		if (!showMultiDialog(allResults, items))
			return;

		if (comp == null)
		{
			comp = new Comparator<PeakResult>()
			{
				public int compare(PeakResult o1, PeakResult o2)
				{
					int result = o1.getId() - o2.getId();
					if (result != 0)
						return result;
					return o1.getFrame() - o2.getFrame();
				}
			};
		}

		exportFormat = getExportFormat();

		for (MemoryPeakResults results : allResults)
			export(results);
	}

	private boolean showDialog()
	{
		if (FORMAT_NAMES == null)
		{
			FORMAT_NAMES = SettingsManager.getNames((Object[]) ExportFormat.values());
		}
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Export traces to a directory");
		gd.addDirectoryField("Directory", directory, 30);
		gd.addSlider("Min_length", 2, 20, minLength);
		gd.addChoice("Format", FORMAT_NAMES, FORMAT_NAMES[format]);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;
		directory = gd.getNextString();
		minLength = (int) Math.abs(gd.getNextNumber());
		format = gd.getNextChoiceIndex();
		return true;
	}

	private boolean showMultiDialog(ArrayList<MemoryPeakResults> allResults, MemoryResultsItems items)
	{
		// Show a list box containing all the results. This should remember the last set of chosen items.
		MultiDialog md = new MultiDialog(TITLE, items);
		md.addSelected(selected);

		md.showDialog();

		if (md.wasCanceled())
			return false;

		selected = md.getSelectedResults();
		if (selected.isEmpty())
		{
			IJ.error(TITLE, "No results were selected");
			return false;
		}

		for (String name : selected)
		{
			MemoryPeakResults r = MemoryPeakResults.getResults(name);
			if (r != null)
				allResults.add(r);
		}

		return !allResults.isEmpty();
	}

	private ExportFormat getExportFormat()
	{
		if (format >= 0 && format < FORMAT_NAMES.length)
			return ExportFormat.values()[format];
		return ExportFormat.SPOT_ON;
	}

	private void export(MemoryPeakResults results)
	{
		// Copy to allow manipulation
		results = results.copy();

		// Strip results with no trace Id
		results.removeIf(new PeakResultPredicate()
		{
			public boolean test(PeakResult t)
			{
				return t.getId() <= 0;
			}
		});

		// Sort by ID then time
		results.sort(comp);

		// Count each ID and remove short traces
		int count = 0, id = 0;
		final TIntHashSet remove = new TIntHashSet();
		for (int i = 0, size = results.size(); i < size; i++)
		{
			if (results.get(i).getId() != id)
			{
				if (count < minLength)
					remove.add(id);
				count = 0;
				id = results.get(i).getId();
			}
			count++;
		}

		if (!remove.isEmpty())
		{
			results.removeIf(new PeakResultPredicate()
			{
				public boolean test(PeakResult t)
				{
					return remove.contains(t.getId());
				}
			});
			results.sort(comp);
		}

		switch (exportFormat)
		{
			case SPOT_ON:
			default:
				exportSpotOn(results);
		}
	}

	private void exportSpotOn(MemoryPeakResults results)
	{
		// Simple Spot-On CSV file format:
		// https://spoton.berkeley.edu/SPTGUI/docs/latest#input-formats
		// frame, t (seconds), trajectory (trace id), x (um), y (um)

		BufferedWriter out = null;
		try
		{
			File file = new File(directory, results.getName() + ".csv");
			FileOutputStream fos = new FileOutputStream(file);
			out = new BufferedWriter(new OutputStreamWriter(fos, "UTF-8"));
			out.write("frame,t,trajectory,x,y");
			out.newLine();

			final TypeConverter<TimeUnit> converter = UnitConverterFactory.createConverter(TimeUnit.FRAME,
					TimeUnit.SECOND, results.getCalibrationReader().getExposureTime());

			final BufferedWriter writer = out;
			results.forEach(DistanceUnit.UM, new XYRResultProcedure()
			{
				public void executeXYR(float x, float y, PeakResult result)
				{
					try
					{
						writer.write(Integer.toString(result.getFrame()));
						writer.write(",");
						writer.write(Float.toString(converter.convert(result.getFrame())));
						writer.write(",");
						writer.write(Integer.toString(result.getId()));
						writer.write(",");
						writer.write(Float.toString(x));
						writer.write(",");
						writer.write(Float.toString(y));
						writer.newLine();
					}
					catch (IOException e)
					{
						// Allow clean-up by passing the exception up
						throw new RuntimeException(e);
					}
				}
			});
		}
		catch (Exception e)
		{
		}
		finally
		{
			if (out != null)
			{
				try
				{
					out.close();
				}
				catch (IOException e)
				{
				}
			}
		}
	}
}
