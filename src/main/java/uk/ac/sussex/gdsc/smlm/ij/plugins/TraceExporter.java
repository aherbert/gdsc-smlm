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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gnu.trove.set.hash.TIntHashSet;
import ij.IJ;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterFactory;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.MultiDialog.MemoryResultsItems;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.predicates.PeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XYRResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Plugin to export traced datasets.
 */
public class TraceExporter implements PlugIn
{
    private enum ExportFormat implements NamedObject
    {
        SPOT_ON("Spot-On");

        private final String name;

        ExportFormat(String name)
        {
            this.name = name;
        }

        @Override
        public String getName()
        {
            return name;
        }

        @Override
        public String getShortName()
        {
            return name;
        }
    }

    private static final String TITLE = "Trace Exporter";
    private static ArrayList<String> selected;
    private static String directory = "";
    private static int minLength = 2;
    private static int maxJump = 1;
    private static double wobble = 0;

    private static String[] FORMAT_NAMES;
    private static int format = 0;

    private ExportFormat exportFormat;

    /*
     * (non-Javadoc)
     *
     * @see ij.plugin.PlugIn#run(java.lang.String)
     */
    @Override
    public void run(String arg)
    {
        SMLMUsageTracker.recordPlugin(this.getClass(), arg);

        final MemoryResultsItems items = new MemoryResultsItems(new MultiDialog.MemoryResultsFilter()
        {
            @Override
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

        final ArrayList<MemoryPeakResults> allResults = new ArrayList<>();

        // Pick multiple input datasets together using a list box.
        if (!showMultiDialog(allResults, items))
            return;

        exportFormat = getExportFormat();

        for (final MemoryPeakResults results : allResults)
            export(results);
    }

    private static boolean showDialog()
    {
        if (FORMAT_NAMES == null)
            FORMAT_NAMES = SettingsManager.getNames((Object[]) ExportFormat.values());
        final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
        gd.addMessage("Export traces to a directory");
        gd.addDirectoryField("Directory", directory, 30);
        gd.addNumericField("Min_length", minLength, 0);
        gd.addMessage("Specify the maximum jump allowed within a trace.\nTraces with larger jumps will be split.");
        gd.addNumericField("Max_jump", maxJump, 0);
        gd.addMessage("Specify localistion precision (wobble) to add");
        gd.addNumericField("Wobble", wobble, 0, 6, "nm");
        gd.addChoice("Format", FORMAT_NAMES, FORMAT_NAMES[format]);
        gd.showDialog();
        if (gd.wasCanceled())
            return false;
        directory = gd.getNextString();
        minLength = (int) Math.abs(gd.getNextNumber());
        maxJump = (int) Math.abs(gd.getNextNumber());
        wobble = Math.abs(gd.getNextNumber());
        format = gd.getNextChoiceIndex();
        return true;
    }

    private static boolean showMultiDialog(ArrayList<MemoryPeakResults> allResults, MemoryResultsItems items)
    {
        // Show a list box containing all the results. This should remember the last set of chosen items.
        final MultiDialog md = new MultiDialog(TITLE, items);
        md.addSelected(selected);

        md.showDialog();

        if (md.wasCancelled())
            return false;

        selected = md.getSelectedResults();
        if (selected.isEmpty())
        {
            IJ.error(TITLE, "No results were selected");
            return false;
        }

        for (final String name : selected)
        {
            final MemoryPeakResults r = MemoryPeakResults.getResults(name);
            if (r != null)
                allResults.add(r);
        }

        return !allResults.isEmpty();
    }

    private static ExportFormat getExportFormat()
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
            @Override
            public boolean test(PeakResult t)
            {
                return t.getId() <= 0;
            }
        });

        // Sort by ID then time
        results.sort(IdFramePeakResultComparator.INSTANCE);

        // Split traces with big jumps
        results = splitTraces(results);

        // Count each ID and remove short traces
        int id = 0;
        int count = 0;
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
        // Final ID
        if (count < minLength)
            remove.add(id);

        if (!remove.isEmpty())
        {
            results.removeIf(new PeakResultPredicate()
            {
                @Override
                public boolean test(PeakResult t)
                {
                    return remove.contains(t.getId());
                }
            });
            results.sort(IdFramePeakResultComparator.INSTANCE);
        }

        if (wobble > 0)
        {
            // Just leave any exceptions to trickle up and kill the plugin
            final TypeConverter<DistanceUnit> c = results.getDistanceConverter(DistanceUnit.NM);
            final double w = c.convertBack(wobble);
            final RandomGenerator r = new Well19937c();
            final boolean is3D = results.is3D();
            results.forEach(new PeakResultProcedure()
            {
                @Override
                public void execute(PeakResult peakResult)
                {
                    peakResult.setXPosition((float) (peakResult.getXPosition() + w * r.nextGaussian()));
                    peakResult.setYPosition((float) (peakResult.getYPosition() + w * r.nextGaussian()));
                    if (is3D)
                        peakResult.setZPosition((float) (peakResult.getZPosition() + w * r.nextGaussian()));
                }
            });

        }

        switch (exportFormat)
        {
            case SPOT_ON:
            default:
                exportSpotOn(results);
        }
    }

    private static MemoryPeakResults splitTraces(MemoryPeakResults results)
    {
        if (maxJump < 1)
            // Disabled
            return results;

        int id = 0;
        final int lastT = 0;
        for (int i = 0, size = results.size(); i < size; i++)
        {
            final PeakResult r = results.get(i);
            if (r.getId() != id)
                id = r.getId();
            else if (r.getFrame() - lastT > maxJump)
                return doSplit(results);
        }
        return results;
    }

    private static MemoryPeakResults doSplit(MemoryPeakResults results)
    {
        final MemoryPeakResults results2 = new MemoryPeakResults(results.size());
        results2.copySettings(results);
        int nextId = results.getLast().getId();
        int id = 0, idOut = 0;
        int lastT = 0;
        for (int i = 0, size = results.size(); i < size; i++)
        {
            final PeakResult r = results.get(i);
            if (r.getId() != id)
            {
                id = r.getId();
                idOut = id;
            }
            else if (r.getFrame() - lastT > maxJump)
                idOut = ++nextId;
            final AttributePeakResult r2 = new AttributePeakResult(r);
            r2.setId(idOut);
            results2.add(r2);
            lastT = r.getEndFrame();
        }
        return results2;
    }

    private static void exportSpotOn(MemoryPeakResults results)
    {
        // Simple Spot-On CSV file format:
        // https://spoton.berkeley.edu/SPTGUI/docs/latest#input-formats
        // frame, t (seconds), trajectory (trace id), x (um), y (um)

        try (BufferedWriter out = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(new File(directory, results.getName() + ".csv")), "UTF-8")))
        {
            out.write("frame,t,trajectory,x,y");
            out.newLine();

            final TypeConverter<TimeUnit> converter = UnitConverterFactory.createConverter(TimeUnit.FRAME,
                    TimeUnit.SECOND, results.getCalibrationReader().getExposureTime());

            @SuppressWarnings("resource")
            final BufferedWriter writer = out;
            results.forEach(DistanceUnit.UM, new XYRResultProcedure()
            {
                @Override
                public void executeXYR(float x, float y, PeakResult result)
                {
                    try
                    {
                        if (result.hasEndFrame())
                        {
                            final String sId = Integer.toString(result.getId());
                            final String sx = Float.toString(x);
                            final String sy = Float.toString(y);
                            for (int t = result.getFrame(); t <= result.getEndFrame(); t++)
                            {
                                writer.write(Integer.toString(t));
                                writer.write(",");
                                writer.write(Float.toString(converter.convert(t)));
                                writer.write(",");
                                writer.write(sId);
                                writer.write(",");
                                writer.write(sx);
                                writer.write(",");
                                writer.write(sy);
                                writer.newLine();
                            }
                        }
                        else
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
                    }
                    catch (final IOException e)
                    {
                        // Allow clean-up by passing the exception up
                        throw new RuntimeException(e);
                    }
                }
            });
        }
        catch (final Exception e)
        {
            e.printStackTrace();
        }
    }
}
