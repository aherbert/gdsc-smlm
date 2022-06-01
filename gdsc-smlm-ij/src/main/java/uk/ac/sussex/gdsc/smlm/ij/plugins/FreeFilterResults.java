/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import ij.IJ;
import ij.plugin.PlugIn;
import java.awt.Checkbox;
import java.awt.event.ItemEvent;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.filter.AndFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.AnrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.CoordinateFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.EShiftFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilter2;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiHysteresisFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiHysteresisFilter2;
import uk.ac.sussex.gdsc.smlm.results.filter.OrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionFilter2;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionHysteresisFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SbrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.ShiftFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SignalFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SnrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SnrHysteresisFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.TraceFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.WidthFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.WidthFilter2;
import uk.ac.sussex.gdsc.smlm.results.filter.XyWidthFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.XyWidthFilter2;

/**
 * Filters PeakFit results that are stored in memory using the configured filters.
 */
public class FreeFilterResults implements PlugIn {
  private static final String TITLE = "Free Filter Results";
  private static final AtomicReference<String> INPUT_OPTION_REF = new AtomicReference<>("");

  private String inputOption;
  private GUIFilterSettings.Builder filterSettings;

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      // Ask user if they want to show the demo filters
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.enableYesNoCancel();
      gd.hideCancelButton();
      gd.addMessage("No results in memory. Show the demo filters?");
      gd.showDialog();
      if (gd.wasOKed()) {
        logDemoFilters(TITLE);
      }
      return;
    }

    if (!showDialog()) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(inputOption, false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }

    // Filter results
    final Filter filter = Filter.fromXml(filterSettings.getFreeFilter());
    if (filter == null) {
      IJ.showStatus("ERROR: Unable to create filter");
    } else {
      final MemoryPeakResults newResults = filter.filter(results);
      if (newResults.size() > 0) {
        newResults.setName(results.getName() + " Free Filtered");
        MemoryPeakResults.addResults(newResults);
      }
      IJ.showStatus(String.format("Filtered %d results to %d", results.size(), newResults.size()));
    }
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("free-filter-results"));

    gd.addMessage("Select a dataset to filter");
    ResultsManager.addInput(gd, INPUT_OPTION_REF.get(), InputSource.MEMORY);

    filterSettings = SettingsManager.readGuiFilterSettings(0).toBuilder();

    String text;
    try {
      text = uk.ac.sussex.gdsc.core.utils.XmlUtils.prettyPrintXml(filterSettings.getFreeFilter());
    } catch (final Exception ex) {
      text = filterSettings.getFreeFilter();
    }
    gd.addTextAreas(text, null, 20, 80);
    gd.addCheckbox("Show_demo_filters", false);

    if (ImageJUtils.isShowGenericDialog()) {
      final Checkbox cb = (Checkbox) gd.getCheckboxes().get(0);
      cb.addItemListener(this::itemStateChanged);
    }

    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    inputOption = ResultsManager.getInputSource(gd);
    filterSettings.setFreeFilter(gd.getNextText());
    final boolean demoFilters = gd.getNextBoolean();

    if (demoFilters) {
      logDemoFilters(TITLE);
      return false;
    }

    INPUT_OPTION_REF.set(inputOption);
    return SettingsManager.writeSettings(filterSettings.build());
  }

  private void itemStateChanged(ItemEvent event) {
    // When the checkbox is clicked, output the list of available filters to the ImageJ log

    final Checkbox cb = (Checkbox) event.getSource();
    if (cb.getState()) {
      cb.setState(false);

      logDemoFilters(TITLE);
    }
  }

  private static void logDemoFilters(String title) {
    comment(title + " example filters");
    IJ.log("");
    comment("Filters are described using XML");
    comment("Multiple filters can be combined using AND/OR filters");
    IJ.log("");

    comment("Single filters");
    IJ.log("");
    demo(new WidthFilter(2));
    demo(new WidthFilter2(0.7, 2));
    demo(new XyWidthFilter(2));
    demo(new XyWidthFilter2(0.7, 2));
    demo(new SbrFilter(15));
    demo(new ShiftFilter(0.7));
    demo(new EShiftFilter(0.8));
    demo(new SignalFilter(1000));
    demo(new SnrFilter(10));
    demo(new AnrFilter(11));
    demo(new PrecisionFilter(30));
    demo(new PrecisionFilter2(30));
    demo(new SnrHysteresisFilter(50, 1, 2, 1, 10, 20));
    demo(new PrecisionHysteresisFilter(2, 0, 1, 0, 20, 30));
    demo(new TraceFilter(0.5, 1));
    demo(new CoordinateFilter(15.5f, 234.5f, 80.99f, 133f));
    demo(new MultiFilter(30, 45f, 0.7, 1.5, 0.5, 0.6, 45, -10, 10));
    demo(new MultiFilter2(30, 45f, 0.7, 1.5, 0.5, 0.6, 45, -10, 10));
    demo(new MultiHysteresisFilter(2, 0, 1, 0, 20, 10, 40f, 20f, 0.8, 0.2, 1.2, 0.4, 0.3, 0.8, 20,
        30));
    demo(new MultiHysteresisFilter2(2, 0, 2, 1, 20, 10, 40f, 20f, 0.8, 0.2, 1.2, 0.4, 0.3, 0.8, 20,
        30));

    comment("Combined filters");
    IJ.log("");
    demo(new AndFilter(new SnrFilter(10), new WidthFilter(2)));
    demo(new OrFilter(new SnrFilter(10), new PrecisionFilter(30)));
    demo(new OrFilter(new AndFilter(new SnrFilter(10), new PrecisionFilter(30)),
        new TraceFilter(0.5, 1)));
  }

  private static void demo(Filter filter) {
    comment(filter.getClass().getSimpleName() + ": " + filter.getDescription());
    IJ.log(filter.toXml());
    IJ.log("");
  }

  private static void comment(String text) {
    IJ.log(TextUtils.wrap("<!-- " + text + " -->", 80));
  }
}
