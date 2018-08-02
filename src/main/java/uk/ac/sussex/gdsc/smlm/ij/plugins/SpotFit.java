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

import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.util.Arrays;
import java.util.regex.Pattern;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.Toolbar;
import ij.plugin.PlugIn;
import ij.plugin.tool.PlugInTool;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
import ij.text.TextPanel;
import ij.text.TextWindow;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.SpotFitSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtosHelper;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.SimplePeakResultValidationData;
import uk.ac.sussex.gdsc.smlm.filters.BlockMeanFilter;
import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.utils.ImageConverter;

/**
 * Plugin to fit spots interactively.
 */
public class SpotFit implements PlugIn
{
    private static final String TITLE = "Spot Fit";

    /**
     * All the work for this plugin is done with the plugin tool.
     * It handles mouse click events from an image.
     */
    private static class SpotFitPluginTool extends PlugInTool
    {
        private static TextWindow resultsWindow = null;

        private final SpotFitSettings.Builder settings;
        private boolean active = true;
        private boolean logging = false;

        private final BlockMeanFilter f = new BlockMeanFilter();
        private final FitConfiguration config;
        private SimplePeakResultValidationData validationData;
        private final Gaussian2DFitter gf;

        private final double[] lower, upper;

        final static Pattern pattern = Pattern.compile("\t");

        private SpotFitPluginTool()
        {
            settings = SettingsManager.readSpotFitSettings(0).toBuilder();
            config = createFitConfiguration();
            gf = new Gaussian2DFitter(config);

            // Support bounded fit on the coordinates
            final int n = Gaussian2DFunction.PARAMETERS_PER_PEAK + 1;
            lower = new double[n];
            upper = new double[n];
            for (int i = 0; i < n; i++)
            {
                lower[i] = Double.NEGATIVE_INFINITY;
                upper[i] = Double.POSITIVE_INFINITY;
            }
        }

        @Override
        public String getToolName()
        {
            return TITLE + " Tool";
        }

        @Override
        public String getToolIcon()
        {
            // A blue dot with red circle outside (a bull's eye target)
            return "C00f0o4466Cf00O11bb";
        }

        @Override
        public void showOptionsDialog()
        {
            final GenericDialog gd = new GenericDialog(TITLE + " Tool Options");
            gd.addMessage(
            //@formatter:off
				TextUtils.wrap(
				"Click on an image and fit a spot in a selected channel. " +
				"The maxima within a search range is used to centre the " +
				"fit window for Gaussian 2D fitting.", 80));
				//@formatter:on
            gd.addNumericField("Channel", settings.getChannel(), 0);
            gd.addSlider("Search_range", 1, 10, settings.getSearchRadius());
            gd.addSlider("Fit_radius", 3, 10, settings.getFitRadius());
            gd.addNumericField("SNR_threshold", settings.getSnrThreshold(), 0);
            gd.addCheckbox("Show_fit_ROI", settings.getShowFitRoi());
            gd.addCheckbox("Show_overlay", settings.getShowOverlay());
            gd.addCheckbox("Attach_to_slice", settings.getAttachToSlice());
            gd.addCheckbox("Log_progress", settings.getLogProgress());

            gd.showDialog();
            if (gd.wasCanceled())
                return;
            synchronized (this)
            {
                settings.setChannel((int) gd.getNextNumber());
                settings.setSearchRadius((int) gd.getNextNumber());
                settings.setFitRadius((int) gd.getNextNumber());
                settings.setSnrThreshold(gd.getNextNumber());
                settings.setShowFitRoi(gd.getNextBoolean());
                settings.setShowOverlay(gd.getNextBoolean());
                settings.setAttachToSlice(gd.getNextBoolean());
                logging = gd.getNextBoolean();
                settings.setLogProgress(logging);
                // Only active if the settings are valid
                active = (settings.getFitRadius() > 1);
                if (!active)
                    IJ.error(TITLE, "Settings are invalid");
                SettingsManager.writeSettings(settings, 0);
            }
        }

        @Override
        public void mouseClicked(ImagePlus imp, MouseEvent e)
        {
            if (!active)
                return;
            final int c = imp.getNChannels();
            if (c < settings.getChannel())
            {
                // Always warn if the channel is incorrect for the image
                //if (logging)
                Utils.log(TITLE + ": Image %s does not contain channel %d", imp.getTitle(), settings.getChannel());
                return;
            }

            // Mark this event as handled
            e.consume();

            // TODO - More control over fitting.

            // Ensure rapid mouse click / new options does not break things
            synchronized (this)
            {
                final ImageCanvas ic = imp.getCanvas();
                int x = ic.offScreenX(e.getX());
                int y = ic.offScreenY(e.getY());

                if (logging)
                    Utils.log("Clicked %d,%d", x, y);

                // Get the data
                final int channel = settings.getChannel();
                final int slice = imp.getZ();
                final int frame = imp.getFrame();

                final int stackIndex = imp.getStackIndex(channel, slice, frame);

                final ImageExtractor ie = new ImageExtractor(null, imp.getWidth(), imp.getHeight());

                if (isRemoveEvent(e))
                {
                    removeSpots(imp, channel, slice, frame, x, y, ie);
                    return;
                }

                final ImageStack stack = imp.getImageStack();
                final ImageProcessor ip = stack.getProcessor(stackIndex);

                // Search for the maxima using the search radius
                final int index = findMaxima(ip, ie, x, y);

                // Fit the maxima
                x = index % ip.getWidth();
                y = index / ip.getWidth();
                if (logging)
                    Utils.log("Fitting %d,%d", x, y);
                final Rectangle bounds = ie.getBoxRegionBounds(x, y, settings.getFitRadius());
                if (settings.getShowFitRoi())
                    imp.setRoi(bounds);
                final FitResult fitResult = fitMaxima(ip, bounds, x, y);

                if (logging)
                    Utils.log("Fit estimate = %s", Arrays.toString(fitResult.getInitialParameters()));
                if (logging)
                {
                    String msg = "Fit status = " + fitResult.getStatus();
                    final Object data = fitResult.getStatusData();
                    if (data != null)
                        msg += " : " + SimpleArrayUtils.toString(data);
                    Utils.log(msg);
                }

                if (fitResult.getStatus() != FitStatus.OK)
                    // Q. Do something?
                    return;

                // Add result
                createResultsWindow();
                addResult(imp, channel, slice, frame, bounds, fitResult);

                if (settings.getShowOverlay())
                    addOverlay(imp, channel, slice, frame, fitResult);
            }
        }

        private static boolean isRemoveEvent(MouseEvent e)
        {
            return e.isAltDown() || e.isShiftDown() || e.isControlDown();
        }

        private int findMaxima(ImageProcessor ip, ImageExtractor ie, int x, int y)
        {
            final int n = settings.getSearchRadius();
            if (n <= 0)
                // No search
                return y * ip.getWidth() + x;
            // Get a region
            final Rectangle bounds = ie.getBoxRegionBounds(x, y, n);
            final float[] data = new ImageConverter().getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds,
                    null);
            // Smooth
            f.blockFilter(data, bounds.width, bounds.height, 1);
            int index = 0;
            // Find maxima
            for (int i = 1; i < data.length; i++)
                if (data[index] < data[i])
                    index = i;
            // Convert back
            x = bounds.x + index % bounds.width;
            y = bounds.y + index / bounds.width;
            return ip.getWidth() * y + x;
        }

        private FitResult fitMaxima(ImageProcessor ip, Rectangle bounds, int x, int y)
        {
            config.setInitialPeakStdDev(0);

            // Get a region
            final double[] data = new ImageConverter().getDoubleData(ip.getPixels(), ip.getWidth(), ip.getHeight(),
                    bounds, null);

            setupPeakFiltering(config, bounds, ip);

            //Utils.display(TITLE + " fit data", data, bounds.width, bounds.height, 0);

            // Find the index of the maxima
            final int ox = x - bounds.x;
            final int oy = y - bounds.y;
            final int index = oy * bounds.width + ox;

            // Limit the range for the XY position
            final int range = settings.getFitRadius() / 2;
            lower[Gaussian2DFunction.X_POSITION] = ox - range;
            lower[Gaussian2DFunction.Y_POSITION] = oy - range;
            upper[Gaussian2DFunction.X_POSITION] = ox + range;
            upper[Gaussian2DFunction.Y_POSITION] = oy + range;
            gf.setBounds(lower, upper);

            // Leave to the fitter to estimate background, width and height
            final FitResult fitResult = gf.fit(data, bounds.width, bounds.height, new int[] { index });
            if (fitResult.getStatus() == FitStatus.OK)
            {
                final double[] params = fitResult.getParameters();
                // Add the pixel offset
                params[Gaussian2DFunction.X_POSITION] += bounds.x + 0.5;
                params[Gaussian2DFunction.Y_POSITION] += bounds.y + 0.5;
            }
            return fitResult;
        }

        private static FitConfiguration createFitConfiguration()
        {
            final FitConfiguration config = new FitConfiguration();
            config.setFitSolver(FitSolver.LVM_LSE);
            config.setPSF(PSFProtosHelper.getDefaultPSF(PSFType.ONE_AXIS_GAUSSIAN_2D));
            //config.setMaxIterations(getMaxIterations());
            //config.setRelativeThreshold(relativeThreshold);
            //config.setAbsoluteThreshold(absoluteThreshold);
            config.setInitialPeakStdDev(0);
            config.setComputeDeviations(false);

            // Set-up peak filtering
            config.setDisableSimpleFilter(false);

            //if (isLogProgress())
            //{
            //	config.setLog(new IJLogger());
            //}

            config.setBackgroundFitting(true);

            return config;
        }

        protected void setupPeakFiltering(FitConfiguration config, Rectangle bounds, ImageProcessor ip)
        {
            config.setCoordinateShift(settings.getSearchRadius());
            config.setSignalStrength(settings.getSnrThreshold());
            config.setMinWidthFactor(0.5);
            //config.setWidthFactor(3);

            validationData = new SimplePeakResultValidationData(
                    new GaussianFunctionFactory(config.getFunctionFlags(), config.getAstigmatismZModel()), bounds.x,
                    bounds.y, ip.getPixels(), ip.getWidth(), ip.getHeight());
            config.setPeakResultValidationData(validationData);
        }

        /**
         * Create the result window (if it is not available)
         */
        private static void createResultsWindow()
        {
            if (resultsWindow == null || !resultsWindow.isShowing())
                resultsWindow = new TextWindow(TITLE + " Results", createHeader(), "", 700, 300);
        }

        private static String createHeader()
        {
            final StringBuilder sb = new StringBuilder();
            sb.append("Image\t");
            sb.append("C\t");
            sb.append("Z\t");
            sb.append("T\t");
            sb.append("Region\t");
            sb.append("Background\t");
            sb.append("Intensity\t");
            sb.append("X\t");
            sb.append("Y\t");
            sb.append("S\t");
            sb.append("Noise\t");
            sb.append("SNR");
            return sb.toString();
        }

        private void addResult(ImagePlus imp, int channel, int slice, int frame, Rectangle bounds, FitResult fitResult)
        {
            createResultsWindow();

            final StringBuilder sb = new StringBuilder();
            sb.append(imp.getTitle()).append('\t');
            sb.append(channel).append('\t');
            sb.append(slice).append('\t');
            sb.append(frame).append('\t');
            sb.append(bounds.x).append(',');
            sb.append(bounds.y).append(' ');
            sb.append(bounds.width).append('x');
            sb.append(bounds.height);
            final double[] params = fitResult.getParameters();
            sb.append('\t').append(Utils.rounded(params[Gaussian2DFunction.BACKGROUND]));
            final double signal = params[Gaussian2DFunction.SIGNAL];
            sb.append('\t').append(Utils.rounded(signal));
            sb.append('\t').append(Utils.rounded(params[Gaussian2DFunction.X_POSITION]));
            sb.append('\t').append(Utils.rounded(params[Gaussian2DFunction.Y_POSITION]));
            final double xsd = params[Gaussian2DFunction.X_SD];
            sb.append('\t').append(Utils.rounded(xsd));
            final double noise = validationData.getNoise();
            final double snr = Gaussian2DPeakResultHelper.getMeanSignalUsingP05(signal, xsd, xsd) / noise;
            sb.append('\t').append(Utils.rounded(noise));
            sb.append('\t').append(Utils.rounded(snr));
            resultsWindow.append(sb.toString());
        }

        private void addOverlay(ImagePlus imp, int channel, int slice, int frame, FitResult fitResult)
        {
            final double[] params = fitResult.getParameters();
            Overlay o = imp.getOverlay();
            if (o == null)
                o = new Overlay();
            final PointRoi roi = new PointRoi(params[Gaussian2DFunction.X_POSITION],
                    params[Gaussian2DFunction.Y_POSITION]);
            roi.setPointType(3);
            if (imp.isDisplayedHyperStack())
                roi.setPosition(channel, (settings.getAttachToSlice()) ? slice : 0, frame);
            else if (settings.getAttachToSlice())
                roi.setPosition(imp.getCurrentSlice());
            o.add(roi);
            imp.setOverlay(o);
        }

        private void removeSpots(ImagePlus imp, int channel, int slice, int frame, int x, int y, ImageExtractor ie)
        {
            // Get a region to search for spots
            final Rectangle bounds = ie.getBoxRegionBounds(x, y, Math.max(0, settings.getSearchRadius()));
            if (bounds.width == 0 || bounds.height == 0)
                return;

            final boolean isDisplayedHyperStack = imp.isDisplayedHyperStack();

            // Remove all the overlay components
            Overlay o = imp.getOverlay();
            if (o != null)
            {
                final Roi[] rois = o.toArray();
                final int size = o.size();
                o = new Overlay();
                for (int i = 0; i < rois.length; i++)
                {
                    if (rois[i] instanceof PointRoi)
                    {
                        final PointRoi roi = (PointRoi) rois[i];
                        boolean boundsCheck = true;
                        if (isDisplayedHyperStack)
                            // Must be on the same channel/slice/frame
                            boundsCheck = roi.getCPosition() == channel && roi.getTPosition() == frame &&
                                    (roi.getZPosition() == 0 || roi.getZPosition() == slice);
                        if (boundsCheck)
                        {
                            final FloatPolygon poly = roi.getFloatPolygon();
                            if (bounds.contains(poly.xpoints[0], poly.ypoints[0]))
                                continue;
                        }
                    }
                    o.add(rois[i]);
                }
                if (o.size() != size)
                    if (o.size() == 0)
                        imp.setOverlay(null);
                    else
                        imp.setOverlay(o);
            }

            if (resultsWindow != null && resultsWindow.isShowing())
            {
                final TextPanel tp = resultsWindow.getTextPanel();
                final String title = imp.getTitle();
                for (int i = 0; i < tp.getLineCount(); i++)
                {
                    final String line = tp.getLine(i);
                    // Check the image name
                    final int startIndex = line.indexOf('\t');
                    if (startIndex == -1)
                        continue;
                    if (!title.equals(line.substring(0, startIndex)))
                        continue;

                    final String[] fields = pattern.split(line, 0);

                    try
                    {
                        if (isDisplayedHyperStack)
                        {
                            final int c = Integer.parseInt(fields[1]);
                            // Ignore z as the user may be on the wrong slice but can still see
                            // the overlay if it is not tied to the slice position
                            //int z = Integer.parseInt(fields[2]);
                            final int t = Integer.parseInt(fields[3]);
                            if (c != channel || t != frame)
                                continue;
                        }
                        final float xp = Float.parseFloat(fields[7]);
                        final float yp = Float.parseFloat(fields[8]);
                        if (bounds.contains(xp, yp))
                        {
                            tp.setSelection(i, i);
                            tp.clearSelection();
                            // Since i will be incremented for the next line,
                            // decrement to check the current line again.
                            i--;
                        }
                    }
                    catch (final NumberFormatException ex)
                    {
                        // Ignore
                    }
                }
            }
        }
    }

    private static SpotFitPluginTool toolInstance = null;

    /**
     * Initialise the spot fit tool. This is to allow support for calling within macro toolsets.
     */
    public static void addPluginTool()
    {
        if (toolInstance == null)
            toolInstance = new SpotFitPluginTool();

        // Add the tool
        Toolbar.addPlugInTool(toolInstance);
        IJ.showStatus("Added " + TITLE + " Tool");
    }

    @Override
    public void run(String arg)
    {
        SMLMUsageTracker.recordPlugin(this.getClass(), arg);

        addPluginTool();

        // Fiji restores the toolbar from the last session.
        // Do not show the options if this is happening.
        final ImageJ ij = IJ.getInstance();
        if (ij == null || !ij.isVisible())
            return;

        toolInstance.showOptionsDialog();
    }
}
