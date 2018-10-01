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

import java.awt.AWTEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import org.apache.commons.math3.random.RandomDataGenerator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.ImageProcessor;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.NonBlockingExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.NucleusMaskSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;

/**
 * This plugin creates a mask image stack using an XY and XZ mask image
 */
public class NucleusMask implements PlugIn, MouseListener, DialogListener
{
    private static final String TITLE = "Nucleus Mask";

    private static final String[] MODE = { "Random", "User Input" };

    private NucleusMaskSettings.Builder settings;

    private ImagePlus imp = null;
    private ImageStack sphere = null;

    /** {@inheritDoc} */
    @Override
    public void run(String arg)
    {
        SMLMUsageTracker.recordPlugin(this.getClass(), arg);

        if (!showDialog())
            return;

        createMask();
    }

    private double diameter, nmPerPixel, nmPerSlice;

    private boolean showDialog()
    {
        ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
        gd.addHelp(About.HELP_URL);

        gd.addMessage("Create a mask stack using spheres");

        settings = SettingsManager.readNucleusMaskSettings(0).toBuilder();
        gd.addChoice("Mode", MODE, settings.getMode());
        gd.addNumericField("Field_width", settings.getFieldWidth(), 2, 6, "px");
        gd.addNumericField("Pixel_width", settings.getNmPerPixel(), 2, 6, "nm");
        gd.addNumericField("Pixel_depth", settings.getNmPerSlice(), 2, 6, "nm");

        gd.showDialog();

        if (gd.wasCanceled())
            return false;

        settings.setMode(gd.getNextChoiceIndex());
        settings.setFieldWidth((int) gd.getNextNumber());
        settings.setNmPerPixel(gd.getNextNumber());
        settings.setNmPerSlice(gd.getNextNumber());

        // Check arguments
        try
        {
            Parameters.isAboveZero("Field width", settings.getFieldWidth());
            Parameters.isAboveZero("Pixel width", settings.getNmPerPixel());
            Parameters.isAboveZero("Pixel depth", settings.getNmPerSlice());
        }
        catch (final IllegalArgumentException e)
        {
            IJ.error(TITLE, e.getMessage());
            return false;
        }

        if (settings.getMode() == 0)
        {
            gd = new ExtendedGenericDialog(TITLE);
            gd.addHelp(About.HELP_URL);

            gd.addMessage("Create a mask stack using uniform random spheres");

            gd.addNumericField("y_dither", settings.getYDither(), 2, 6, "um");
            gd.addNumericField("z_dither", settings.getZDither(), 2, 6, "um");
            gd.addNumericField("Diameter", settings.getDiameter(), 2, 6, "um");

            gd.showDialog();

            if (gd.wasCanceled())
                return false;

            settings.setYDither(gd.getNextNumber());
            settings.setZDither(gd.getNextNumber());
            settings.setDiameter(gd.getNextNumber());

            try
            {
                Parameters.isAboveZero("Diameter", settings.getDiameter());
            }
            catch (final IllegalArgumentException e)
            {
                IJ.error(TITLE, e.getMessage());
                return false;
            }
        }

        SettingsManager.writeSettings(settings);

        return true;
    }

    private void createMask()
    {
        diameter = settings.getDiameter();
        nmPerPixel = settings.getNmPerPixel();
        nmPerSlice = settings.getNmPerSlice();

        // Create the dimensions using the scale.
        // Scale diameter in um to nm
        final int radius = (int) Math.ceil(diameter * 500 / nmPerPixel);
        final int radiusz = (int) Math.ceil(diameter * 500 / nmPerSlice);

        final int inc = 2 * radius + 1;
        final int incz = 2 * radiusz + 1;
        final int maxx = settings.getFieldWidth();
        final int maxy = maxx;
        final int ditherHeight = (settings.getYDither() > 0)
                ? (int) Math.ceil(settings.getYDither() * 1000 / nmPerPixel)
                : 0;
        final int ditherDepth = (settings.getZDither() > 0) ? (int) Math.ceil(settings.getZDither() * 1000 / nmPerSlice)
                : 0;
        final int maxz = ditherDepth + incz;
        final ImageStack stack = new ImageStack(maxx, maxy, maxz);
        byte[] mask = new byte[maxx * maxy];
        for (int z = 0; z < maxz; z++)
        {
            mask = (z == 0) ? mask : mask.clone();
            stack.setPixels(mask, z + 1);
        }

        if (settings.getMode() == 0)
        {
            final ImageStack stack2 = createEllipsoid(inc, inc, incz);

            //Utils.display(TITLE + " nucleus", stack2, 0);

            // Dither
            int cx = radius;
            final int lowerz = (maxz - ditherDepth) / 2;
            final int upperz = (maxz + ditherDepth) / 2;
            final RandomDataGenerator r = new RandomDataGenerator();

            while (cx < maxx)
            {
                final int xloc = cx - radius;
                int cy = radius + r.nextInt(0, ditherHeight);
                while (cy < maxy)
                {
                    final int yloc = cy - radius;
                    final int offset = r.nextInt(lowerz, upperz) - radiusz;
                    for (int slice = 1; slice <= stack2.getSize(); slice++)
                    {
                        final int i = slice + offset;
                        if (i < 1 || i > maxz)
                            continue;
                        final ImageProcessor ip = stack.getProcessor(i);
                        final ImageProcessor ip2 = stack2.getProcessor(slice);
                        ip.copyBits(ip2, xloc, yloc, Blitter.MAX);
                    }
                    cy += inc + 1 + r.nextInt(0, ditherHeight);
                }
                cx += inc + 1;
            }
        }

        // The final image will have a scale added to it.
        imp = Utils.display(TITLE, stack);
        calibrate(imp);

        imp.setSlice(maxz / 2);

        if (settings.getMode() == 1)
        {
            // Allow mouse click to draw spheres
            imp.getCanvas().addMouseListener(this);

            final NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
            gd.addHelp(About.HELP_URL);

            gd.addMessage("Click the image to add a sphere");
            gd.addNumericField("Diameter", diameter, 2, 6, "um");
            gd.addDialogListener(this);
            gd.hideCancelButton();
            gd.setOKLabel("Close");
            gd.showDialog();

            imp.getCanvas().removeMouseListener(this);

            if (diameter != settings.getDiameter())
            {
                settings.setDiameter(diameter);
                SettingsManager.writeSettings(settings);
            }
        }
    }

    private void calibrate(ImagePlus imp)
    {
        // Calibrate
        final Calibration cal = new Calibration();
        cal.setUnit("um");
        cal.pixelWidth = cal.pixelHeight = nmPerPixel / 1000;
        cal.pixelDepth = nmPerSlice / 1000;
        imp.setCalibration(cal);
    }

    /**
     * Create an ellipsoid using the given dimensions.
     *
     * @param widthX
     *            the width X
     * @param heightX
     *            the height X
     * @param depthX
     *            the depth X
     * @return An ellipsoid
     */
    public static ImageStack createEllipsoid(double widthX, double heightX, double depthX)
    {
        // Precompute squares for the distances
        final double[] sx = getSquareDistances(widthX, depthX);
        final double[] sy = getSquareDistances(heightX, depthX);
        final double[] sz = getSquareDistances(depthX, depthX);

        final int maxx = sx.length;
        final int maxy = sy.length;
        final int maxz = sz.length;

        final ImageStack stack = new ImageStack(maxx, maxy, maxz);
        final byte on = (byte) 255;
        // Squared distances are relative to the radius of the depth
        final double r2 = Maths.pow2(depthX * 0.5);
        for (int iz = 0; iz < maxz; iz++)
        {
            final byte[] mask = new byte[maxx * maxy];
            final double z2 = sz[iz];
            for (int iy = 0, i = 0; iy < maxy; iy++)
            {
                final double y2z2 = sy[iy] + z2;
                for (int ix = 0; ix < maxx; ix++, i++)
                {
                    final double d2 = sx[ix] + y2z2;
                    if (d2 < r2)
                        mask[i] = on;
                }
            }
            stack.setPixels(mask, iz + 1);
        }

        return stack;
    }

    /**
     * Gets the scaled square distances.
     *
     * @param size
     *            the size of the axis
     * @param reference
     *            the size of the reference axis
     * @return the scaled square distances
     */
    private static double[] getSquareDistances(double size, double reference)
    {
        final int n = (int) Math.ceil(size);
        final double centre = n * 0.5;
        final double[] s = new double[n];
        final double scale = reference / size;
        for (int i = 0; i < n; i++)
            s[i] = Maths.pow2((centre - (i + 0.5)) * scale);
        return s;
    }

    @Override
    public void mouseClicked(MouseEvent e)
    {
        if (e == null)
            return;
        final ImageCanvas ic = imp.getCanvas();
        final int cx = ic.offScreenX(e.getX());
        final int cy = ic.offScreenY(e.getY());
        final int radius = (int) Math.ceil(diameter * 500 / nmPerPixel);
        final int radiusz = (int) Math.ceil(diameter * 500 / nmPerSlice);
        if (sphere == null)
        {
            final int inc = 2 * radius + 1;
            final int incz = 2 * radiusz + 1;
            sphere = createEllipsoid(inc, inc, incz);
        }
        final int xloc = cx - radius;
        final int yloc = cy - radius;
        final int offset = imp.getCurrentSlice() - radiusz;
        final ImageStack stack = imp.getImageStack();
        for (int slice = 1; slice <= sphere.getSize(); slice++)
        {
            final int i = slice + offset;
            if (i < 1 || i > stack.getSize())
                continue;
            final ImageProcessor ip = stack.getProcessor(i);
            final ImageProcessor ip2 = sphere.getProcessor(slice);
            ip.copyBits(ip2, xloc, yloc, Blitter.MAX);
        }
        imp.updateAndDraw();
    }

    @Override
    public void mousePressed(MouseEvent e)
    {
        // Ignore
    }

    @Override
    public void mouseReleased(MouseEvent e)
    {
        // Ignore
    }

    @Override
    public void mouseEntered(MouseEvent e)
    {
        // Ignore
    }

    @Override
    public void mouseExited(MouseEvent e)
    {
        // Ignore
    }

    /** {@inheritDoc} */
    @Override
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
    {
        final double old = diameter;
        diameter = gd.getNextNumber();
        if (diameter != old)
            sphere = null;
        return !gd.invalidNumber();
    }
}
