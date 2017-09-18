package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import org.apache.commons.math3.random.RandomDataGenerator;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.core.ij.Utils;
import gdsc.core.utils.Maths;
import gdsc.smlm.data.config.GUIProtos.NucleusMaskSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.ImageProcessor;

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

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (!showDialog())
			return;

		createMask();
	}

	double diameter, nmPerPixel, nmPerSlice;

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
		catch (IllegalArgumentException e)
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
			catch (IllegalArgumentException e)
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

		int inc = 2 * radius + 1;
		int incz = 2 * radiusz + 1;
		int maxx = settings.getFieldWidth();
		int maxy = maxx;
		int ditherHeight = (settings.getYDither() > 0) ? (int) Math.ceil(settings.getYDither() * 1000 / nmPerPixel) : 0;
		int ditherDepth = (settings.getZDither() > 0) ? (int) Math.ceil(settings.getZDither() * 1000 / nmPerSlice) : 0;
		int maxz = ditherDepth + incz;
		ImageStack stack = new ImageStack(maxx, maxy, maxz);
		byte[] mask = new byte[maxx * maxy];
		for (int z = 0; z < maxz; z++)
		{
			mask = (z == 0) ? mask : mask.clone();
			stack.setPixels(mask, z + 1);
		}

		if (settings.getMode() == 0)
		{
			ImageStack stack2 = createEllipsoid(inc, inc, incz);

			//Utils.display(TITLE + " nucleus", stack2, 0);

			// Dither 
			int cx = radius;
			int lowerz = (maxz - ditherDepth) / 2;
			int upperz = (maxz + ditherDepth) / 2;
			RandomDataGenerator r = new RandomDataGenerator();

			while (cx < maxx)
			{
				int xloc = cx - radius;
				int cy = radius + r.nextInt(0, ditherHeight);
				while (cy < maxy)
				{
					int yloc = cy - radius;
					int offset = r.nextInt(lowerz, upperz) - radiusz;
					for (int slice = 1; slice <= stack2.getSize(); slice++)
					{
						int i = slice + offset;
						if (i < 1 || i > maxz)
							continue;
						ImageProcessor ip = stack.getProcessor(i);
						ImageProcessor ip2 = stack2.getProcessor(slice);
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

			NonBlockingExtendedGenericDialog gd = new NonBlockingExtendedGenericDialog(TITLE);
			gd.addHelp(About.HELP_URL);

			gd.addMessage("Click the image to add a sphere");
			gd.addNumericField("Diameter", diameter, 2, 6, "um");
			gd.addDialogListener(this);
			gd.hideCancelButton();
			gd.setOKLabel("Close");
			gd.showDialog();

			imp.getWindow().removeMouseListener(this);

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
		Calibration cal = new Calibration();
		cal.setUnit("um");
		cal.pixelWidth = cal.pixelHeight = nmPerPixel / 1000;
		cal.pixelDepth = nmPerSlice / 1000;
		imp.setCalibration(cal);
	}

	/**
	 * Create an ellipsoid using the given dimensions.
	 *
	 * @param width
	 *            the width
	 * @param height
	 *            the height
	 * @param depth
	 *            the depth
	 * @return An ellipsoid
	 */
	public static ImageStack createEllipsoid(double widthX, double heightX, double depthX)
	{
		// Precompute squares for the distances
		double[] sx = getSquareDistances(widthX, depthX);
		double[] sy = getSquareDistances(heightX, depthX);
		double[] sz = getSquareDistances(depthX, depthX);

		int maxx = sx.length;
		int maxy = sy.length;
		int maxz = sz.length;

		ImageStack stack = new ImageStack(maxx, maxy, maxz);
		final byte on = (byte) 255;
		// Squared distances are relative to the radius of the depth
		final double r2 = Maths.pow2(depthX * 0.5);
		for (int iz = 0; iz < maxz; iz++)
		{
			byte[] mask = new byte[maxx * maxy];
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
	 * @param norm
	 *            the size of the reference axis
	 * @return the scaled square distances
	 */
	private static double[] getSquareDistances(double size, double reference)
	{
		int n = (int) Math.ceil(size);
		double centre = n * 0.5;
		double[] s = new double[n];
		double scale = reference / size;
		for (int i = 0; i < n; i++)
		{
			s[i] = Maths.pow2((centre - (i + 0.5)) * scale);
		}
		return s;
	}

	public void mouseClicked(MouseEvent e)
	{
		if (e == null)
			return;
		ImageCanvas ic = imp.getCanvas();
		int cx = ic.offScreenX(e.getX());
		int cy = ic.offScreenY(e.getY());
		final int radius = (int) Math.ceil(diameter * 500 / nmPerPixel);
		final int radiusz = (int) Math.ceil(diameter * 500 / nmPerSlice);
		if (sphere == null)
		{
			int inc = 2 * radius + 1;
			int incz = 2 * radiusz + 1;
			sphere = createEllipsoid(inc, inc, incz);
		}
		int xloc = cx - radius;
		int yloc = cy - radius;
		int offset = imp.getCurrentSlice() - radiusz;
		ImageStack stack = imp.getImageStack();
		for (int slice = 1; slice <= sphere.getSize(); slice++)
		{
			int i = slice + offset;
			if (i < 1 || i > stack.getSize())
				continue;
			ImageProcessor ip = stack.getProcessor(i);
			ImageProcessor ip2 = sphere.getProcessor(slice);
			ip.copyBits(ip2, xloc, yloc, Blitter.MAX);
		}

	}

	public void mousePressed(MouseEvent e)
	{
	}

	public void mouseReleased(MouseEvent e)
	{
	}

	public void mouseEntered(MouseEvent e)
	{
	}

	public void mouseExited(MouseEvent e)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		double old = diameter;
		diameter = gd.getNextNumber();
		if (diameter != old)
			sphere = null;
		return !gd.invalidNumber();
	}
}
