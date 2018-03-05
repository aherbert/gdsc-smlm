package gdsc.smlm.ij.ij3d;

import java.util.Arrays;

import org.scijava.java3d.GeometryArray;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import customnode.CustomMesh;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.ImageProcessor;
import vib.InterpolatedImage;

/**
 * Provide helper functionality for dealing with the CustomMesh class.
 */
public class CustomMeshHelper
{
	CustomMesh mesh;

	/**
	 * Instantiates a new custom mesh helper.
	 *
	 * @param mesh
	 *            the mesh
	 */
	public CustomMeshHelper(CustomMesh mesh)
	{
		this.mesh = mesh;
	}

	/**
	 * Load surface colors from a 2D image. This only works if the xy coordinates from the mesh can be mapped directly
	 * onto the image pixel coordinates by dividing by the input image pixel width / height (i.e. the input image must
	 * be calibrated appropriately).
	 *
	 * @param imp
	 *            the imp
	 */
	public void loadSurfaceColorsFromImage2D(ImagePlus imp)
	{
		final GeometryArray ga = (GeometryArray) mesh.getGeometry();
		if (ga == null)
			return;

		ImageProcessor ip = imp.getProcessor();
		if (imp.getType() != ImagePlus.COLOR_RGB)
		{
			ip = ip.duplicate().convertToRGB();
		}
		// Create a stack 
		ImageStack stack = new ImageStack(ip.getWidth(), ip.getHeight());
		stack.addSlice(ip);
		stack.addSlice(ip);
		final Calibration cal = imp.getCalibration();
		imp = new ImagePlus(null, stack);
		
		final InterpolatedImage ii = new InterpolatedImage(imp);

		final int N = ga.getValidVertexCount();
		final Color3f[] colors = new Color3f[N];
		final double pw = cal.pixelWidth;
		final double ph = cal.pixelHeight;
		final Point3f coord = new Point3f();
		for (int i = 0; i < N; i++)
		{
			ga.getCoordinate(i, coord);
			// Ignore the z-coordinate
			final int v = (int) Math.round(ii.interpol.get(coord.x / pw, coord.y / ph, 0));
			colors[i] = new Color3f(((v & 0xff0000) >> 16) / 255f, ((v & 0xff00) >> 8) / 255f, (v & 0xff) / 255f);
		}
		mesh.setColor(Arrays.asList(colors));
	}

}
