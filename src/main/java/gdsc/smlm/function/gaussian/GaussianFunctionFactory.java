package gdsc.smlm.function.gaussian;

import gdsc.smlm.function.gaussian.erf.SingleCircularErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFixedErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.MultiAstigmatismErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.MultiCircularErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.MultiFixedErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.MultiFreeCircularErfGaussian2DFunction;
import gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Creates the appropriate Gaussian function
 */
public class GaussianFunctionFactory
{
	public static final int FIT_BACKGROUND = 0x00000001;
	public static final int FIT_ANGLE = 0x00000002;
	public static final int FIT_X_WIDTH = 0x00000004;
	public static final int FIT_Y_WIDTH = 0x00000008;
	public static final int FIT_SIGNAL = 0x00000010;
	public static final int FIT_Z = 0x00000020;

	public static final int FIT_ELLIPTICAL = FIT_BACKGROUND | FIT_ANGLE | FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL;
	public static final int FIT_FREE_CIRCLE = FIT_BACKGROUND | FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL;
	public static final int FIT_CIRCLE = FIT_BACKGROUND | FIT_X_WIDTH | FIT_SIGNAL;
	public static final int FIT_FIXED = FIT_BACKGROUND | FIT_SIGNAL;

	public static final int FIT_NB_ELLIPTICAL = FIT_ANGLE | FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL;
	public static final int FIT_NB_FREE_CIRCLE = FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL;
	public static final int FIT_NB_CIRCLE = FIT_X_WIDTH | FIT_SIGNAL;
	public static final int FIT_NB_FIXED = FIT_SIGNAL;

	// Special options for fixed function to be able to fix the signal
	public static final int FIT_NS_FIXED = FIT_BACKGROUND;
	public static final int FIT_NS_NB_FIXED = 0;

	// Flags for ERF Gaussian functions
	public static final int FIT_ERF = 0x00000100;
	public static final int FIT_ERF_FREE_CIRCLE = FIT_FREE_CIRCLE | FIT_ERF;
	public static final int FIT_ERF_CIRCLE = FIT_CIRCLE | FIT_ERF;
	public static final int FIT_ERF_FIXED = FIT_FIXED | FIT_ERF;
	public static final int FIT_ERF_ASTIGMATISM = FIT_BACKGROUND | FIT_Z | FIT_SIGNAL | FIT_ERF;

	/**
	 * Create the correct 2D Gaussian function for the specified parameters.
	 *
	 * @param nPeaks
	 *            The number of peaks (N)
	 * @param maxx
	 *            The maximum X-dimension
	 * @param maxy
	 *            The maximum Y-dimension
	 * @param flags
	 *            Enable all the parameters that should evaluate gradient
	 * @param zModel
	 *            the z model
	 * @return The function
	 */
	public static Gaussian2DFunction create2D(int nPeaks, int maxx, int maxy, int flags, AstimatismZModel zModel)
	{
		if ((flags & FIT_ERF) == FIT_ERF)
		{
			if (nPeaks == 1)
			{
				if ((flags & FIT_BACKGROUND) == FIT_BACKGROUND)
				{
					// Independent X/Y width
					if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
						return new SingleFreeCircularErfGaussian2DFunction(maxx, maxy);
					// Combined X/Y width
					if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
						return new SingleCircularErfGaussian2DFunction(maxx, maxy);
					// Z-depth function
					if ((flags & FIT_Z) == FIT_Z)
						return new SingleAstigmatismErfGaussian2DFunction(maxx, maxy, zModel);
					// Fixed width
					if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
						return new SingleFixedErfGaussian2DFunction(maxx, maxy);
				}
			}
			else
			{
				if ((flags & FIT_BACKGROUND) == FIT_BACKGROUND)
				{
					// Independent X/Y width
					if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
						return new MultiFreeCircularErfGaussian2DFunction(nPeaks, maxx, maxy);
					// Combined X/Y width
					if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
						return new MultiCircularErfGaussian2DFunction(nPeaks, maxx, maxy);
					// Z-depth function
					if ((flags & FIT_Z) == FIT_Z)
						return new MultiAstigmatismErfGaussian2DFunction(nPeaks, maxx, maxy, zModel);
					// Fixed width
					if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
						return new MultiFixedErfGaussian2DFunction(nPeaks, maxx, maxy);
				}
			}
		}

		if (nPeaks == 1)
		{
			if ((flags & FIT_BACKGROUND) == FIT_BACKGROUND)
			{
				if ((flags & FIT_ANGLE) == FIT_ANGLE)
					return new SingleEllipticalGaussian2DFunction(maxx, maxy);
				if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
					return new SingleFreeCircularGaussian2DFunction(maxx, maxy);
				if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
					return new SingleCircularGaussian2DFunction(maxx, maxy);

				// Fixed function
				if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
					return new SingleFixedGaussian2DFunction(maxx, maxy);

				return new SingleNSFixedGaussian2DFunction(maxx, maxy);
			}

			if ((flags & FIT_ANGLE) == FIT_ANGLE)
				return new SingleNBEllipticalGaussian2DFunction(maxx, maxy);
			if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
				return new SingleNBFreeCircularGaussian2DFunction(maxx, maxy);
			if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
				return new SingleNBCircularGaussian2DFunction(maxx, maxy);

			// Fixed function
			if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
				return new SingleNBFixedGaussian2DFunction(maxx, maxy);

			return new SingleNSNBFixedGaussian2DFunction(maxx, maxy);
		}
		else
		{
			if ((flags & FIT_BACKGROUND) == FIT_BACKGROUND)
			{
				if ((flags & FIT_ANGLE) == FIT_ANGLE)
					return new EllipticalGaussian2DFunction(nPeaks, maxx, maxy);
				if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
					return new FreeCircularGaussian2DFunction(nPeaks, maxx, maxy);
				if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
					return new CircularGaussian2DFunction(nPeaks, maxx, maxy);

				// Fixed function
				if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
					return new FixedGaussian2DFunction(nPeaks, maxx, maxy);

				return new NSFixedGaussian2DFunction(nPeaks, maxx, maxy);
			}

			if ((flags & FIT_ANGLE) == FIT_ANGLE)
				return new NBEllipticalGaussian2DFunction(nPeaks, maxx, maxy);
			if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
				return new NBFreeCircularGaussian2DFunction(nPeaks, maxx, maxy);
			if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
				return new NBCircularGaussian2DFunction(nPeaks, maxx, maxy);

			// Fixed function
			if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
				return new NBFixedGaussian2DFunction(nPeaks, maxx, maxy);

			return new NSNBFixedGaussian2DFunction(nPeaks, maxx, maxy);
		}
	}
}
