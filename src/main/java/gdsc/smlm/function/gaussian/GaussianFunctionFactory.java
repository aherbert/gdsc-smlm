package gdsc.smlm.function.gaussian;

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
	public static final int FIT_BACKGROUND = 1;
	public static final int FIT_ANGLE = 2;
	public static final int FIT_X_WIDTH = 4;
	public static final int FIT_Y_WIDTH = 8;
	public static final int FIT_SIGNAL = 16;

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
	
	/**
	 * Create the correct 2D Gaussian function for the specified parameters
	 * 
	 * @param nPeaks
	 *            The number of peaks (N)
	 * @param maxx
	 *            The maximum X-dimension
	 * @param flags
	 *            Enable all the parameters that should evaluate gradient
	 * @return The function
	 */
	public static Gaussian2DFunction create2D(int nPeaks, int maxx, int flags)
	{
		if (nPeaks == 1)
		{
			if ((flags & FIT_BACKGROUND) == FIT_BACKGROUND)
			{
				if ((flags & FIT_ANGLE) == FIT_ANGLE)
					return new SingleEllipticalGaussian2DFunction(maxx);
				if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
					return new SingleFreeCircularGaussian2DFunction(maxx);
				if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
					return new SingleCircularGaussian2DFunction(maxx);
				
				// Fixed function
				if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
					return new SingleFixedGaussian2DFunction(maxx);

				return new SingleNSFixedGaussian2DFunction(maxx);
			}

			if ((flags & FIT_ANGLE) == FIT_ANGLE)
				return new SingleNBEllipticalGaussian2DFunction(maxx);
			if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
				return new SingleNBFreeCircularGaussian2DFunction(maxx);
			if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
				return new SingleNBCircularGaussian2DFunction(maxx);

			// Fixed function
			if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
				return new SingleNBFixedGaussian2DFunction(maxx);

			return new SingleNSNBFixedGaussian2DFunction(maxx);
		}
		else
		{
			if ((flags & FIT_BACKGROUND) == FIT_BACKGROUND)
			{
				if ((flags & FIT_ANGLE) == FIT_ANGLE)
					return new EllipticalGaussian2DFunction(nPeaks, maxx);
				if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
					return new FreeCircularGaussian2DFunction(nPeaks, maxx);
				if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
					return new CircularGaussian2DFunction(nPeaks, maxx);

				// Fixed function
				if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
					return new FixedGaussian2DFunction(nPeaks, maxx);

				return new NSFixedGaussian2DFunction(nPeaks, maxx);
			}

			if ((flags & FIT_ANGLE) == FIT_ANGLE)
				return new NBEllipticalGaussian2DFunction(nPeaks, maxx);
			if ((flags & FIT_Y_WIDTH) == FIT_Y_WIDTH)
				return new NBFreeCircularGaussian2DFunction(nPeaks, maxx);
			if ((flags & FIT_X_WIDTH) == FIT_X_WIDTH)
				return new NBCircularGaussian2DFunction(nPeaks, maxx);

			// Fixed function
			if ((flags & FIT_SIGNAL) == FIT_SIGNAL)
				return new NBFixedGaussian2DFunction(nPeaks, maxx);

			return new NSNBFixedGaussian2DFunction(nPeaks, maxx);
		}
	}
}
