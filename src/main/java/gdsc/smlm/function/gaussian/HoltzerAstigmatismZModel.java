package gdsc.smlm.function.gaussian;

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

/**
 * Implements a astigmatism model of a 2D Gaussian function, where z-depth determines the x and y width.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note).
 * <p>
 * Ref: Holtzer, L., Meckel, T. & Schmidt, T. Nanometric three-dimensional tracking of individual quantum dots in cells.
 * Applied Physics Letters 90, 1–3 (2007).
 */
public class HoltzerAstigmatismZModel implements AstigmatismZModel
{
	public final double gamma, one_d2, Ax, Bx, Ay, By;

	/**
	 * Static constructor.
	 * <p>
	 * Note that a positive gamma puts the focal plane for the X-dimension above the z-centre (positive Z) and the focal
	 * plane for the Y-dimension below the z-centre (negative Z). If gamma is negative then the orientation of the focal 
	 * planes of X and Y are reversed.
	 *
	 * @param gamma
	 *            the gamma parameter (half the distance between the focal planes)
	 * @param d
	 *            the depth of focus
	 * @param Ax
	 *            Empirical constant A for the x-astigmatism of the PSF
	 * @param Bx
	 *            Empirical constant B for the x-astigmatism of the PSF
	 * @param Ay
	 *            Empirical constant A for the y-astigmatism of the PSF
	 * @param By
	 *            Empirical constant B for the y-astigmatism of the PSF
	 * @return the holtzer astimatism Z model
	 */
	public static HoltzerAstigmatismZModel create(double gamma, double d, double Ax, double Bx, double Ay, double By)
	{
		final double d2 = d * d;
		return new HoltzerAstigmatismZModel(gamma, 1.0 / d2, Ax, Bx, Ay, By);
	}

	/**
	 * Constructor.
	 * <p>
	 * Note that a positive gamma puts the focal plane for the X-dimension above the z-centre (positive Z) and the focal
	 * plane for the Y-dimension below the z-centre (negative Z). If gamma is negative then the orientation of the focal 
	 * planes of X and Y are reversed.
	 *
	 * @param gamma
	 *            the gamma parameter (half the distance between the focal planes)
	 * @param one_d2
	 *            one over the depth of focus squared (1/d^2)
	 * @param Ax
	 *            Empirical constant A for the x-astigmatism of the PSF
	 * @param Bx
	 *            Empirical constant B for the x-astigmatism of the PSF
	 * @param Ay
	 *            Empirical constant A for the y-astigmatism of the PSF
	 * @param By
	 *            Empirical constant B for the y-astigmatism of the PSF
	 */
	public HoltzerAstigmatismZModel(double gamma, double one_d2, double Ax, double Bx, double Ay, double By)
	{
		this.gamma = gamma;
		this.one_d2 = one_d2;
		this.Ax = Ax;
		this.Bx = Bx;
		this.Ay = Ay;
		this.By = By;
	}

	/**
	 * Gets the standard deviation, first and second derivatives for the z-depth.
	 *
	 * @param z
	 *            the z
	 * @param one_d2
	 *            one over the depth of focus squared (1/d^2)
	 * @param A
	 *            Empirical constant A for the astigmatism of the PSF
	 * @param B
	 *            Empirical constant B for the astigmatism of the PSF
	 * @param ds_dz
	 *            the first and second derivative of s given z
	 * @return the standard deviation
	 */
	public static double getS2(double z, double one_d2, double A, double B, double[] ds_dz)
	{
		final double z2 = z * z;
		final double z3 = z2 * z;
		final double z4 = z2 * z2;
		// Eq. 17a
		final double s = Math.sqrt(1 + one_d2 * (z2 + A * z3 + B * z4));
		// Eq. 19a
		ds_dz[0] = (one_d2 * (2 * z + A * 3 * z2 + B * 4 * z3)) / (2 * s);
		// Eq. 19b 
		ds_dz[1] = (one_d2 * (2 + A * 6 * z + B * 12 * z2)) / (2 * s) -
				pow2(one_d2 * (2 * z + A * 3 * z2 + B * 4 * z3)) / (4 * s * s * s);
		return s;
	}

	private static double pow2(final double d)
	{
		return d * d;
	}

	/**
	 * Gets the standard deviation and first derivative for the z-depth.
	 *
	 * @param z
	 *            the z
	 * @param one_d2
	 *            one over the depth of focus squared (1/d^2)
	 * @param A
	 *            Empirical constant A for the astigmatism of the PSF
	 * @param B
	 *            Empirical constant B for the astigmatism of the PSF
	 * @param ds_dz
	 *            the first derivative of s given z
	 * @return the standard deviation
	 */
	public static double getS1(double z, double one_d2, double A, double B, double[] ds_dz)
	{
		final double z2 = z * z;
		final double z3 = z2 * z;
		final double z4 = z2 * z2;
		// Eq. 17a
		final double s = Math.sqrt(1 + one_d2 * (z2 + A * z3 + B * z4));
		// Eq. 19a
		ds_dz[0] = (one_d2 * (2 * z + A * 3 * z2 + B * 4 * z3)) / (2 * s);
		return s;
	}

	/**
	 * Gets the standard deviation for the z-depth.
	 *
	 * @param z
	 *            the z
	 * @param one_d2
	 *            one over the depth of focus squared (1/d^2)
	 * @param A
	 *            Empirical constant A for the astigmatism of the PSF
	 * @param B
	 *            Empirical constant B for the astigmatism of the PSF
	 * @return the standard deviation
	 */
	public static double getS(double z, double one_d2, double A, double B)
	{
		final double z2 = z * z;
		final double z3 = z2 * z;
		final double z4 = z2 * z2;
		// Eq. 17a
		return Math.sqrt(1 + one_d2 * (z2 + A * z3 + B * z4));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.AstimatismZModel#getSx(double)
	 */
	public double getSx(double z)
	{
		return getS(z - gamma, one_d2, Ax, Bx);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.AstimatismZModel#getSx(double, double[])
	 */
	public double getSx(double z, double[] ds_dz)
	{
		return getS1(z - gamma, one_d2, Ax, Bx, ds_dz);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.AstimatismZModel#getSx2(double, double[])
	 */
	public double getSx2(double z, double[] ds_dz)
	{
		return getS2(z - gamma, one_d2, Ax, Bx, ds_dz);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.AstimatismZModel#getSy(double)
	 */
	public double getSy(double z)
	{
		return getS(z + gamma, one_d2, Ay, By);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.AstimatismZModel#getSy(double, double[])
	 */
	public double getSy(double z, double[] ds_dz)
	{
		return getS1(z + gamma, one_d2, Ay, By, ds_dz);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.gaussian.AstimatismZModel#getSy2(double, double[])
	 */
	public double getSy2(double z, double[] ds_dz)
	{
		return getS2(z + gamma, one_d2, Ay, By, ds_dz);
	}
}
