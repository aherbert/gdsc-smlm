package gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.Pair;

import gdsc.smlm.function.ExtendedNonLinearFunction;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.GradientFunction;
import gdsc.smlm.function.NoiseModel;
import gdsc.smlm.function.ValueProcedure;

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
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The function will calculate the value of the Gaussian and evaluate the gradient of a set of parameters. The class can
 * specify which of the following parameters the function will evaluate:<br/>
 * background, signal, angle, position0, position1, sd0, sd1
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class Gaussian2DFunction implements ExtendedNonLinearFunction, GradientFunction
{
	/**
	 * The factor for converting a Gaussian standard deviation to Full Width at Half Maxima (FWHM)
	 */
	public static final double SD_TO_FWHM_FACTOR = (2.0 * Math.sqrt(2.0 * Math.log(2.0)));

	/**
	 * The factor for converting a Gaussian standard deviation to Half Width at Half Maxima (FWHM)
	 */
	public static final double SD_TO_HWHM_FACTOR = (Math.sqrt(2.0 * Math.log(2.0)));

	private NoiseModel noiseModel = null;

	public static final double ONE_OVER_TWO_PI = 0.5 / Math.PI;

	public static final int BACKGROUND = 0;
	public static final int SIGNAL = 1;
	public static final int SHAPE = 2;
	public static final int X_POSITION = 3;
	public static final int Y_POSITION = 4;
	public static final int X_SD = 5;
	public static final int Y_SD = 6;

	/**
	 * Gets the name of the parameter assuming a 2D Gaussian function packed as: background + n * [signal, shape,
	 * position0, position1, sd0, sd1].
	 *
	 * @param index
	 *            the index (zero or above)
	 * @return the name
	 */
	public String getName(int index)
	{
		final int i = 1 + (index - 1) % 6;
		switch (i)
		{
			//@formatter:off
			case 0: return "Background";
			case 1: return "Signal";
			case 2: return getShapeName();
			case 3: return "X";
			case 4: return "Y";
			case 5: return "X SD";
			case 6: return "Y SD";
			default: return "Unknown: "+index;
			//@formatter:on
		}
	}

	protected String getShapeName()
	{
		// This method provides a link between the simple Gaussian functions in this package 
		// which evaluate an elliptical Gaussian with an angle of rotation and the functions 
		// in the erf sub-package which support z-depth.
		return "Angle";
	}

	protected final int maxx, maxy;

	/**
	 * Instantiates a new gaussian 2 D function.
	 *
	 * @param maxx
	 *            The maximum x value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 * @param maxy
	 *            The maximum y value of the 2-dimensional data (used to unpack a linear index into coordinates)
	 */
	public Gaussian2DFunction(int maxx, int maxy)
	{
		this.maxx = (maxx < 1) ? 1 : maxx;
		this.maxy = (maxy < 1) ? 1 : maxy;
	}

	/**
	 * @return the dimensions
	 */
	public int[] getDimensions()
	{
		return new int[] { maxx, maxy };
	}

	/**
	 * @return the maximum size in the first dimension
	 */
	public int getMaxX()
	{
		return maxx;
	}

	/**
	 * @return the maximum size in the second dimension
	 */
	public int getMaxY()
	{
		return maxy;
	}

	/**
	 * Gets the order of partial derivatives that can be computed.
	 *
	 * @return the derivative order
	 */
	public int getDerivativeOrder()
	{
		return 1;
	}

	/**
	 * Checks if there is an overhead to computing derivatives of the given order. The input order is expected to be
	 * less than the value returned from {@link #getDerivativeOrder()}.
	 *
	 * @param derivativeOrder
	 *            the derivative order
	 * @return true, if is overhead
	 */
	public boolean isOverhead(int derivativeOrder)
	{
		return false;
	}

	/**
	 * Creates a new function with the ability to compute derivatives of the chosen order. Note that computation of
	 * higher order derivatives may incur an overhead and so a function can be created with only the derivatives
	 * required.
	 * <p>
	 * If the derivative order matches the current derivative order for this function then the same object may be
	 * returned.
	 *
	 * @param derivativeOrder
	 *            the derivative order
	 * @return the gaussian 2D function
	 */
	public Gaussian2DFunction create(int derivativeOrder)
	{
		return this;
	}

	/**
	 * Copy the function.
	 *
	 * @return a copy
	 */
	abstract public Gaussian2DFunction copy();

	/**
	 * @return the number of peaks
	 */
	public abstract int getNPeaks();

	/**
	 * @return True if the function can evaluate the background gradient
	 */
	public abstract boolean evaluatesBackground();

	/**
	 * @return True if the function can evaluate the signal gradient
	 */
	public abstract boolean evaluatesSignal();

	/**
	 * @return True if the function can evaluate the shape gradient
	 */
	public abstract boolean evaluatesShape();

	/**
	 * @return True if the function can evaluate the position gradient
	 */
	public abstract boolean evaluatesPosition();

	/**
	 * @return True if the function can evaluate the standard deviation gradient for the 1st dimension
	 */
	public abstract boolean evaluatesSD0();

	/**
	 * @return True if the function can evaluate the standard deviation gradient for the 2nd dimension
	 */
	public abstract boolean evaluatesSD1();

	/**
	 * @return The number of parameters per peak
	 */
	public abstract int getParametersPerPeak();

	/**
	 * Execute the {@link #eval(int, float[])} method and set the expected variance using the noise model
	 * 
	 * @throws NullPointerException
	 *             if the noise model is null
	 * @see gdsc.smlm.function.NonLinearFunction#eval(int, float[], float[])
	 */
	public double eval(final int x, final double[] dyda, final double[] w) throws NullPointerException
	{
		final double value = eval(x, dyda);
		//w[0] = (noiseModel == null) ? 1 : noiseModel.variance(value);
		// Just throw a null pointer exception if noiseModel is null
		w[0] = noiseModel.variance(value);
		return value;
	}

	/**
	 * Execute the {@link #eval(int)} method and set the expected variance using the noise model
	 * 
	 * @throws NullPointerException
	 *             if the noise model is null
	 * @see gdsc.smlm.function.NonLinearFunction#evalw(int, double[])
	 */
	public double evalw(int x, double[] w)
	{
		final double value = eval(x);
		//w[0] = (noiseModel == null) ? 1 : noiseModel.variance(value);
		// Just throw a null pointer exception if noiseModel is null
		w[0] = noiseModel.variance(value);
		return value;
	}

	/**
	 * @return the noise model
	 */
	public NoiseModel getNoiseModel()
	{
		return noiseModel;
	}

	/**
	 * Set the noise model used in {@link #eval(int, float[], float[])}.
	 * 
	 * @param noiseModel
	 *            the noise model to set
	 */
	public void setNoiseModel(NoiseModel noiseModel)
	{
		this.noiseModel = noiseModel;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.function.NonLinearFunction#canComputeWeights()
	 */
	public boolean canComputeWeights()
	{
		return (noiseModel != null);
	}

	/**
	 * Build the index array that maps the gradient index back to the original parameter index so that:<br/>
	 * a[indices[i]] += dy_da[i]
	 * 
	 * @param nPeaks
	 * @return The indices
	 */
	protected int[] createGradientIndices(int nPeaks)
	{
		return createGradientIndices(nPeaks, this);
	}

	protected static int[] createGradientIndices(int nPeaks, Gaussian2DFunction gf)
	{
		// Parameters are: 
		// Background + n * { Signal, Angle, Xpos, Ypos, Xsd, Ysd }
		int nparams = (gf.evaluatesBackground() ? 1 : 0) + nPeaks * gf.getParametersPerPeak();
		int[] indices = new int[nparams];

		int p = 0;
		if (gf.evaluatesBackground())
			indices[p++] = 0;
		for (int n = 0, i = 0; n < nPeaks; n++, i += 6)
		{
			if (gf.evaluatesSignal())
				indices[p++] = i + SIGNAL;
			if (gf.evaluatesShape())
				indices[p++] = i + SHAPE;
			// All functions evaluate the position gradient
			indices[p++] = i + X_POSITION;
			indices[p++] = i + Y_POSITION;
			if (gf.evaluatesSD0())
				indices[p++] = i + X_SD;
			if (gf.evaluatesSD1())
				indices[p++] = i + Y_SD;
		}

		return indices;
	}

	/**
	 * Locate the index within the gradient indices for the specified parameter
	 * 
	 * @param parameterIndex
	 * @return the gradient index (or -1 if not present)
	 */
	public int findGradientIndex(int parameterIndex)
	{
		int[] gradientIndices = gradientIndices();
		for (int i = 0; i < gradientIndices.length; i++)
			if (gradientIndices[i] == parameterIndex)
				return i;
		return -1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeValues(double[])
	 */
	public double[] computeValues(double[] variables)
	{
		initialise(variables);
		final double[] values = new double[size()];
		forEach(new ValueProcedure()
		{
			int i = 0;

			public void execute(double value)
			{
				values[i++] = value;
			}
		});
		return values;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeJacobian(double[])
	 */
	public double[][] computeJacobian(double[] variables)
	{
		return computeValuesAndJacobian(variables).getSecond();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#canComputeValuesAndJacobian()
	 */
	public boolean canComputeValuesAndJacobian()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeValuesAndJacobian(double[])
	 */
	public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables)
	{
		initialise(variables);
		final int n = size();
		final double[][] jacobian = new double[n][];
		final double[] values = new double[n];
		forEach(new Gradient1Procedure()
		{
			int i = 0;

			public void execute(double value, double[] dy_da)
			{
				values[i] = value;
				jacobian[i++] = dy_da.clone();
			}
		});
		return new Pair<double[], double[][]>(values, jacobian);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#size()
	 */
	public int size()
	{
		return maxx * maxy;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#getNumberOfGradients()
	 */
	public int getNumberOfGradients()
	{
		return gradientIndices().length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.ValueProcedure)
	 */
	public void forEach(ValueProcedure procedure)
	{
		for (int i = 0, n = size(); i < n; i++)
			procedure.execute(eval(i));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientFunction#forEach(gdsc.smlm.function.Gradient1Procedure)
	 */
	public void forEach(Gradient1Procedure procedure)
	{
		final double[] duda = new double[getNumberOfGradients()];
		for (int i = 0, n = size(); i < n; i++)
		{
			final double value = eval(i, duda);
			procedure.execute(value, duda);
		}
	}
}
