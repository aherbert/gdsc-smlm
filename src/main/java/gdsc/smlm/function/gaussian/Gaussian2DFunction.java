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
package gdsc.smlm.function.gaussian;

import java.util.Arrays;

import gdsc.smlm.function.ExtendedNonLinearFunction;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.Gradient1Procedure;
import gdsc.smlm.function.IntegralValueProcedure;
import gdsc.smlm.function.NamedFunction;
import gdsc.smlm.function.NoiseModel;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.utils.Pair;


/**
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
 * <p>
 * The function will calculate the value of the Gaussian and evaluate the gradient of a set of parameters. The class can
 * specify which of the following parameters the function will evaluate:<br/>
 * background, signal, position0, position1, sd0, sd1, angle. A parameter is provided for position2 (z-depth) to support
 * 3D function using astimatism.
 * <p>
 * The class provides an index of the position in the parameter array where the parameter is expected.
 */
public abstract class Gaussian2DFunction implements ExtendedNonLinearFunction, Gradient1Function, NamedFunction
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

	/** Index of the background in the parameters array */
	public static final int BACKGROUND = 0;
	/** Index of the signal intensity in the parameters array */
	public static final int SIGNAL = 1;
	/** Index of the x-position in the parameters array */
	public static final int X_POSITION = 2;
	/** Index of the y-position in the parameters array */
	public static final int Y_POSITION = 3;
	/** Index of the z-position in the parameters array */
	public static final int Z_POSITION = 4;
	/** Index of the x-standard deviation in the parameters array */
	public static final int X_SD = 5;
	/** Index of the y-standard deviation in the parameters array */
	public static final int Y_SD = 6;
	/** Index of the angle in the parameters array */
	public static final int ANGLE = 7;

	/** The number of parameters per Gaussian peak */
	public static final int PARAMETERS_PER_PEAK = 7;

	/**
	 * Gets the name of the parameter assuming a 2D Gaussian function.
	 *
	 * @param index
	 *            the index (zero or above)
	 * @return the name
	 */
	public static String getName(int index)
	{
		final int i = 1 + (index - 1) % PARAMETERS_PER_PEAK;
		switch (i)
		{
			//@formatter:off
			case BACKGROUND: return "Background";
			case SIGNAL: return "Signal";
			case X_POSITION: return "X";
			case Y_POSITION: return "Y";
			case Z_POSITION: return "Z";
			case X_SD: return "X SD";
			case Y_SD: return "Y SD";
			case ANGLE: return "Angle";
			default: return "Unknown: "+index;
			//@formatter:on
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.NamedFunction#getParameterName(int)
	 */
	public String getParameterName(int i)
	{
		return getName(i);
	}

	/**
	 * Gets the peak number (zero-based index) of the parameter assuming a 2D Gaussian function.
	 *
	 * @param index
	 *            the index (zero or above)
	 * @return the peak number
	 */
	public static int getPeak(int index)
	{
		if (index < 1)
			return 0;
		return (index - 1) / PARAMETERS_PER_PEAK;
	}

	/**
	 * Gets the index of the parameter in a multi-peak parameter array assuming a 2D Gaussian function.
	 *
	 * @param peak
	 *            the peak number (zero-based index)
	 * @param parameterIndex
	 *            the parameter index for a single peak (this can use the class constants, e.g.
	 *            {@link Gaussian2DFunction#SIGNAL})
	 * @return the index
	 */
	public static int getIndex(int peak, int parameterIndex)
	{
		if (parameterIndex < 1)
			return 0;
		return peak * PARAMETERS_PER_PEAK + parameterIndex;
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
	 * @return True if the function can evaluate the XY-position gradient
	 */
	public abstract boolean evaluatesPosition();

	/**
	 * @return True if the function can evaluate the Z-position gradient
	 */
	public boolean evaluatesZ()
	{
		// No standard Gaussian 2D function evaluates the z-position
		return false;
	}

	/**
	 * @return True if the function can evaluate the standard deviation gradient for the 1st dimension
	 */
	public abstract boolean evaluatesSD0();

	/**
	 * @return True if the function can evaluate the standard deviation gradient for the 2nd dimension
	 */
	public abstract boolean evaluatesSD1();

	/**
	 * @return True if the function can evaluate the angle gradient
	 */
	public abstract boolean evaluatesAngle();

	/**
	 * @return The number of gradient parameters per peak
	 */
	public abstract int getGradientParametersPerPeak();

	/**
	 * Produce an output predicted value for a given set of input
	 * predictors (x) and coefficients (a).
	 * <p>
	 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
	 * <p>
	 * The first coefficient is the Gaussian background level. The coefficients are then packed for each peak
	 * using the indices specified in the Gaussian2DFunction class.
	 * 
	 * @param x
	 *            Input predictor
	 * @return The predicted value
	 * 
	 * @see gdsc.smlm.function.NonLinearFunction#eval(int)
	 */
	public abstract double eval(final int x);

	/**
	 * Produce an output predicted value for a given set of input
	 * predictors (x) and coefficients (a).
	 * <p>
	 * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
	 * <p>
	 * The first coefficient is the Gaussian background level. The coefficients are then packed for each peak
	 * using the indices specified in the Gaussian2DFunction class.
	 * 
	 * @param x
	 *            Input predictor
	 * @param dyda
	 *            Partial gradient of function with respect to each coefficient
	 * @return The predicted value
	 * 
	 * @see gdsc.smlm.function.NonLinearFunction#eval(int, double[])
	 */
	public abstract double eval(final int x, final double[] dyda);

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
		// Background + n * { Signal, Shape, Xpos, Ypos, Xsd, Ysd }
		int nparams = (gf.evaluatesBackground() ? 1 : 0) + nPeaks * gf.getGradientParametersPerPeak();
		int[] indices = new int[nparams];

		int p = 0;
		if (gf.evaluatesBackground())
			indices[p++] = 0;
		for (int n = 0, i = 0; n < nPeaks; n++, i += PARAMETERS_PER_PEAK)
		{
			if (gf.evaluatesSignal())
				indices[p++] = i + SIGNAL;
			// All functions evaluate the position gradient
			indices[p++] = i + X_POSITION;
			indices[p++] = i + Y_POSITION;
			if (gf.evaluatesZ())
				indices[p++] = i + Z_POSITION;
			if (gf.evaluatesSD0())
				indices[p++] = i + X_SD;
			if (gf.evaluatesSD1())
				indices[p++] = i + Y_SD;
			if (gf.evaluatesAngle())
				indices[p++] = i + ANGLE;
		}

		return indices;
	}

	/**
	 * Gets the name of the gradient parameter.
	 *
	 * @param index
	 *            the index (must be within the array returned from {@link #gradientIndices()})
	 * @return the name
	 */
	public String getGradientParameterName(int index)
	{
		return getName(gradientIndices()[index]);
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
		initialise0(variables);
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

	/**
	 * Compute the integral. This is the sum of the values.
	 *
	 * @param a
	 *            an array of coefficients
	 * @return the integral
	 */
	public double integral(double[] a)
	{
		return new IntegralValueProcedure().getIntegral(this, a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ExtendedNonLinearFunction#computeJacobian(double[])
	 */
	public double[][] computeJacobian(double[] variables)
	{
		return computeValuesAndJacobian(variables).b;
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
		initialise1(variables);
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueFunction#initialise0(double[])
	 */
	public void initialise0(double[] a)
	{
		// TODO - Update these functions to support initialisation
		// for computing the value only
		initialise(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Function#initialise1(double[])
	 */
	public void initialise1(double[] a)
	{
		initialise(a);
	}

	/**
	 * Check if NaN (invalid) gradients.
	 *
	 * @param a
	 *            the gradients
	 * @return true, if successful
	 */
	protected static boolean invalidGradients(double[] a)
	{
		for (int i = 0; i < a.length; i++)
			if (Double.isNaN(a[i]))
			{
				System.out.println(Arrays.toString(a));
				return true;
			}
		return false;
	}
}
