/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Abstract base class for a 2-dimensional Gaussian function for a configured number of peaks.
 *
 * <p>The function will calculate the value of the Gaussian and evaluate the gradient of a set of
 * parameters. The class can specify which of the following parameters the function will
 * evaluate:<br> background, signal, z-depth, position0, position1, sd0, sd1
 *
 * <p>The class provides an index of the position in the parameter array where the parameter is
 * expected.
 */
public abstract class ErfGaussian2DFunction extends Gaussian2DFunction
    implements Gradient2Function, ExtendedGradient2Function {
  /** The constant <code>1.0 / Math.sqrt(2)</code>. */
  protected static final double ONE_OVER_ROOT2 = 1.0 / Math.sqrt(2);
  /** The constant <code>1.0 / Math.sqrt(2 * Math.PI)</code>. */
  protected static final double ONE_OVER_ROOT2PI = 1.0 / Math.sqrt(2 * Math.PI);

  // Required for the PSF

  /** The delta Ex array (difference in Gaussian integral in the x-dimension for each pixel). */
  protected double[] deltaEx;
  /** The delta Ey array (difference in Gaussian integral in the y-dimension for each pixel). */
  protected double[] deltaEy;

  /** The background. */
  protected double tb;

  // Required for the first gradients

  /** The x position pre-factors for first-order partial derivatives. */
  protected double[] duDtx;
  /** The y position pre-factors for first-order partial derivatives. */
  protected double[] duDty;
  /** The sx pre-factors for first-order partial derivatives. */
  protected double[] duDtsx;
  /** The sy pre-factors for first-order partial derivatives. */
  protected double[] duDtsy;

  // Required for the second gradients

  /** The x position pre-factors for second-order partial derivatives. */
  protected double[] d2uDtx2;
  /** The y position pre-factors for second-order partial derivatives. */
  protected double[] d2uDty2;
  /** The sx pre-factors for second-order partial derivatives. */
  protected double[] d2uDtsx2;
  /** The sy pre-factors for second-order partial derivatives. */
  protected double[] d2uDtsy2;

  // Required for the extended second gradients

  /** The x pre-factors for extended second-order partial derivatives. */
  protected double[] d2deltaExDtsxDx;
  /** The y pre-factors for extended second-order partial derivatives. */
  protected double[] d2deltaEyDtsyDy;

  /** The erf function. */
  private ErfFunction erfFunction = ErfFunction.FAST;

  /**
   * The instance of the error function. This should not be null. It is updated when the Erf
   * function property is changed.
   */
  private ErrorFunction errorFunction = FastErrorFunction.INSTANCE;

  /**
   * Define the error function used.
   */
  public enum ErfFunction {
    /** Use a fast approximation for the Error function. */
    FAST,

    /** Use the Apache commons math library Error function. */
    COMMONS_MATH;
  }

  @FunctionalInterface
  private interface ErrorFunction {
    public double erf(double x);
  }

  private static class FastErrorFunction implements ErrorFunction {
    static final FastErrorFunction INSTANCE = new FastErrorFunction();

    @Override
    public double erf(double x) {
      return uk.ac.sussex.gdsc.smlm.function.Erf.erf(x);
    }
  }

  private static class CommontsMathErrorFunction implements ErrorFunction {
    static final CommontsMathErrorFunction INSTANCE = new CommontsMathErrorFunction();

    @Override
    public double erf(double x) {
      return org.apache.commons.math3.special.Erf.erf(x);
    }
  }

  /**
   * Instantiates a new erf gaussian 2D function.
   *
   * @param numberOfPeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public ErfGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(maxx, maxy);
    deltaEx = new double[numberOfPeaks * this.maxx];
    deltaEy = new double[numberOfPeaks * this.maxy];
  }

  /**
   * Creates the arrays needed to compute the first-order partial derivatives.
   */
  protected void create1Arrays() {
    if (duDtx != null) {
      return;
    }
    duDtx = new double[deltaEx.length];
    duDty = new double[deltaEy.length];
    duDtsx = new double[deltaEx.length];
    duDtsy = new double[deltaEy.length];
  }

  /**
   * Creates the arrays needed to compute the first and second order partial derivatives.
   */
  protected void create2Arrays() {
    if (d2uDtx2 != null) {
      return;
    }
    d2uDtx2 = new double[deltaEx.length];
    d2uDty2 = new double[deltaEy.length];
    d2uDtsx2 = new double[deltaEx.length];
    d2uDtsy2 = new double[deltaEy.length];
    create1Arrays();
  }

  /**
   * Creates the arrays needed to compute the first and extended second order partial derivatives.
   */
  protected void createEx2Arrays() {
    if (d2deltaExDtsxDx != null) {
      return;
    }
    d2deltaExDtsxDx = new double[deltaEx.length];
    d2deltaEyDtsyDy = new double[deltaEy.length];
    create2Arrays();
  }

  /**
   * Copy the function.
   *
   * @return a copy
   */
  @Override
  public abstract ErfGaussian2DFunction copy();

  @Override
  public boolean evaluatesAngle() {
    // None of the ERF functions support rotation due to the strict XY separation
    return false;
  }

  // Force new implementation from the base Gaussian2DFunction
  @Override
  public abstract void forEach(Gradient1Procedure procedure);

  @Override
  public void initialise(double[] a) {
    // The base Gaussian2DFunction does all the work in NonLinearFunction#initialise(double[]).
    // The ERF Gaussian2DFunction does all the work in Gradient1Function.initialise1(double[])
    initialise1(a);
  }

  // Force new implementation from the base Gaussian2DFunction
  // (which just calls NonLinearFunction#initialise(double[]))
  @Override
  public abstract void initialise0(double[] a);

  // Force new implementation from the base Gaussian2DFunction
  // (which just calls NonLinearFunction#initialise(double[]))
  @Override
  public abstract void initialise1(double[] a);

  /**
   * Evaluates a 2-dimensional Gaussian function for a given input predictor (x) and partial
   * gradients for each of the coefficients (a).
   *
   * <p>Note: This is a logical extension of the support for the {@link NonLinearFunction} interface
   * to allow access to the second order derivatives using an input predictor. The function should
   * first be initialised using {@link Gradient2Function#initialise2(double[])}.
   *
   * @param x Predictor
   * @param dyda First order partial gradient of function with respect to each coefficient
   *        identified by {@link #gradientIndices()}
   * @param d2yda2 Second order partial gradient of function with respect to each coefficient
   *        identified by {@link #gradientIndices()}
   * @return The predicted value y
   */
  public abstract double eval2(final int x, final double[] dyda, final double[] d2yda2);

  /**
   * Gets the erf function.
   *
   * @return the erf function
   */
  public ErfFunction getErfFunction() {
    return erfFunction;
  }

  /**
   * Sets the erf function.
   *
   * @param erfFunction the new erf function
   * @throws IllegalArgumentException If the error function is unknown
   */
  public void setErfFunction(ErfFunction erfFunction) {
    switch (erfFunction) {
      case COMMONS_MATH:
        errorFunction = CommontsMathErrorFunction.INSTANCE;
        break;
      case FAST:
        errorFunction = FastErrorFunction.INSTANCE;
        break;
      default:
        throw new IllegalArgumentException("Unknown error function: " + erfFunction);
    }
    this.erfFunction = erfFunction;
  }

  /**
   * Returns the error function.
   *
   * <p>erf(x) = 2/&radic;&pi; <sub>0</sub>&int;<sup>x</sup> e<sup>-t*t</sup>dt </p>
   *
   * <p>Uses the configured implementation (see {@link #getErfFunction()}). </p>
   *
   * <p>The value returned is always between -1 and 1 (inclusive). If {@code Math.abs(x) > 40}, then
   * {@code erf(x)} is indistinguishable from either 1 or -1 as a double, so the appropriate extreme
   * value is returned. </p>
   *
   * @param x the value.
   * @return the error function erf(x)
   */
  public double erf(double x) {
    return errorFunction.erf(x);
  }

  /**
   * Compute the 1D integral from 0 to n. This is the sum of the Gaussian function using the error
   * function for all of the pixels from 0 to n.
   *
   * @param oneOverSsqrt2 one over (s times sqrt(2))
   * @param n the n
   * @param u the mean of the Gaussian
   * @return the integral
   */
  protected double compute1DIntegral(double oneOverSsqrt2, int n, double u) {
    return 0.5 * (erf((n - u) * oneOverSsqrt2) - erf(-u * oneOverSsqrt2));
  }
}
