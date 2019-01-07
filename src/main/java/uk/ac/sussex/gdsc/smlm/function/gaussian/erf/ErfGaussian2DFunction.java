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
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
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
  protected double tB;

  // Required for the first gradients

  /** The x position pre-factors for first-order partial derivatives. */
  protected double[] du_dtx;
  /** The y position pre-factors for first-order partial derivatives. */
  protected double[] du_dty;
  /** The sx pre-factors for first-order partial derivatives. */
  protected double[] du_dtsx;
  /** The sy pre-factors for first-order partial derivatives. */
  protected double[] du_dtsy;

  // Required for the second gradients

  /** The x position pre-factors for second-order partial derivatives. */
  protected double[] d2u_dtx2;
  /** The y position pre-factors for second-order partial derivatives. */
  protected double[] d2u_dty2;
  /** The sx pre-factors for second-order partial derivatives. */
  protected double[] d2u_dtsx2;
  /** The sy pre-factors for second-order partial derivatives. */
  protected double[] d2u_dtsy2;

  // Required for the extended second gradients

  /** The x pre-factors for extended second-order partial derivatives. */
  protected double[] d2deltaEx_dtsxdx;
  /** The y pre-factors for extended second-order partial derivatives. */
  protected double[] d2deltaEy_dtsydy;

  /**
   * Define the error function used.
   */
  public enum ErfFunction {
    /** Use a fast approximation for the Error function. */
    FAST,

    /** Use the Apache commons math library Error function. */
    COMMONS_MATH;
  }

  private ErfFunction erfFunction = ErfFunction.FAST;

  private static interface ErrorFunction {
    public double erf(double x);
  }

  private static class FastErrorFunction implements ErrorFunction {
    @Override
    public double erf(double x) {
      return uk.ac.sussex.gdsc.smlm.function.Erf.erf(x);
    }
  }

  private static class CommontsMathErrorFunction implements ErrorFunction {
    @Override
    public double erf(double x) {
      return org.apache.commons.math3.special.Erf.erf(x);
    }
  }

  private static final FastErrorFunction fastErrorFunction = new FastErrorFunction();
  private static final CommontsMathErrorFunction commontsMathErrorFunction =
      new CommontsMathErrorFunction();
  private ErrorFunction errorFunction = fastErrorFunction;

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
    if (du_dtx != null) {
      return;
    }
    du_dtx = new double[deltaEx.length];
    du_dty = new double[deltaEy.length];
    du_dtsx = new double[deltaEx.length];
    du_dtsy = new double[deltaEy.length];
  }

  /**
   * Creates the arrays needed to compute the first and second order partial derivatives.
   */
  protected void create2Arrays() {
    if (d2u_dtx2 != null) {
      return;
    }
    d2u_dtx2 = new double[deltaEx.length];
    d2u_dty2 = new double[deltaEy.length];
    d2u_dtsx2 = new double[deltaEx.length];
    d2u_dtsy2 = new double[deltaEy.length];
    create1Arrays();
  }

  /**
   * Creates the arrays needed to compute the first and extended second order partial derivatives.
   */
  protected void createEx2Arrays() {
    if (d2deltaEx_dtsxdx != null) {
      return;
    }
    d2deltaEx_dtsxdx = new double[deltaEx.length];
    d2deltaEy_dtsydy = new double[deltaEy.length];
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

  /**
   * Evaluates an 2-dimensional Gaussian function for a single peak.
   *
   * @param i Input predictor
   * @param duda Partial first gradient of function with respect to each coefficient
   * @param d2uda2 Partial second gradient of function with respect to each coefficient
   * @return The predicted value
   */
  @Override
  public abstract double eval(final int i, final double[] duda, final double[] d2uda2);

  // Force new implementation from the base Gaussian2DFunction
  @Override
  public abstract void forEach(Gradient1Procedure procedure);

  /** {@inheritDoc} */
  @Override
  public void initialise(double[] a) {
    // The base Gaussian2DFunction does all the work in NonLinearFunction#initialise(double[]).
    // The ERF Gaussian2DFunction does all the work in Gradient1Function.initialise1(double[])
    initialise1(a);
  }

  // Force new implementation from the base Gaussian2DFunction
  @Override
  public abstract void initialise0(double[] a);

  // Force new implementation from the base Gaussian2DFunction
  @Override
  public abstract void initialise1(double[] a);

  /**
   * Get the absolute of a value.
   *
   * @param value the value
   * @return the absolute value
   */
  protected static double abs(double value) {
    // return Math.abs(d);
    return (value <= 0.0D) ? 0.0d - value : value;
  }

  /**
   * Safe divide the numerator by the denominator. If the numerator is zero then zero is returned
   * avoiding a potential divide by zero producing a NaN. This can happen when computing gradients
   * if the numerator and denominator are both zero.
   *
   * @param numerator the numerator
   * @param denominator the denominator
   * @return the double
   */
  protected static double safeDivide(double numerator, double denominator) {
    // Note: This could be used when computing gradients:
    // e.g. final double du_dty_tI = du_dty / tI;
    // => final double du_dty_tI = safeDivide(du_dty, tI);
    // Currently this is not performed as the ERF functions are used in the context
    // of bounded parameters so avoiding bad parameters, e.g. tI being zero.
    return (numerator == 0) ? 0 : numerator / denominator;
  }

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
        errorFunction = commontsMathErrorFunction;
        break;
      case FAST:
        errorFunction = fastErrorFunction;
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
   * <p>The value returned is always between -1 and 1 (inclusive). If {@code abs(x) > 40}, then
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
   * @param one_sSqrt2 one over (s times sqrt(2))
   * @param n the n
   * @param u the mean of the Gaussian
   * @return the integral
   */
  protected double compute1DIntegral(double one_sSqrt2, int n, double u) {
    return 0.5 * (erf((n - u) * one_sSqrt2) - erf(-u * one_sSqrt2));
  }
}
