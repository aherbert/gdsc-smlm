/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

/**
 * This is a wrapper for any function to compute the negative log-likelihood assuming a per-pixel
 * Poisson distribution with Gaussian noise (i.e. the noise from a sCMOS camera). This uses the
 * MLE-sCMOS formula from Huang, et al (2013), Supplementary Notes Eq 3.3:<br> P_sCMOS
 * (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi) = e^-(ui+vari/gi^2) (ui+vari/gi^2)^x / gamma(x+1) <br>
 * Where:<br> i = the pixel index <br> vari = the variance of the pixel <br> gi = the gain of the
 * pixel <br> oi = the offset of the pixel <br> ui = the function value (expected number of photons)
 * <br> Di = the observed value at the pixel x = the observed random variable (observed number of
 * photons adjusted by a pixel dependent constant) <br>
 *
 * <p>The negative log-likelihood function is: <br> -LL(P_sCMOS (x=[(Di-oi)/gi +
 * vari/gi^2]|ui,vari,gi,oi)) <br> = (ui+vari/gi^2) - x * ln(ui+vari/gi^2) + ln(gamma(x+1)) <br>
 *
 * <p>The negative log-likelihood (and gradient) can be evaluated over the entire set of observed
 * values or for a chosen observed value.
 *
 * <p>To allow a likelihood to be computed: (a) when the function predicts negative photon count
 * data the function prediction is set to zero; (b) if the observed random variable (x) is negative
 * it is also set to zero. This occurs when true signal readout from the sCMOS camera is low enough
 * to be negated by readout noise. In this case the noise can be ignored.
 *
 * <p>See: Hunag, et al (2013) Video-rate nanoscopy using sCMOS camera–specific single-molecule
 * localization algorithms. Nature Methods 10, 653–658.
 */
public class ScmosLikelihoodWrapper extends LikelihoodWrapper {
  private final double logNormalisation;
  private final double[] varG2;
  private final double[] x;
  private final double[] logG;

  private double observedLikelihood = Double.NaN;

  /**
   * Initialise the function.
   *
   * <p>The input parameters must be the full parameters for the non-linear function. Only those
   * parameters with gradient indices should be passed in to the functions to obtain the value (and
   * gradient).
   *
   * @param func The function to be used to calculated the expected values (Note that the expected
   *        value is the number of photons)
   * @param a The initial parameters for the function
   * @param data The observed values from the sCMOS camera
   * @param n The number of observed values
   * @param var the variance of each pixel
   * @param gain the gain of each pixel
   * @param offset the offset of each pixel
   * @throws IllegalArgumentException if the input observed values are not integers
   */
  public ScmosLikelihoodWrapper(NonLinearFunction func, double[] a, double[] data, int n,
      float[] var, float[] gain, float[] offset) {
    super(func, a, data, n);

    varG2 = new double[n];
    x = new double[n];
    logG = new double[n];

    // Pre-compute the sum over the data
    double sum = 0;
    for (int i = 0; i < n; i++) {
      varG2[i] = var[i] / (gain[i] * gain[i]);
      x[i] = Math.max(0, (data[i] - offset[i]) / gain[i] + varG2[i]);
      logG[i] = Math.log(gain[i]);

      sum += logGamma1(x[i]) + logG[i];
    }

    logNormalisation = sum;
  }

  /**
   * Initialise the function using pre-computed per pixel working variables. This allows the
   * pre-computation to be performed once for the sCMOS pixels for all likelihood computations.
   *
   * <p>The input parameters must be the full parameters for the non-linear function. Only those
   * parameters with gradient indices should be passed in to the functions to obtain the value (and
   * gradient).
   *
   * @param func The function to be used to calculated the expected values (Note that the expected
   *        value is the number of photons)
   * @param a The initial parameters for the function
   * @param x The observed values from the sCMOS camera mapped using [x=max(0, (k-o)/g + var/g^2)]
   *        per pixel
   * @param n The number of observed values
   * @param varG2 the variance of each pixel divided by the gain squared
   * @param logG the log of the gain of each pixel
   * @throws IllegalArgumentException if the input observed values are not integers
   */
  public ScmosLikelihoodWrapper(NonLinearFunction func, double[] a, double[] x, int n,
      double[] varG2, double[] logG) {
    super(func, a, x, n);

    this.varG2 = varG2;
    this.x = x;
    this.logG = logG;

    // Pre-compute the sum over the data
    double sum = 0;
    for (int i = 0; i < n; i++) {
      sum += logGamma1(x[i]) + logG[i];
    }

    logNormalisation = sum;
  }

  /**
   * Copy constructor.
   *
   * @param func The function to be used to calculated the expected values (Note that the expected
   *        value is the number of photons)
   * @param a The initial parameters for the function
   * @param x The observed values from the sCMOS camera mapped using [x=max(0, (k-o)/g + var/g^2)]
   *        per pixel
   * @param n The number of observed values
   * @param varG2 the variance of each pixel divided by the gain squared
   * @param logG the log of the gain of each pixel
   * @param logNormalisation the log normalisation
   * @throws IllegalArgumentException if the input observed values are not integers
   */
  private ScmosLikelihoodWrapper(NonLinearFunction func, double[] a, double[] x, int n,
      double[] varG2, double[] logG, double logNormalisation) {
    super(func, a, x, n);
    this.varG2 = varG2;
    this.x = x;
    this.logG = logG;
    this.logNormalisation = logNormalisation;
  }

  /**
   * Builds a new instance with a new function. All pre-computation on the data is maintained.
   *
   * @param func The function to be used to calculated the expected values (Note that the expected
   *        value is the number of photons)
   * @param a The initial parameters for the function
   * @return the SCMOS likelihood wrapper
   */
  public ScmosLikelihoodWrapper build(NonLinearFunction func, double[] a) {
    return new ScmosLikelihoodWrapper(func, a, x, dataSize, varG2, logG, logNormalisation);
  }

  /**
   * Compute variance divided by the gain squared. This can be used in the
   * {@link #ScmosLikelihoodWrapper(NonLinearFunction, double[], double[], int, double[], double[])}
   * constructor.
   *
   * @param var the variance of each pixel
   * @param gain the gain of each pixel
   * @return the variance divided by the gain squared
   */
  public static double[] computeVarG2(float[] var, float[] gain) {
    final int n = Math.min(var.length, gain.length);
    final double[] varG2 = new double[n];
    for (int i = 0; i < n; i++) {
      varG2[i] = var[i] / (gain[i] * gain[i]);
    }
    return varG2;
  }

  /**
   * Compute log of the gain.
   *
   * @param gain the gain of each pixel
   * @return the log of the gain
   */
  public static double[] computeLogG(float[] gain) {
    final int n = gain.length;
    final double[] logG = new double[n];
    for (int i = 0; i < n; i++) {
      logG[i] = Math.log(gain[i]);
    }
    return logG;
  }

  /**
   * Compute X (the mapped observed values from the sCMOS camera).
   *
   * @param data The observed values from the sCMOS camera
   * @param varG2 the variance divided by the gain squared
   * @param gain the gain of each pixel
   * @param offset the offset of each pixel
   * @return The observed values from the sCMOS camera mapped using [x=max(0, (k-o)/g + var/g^2)]
   *         per pixel
   */
  public static double[] computeX(double[] data, float[] varG2, float[] gain, float[] offset) {
    final int n = data.length;
    final double[] x = new double[n];
    for (int i = 0; i < n; i++) {
      x[i] = Math.max(0, (data[i] - offset[i]) / gain[i] + varG2[i]);
    }
    return x;
  }

  @Override
  public double computeLikelihood() {
    // Compute the negative log-likelihood to be minimised:
    // (ui+vari/gi^2) - x * ln(ui+vari/gi^2) + ln(gamma(x+1))
    double ll = 0;
    for (int i = 0; i < dataSize; i++) {
      double ui = function.eval(i);

      if (ui < 0) {
        ui = 0;
      }

      final double l = ui + varG2[i];

      ll += l;
      if (x[i] != 0) {
        ll -= x[i] * Math.log(l);
      }
    }
    return ll + logNormalisation;
  }

  @Override
  public double computeLikelihood(double[] gradient) {
    // Compute the negative log-likelihood to be minimised
    // (ui+vari/gi^2) - x * ln(ui+vari/gi^2) + ln(gamma(x+1))
    // with x as the mapped observed value: x = (k-o)/g + var/g^2

    // To compute the gradient we do the same as for a Poisson distribution:
    // f(x) = l(x) - k * ln(l(x)) + log(gamma(k+1))
    // with l(x) as the Poisson mean (the output dependent on the function variables x)
    // and k the observed value.
    //
    // Since (k * ln(l(x)))' = (k * ln(l(x))') * l'(x)
    // = (k / l(x)) * l'(x)

    // f'(x) = l'(x) - (k/l(x) * l'(x))
    // f'(x) = l'(x) * (1 - k/l(x))

    double ll = 0;
    for (int j = 0; j < numberOfVariables; j++) {
      gradient[j] = 0;
    }
    final double[] dlda = new double[numberOfVariables];
    for (int i = 0; i < dataSize; i++) {
      double ui = function.eval(i, dlda);

      if (ui < 0) {
        ui = 0;
      }

      final double l = ui + varG2[i];

      ll += l;
      if (x[i] != 0) {
        ll -= x[i] * Math.log(l);
      }

      // Note: if l==0 then we get divide by zero and a NaN value
      final double factor = (1 - x[i] / l);
      for (int j = 0; j < gradient.length; j++) {
        gradient[j] += dlda[j] * factor;
      }
    }
    return ll + logNormalisation;
  }

  @Override
  public double computeLikelihood(int index) {
    double ui = function.eval(index);

    if (ui < 0) {
      ui = 0;
    }

    final double l = ui + varG2[index];

    double ll = l + logG[index];
    if (x[index] != 0) {
      ll += logGamma1(x[index]) - x[index] * Math.log(l);
    }

    return ll;
  }

  @Override
  public double computeLikelihood(double[] gradient, int index) {
    for (int j = 0; j < numberOfVariables; j++) {
      gradient[j] = 0;
    }
    final double[] dlda = new double[numberOfVariables];

    double ui = function.eval(index, dlda);

    if (ui < 0) {
      ui = 0;
    }

    final double l = ui + varG2[index];

    final double factor = (1 - x[index] / l);
    for (int j = 0; j < gradient.length; j++) {
      gradient[j] = dlda[j] * factor;
    }

    double ll = l + logG[index];
    if (x[index] != 0) {
      ll += logGamma1(x[index]) - x[index] * Math.log(l);
    }

    return ll;
  }

  private static double logGamma1(double value) {
    if (value <= 1) {
      return 0;
    }
    return Gamma.logGamma(value + 1);
  }

  /**
   * Compute the observed negative log likelihood. This is the value of {@link #computeLikelihood()}
   * if the function were to return the observed values for each point.
   *
   * @return the observed negative log likelihood
   */
  public double computeObservedLikelihood() {
    if (Double.isNaN(observedLikelihood)) {
      double ll = 0;
      for (int i = 0; i < dataSize; i++) {
        // We need to input the observed value as the expected value.
        // So we need (k-o)/g as the expected value. We did not store this so
        // compute it by subtracting varG2 from x.
        // Then perform the same likelihood computation.

        // double u = x[i] - varG2[i];
        //
        // if (u < 0)
        // u = 0;
        //
        // double l = u + varG2[i];

        // We can do this in one step ...
        final double l = (x[i] < varG2[i]) ? varG2[i] : x[i];

        ll += l;
        if (x[i] != 0) {
          ll -= x[i] * Math.log(l);
        }
      }
      observedLikelihood = ll + logNormalisation;
    }
    return observedLikelihood;
  }

  /**
   * Compute log likelihood ratio.
   *
   * @param ll the negative log likelihood of the function
   * @return the log likelihood ratio
   */
  public double computeLogLikelihoodRatio(double ll) {
    // From https://en.wikipedia.org/wiki/Likelihood-ratio_test#Use:
    // LLR = -2 * [ ln(likelihood for alternative model) - ln(likelihood for null model)]
    // The model with more parameters (here alternative) will always fit at least as well—
    // i.e., have the same or greater log-likelihood—than the model with fewer parameters
    // (here null)

    final double llAlternative = computeObservedLikelihood();
    double llNull = ll;

    // The alternative should always fit better than the null model
    if (llNull < llAlternative) {
      llNull = llAlternative;
    }

    // Since we have negative log likelihood we reverse the sign
    // return 2 * (-llAlternative - -llNull);
    return -2 * (llAlternative - llNull);
  }

  /**
   * Compute the q-value of the log-likelihood ratio. This is the probability that a value of LLR as
   * poor as the value should occur by chance.
   *
   * @param ll the minimum negative log likelihood of the function (the null model)
   * @return the p-value
   */
  public double computeQValue(double ll) {
    final double llr = computeLogLikelihoodRatio(ll);
    final int degreesOfFreedom = x.length - numberOfVariables;
    return ChiSquaredDistributionTable.computeQValue(llr, degreesOfFreedom);
  }

  /**
   * Compute the negative log likelihood.
   *
   * @param ui the expected value (number of photoelectrons)
   * @param var the variance of the pixel
   * @param gain the gain of the pixel
   * @param offset the offset of the pixel
   * @param data The observed value (count from the sCMOS pixel)
   * @return the negative log likelihood
   */
  public static double negativeLogLikelihood(double ui, float var, float gain, float offset,
      double data) {
    final double varG2 = var / (gain * gain);
    final double x = Math.max(0, (data - offset) / gain + varG2);
    if (ui < 0) {
      ui = 0;
    }
    final double l = ui + varG2;
    // Note we need the Math.log(g) to normalise the Poisson distribution to 1
    // since the observed values (k) are scaled by the gain
    double ll = l + Math.log(gain);
    if (x != 0) {
      ll += logGamma1(x) - x * Math.log(l);
    }

    return ll;
  }

  /**
   * Compute the likelihood.
   *
   * @param ui the expected value (number of photoelectrons)
   * @param var the variance of the pixel
   * @param gain the gain of the pixel
   * @param offset the offset of the pixel
   * @param data The observed value (count from the sCMOS pixel)
   * @return the likelihood
   */
  public static double likelihood(double ui, float var, float gain, float offset, double data) {
    // double varG2 = var / (g * g);
    // double x = Math.max(0, (k - o) / g + varG2);
    // double l = u + varG2;
    // double v = FastMath.exp(-l) * Math.pow(l, x) / gamma1(x);
    // if (v != v)
    // throw new RuntimeException("Failed computation");
    // return v;

    final double nll = negativeLogLikelihood(ui, var, gain, offset, data);
    return FastMath.exp(-nll);
  }

  @Override
  public boolean canComputeGradient() {
    return true;
  }

  /**
   * Compute the Fisher's Information Matrix (I) for fitted variables.
   *
   * <pre>
   * Iab = sum(k) 1/(uk+vark/gk^2)  (duk da) * (duk db)
   * </pre>
   *
   * @param variables The variables of the function
   * @return Fisher's Information Matrix (I)
   */
  @Override
  public double[][] fisherInformation(final double[] variables) {
    initialiseFunction(variables);

    final double[] du_da = new double[numberOfVariables];

    final double[][] I = new double[numberOfVariables][numberOfVariables];

    for (int k = 0; k < dataSize; k++) {
      final double uk = function.eval(k, du_da);
      final double yk = 1 / (uk + varG2[k]);
      for (int i = 0; i < numberOfVariables; i++) {
        final double du_dai = yk * du_da[i];
        for (int j = 0; j <= i; j++) {
          I[i][j] += du_dai * du_da[j];
        }
      }
    }

    // Generate symmetric matrix
    for (int i = 0; i < numberOfVariables - 1; i++) {
      for (int j = i + 1; j < numberOfVariables; j++) {
        I[i][j] = I[j][i];
      }
    }

    return I;
  }
}
