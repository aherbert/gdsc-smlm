package gdsc.smlm.function;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

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
 * This is a wrapper for any function to compute the negative log-likelihood assuming a per-pixel Poisson distribution
 * with Gaussian noise (i.e. the noise from a sCMOS camera). This uses the MLE-sCMOS formula from Huang, et al
 * (2013), Supplementary Notes Eq 3.3:<br/>
 * P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi) = e^-(ui+vari/gi^2) (ui+vari/gi^2)^x / gamma(x+1) <br/>
 * Where:<br/>
 * i = the pixel index <br/>
 * vari = the variance of the pixel <br/>
 * gi = the gain of the pixel <br/>
 * oi = the offset of the pixel <br/>
 * ui = the function value (expected number of photons) <br/>
 * Di = the observed value at the pixel
 * x = the observed random variable (observed number of photons adjusted by a pixel dependent constant) <br/>
 * <p>
 * The negative log-likelihood function is: <br/>
 * -LL(P_sCMOS (x=[(Di-oi)/gi + vari/gi^2]|ui,vari,gi,oi)) <br/>
 * = (ui+vari/gi^2) - x * ln(ui+vari/gi^2) + ln(gamma(x+1)) <br/>
 * <p>
 * The negative log-likelihood (and gradient) can be evaluated over the entire set of observed values or for a chosen
 * observed value.
 * <p>
 * To allow a likelihood to be computed: (a) when the function predicts negative photon count data the function
 * prediction is set to zero; (b) if the observed random variable (x) is negative it is also set to zero. This occurs
 * when true signal readout from the sCMOS camera is low enough to be negated by readout noise. In this case the noise
 * can be ignored.
 * 
 * @see Hunag, et al (2013) Video-rate nanoscopy using sCMOS camera–specific single-molecule localization algorithms.
 *      Nature Methods 10, 653–658.
 */
public class SCMOSLikelihoodWrapper extends LikelihoodWrapper
{
	private final double logNormalisation;
	//private final float[] var, g, o;
	private final double[] var_g2, x, logG;

	/**
	 * Initialise the function.
	 * <p>
	 * The input parameters must be the full parameters for the non-linear function. Only those parameters with gradient
	 * indices should be passed in to the functions to obtain the value (and gradient).
	 *
	 * @param f
	 *            The function to be used to calculated the expected values (Note that the expected value is the number
	 *            of photons)
	 * @param a
	 *            The initial parameters for the function
	 * @param k
	 *            The observed values from the sCMOS camera
	 * @param n
	 *            The number of observed values
	 * @param var
	 *            the variance of each pixel
	 * @param g
	 *            the gain of each pixel
	 * @param o
	 *            the offset of each pixel
	 * @throws IllegalArgumentException
	 *             if the input observed values are not integers
	 */
	public SCMOSLikelihoodWrapper(NonLinearFunction f, double[] a, double[] k, int n, float[] var, float[] g, float[] o)
	{
		super(f, a, k, n);

		//this.var = var;
		//this.g = g;
		//this.o = o;

		var_g2 = new double[n];
		x = new double[n];
		logG = new double[n];

		// Pre-compute the sum over the data
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			var_g2[i] = var[i] / (g[i] * g[i]);
			x[i] = FastMath.max(0, (k[i] - o[i]) / g[i] + var_g2[i]);
			logG[i] = Math.log(g[i]);

			sum += logGamma1(x[i]) + logG[i];
		}

		logNormalisation = sum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood()
	 */
	public double computeLikelihood()
	{
		// Compute the negative log-likelihood to be minimised:
		// (ui+vari/gi^2) - x * ln(ui+vari/gi^2) + ln(gamma(x+1))
		double ll = 0;
		for (int i = 0; i < n; i++)
		{
			double u = f.eval(i);

			if (u < 0)
				u = 0;

			double l = u + var_g2[i];

			ll += l - x[i] * Math.log(l);
		}
		return ll + logNormalisation;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood(double[])
	 */
	public double computeLikelihood(double[] gradient)
	{
		// Compute the negative log-likelihood to be minimised
		// (ui+vari/gi^2) - x * ln(ui+vari/gi^2) + ln(gamma(x+1))
		// 

		// TODO - check this is still valid

		// f(x) = l(x) - k * ln(l(x)) + log(gamma(k+1))
		// 
		// Since (k * ln(l(x)))' = (k * ln(l(x))') * l'(x)
		//                       = (k / l(x)) * l'(x)

		// f'(x) = l'(x) - (k/l(x) * l'(x))
		// f'(x) = l'(x) * (1 - k/l(x))

		double ll = 0;
		for (int j = 0; j < nVariables; j++)
			gradient[j] = 0;
		double[] dl_da = new double[nVariables];
		for (int i = 0; i < n; i++)
		{
			double u = f.eval(i, dl_da);

			if (u < 0)
				u = 0;

			double l = u + var_g2[i];

			ll += l - x[i] * Math.log(l);

			// Note: if l==0 then we get divide by zero and a NaN value
			final double factor = (1 - x[i] / l);
			for (int j = 0; j < gradient.length; j++)
			{
				gradient[j] += dl_da[j] * factor;
			}
		}
		return ll + logNormalisation;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood(int)
	 */
	public double computeLikelihood(int i)
	{
		double u = f.eval(i);

		if (u < 0)
			u = 0;

		double l = u + var_g2[i];

		return l - x[i] * Math.log(l) + logGamma1(x[i]) + logG[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#computeLikelihood(double[], int)
	 */
	public double computeLikelihood(double[] gradient, int i)
	{
		for (int j = 0; j < nVariables; j++)
			gradient[j] = 0;
		double[] dl_da = new double[nVariables];

		double u = f.eval(i, dl_da);

		if (u < 0)
			u = 0;

		double l = u + var_g2[i];

		final double factor = (1 - x[i] / l);
		for (int j = 0; j < gradient.length; j++)
		{
			gradient[j] = dl_da[j] * factor;
		}

		return l - x[i] * Math.log(l) + logGamma1(x[i]) + logG[i];
	}

	private static double logGamma1(double k)
	{
		if (k <= 1)
			return 0;
		return Gamma.logGamma(k + 1);
	}

	/**
	 * Compute the negative log likelihood.
	 *
	 * @param u
	 *            the expected value (number of photoelectrons)
	 * @param var
	 *            the variance of the pixel
	 * @param g
	 *            the gain of the pixel
	 * @param o
	 *            the offset of the pixel
	 * @param k
	 *            The observed value (count from the sCMOS pixel)
	 * @return the negative log likelihood
	 */
	public static double negativeLogLikelihood(double u, float var, float g, float o, double k)
	{
		double var_g2 = var / (g * g);
		double x = Math.max(0, (k - o) / g + var_g2);
		if (u < 0)
			u = 0;
		double l = u + var_g2;
		// Note we need the Math.log(g) to normalise the Poisson distribution to 1 
		// since the observed values (k) are scaled by the gain
		return l - x * Math.log(l) + logGamma1(x) + Math.log(g);
	}

	/**
	 * Compute the likelihood.
	 *
	 * @param u
	 *            the expected value (number of photoelectrons)
	 * @param var
	 *            the variance of the pixel
	 * @param g
	 *            the gain of the pixel
	 * @param o
	 *            the offset of the pixel
	 * @param k
	 *            The observed value (count from the sCMOS pixel)
	 * @return the likelihood
	 */
	public static double likelihood(double u, float var, float g, float o, double k)
	{
		//double var_g2 = var / (g * g);
		//double x = FastMath.max(0, (k - o) / g + var_g2);
		//double l = u + var_g2;
		//double v = FastMath.exp(-l) * Math.pow(l, x) / gamma1(x);
		//if (v != v)
		//	throw new RuntimeException("Failed computation");
		//return v;

		final double nll = negativeLogLikelihood(u, var, g, o, k);
		return FastMath.exp(-nll);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodWrapper#canComputeGradient()
	 */
	@Override
	public boolean canComputeGradient()
	{
		return true;
	}
}