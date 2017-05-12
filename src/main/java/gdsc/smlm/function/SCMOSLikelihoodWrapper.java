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

	/**
	 * Compute variance divided by the gain squared. This can be used in the
	 * {@link #SCMOSLikelihoodWrapper(NonLinearFunction, double[], double[], int, double[], double[])} constructor.
	 *
	 * @param var
	 *            the variance of each pixel
	 * @param g
	 *            the gain of each pixel
	 * @return the variance divided by the gain squared
	 */
	public static double[] compute_var_g2(float[] var, float[] g)
	{
		int n = Math.min(var.length, g.length);
		double[] var_g2 = new double[n];
		for (int i = 0; i < n; i++)
			var_g2[i] = var[i] / (g[i] * g[i]);
		return var_g2;
	}

	/**
	 * Compute log of the gain.
	 *
	 * @param g
	 *            the gain of each pixel
	 * @return the log of the gain
	 */
	public static double[] compute_logG(float[] g)
	{
		int n = g.length;
		double[] logG = new double[n];
		for (int i = 0; i < n; i++)
			logG[i] = Math.log(g[i]);
		return logG;
	}

	/**
	 * Compute X (the mapped observed values from the sCMOS camera)
	 *
	 * @param k
	 *            The observed values from the sCMOS camera
	 * @param var_g2
	 *            the variance divided by the gain squared
	 * @param g
	 *            the gain of each pixel
	 * @param o
	 *            the offset of each pixel
	 * @return The observed values from the sCMOS camera mapped using [x=max(0, (k-o)/g + var/g^2)] per pixel
	 */
	public static double[] computeX(double[] k, float[] var_g2, float[] g, float[] o)
	{
		int n = k.length;
		double[] x = new double[n];
		for (int i = 0; i < n; i++)
		{
			x[i] = FastMath.max(0, (k[i] - o[i]) / g[i] + var_g2[i]);
		}
		return x;
	}

	/**
	 * Initialise the function using pre-computed per pixel working variables. This allows the pre-computation to be
	 * performed once for the sCMOS pixels for all likelihood computations.
	 * <p>
	 * The input parameters must be the full parameters for the non-linear function. Only those parameters with gradient
	 * indices should be passed in to the functions to obtain the value (and gradient).
	 *
	 * @param f
	 *            The function to be used to calculated the expected values (Note that the expected value is the number
	 *            of photons)
	 * @param a
	 *            The initial parameters for the function
	 * @param x
	 *            The observed values from the sCMOS camera mapped using [x=max(0, (k-o)/g + var/g^2)] per pixel
	 * @param n
	 *            The number of observed values
	 * @param var_g2
	 *            the variance of each pixel divided by the gain squared
	 * @param logG
	 *            the log of the gain of each pixel
	 * @throws IllegalArgumentException
	 *             if the input observed values are not integers
	 */
	public SCMOSLikelihoodWrapper(NonLinearFunction f, double[] a, double[] x, int n, double[] var_g2, double[] logG)
	{
		super(f, a, x, n);

		this.var_g2 = var_g2;
		this.x = x;
		this.logG = logG;

		// Pre-compute the sum over the data
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += logGamma1(x[i]) + logG[i];
		}

		logNormalisation = sum;
	}

	/**
	 * Copy constructor
	 *
	 * @param f
	 *            The function to be used to calculated the expected values (Note that the expected value is the number
	 *            of photons)
	 * @param a
	 *            The initial parameters for the function
	 * @param x
	 *            The observed values from the sCMOS camera mapped using [x=max(0, (k-o)/g + var/g^2)] per pixel
	 * @param n
	 *            The number of observed values
	 * @param var_g2
	 *            the variance of each pixel divided by the gain squared
	 * @param logG
	 *            the log of the gain of each pixel
	 * @param logNormalisation
	 *            the log normalisation
	 * @throws IllegalArgumentException
	 *             if the input observed values are not integers
	 */
	private SCMOSLikelihoodWrapper(NonLinearFunction f, double[] a, double[] x, int n, double[] var_g2, double[] logG,
			double logNormalisation)
	{
		super(f, a, x, n);
		this.var_g2 = var_g2;
		this.x = x;
		this.logG = logG;
		this.logNormalisation = logNormalisation;
	}

	/**
	 * Builds a new instance with a new function. All pre-computation on the data is maintained.
	 *
	 * @param f
	 *            The function to be used to calculated the expected values (Note that the expected value is the number
	 *            of photons)
	 * @param a
	 *            The initial parameters for the function
	 * @return the SCMOS likelihood wrapper
	 */
	public SCMOSLikelihoodWrapper build(NonLinearFunction f, double[] a)
	{
		return new SCMOSLikelihoodWrapper(f, a, x, n, var_g2, logG, logNormalisation);
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

			ll += l;
			if (x[i] != 0)
				ll -= x[i] * Math.log(l);
		}
		return ll + logNormalisation;
	}

	private double observedLikelihood = Double.NaN;

	/**
	 * Compute the observed negative log likelihood. This is the value of {@link #computeLikelihood()} if the function
	 * were to return the observed values for each point.
	 *
	 * @return the observed negative log likelihood
	 */
	public double computeObservedLikelihood()
	{
		if (Double.isNaN(observedLikelihood))
		{
			double ll = 0;
			for (int i = 0; i < n; i++)
			{
				// We need to input the observed value as the expected value.
				// So we need (k-o)/g as the expected value. We did not store this so
				// compute it by subtracting var_g2 from x.
				// Then perform the same likelihood computation.

				//double u = x[i] - var_g2[i];
				//
				//if (u < 0)
				//	u = 0;
				//
				//double l = u + var_g2[i];

				// We can do this in one step ...
				double l = (x[i] < var_g2[i]) ? var_g2[i] : x[i];

				ll += l;
				if (x[i] != 0)
					ll -= x[i] * Math.log(l);
			}
			observedLikelihood = ll + logNormalisation;
		}
		return observedLikelihood;
	}

	/**
	 * Compute log likelihood ratio.
	 *
	 * @param ll
	 *            the negative log likelihood of the function
	 * @return the log likelihood ratio
	 */
	public double computeLogLikelihoodRatio(double ll)
	{
		// From https://en.wikipedia.org/wiki/Likelihood-ratio_test#Use:
		// LLR = -2 * [ ln(likelihood for alternative model) - ln(likelihood for null model)]
		// The model with more parameters (here alternative) will always fit at least as well—
		// i.e., have the same or greater log-likelihood—than the model with fewer parameters 
		// (here null)

		double llAlternative = computeObservedLikelihood();
		double llNull = ll;

		// The alternative should always fit better than the null model 
		if (llNull < llAlternative)
			llNull = llAlternative;

		// Since we have negative log likelihood we reverse the sign
		//return 2 * (-llAlternative - -llNull);
		return -2 * (llAlternative - llNull);
	}

	/**
	 * Compute the q-value of the log-likelihood ratio. This is the probability that a value of LLR as poor as the value
	 * should occur by chance.
	 *
	 * @param ll
	 *            the minimum negative log likelihood of the function (the null model)
	 * @return the p-value
	 */
	public double computeQValue(double ll)
	{
		double llr = computeLogLikelihoodRatio(ll);
		int degreesOfFreedom = x.length - nVariables;
		return ChiSquaredDistributionTable.computeQValue(llr, degreesOfFreedom);
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
		// with x as the mapped observed value: x = (k-o)/g + var/g^2

		// To compute the gradient we do the same as for a Poisson distribution:
		// f(x) = l(x) - k * ln(l(x)) + log(gamma(k+1))
		// with l(x) as the Poisson mean (the output dependent on the function variables x)
		// and k the observed value.
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

			ll += l;
			if (x[i] != 0)
				ll -= x[i] * Math.log(l);

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

		double ll = l + logG[i];
		if (x[i] != 0)
			ll += logGamma1(x[i]) - x[i] * Math.log(l);

		return ll;
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

		double ll = l + logG[i];
		if (x[i] != 0)
			ll += logGamma1(x[i]) - x[i] * Math.log(l);

		return ll;
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
		double ll = l + Math.log(g);
		if (x != 0)
			ll += logGamma1(x) - x * Math.log(l);

		return ll;
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

	/**
	 * Compute the Fisher's Information Matrix (I) for fitted variables:
	 * 
	 * <pre>
	 * Iab = sum(k) 1/(uk+vark/gk^2)  (duk da) * (duk db)
	 * </pre>
	 * 
	 * @param variables
	 *            The variables of the function
	 * @return Fisher's Information Matrix (I)
	 */
	@Override
	public double[][] fisherInformation(final double[] variables)
	{
		initialiseFunction(variables);

		double[] du_da = new double[nVariables];

		final double[][] I = new double[nVariables][nVariables];

		for (int k = 0; k < n; k++)
		{
			final double uk = f.eval(k, du_da);
			final double yk = 1 / (uk + var_g2[k]);
			for (int i = 0; i < nVariables; i++)
			{
				double du_dai = yk * du_da[i];
				for (int j = 0; j <= i; j++)
				{
					I[i][j] += du_dai * du_da[j];
				}
			}
		}

		// Generate symmetric matrix
		for (int i = 0; i < nVariables - 1; i++)
			for (int j = i + 1; j < nVariables; j++)
				I[i][j] = I[j][i];

		return I;
	}
}