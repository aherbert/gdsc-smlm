package gdsc.smlm.function;

import java.util.Arrays;

import uk.ac.sussex.gdsc.core.utils.NotImplementedException;

/**
 * This is a wrapper for any function to compute the negative log-likelihood
 */
public abstract class LikelihoodWrapper
{
	final protected NonLinearFunction f;
	final protected double[] a, data;
	final protected int n;
	final protected int nVariables;

	private double lastScore;
	private double[] lastVariables;

	/**
	 * Initialise the function.
	 * <p>
	 * The input parameters must be the full parameters for the non-linear function. Only those parameters with gradient
	 * indices should be passed in to the functions to obtain the value (and gradient).
	 * 
	 * @param f
	 *            The function to be used to calculated the expected values
	 * @param a
	 *            The initial parameters for the function
	 * @param k
	 *            The observed values
	 * @param n
	 *            The number of observed values
	 */
	public LikelihoodWrapper(NonLinearFunction f, double[] a, double[] k, int n)
	{
		this.f = f;
		this.a = Arrays.copyOf(a, a.length);
		this.data = k;
		this.n = n;
		nVariables = f.gradientIndices().length;
	}

	/**
	 * Copy the variables into the appropriate parameter positions for the NonLinearFunction
	 * 
	 * @param variables
	 */
	private void initialiseFunction(double[] variables)
	{
		final int[] gradientIndices = f.gradientIndices();
		for (int i = 0; i < gradientIndices.length; i++)
			a[gradientIndices[i]] = variables[i];
		f.initialise(a);
	}

	/**
	 * Check if the variable match those last used for computation of the value
	 * 
	 * @param variables
	 * @return True if the variables are the same
	 */
	private boolean sameVariables(double[] variables)
	{
		if (lastVariables != null)
		{
			for (int i = 0; i < variables.length; i++)
				if (variables[i] != lastVariables[i])
					return false;
			return true;
		}
		return false;
	}

	/**
	 * Compute the likelihood score. Returns positive infinity if the likelihood is zero at any point in the observed
	 * values.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @return The negative log likelihood
	 */
	public double likelihood(double[] variables)
	{
		// Check if we have a cached score
		if (sameVariables(variables))
			return lastScore;

		initialiseFunction(variables);
		lastVariables = variables.clone();
		return lastScore = computeLikelihood();
	}

	/**
	 * Compute the likelihood score. Returns positive infinity if the likelihood is zero at any point in the observed
	 * values.
	 * <p>
	 * The wrapped NonLinearFunction will be correctly initialised before this function is called.
	 * 
	 * @return The negative log likelihood
	 */
	protected abstract double computeLikelihood();

	/**
	 * Compute the likelihood score and the gradient. Returns positive infinity if the likelihood is zero
	 * at any point in the observed values. In this case the gradient computed is invalid.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @param gradient
	 *            The gradient (must be equal length to the variables array)
	 * @return The negative log likelihood
	 * @throws NotImplementedException
	 *             If the sub-class cannot provide gradients
	 */
	public double likelihood(double[] variables, double[] gradient) throws NotImplementedException
	{
		initialiseFunction(variables);
		lastVariables = variables.clone();
		return lastScore = computeLikelihood(gradient);
	}

	/**
	 * Compute the likelihood score and the gradient. Returns positive infinity if the likelihood is zero
	 * at any point in the observed values. In this case the gradient computed is invalid.
	 * <p>
	 * The wrapped NonLinearFunction will be correctly initialised before this function is called
	 * 
	 * @param gradient
	 *            The gradient (must be equal length to the variables array)
	 * @return The negative log likelihood
	 * @throws NotImplementedException
	 *             If the sub-class cannot provide gradients
	 */
	protected double computeLikelihood(double[] gradient) throws NotImplementedException
	{
		throw new NotImplementedException();
	}

	/**
	 * Compute the likelihood score at observed value i. Returns positive infinity if the likelihood is zero at the
	 * observed value.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @param i
	 *            Observed value i
	 * @return The negative log likelihood
	 */
	public double likelihood(double[] variables, int i)
	{
		initialiseFunction(variables);
		return computeLikelihood(i);
	}

	/**
	 * Compute the likelihood score at observed value i. Returns positive infinity if the likelihood is zero at the
	 * observed value.
	 * <p>
	 * The wrapped NonLinearFunction will be correctly initialised before this function is called
	 * 
	 * @param i
	 *            Observed value i
	 * @return The negative log likelihood
	 */
	protected abstract double computeLikelihood(int i);

	/**
	 * Compute the likelihood score and gradient of the function at observed value i. Returns positive infinity if the
	 * likelihood is zero at the observed value. In this case the gradient computed will be invalid.
	 * 
	 * @param variables
	 *            The variables of the function
	 * @param gradient
	 *            The gradient (must be equal length to the variables array)
	 * @param i
	 *            Observed value i
	 * @return The negative log likelihood
	 * @throws NotImplementedException
	 *             If the sub-class cannot provide gradients
	 */
	public double likelihood(double[] variables, double[] gradient, int i) throws NotImplementedException
	{
		initialiseFunction(variables);
		return computeLikelihood(gradient, i);
	}

	/**
	 * Compute the likelihood score and gradient of the function at observed value i. Returns positive infinity if the
	 * likelihood is zero at the observed value. In this case the gradient computed will be invalid.
	 * <p>
	 * The wrapped NonLinearFunction will be correctly initialised before this function is called
	 * 
	 * @param gradient
	 *            The gradient (must be equal length to the variables array)
	 * @param i
	 *            Observed value i
	 * @return The negative log likelihood
	 * @throws NotImplementedException
	 *             If the sub-class cannot provide gradients
	 */
	protected double computeLikelihood(double[] gradient, int i) throws NotImplementedException
	{
		throw new NotImplementedException();
	}

	/**
	 * Specify if the likelihood function can compute gradients. If false then the calls to the likelihood functions to
	 * compute the gradient will throw a {@link gdsc.core.utils.NotImplementedException }
	 * 
	 * @return True if the likelihood function can compute gradients
	 */
	public abstract boolean canComputeGradient();
}