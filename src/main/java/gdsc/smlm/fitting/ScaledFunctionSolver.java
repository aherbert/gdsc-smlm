package gdsc.smlm.fitting;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Wrap a function solver to scale the function value. Parameters that are scaled must be provided in the constructor.
 * It is assumed that a linear scale can be applied to all these parameters with the effect that the output function
 * value is reduced by the scale factor. This will be the case for offset parameters (if the offset is from zero) and
 * magnitude parameters. An example is y=m*x+c can be scaled to a*y=a*m*x+a*c. Interface methods that compute the
 * function value are rescaled after computation.
 */
public class ScaledFunctionSolver extends WrappedFunctionSolver
{
	protected final double upScale;
	protected final double downScale;
	protected final int[] indices;

	/**
	 * Instantiates a new scaled function solver.
	 * <p>
	 * Indexed parameters are up-scaled prior to calling the inner function solver. Output parameters, deviations and
	 * the function value are are down-scaled upon completion.
	 *
	 * @param solver
	 *            the solver
	 * @param scale
	 *            the scale
	 * @param indices
	 *            the indices of the parameters to scale
	 */
	public ScaledFunctionSolver(FunctionSolver solver, double scale, int[] indices)
	{
		super(solver);
		this.upScale = scale;
		this.downScale = 1.0 / scale;
		this.indices = indices;
	}

	public FitStatus fit(double[] y, double[] f, double[] a, double[] aDev)
	{
		// Do not break the view that the solver has on the data
		double[] f2 = (f == null) ? null : new double[f.length];
		double[] a2 = cloneAndScaleParameters(a, upScale);
		double[] aDev2 = (aDev == null) ? null : new double[aDev.length];
		FitStatus result = solver.fit(y, f2, a2, aDev2);
		scaleFunctionValue(f2, f, downScale);
		scaleParameters(a2, a, downScale);
		scaleDeviations(aDev2, aDev, downScale);
		return result;
	}

	public void setBounds(double[] lower, double[] upper)
	{
		if (lower != null)
			lower = cloneAndScaleParameters(lower, upScale);
		if (upper != null)
			upper = cloneAndScaleParameters(upper, upScale);
		solver.setBounds(lower, upper);
	}

	public void setConstraints(double[] lower, double[] upper)
	{
		if (lower != null)
			lower = cloneAndScaleParameters(lower, upScale);
		if (upper != null)
			upper = cloneAndScaleParameters(upper, upScale);
		solver.setConstraints(lower, upper);
	}

	public boolean evaluate(double[] y, double[] f, double[] a)
	{
		// Do not break the view that the solver has on the data
		double[] f2 = (f == null) ? null : new double[f.length];
		double[] a2 = cloneAndScaleParameters(a, upScale);
		boolean result = solver.evaluate(y, f2, a2);
		if (result)
		{
			scaleFunctionValue(f2, f, downScale);
			scaleParameters(a2, a, downScale);
		}
		return result;
	}

	public boolean computeDeviations(double[] y, double[] a, double[] aDev)
	{
		// Do not break the view that the solver has on the data
		double[] a2 = cloneAndScaleParameters(a, upScale);
		double[] aDev2 = new double[aDev.length]; // Assume aDev is not null
		boolean result = solver.computeDeviations(y, a2, aDev2);
		if (result)
			scaleDeviations(aDev2, aDev, downScale);
		return result;
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * This value is NOT scaled. It is assumed that the caller requires the value of the solution unmodified.
	 * 
	 * @see gdsc.smlm.fitting.WrappedFunctionSolver#getValue()
	 */
	@Override
	public double getValue()
	{
		return super.getValue();
	}

	private void scaleParameters(double[] in, double[] out, double scale)
	{
		for (int i = indices.length; i-- > 0;)
			out[indices[i]] = in[indices[i]] * scale;
	}

	private void scaleDeviations(double[] in, double[] out, double scale)
	{
		// Only on output so check if null
		if (in == null)
			return;
		// Deviations are variances so multiply by the square of the scale
		scale *= scale;
		for (int i = indices.length; i-- > 0;)
			out[indices[i]] = in[indices[i]] * scale;
	}

	private double[] cloneAndScaleParameters(double[] a, double scale)
	{
		double[] out = a.clone();
		scaleParameters(a, out, scale);
		return out;
	}

	private void scaleFunctionValue(double[] in, double[] out, double scale)
	{
		// Only on output so check if null
		if (in == null)
			return;
		for (int i = in.length; i-- > 0;)
			out[i] = in[i] * scale;
	}
}
