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

package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Wrap a function solver to scale the function value. Parameters that are scaled must be provided
 * in the constructor. It is assumed that a linear scale can be applied to all these parameters with
 * the effect that the output function value is reduced by the scale factor. This will be the case
 * for offset parameters (if the offset is from zero) and magnitude parameters. An example is
 * y=m*x+c can be scaled to a*y=a*m*x+a*c. Interface methods that compute the function value are
 * rescaled after computation.
 */
public class ScaledFunctionSolver extends WrappedFunctionSolver {
  /** The up scale. */
  protected final double upScale;

  /** The down scale. */
  protected final double downScale;

  /** The indices. */
  protected final int[] indices;

  /**
   * Instantiates a new scaled function solver.
   *
   * <p>Indexed parameters are up-scaled prior to calling the inner function solver. Output
   * parameters, deviations and the function value are are down-scaled upon completion.
   *
   * @param solver the solver
   * @param scale the scale
   * @param indices the indices of the parameters to scale
   */
  public ScaledFunctionSolver(FunctionSolver solver, double scale, int[] indices) {
    super(solver);
    this.upScale = scale;
    this.downScale = 1.0 / scale;
    this.indices = indices;
  }

  @Override
  public FitStatus fit(double[] y, double[] f, double[] a, double[] aDev) {
    // Do not break the view that the solver has on the data
    final double[] f2 = (f == null) ? null : new double[f.length];
    final double[] a2 = cloneAndScaleParameters(a, upScale);
    final double[] aDev2 = (aDev == null) ? null : new double[aDev.length];
    final FitStatus result = solver.fit(y, f2, a2, aDev2);
    scaleFunctionValue(f2, f, downScale);
    scaleParameters(a2, a, downScale);
    scaleDeviations(aDev2, aDev, downScale);
    return result;
  }

  @Override
  public void setBounds(double[] lower, double[] upper) {
    if (lower != null) {
      lower = cloneAndScaleParameters(lower, upScale);
    }
    if (upper != null) {
      upper = cloneAndScaleParameters(upper, upScale);
    }
    solver.setBounds(lower, upper);
  }

  @Override
  public void setConstraints(double[] lower, double[] upper) {
    if (lower != null) {
      lower = cloneAndScaleParameters(lower, upScale);
    }
    if (upper != null) {
      upper = cloneAndScaleParameters(upper, upScale);
    }
    solver.setConstraints(lower, upper);
  }

  @Override
  public boolean evaluate(double[] y, double[] f, double[] a) {
    // Do not break the view that the solver has on the data
    final double[] f2 = (f == null) ? null : new double[f.length];
    final double[] a2 = cloneAndScaleParameters(a, upScale);
    final boolean result = solver.evaluate(y, f2, a2);
    if (result) {
      scaleFunctionValue(f2, f, downScale);
      scaleParameters(a2, a, downScale);
    }
    return result;
  }

  @Override
  public boolean computeDeviations(double[] y, double[] a, double[] aDev) {
    // Do not break the view that the solver has on the data
    final double[] a2 = cloneAndScaleParameters(a, upScale);
    final double[] aDev2 = new double[aDev.length]; // Assume aDev is not null
    final boolean result = solver.computeDeviations(y, a2, aDev2);
    if (result) {
      scaleDeviations(aDev2, aDev, downScale);
    }
    return result;
  }

  /**
   * {@inheritDoc}
   *
   * <p>This value is NOT scaled. It is assumed that the caller requires the value of the solution
   * unmodified.
   *
   * @see uk.ac.sussex.gdsc.smlm.fitting.WrappedFunctionSolver#getValue()
   */
  @Override
  public double getValue() {
    return super.getValue();
  }

  private void scaleParameters(double[] in, double[] out, double scale) {
    // Copy back all fitted parameters
    System.arraycopy(in, 0, out, 0, in.length);
    for (int i = indices.length; i-- > 0;) {
      out[indices[i]] = in[indices[i]] * scale;
    }
  }

  private void scaleDeviations(double[] in, double[] out, double scale) {
    // Only on output so check if null
    if (in == null) {
      return;
    }
    // Copy back all fitted parameters
    System.arraycopy(in, 0, out, 0, in.length);
    // Deviations are variances so multiply by the square of the scale
    scale *= scale;
    for (int i = indices.length; i-- > 0;) {
      out[indices[i]] = in[indices[i]] * scale;
    }
  }

  private double[] cloneAndScaleParameters(double[] a, double scale) {
    final double[] out = a.clone();
    scaleParameters(a, out, scale);
    return out;
  }

  private static void scaleFunctionValue(double[] in, double[] out, double scale) {
    // Only on output so check if null
    if (in == null) {
      return;
    }
    for (int i = in.length; i-- > 0;) {
      out[i] = in[i] * scale;
    }
  }
}
