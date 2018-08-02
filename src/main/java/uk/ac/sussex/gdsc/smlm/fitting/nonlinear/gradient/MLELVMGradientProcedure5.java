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
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Calculates the scaled Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the scaled gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * This procedure computes a modified Chi-squared expression to perform Maximum Likelihood Estimation assuming Poisson
 * model. See Laurence &amp; Chromy (2010) Efficient maximum likelihood estimator. Nature Methods 7, 338-339. The input
 * data
 * must be Poisson distributed for this to be relevant.
 */
public class MLELVMGradientProcedure5 extends MLELVMGradientProcedure
{
    /**
     * @param y
     *            Data to fit (must be positive)
     * @param func
     *            Gradient function
     */
    public MLELVMGradientProcedure5(final double[] y, final Gradient1Function func)
    {
        super(y, func);
        if (n != 5)
            throw new IllegalArgumentException("Function must compute 5 gradients");
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
     */
    @Override
    public void execute(double fi, double[] dfi_da)
    {
        ++yi;
        if (fi > 0.0)
        {
            final double xi = y[yi];

            // We assume y[i] is positive but must handle zero
            if (xi > 0.0)
            {
                value += (fi - xi - xi * Math.log(fi / xi));

                final double xi_fi2 = xi / fi / fi;
                final double e = 1 - (xi / fi);

                beta[0] -= e * dfi_da[0];
                beta[1] -= e * dfi_da[1];
                beta[2] -= e * dfi_da[2];
                beta[3] -= e * dfi_da[3];
                beta[4] -= e * dfi_da[4];

                alpha[0] += dfi_da[0] * xi_fi2 * dfi_da[0];
                double w;
                w = dfi_da[1] * xi_fi2;
                alpha[1] += w * dfi_da[0];
                alpha[2] += w * dfi_da[1];
                w = dfi_da[2] * xi_fi2;
                alpha[3] += w * dfi_da[0];
                alpha[4] += w * dfi_da[1];
                alpha[5] += w * dfi_da[2];
                w = dfi_da[3] * xi_fi2;
                alpha[6] += w * dfi_da[0];
                alpha[7] += w * dfi_da[1];
                alpha[8] += w * dfi_da[2];
                alpha[9] += w * dfi_da[3];
                w = dfi_da[4] * xi_fi2;
                alpha[10] += w * dfi_da[0];
                alpha[11] += w * dfi_da[1];
                alpha[12] += w * dfi_da[2];
                alpha[13] += w * dfi_da[3];
                alpha[14] += w * dfi_da[4];
            }
            else
            {
                value += fi;
                beta[0] -= dfi_da[0];
                beta[1] -= dfi_da[1];
                beta[2] -= dfi_da[2];
                beta[3] -= dfi_da[3];
                beta[4] -= dfi_da[4];
            }
        }
    }

    @Override
    protected void initialiseGradient()
    {
        GradientProcedureHelper.initialiseWorkingMatrix5(alpha);
        beta[0] = 0;
        beta[1] = 0;
        beta[2] = 0;
        beta[3] = 0;
        beta[4] = 0;
    }

    @Override
    public void getAlphaMatrix(double[][] alpha)
    {
        GradientProcedureHelper.getMatrix5(this.alpha, alpha);
    }

    @Override
    public void getAlphaLinear(double[] alpha)
    {
        GradientProcedureHelper.getMatrix5(this.alpha, alpha);
    }
}
