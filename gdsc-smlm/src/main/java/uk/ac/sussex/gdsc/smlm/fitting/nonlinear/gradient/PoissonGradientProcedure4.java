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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Calculates the Fisher information matrix for a Poisson process.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class PoissonGradientProcedure4 extends PoissonGradientProcedure {

  /**
   * Instantiates a new procedure.
   *
   * @param func Gradient function
   */
  public PoissonGradientProcedure4(final Gradient1Function func) {
    super(func);
    ValidationUtils.checkArgument(numberOfGradients == 4, "Function must compute 4 gradients");
  }

  @Override
  public void execute(double value, double[] dyDa) {
    if (value > 0) {
      final double function = 1.0 / value;

      data[0] += dyDa[0] * function * dyDa[0];
      double wgt;
      wgt = dyDa[1] * function;
      data[1] += wgt * dyDa[0];
      data[2] += wgt * dyDa[1];
      wgt = dyDa[2] * function;
      data[3] += wgt * dyDa[0];
      data[4] += wgt * dyDa[1];
      data[5] += wgt * dyDa[2];
      wgt = dyDa[3] * function;
      data[6] += wgt * dyDa[0];
      data[7] += wgt * dyDa[1];
      data[8] += wgt * dyDa[2];
      data[9] += wgt * dyDa[3];
    }
  }

  @Override
  protected void initialiseWorkingMatrix() {
    GradientProcedureHelper.initialiseWorkingMatrix4(data);
  }

  @Override
  public void getMatrix(double[][] matrix) {
    GradientProcedureHelper.getMatrix4(data, matrix);
  }

  @Override
  public void getLinear(double[] matrix) {
    GradientProcedureHelper.getMatrix4(data, matrix);
  }
}
