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

package uk.ac.sussex.gdsc.smlm.function;

/**
 * Class for evaluating a function.
 */
public class StandardGradient1Procedure implements Gradient1Procedure {
  private int index;

  /**
   * The values from the last call to {@link #getValues(Gradient1Function, double[])}.
   */
  public double[] values;
  /**
   * The gradients from the last call to {@link #getValues(Gradient1Function, double[])}.
   */
  public double[][] gradients;

  /**
   * Gets the values.
   *
   * @param function the function
   * @param parameters the function coefficients
   * @return the values
   */
  public double[] getValues(Gradient1Function function, double[] parameters) {
    values = new double[function.size()];
    gradients = new double[values.length][];
    index = 0;
    function.initialise1(parameters);
    function.forEach(this);
    return values;
  }

  @Override
  public void execute(double value, double[] gradient) {
    values[index] = value;
    gradients[index] = gradient.clone();
    index++;
  }
}
