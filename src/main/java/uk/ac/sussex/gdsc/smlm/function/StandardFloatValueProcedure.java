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
public class StandardFloatValueProcedure implements ValueProcedure {
  private int index;
  /**
   * The values from the last call to {@link #getValues(ValueFunction, double[])}.
   */
  public float[] values;

  /**
   * Gets the values.
   *
   * @param function the function
   * @param parameters the function coefficients
   * @return the values
   */
  public float[] getValues(ValueFunction function, double[] parameters) {
    values = new float[function.size()];
    index = 0;
    function.initialise0(parameters);
    function.forEach(this);
    return values;
  }

  /**
   * Gets the values into the buffer.
   *
   * @param function the function
   * @param parameters the function coefficients
   * @param buffer the buffer
   * @param offset the offset
   */
  public void getValues(ValueFunction function, double[] parameters, float[] buffer, int offset) {
    if (buffer == null || buffer.length < offset + function.size()) {
      throw new IllegalArgumentException("Buffer is not large enough for the function values");
    }
    values = buffer;
    index = offset;
    function.initialise0(parameters);
    function.forEach(this);
  }

  /** {@inheritDoc} */
  @Override
  public void execute(double value) {
    values[index++] = (float) value;
  }
}
