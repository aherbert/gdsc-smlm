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

package uk.ac.sussex.gdsc.smlm.ij.utils;

import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

@SuppressWarnings({"javadoc"})
public class DoubleImage2DTest extends Image2DTest {
  @Override
  protected DoubleImage2D createData(int w, int h) {
    final double[] data = SimpleArrayUtils.newArray(w * h, 1.0, 1.0);
    return new DoubleImage2D(w, h, data);
  }

  @Override
  protected DoubleImage2D createEmptyData(int w, int h) {
    return new DoubleImage2D(w, h);
  }
}
