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

package uk.ac.sussex.gdsc.smlm.utils;

/**
 * A generic triplet.
 *
 * @author Alex Herbert
 * @param <A> the generic type
 * @param <B> the generic type
 * @param <C> the generic type
 */
public class Triplet<A, B, C> {
  /** First item of the triplet. */
  public final A a;
  /** Second item of the triplet. */
  public final B b;
  /** Third item of the triplet. */
  public final C c;

  /**
   * Instantiates a new triplet.
   *
   * @param a the first item of the triplet
   * @param b the second item of the triplet
   * @param c the third item of the triplet
   */
  public Triplet(A a, B b, C c) {
    this.a = a;
    this.b = b;
    this.c = c;
  }
}
