/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
 * Implement the {@link FastLog} methods using {@link Math#log(double)}.
 */
public class NonFastLog extends FastLog {
  /** An instance of the class. */
  public static final NonFastLog INSTANCE = new NonFastLog();

  @Override
  public double getBase() {
    return Math.E;
  }

  @Override
  public double getScale() {
    return LN2;
  }

  @Override
  public int getN() {
    return 52;
  }

  @Override
  public float log2(float x) {
    return (float) (Math.log(x) / LN2);
  }

  @Override
  public float log2(double x) {
    return (float) (Math.log(x) / LN2);
  }

  @Override
  public float fastLog2(float x) {
    return (float) (Math.log(x) / LN2);
  }

  @Override
  public float fastLog2(double x) {
    return (float) (Math.log(x) / LN2);
  }

  @Override
  public float log(float x) {
    return (float) Math.log(x);
  }

  @Override
  public float log(double x) {
    return (float) Math.log(x);
  }

  @Override
  public float fastLog(float x) {
    return (float) Math.log(x);
  }

  @Override
  public float fastLog(double x) {
    return (float) Math.log(x);
  }

  @Override
  public double log2D(double x) {
    return Math.log(x) / LN2;
  }

  @Override
  public double fastLog2D(double x) {
    return Math.log(x) / LN2;
  }

  @Override
  public double logD(double x) {
    return Math.log(x);
  }

  @Override
  public double fastLogD(double x) {
    return Math.log(x);
  }
}
