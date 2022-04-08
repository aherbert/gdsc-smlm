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

package uk.ac.sussex.gdsc.smlm.function;

/**
 * Implementation of the ICSILog algorithm as described in O. Vinyals, G. Friedland, N. Mirghafori
 * "Revisiting a basic function on current CPUs: A fast logarithm implementation with adjustable
 * accuracy" (2007).
 *
 * <p>This class is based on the original algorithm description. It can have large errors when the
 * unbiased exponent is zero.
 *
 * @see <a href=
 *      "http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf">http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf</a>
 */
public final class IcsiFastLog extends FastLog {
  /**
   * Specify the data type that will be used in the call to the log functions.
   */
  public enum DataType {
    /** float. */
    FLOAT,
    /** double. */
    DOUBLE,
    /** float and double. */
    BOTH
  }

  /** The number of bits to remove from a float mantissa. */
  // CHECKSTYLE.OFF: MemberName
  private final int q;
  // CHECKSTYLE.ON: MemberName
  /** The number of bits to remove from a double mantissa. */
  private final int qd;
  /**
   * The table of the log2 value of binary number 1.0000... to 1.1111...., depending on the
   * precision. The table has had the float bias (-127) pre-subtracted.
   */
  private final float[] data;
  /**
   * The table of the log2 value of binary number 1.0000... to 1.1111...., depending on the
   * precision. The table has had the double bias (-1023) pre-subtracted.
   */
  private final float[] ddata;

  /**
   * Create a new natural logarithm calculation instance using the default table size.
   *
   * <p>If the table is not initialised then a call to a log function with the particular data type
   * will throw a null pointer exception.
   *
   * @param dataType the data type
   * @return the fast log instance
   * @throws IllegalArgumentException if n is not in the range 0-23, or no datatype is specified
   */
  public static IcsiFastLog create(DataType dataType) {
    return create(N, dataType);
  }

  /**
   * Create a new natural logarithm calculation instance. This will hold the pre-calculated log
   * values for base 2 and a table size depending on a given mantissa precision.
   *
   * <p>If the table is not initialised then a call to a log function with the particular data type
   * will throw a null pointer exception.
   *
   * @param n The number of bits to keep from the mantissa. Table storage = 2^n * 4 bytes, e.g. 32Kb
   *        for n=13.
   * @param dataType the data type
   * @return the fast log instance
   * @throws IllegalArgumentException if n is not in the range 0-23, or no datatype is specified
   */
  public static IcsiFastLog create(int n, DataType dataType) {
    if (n < 0 || n > 23) {
      throw new IllegalArgumentException("N must be in the range 0<=n<=23");
    }
    if (dataType == null) {
      throw new IllegalArgumentException("Must create for float/double computation");
    }
    return new IcsiFastLog(n, dataType);
  }

  /**
   * Create a new natural logarithm calculation instance. This will hold the pre-calculated log
   * values for base 2 and a table size depending on a given mantissa precision.
   *
   * @param n The number of bits to keep from the mantissa. Table storage = 2^n * 4 bytes, e.g. 32Kb
   *        for n=13.
   * @param dataType the data type
   */
  private IcsiFastLog(int n, DataType dataType) {
    // Store log2 value of a range of floating point numbers using a limited
    // precision mantissa (m). The purpose of this code is to enumerate all
    // possible mantissas of a float with limited precision (23-q). Note the
    // mantissa represents the digits of a binary number after the binary-point: .10101010101.
    // It is assumed that the digit before the point is a 1 if the exponent
    // is non-zero. Otherwise the binary point is moved to the right of the first
    // digit (i.e. a bit shift left).

    // See Float.intBitsToFloat(int):
    // int s = ((bits >> 31) == 0) ? 1 : -1;
    // int e = ((bits >> 23) & 0xff); // Unsigned exponent
    // int m = (e == 0) ?
    // (bits & 0x7fffff) << 1 :
    // (bits & 0x7fffff) | 0x800000;
    //
    // Then the floating-point result equals the value of the mathematical
    // expression s x m x 2^(e-150):
    // e-127 is the unbiased exponent. 23 is the mantissa precision
    // = s x m x 2^(e-127-23)

    // E.g. For a precision of n=(23-q)=6
    // We enumerate:
    // (1.000000 to 1.111111)

    // The mantissa is incremented using an integer representation to allow
    // exact enumeration. This is then converted to a float for the call to
    // exactLog2(double).

    q = 23 - n;
    qd = 52 - n;
    int x = 0x3F800000; // Set the exponent to 0 so the float value=1.0
    // assert Float.intBitsToFloat(x) == 1.0f : "value is not 1.0f";
    final int inc = 1 << q; // Amount to increase the mantissa

    // Allow computation of float/double datatypes
    final int size = 1 << n;
    switch (dataType) {
      case BOTH:
        data = new float[size];
        ddata = new float[size];
        break;
      case DOUBLE:
        data = null;
        ddata = new float[size];
        break;
      case FLOAT:
        data = new float[size];
        ddata = null;
        break;
      default:
        throw new IllegalArgumentException("Unknown datatype: " + dataType);
    }

    for (int i = 0; i < size; i++) {
      // if (i<50 || i>size-50)
      // System.out.println(Integer.toBinaryString(x));

      final float value = Float.intBitsToFloat(x);
      // Subtract the bias here so we don't need to do it in fastLog()
      final double log2 = exactLog2(value);

      // Note: Pre-subtract the exponent bias
      if (data != null) {
        data[i] = (float) (log2 - 127);
      }
      if (ddata != null) {
        ddata[i] = (float) (log2 - 1023);
      }
      x += inc;

      // float log2 = data[i] + 127f;
      // assert uk.ac.sussex.gdsc.core.utils.FloatEquality.almostEqualRelativeOrAbsolute(log2,
      // fastLog2(value), 1e-6f,
      // 1e-16f) : String.format("[%d] data[i](%g) %g != %g %g", i, value, log2, fastLog2(value),
      // uk.ac.sussex.gdsc.core.utils.FloatEquality.relativeError(log2, fastLog2(value)));
    }
  }

  @Override
  public int getN() {
    return 23 - q;
  }

  @Override
  public double getScale() {
    return LN2;
  }

  @Override
  public double getBase() {
    return Math.E;
  }

  @Override
  public float log2(float x) {
    final int bits = Float.floatToRawIntBits(x);

    // Note the documentation from Float.intBitsToFloat(int):
    // int s = ((bits >> 31) == 0) ? 1 : -1;
    // int e = ((bits >> 23) & 0xff);
    // int m = (e == 0) ?
    // (bits & 0x7fffff) << 1 :
    // (bits & 0x7fffff) | 0x800000;
    // expression s x m x 2^(e-150):
    // e-127 is the unbiased exponent. 23 is the mantissa precision
    // = s x m x 2^(e-127-23)

    // Get the biased exponent
    final int e = (bits >> 23) & 0xff;
    // Mantissa
    final int m = (bits & 0x7fffff);

    if (e == 255) {
      // All bits set is a special case
      if (m != 0) {
        return Float.NaN;
      }
      // +/- Infinity
      return ((bits >> 31) != 0) ? Float.NaN : Float.POSITIVE_INFINITY;
    }

    // Check for negatives
    if ((bits >> 31) != 0) {
      // Only -0 is allowed
      return (e == 0 && m == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;
    }

    if (e == 0) {
      return (m == 0)
          // Special case for +0
          ? Float.NEGATIVE_INFINITY
          : (data[m >>> q]);
    }

    return (e + data[m >>> q]);
  }

  @Override
  public float log2(double x) {
    final long bits = Double.doubleToRawLongBits(x);

    // Note the documentation from Double.longBitsToDouble(int):
    // int s = ((bits >> 63) == 0) ? 1 : -1;
    // int e = (int)((bits >>> 52) & 0x7ffL);
    // long m = (e == 0) ?
    // (bits & 0xfffffffffffffL) << 1 :
    // (bits & 0xfffffffffffffL) | 0x10000000000000L;
    // Then the floating-point result equals the value of the mathematical
    // expression s x m x 2^(e-1075):
    // e-1023 is the unbiased exponent. 52 is the mantissa precision
    // = s x m x 2^(e-1023-52)

    // Get the biased exponent
    final int e = (int) ((bits >>> 52) & 0x7ffL);
    // Mantissa
    final long m = (bits & 0xfffffffffffffL);

    if (e == 2047) {
      // All bits set is a special case
      if (m != 0) {
        return Float.NaN;
      }
      // +/- Infinity
      return ((bits >> 63) != 0L) ? Float.NaN : Float.POSITIVE_INFINITY;
    }

    // Check for negatives
    if ((bits >> 63) != 0L) {
      // Only -0 is allowed
      return (e == 0 && m == 0L) ? Float.NEGATIVE_INFINITY : Float.NaN;
    }

    if (e == 0) {
      return (m == 0L)
          // Special case for +0
          ? Float.NEGATIVE_INFINITY
          : (ddata[(int) (m >>> qd)]);
    }

    return (e + ddata[(int) (m >>> qd)]);
  }

  /**
   * Calculate the logarithm using base 2. Requires the argument be finite and positive.
   *
   * <p>Special cases: <ul> <li>If the argument is NaN, then the result is incorrect
   * ({@code >fastLog2(Float.MAX_VALUE)}). <li>If the argument is negative, then the result is
   * incorrect ({@code fastLog2(-x)}). <li>If the argument is positive infinity, then the result is
   * incorrect ({@code fastLog2(Float.MAX_VALUE)}). <li>If the argument is positive zero or negative
   * zero, then the result is incorrect ({@code fastLog2(Float.MIN_VALUE)}). </ul>
   *
   * @param x the argument (must be strictly positive)
   * @return log2(x)
   */
  @Override
  public float fastLog2(float x) {
    final int bits = Float.floatToRawIntBits(x);
    final int e = ((bits >> 23) & 0xff);
    final int m = (bits & 0x7fffff);
    return (e + data[m >>> q]);
  }

  /**
   * Calculate the logarithm using base 2. Requires the argument be finite and positive.
   *
   * <p>Special cases: <ul> <li>If the argument is NaN, then the result is incorrect
   * ({@code >fastLog2(Float.MAX_VALUE)}). <li>If the argument is negative, then the result is
   * incorrect ({@code fastLog2(-x)}). <li>If the argument is positive infinity, then the result is
   * incorrect ({@code fastLog2(Float.MAX_VALUE)}). <li>If the argument is positive zero or negative
   * zero, then the result is incorrect ({@code fastLog2(Float.MIN_VALUE)}). </ul>
   *
   * @param x the argument (must be strictly positive)
   * @return log(x)
   */
  @Override
  public float fastLog2(double x) {
    final long bits = Double.doubleToRawLongBits(x);
    final int e = (int) ((bits >>> 52) & 0x7ffL);
    final long m = (bits & 0xfffffffffffffL);
    return (e + ddata[(int) (m >>> qd)]);
  }

  @Override
  public float log(float x) {
    // Re-implement to avoid multiplication for all the edge cases
    final int bits = Float.floatToRawIntBits(x);
    final int e = (bits >> 23) & 0xff;
    final int m = (bits & 0x7fffff);

    if (e == 255) {
      if (m != 0) {
        return Float.NaN;
      }
      return ((bits >> 31) != 0) ? Float.NaN : Float.POSITIVE_INFINITY;
    }

    if ((bits >> 31) != 0) {
      return (e == 0 && m == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;
    }

    if (e == 0) {
      return (m == 0) ? Float.NEGATIVE_INFINITY : (data[m >>> q]) * LN2F;
    }

    return (e + data[m >>> q]) * LN2F;
  }

  @Override
  public float log(double x) {
    final long bits = Double.doubleToRawLongBits(x);
    final int e = (int) ((bits >>> 52) & 0x7ffL);
    final long m = (bits & 0xfffffffffffffL);

    if (e == 2047) {
      if (m != 0L) {
        return Float.NaN;
      }
      return ((bits >> 63) != 0L) ? Float.NaN : Float.POSITIVE_INFINITY;
    }

    if ((bits >> 63) != 0L) {
      return (e == 0 && m == 0L) ? Float.NEGATIVE_INFINITY : Float.NaN;
    }

    if (e == 0) {
      return (m == 0L) ? Float.NEGATIVE_INFINITY : (ddata[(int) (m >>> qd)]) * LN2F;
    }

    return (e + ddata[(int) (m >>> qd)]) * LN2F;
  }

  /**
   * Calculate the natural logarithm. Requires the argument be finite and positive.
   *
   * <p>Special cases: <ul> <li>If the argument is NaN, then the result is incorrect
   * ({@code >fastLog(Float.MAX_VALUE)}). <li>If the argument is negative, then the result is
   * incorrect ({@code fastLog(-x)}). <li>If the argument is positive infinity, then the result is
   * incorrect ({@code fastLog(Float.MAX_VALUE)}). <li>If the argument is positive zero or negative
   * zero, then the result is incorrect ({@code fastLog(Float.MIN_VALUE)}). </ul>
   *
   * @param x the argument (must be strictly positive)
   * @return log(x)
   */
  @Override
  public float fastLog(float x) {
    final int bits = Float.floatToRawIntBits(x);
    final int e = ((bits >> 23) & 0xff);
    final int m = (bits & 0x7fffff);
    return (e + data[m >>> q]) * LN2F;
  }

  /**
   * Calculate the natural logarithm. Requires the argument be finite and positive.
   *
   * <p>Special cases: <ul> <li>If the argument is NaN, then the result is incorrect
   * ({@code >fastLog(Float.MAX_VALUE)}). <li>If the argument is negative, then the result is
   * incorrect ({@code fastLog(-x)}). <li>If the argument is positive infinity, then the result is
   * incorrect ({@code fastLog(Float.MAX_VALUE)}). <li>If the argument is positive zero or negative
   * zero, then the result is incorrect ({@code fastLog(Float.MIN_VALUE)}). </ul>
   *
   * @param x the argument (must be strictly positive)
   * @return log(x)
   */
  @Override
  public float fastLog(double x) {
    final long bits = Double.doubleToRawLongBits(x);
    final int e = (int) ((bits >>> 52) & 0x7ffL);
    final long m = (bits & 0xfffffffffffffL);
    return (e + ddata[(int) (m >>> qd)]) * LN2F;
  }
}
