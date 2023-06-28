/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import ij.ImageStack;
import org.jtransforms.dht.FloatDHT_3D;
import pl.edu.icm.jlargearrays.LargeArray;
import uk.ac.sussex.gdsc.core.ij.process.Fht;

/**
 * Wrapper to compute the discrete Hartley transform on 3D data. This uses the JTransforms library.
 */
public class FloatDht3D extends FloatImage3D {
  private boolean isFrequencyDomain;
  private final FloatDHT_3D dht;
  // Used for fast multiply operations
  private double[] h2e;
  private double[] h2o;
  private double[] mag;
  private int[] jj;

  /**
   * Instantiates a new 3D discrete Hartley transform.
   *
   * @param stack the stack
   * @throws IllegalArgumentException If any dimension is less than 2, or if the combined dimensions
   *         is too large for an array
   */
  public FloatDht3D(ImageStack stack) {
    super(stack);
    LargeArray.setMaxSizeOf32bitArray(MAX_SIZE_OF_32_BIT_ARRAY);
    dht = new FloatDHT_3D(ns, nr, nc);
  }

  /**
   * Instantiates a new 3D discrete Hartley transform.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @param data the data
   * @param isFrequencyDomain Set to true if in the frequency domain
   * @throws IllegalArgumentException If any dimension is less than 2, or if the data is not the
   *         correct length
   */
  public FloatDht3D(int nc, int nr, int ns, float[] data, boolean isFrequencyDomain) {
    super(nc, nr, ns, data);
    LargeArray.setMaxSizeOf32bitArray(MAX_SIZE_OF_32_BIT_ARRAY);
    dht = new FloatDHT_3D(ns, nr, nc);
    this.isFrequencyDomain = isFrequencyDomain;
  }

  /**
   * Instantiates a new 3D discrete Hartley transform.
   *
   * @param nc the number of columns
   * @param nr the number of rows
   * @param ns the number of slices
   * @param nrByNc the number of rows multiplied by the number of columns
   * @param data the data
   * @param isFrequencyDomain the is frequency domain
   * @param dht the dht
   */
  private FloatDht3D(int nc, int nr, int ns, int nrByNc, float[] data, boolean isFrequencyDomain,
      FloatDHT_3D dht) {
    super(nc, nr, ns, nrByNc, data);
    this.isFrequencyDomain = isFrequencyDomain;
    this.dht = dht; // This can be reused across objects
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected FloatDht3D(FloatDht3D source) {
    super(source);
    isFrequencyDomain = source.isFrequencyDomain;
    dht = source.dht;
    h2e = source.h2e;
    h2o = source.h2o;
    mag = source.mag;
    jj = source.jj;
  }

  @Override
  public FloatDht3D copy() {
    return new FloatDht3D(this);
  }

  /**
   * Performs a forward transform, converting this image into the frequency domain.
   *
   * @throws IllegalArgumentException If already in the frequency domain
   */
  public void transform() {
    if (isFrequencyDomain) {
      throw new IllegalArgumentException("Already frequency domain DHT");
    }
    dht.forward(data);
    isFrequencyDomain = true;
    resetFastOperations();
  }

  /**
   * Performs an inverse transform, converting this image into the space domain.
   *
   * @throws IllegalArgumentException If already in the space domain
   */
  public void inverseTransform() {
    if (!isFrequencyDomain) {
      throw new IllegalArgumentException("Already space domain DHT");
    }
    dht.inverse(data, true);
    isFrequencyDomain = false;
    resetFastOperations();
  }

  /**
   * Checks if is frequency domain.
   *
   * @return true, if is frequency domain
   */
  public boolean isFrequencyDomain() {
    return isFrequencyDomain;
  }

  /**
   * Initialise fast operations for {@link #multiply(FloatDht3D)} and
   * {@link #conjugateMultiply(FloatDht3D)}. This pre-computes the values needed for the operations.
   *
   * <p>Note: This initialises the FHT object for use as the argument to the operation, for example
   * if a convolution kernel is to be applied to many FHT objects.
   */
  public void initialiseFastMultiply() {
    if (h2e == null) {
      // Do this on new arrays for thread safety (i.e. concurrent initialisation)
      final float[] h2 = getData();
      final double[] lh2e = new double[h2.length];
      final double[] lh2o = new double[lh2e.length];
      final int[] ljj = new int[lh2e.length];
      for (int s = 0, nsMinusS = 0, i = 0; s < ns; s++, nsMinusS = ns - s) {
        for (int r = 0, nrMinusR = 0; r < nr; r++, nrMinusR = nr - r) {
          for (int c = 0, ncMinusC = 0; c < nc; c++, ncMinusC = nc - c, i++) {
            final int j = nsMinusS * nrByNc + nrMinusR * nc + ncMinusC;
            lh2e[i] = ((double) h2[i] + (double) h2[j]) / 2.0;
            lh2o[i] = ((double) h2[i] - (double) h2[j]) / 2.0;
            ljj[i] = j;
          }
        }
      }
      this.h2o = lh2o;
      this.jj = ljj;
      // Assign at the end for thread safety (i.e. concurrent initialisation)
      this.h2e = lh2e;
    }
  }

  /**
   * Initialise fast operations for {@link #multiply(FloatDht3D)},
   * {@link #conjugateMultiply(FloatDht3D)} and {@link #divide(FloatDht3D)}. This pre-computes the
   * values needed for the operations.
   *
   * <p>Note: This initialises the FHT object for use as the argument to the operation, for example
   * if a deconvolution kernel is to be applied to many FHT objects.
   */
  public void initialiseFastOperations() {
    initialiseFastMultiply();
    if (mag == null) {
      // Do this on new arrays for thread safety (i.e. concurrent initialisation)
      final double[] lmag = new double[h2e.length];
      final float[] h2 = getData();
      for (int i = 0; i < h2.length; i++) {
        // Note that pre-computed h2e and h2o are divided by 2 so we also
        // divide the magnitude by 2 to allow reuse of the pre-computed values
        // in the divide operation (which does not require h2e/2 and h2o/2)
        lmag[i] = Math.max(1e-20, h2[i] * h2[i] + h2[jj[i]] * h2[jj[i]]) / 2;
      }
      this.mag = lmag;
    }
  }

  /**
   * Checks if is initialised for fast multiply.
   *
   * @return true, if is fast multiply
   */
  public boolean isFastMultiply() {
    return h2e != null;
  }

  /**
   * Checks if is initialised for fast operations.
   *
   * @return true, if is fast operations
   */
  public boolean isFastOperations() {
    return mag != null;
  }

  private void resetFastOperations() {
    h2e = null;
    h2o = null;
    jj = null;
    mag = null;
  }

  /**
   * Returns the image resulting from the point by point Hartley multiplication of this image and
   * the specified image. Both images are assumed to be in the frequency domain. Multiplication in
   * the frequency domain is equivalent to convolution in the space domain.
   *
   * @param dht the dht
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions
   */
  public FloatDht3D multiply(FloatDht3D dht) {
    return multiply(dht, null);
  }

  /**
   * Returns the image resulting from the point by point Hartley multiplication of this image and
   * the specified image. Both images are assumed to be in the frequency domain. Multiplication in
   * the frequency domain is equivalent to convolution in the space domain.
   *
   * @param dht the dht
   * @param tmp the tmp buffer to use for the result
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions
   */
  public FloatDht3D multiply(FloatDht3D dht, float[] tmp) {
    checkDht(dht);
    return (dht.isFastMultiply()) ? multiply(dht.h2e, dht.h2o, dht.jj, tmp)
        : multiply(dht.getData(), tmp);
  }

  /**
   * Returns the image resulting from the point by point Hartley multiplication of this image and
   * the specified image. Both images are assumed to be in the frequency domain. Multiplication in
   * the frequency domain is equivalent to convolution in the space domain.
   *
   * @param h2 the h 2
   * @param tmp the tmp buffer to use for the result
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions
   */
  private FloatDht3D multiply(float[] h2, float[] tmp) {
    final float[] h1 = this.data;
    if (tmp == null || tmp.length != h1.length) {
      tmp = new float[h1.length];
    }

    for (int s = 0, nsMinusS = 0, i = 0; s < ns; s++, nsMinusS = ns - s) {
      for (int r = 0, nrMinusR = 0; r < nr; r++, nrMinusR = nr - r) {
        for (int c = 0, ncMinusC = 0; c < nc; c++, ncMinusC = nc - c, i++) {
          // This is actually doing for 3D data stored as x[slices][rows][columns]
          // https://en.wikipedia.org/wiki/Discrete_Hartley_transform
          // h2e = (h2[s][r][c] + h2[Ns-s][Nr-r][Nr-c]) / 2;
          // h2o = (h2[s][r][c] - h2[Ns-s][Nr-r][Nr-c]) / 2;
          // tmp[s][r][c] = (float) (h1[s][r][c] * h2e + h1[Ns-s][Nr-r][Nc-c] * h2o);
          final int j = nsMinusS * nrByNc + nrMinusR * nc + ncMinusC;
          final double lh2e = ((double) h2[i] + (double) h2[j]) / 2.0;
          final double lh2o = ((double) h2[i] - (double) h2[j]) / 2.0;
          tmp[i] = (float) (h1[i] * lh2e + h1[j] * lh2o);
        }
      }
    }

    return new FloatDht3D(nc, nr, ns, nrByNc, tmp, true, this.dht);
  }

  /**
   * Returns the image resulting from the point by point Hartley multiplication of this image and
   * the specified image. Both images are assumed to be in the frequency domain. Multiplication in
   * the frequency domain is equivalent to convolution in the space domain.
   *
   * @param h2e the pre-initialised h2e value
   * @param h2o the pre-initialised h2o value
   * @param jj the pre-initialised j index
   * @param tmp the buffer for the result (can be null)
   * @return the result
   */
  private FloatDht3D multiply(double[] h2e, double[] h2o, int[] jj, float[] tmp) {
    final float[] h1 = getData();
    if (tmp == null || tmp.length != h1.length) {
      tmp = new float[h1.length];
    }
    for (int i = 0; i < h1.length; i++) {
      tmp[i] = (float) (h1[i] * h2e[i] + h1[jj[i]] * h2o[i]);
    }
    return new FloatDht3D(nc, nr, ns, nrByNc, tmp, true, this.dht);
  }

  /**
   * Returns the image resulting from the point by point Hartley conjugate multiplication of this
   * image and the specified image. Both images are assumed to be in the frequency domain. Conjugate
   * multiplication in the frequency domain is equivalent to correlation in the space domain.
   *
   * @param dht the dht
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions
   */
  public FloatDht3D conjugateMultiply(FloatDht3D dht) {
    return conjugateMultiply(dht, null);
  }

  /**
   * Returns the image resulting from the point by point Hartley conjugate multiplication of this
   * image and the specified image. Both images are assumed to be in the frequency domain. Conjugate
   * multiplication in the frequency domain is equivalent to correlation in the space domain.
   *
   * @param dht the dht
   * @param tmp the tmp buffer to use for the result
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions
   */
  public FloatDht3D conjugateMultiply(FloatDht3D dht, float[] tmp) {
    checkDht(dht);
    return (dht.isFastMultiply()) ? conjugateMultiply(dht.h2e, dht.h2o, dht.jj, tmp)
        : conjugateMultiply(dht.getData(), tmp);
  }

  /**
   * Returns the image resulting from the point by point Hartley conjugate multiplication of this
   * image and the specified image. Both images are assumed to be in the frequency domain. Conjugate
   * multiplication in the frequency domain is equivalent to correlation in the space domain.
   *
   * @param h2 the h 2
   * @param tmp the tmp buffer to use for the result
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions
   */
  private FloatDht3D conjugateMultiply(float[] h2, float[] tmp) {
    final float[] h1 = this.data;
    if (tmp == null || tmp.length != h1.length) {
      tmp = new float[h1.length];
    }

    for (int s = 0, nsMinusS = 0, i = 0; s < ns; s++, nsMinusS = ns - s) {
      for (int r = 0, nrMinusR = 0; r < nr; r++, nrMinusR = nr - r) {
        for (int c = 0, ncMinusC = 0; c < nc; c++, ncMinusC = nc - c, i++) {
          final int j = nsMinusS * nrByNc + nrMinusR * nc + ncMinusC;
          final double lh2e = ((double) h2[i] + (double) h2[j]) / 2.0;
          final double lh2o = ((double) h2[i] - (double) h2[j]) / 2.0;
          // As per multiply but reverse the addition sign for the conjugate
          tmp[i] = (float) (h1[i] * lh2e - h1[j] * lh2o);
        }
      }
    }

    return new FloatDht3D(nc, nr, ns, nrByNc, tmp, true, this.dht);
  }

  /**
   * Returns the image resulting from the point by point Hartley conjugate multiplication of this
   * image and the specified image. Both images are assumed to be in the frequency domain. Conjugate
   * multiplication in the frequency domain is equivalent to correlation in the space domain.
   *
   * @param h2e the pre-initialised h2e value
   * @param h2o the pre-initialised h2o value
   * @param jj the pre-initialised j index
   * @param tmp the buffer for the result (can be null)
   * @return the fht2
   */
  private FloatDht3D conjugateMultiply(double[] h2e, double[] h2o, int[] jj, float[] tmp) {
    final float[] h1 = getData();
    if (tmp == null || tmp.length != h1.length) {
      tmp = new float[h1.length];
    }
    for (int i = 0; i < h1.length; i++) {
      tmp[i] = (float) (h1[i] * h2e[i] - h1[jj[i]] * h2o[i]);
    }
    return new FloatDht3D(nc, nr, ns, nrByNc, tmp, true, this.dht);
  }

  /**
   * Returns the image resulting from the point by point Hartley division of this image by the
   * specified image. Both images are assumed to be in the frequency domain. Division in the
   * frequency domain is equivalent to deconvolution in the space domain.
   *
   * @param dht the dht
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions or in the frequency
   *         domain
   */
  public FloatDht3D divide(FloatDht3D dht) {
    return divide(dht, null);
  }

  /**
   * Returns the image resulting from the point by point Hartley division of this image by the
   * specified image. Both images are assumed to be in the frequency domain. Division in the
   * frequency domain is equivalent to deconvolution in the space domain.
   *
   * @param dht the dht
   * @param tmp the tmp buffer to use for the result
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions or in the frequency
   *         domain
   */
  public FloatDht3D divide(FloatDht3D dht, float[] tmp) {
    checkDht(dht);
    return (dht.isFastOperations()) ? divide(dht.h2e, dht.h2o, dht.jj, dht.mag, tmp)
        : divide(dht.getData(), tmp);
  }

  /**
   * Returns the image resulting from the point by point Hartley division of this image by the
   * specified image. Both images are assumed to be in the frequency domain. Division in the
   * frequency domain is equivalent to deconvolution in the space domain.
   *
   * @param h2 the h 2
   * @param tmp the tmp buffer to use for the result
   * @return the result
   * @throws IllegalArgumentException if the dht is not the same dimensions or in the frequency
   *         domain
   */
  private FloatDht3D divide(float[] h2, float[] tmp) {
    final float[] h1 = this.data;
    if (tmp == null || tmp.length != h1.length) {
      tmp = new float[h1.length];
    }

    for (int s = 0, nsMinusS = 0, i = 0; s < ns; s++, nsMinusS = ns - s) {
      for (int r = 0, nrMinusR = 0; r < nr; r++, nrMinusR = nr - r) {
        for (int c = 0, ncMinusC = 0; c < nc; c++, ncMinusC = nc - c, i++) {
          // This is a copy of the divide operation in ij.process.FHT
          final int j = nsMinusS * nrByNc + nrMinusR * nc + ncMinusC;
          final double h2i = h2[i];
          final double h2j = h2[j];
          double lmag = h2i * h2i + h2j * h2j;
          if (lmag < 1e-20) {
            lmag = 1e-20;
          }
          final double lh2e = (h2i + h2j);
          final double lh2o = (h2i - h2j);
          tmp[i] = (float) ((h1[i] * lh2e - h1[j] * lh2o) / lmag);
        }
      }
    }

    return new FloatDht3D(nc, nr, ns, nrByNc, tmp, true, this.dht);
  }

  /**
   * Returns the image resulting from the point by point Hartley division of this image by the
   * specified image. Both images are assumed to be in the frequency domain. Division in the
   * frequency domain is equivalent to deconvolution in the space domain.
   *
   * @param h2e the pre-initialised h2e value
   * @param h2o the pre-initialised h2o value
   * @param jj the pre-initialised j index
   * @param mag the pre-initialised magnitude value
   * @param tmp the buffer for the result (can be null)
   */
  private FloatDht3D divide(double[] h2e, double[] h2o, int[] jj, double[] mag, float[] tmp) {
    final float[] h1 = getData();
    if (tmp == null || tmp.length != h1.length) {
      tmp = new float[h1.length];
    }
    for (int i = 0; i < h1.length; i++) {
      tmp[i] = (float) ((h1[i] * h2e[i] - h1[jj[i]] * h2o[i]) / mag[i]);
    }
    return new FloatDht3D(nc, nr, ns, nrByNc, tmp, true, this.dht);
  }

  /**
   * Check the DHT matches the dimensions of this DHT. Check both are in the frequency domain.
   *
   * @param dht the dht
   * @throws IllegalArgumentException If multiplication is not possible
   */
  private void checkDht(FloatDht3D dht) {
    if (dht.ns != ns || dht.nr != nr || dht.nc != nc) {
      throw new IllegalArgumentException("Dimension mismatch");
    }
    DhtHelper.checkFrequencyDomain(dht.isFrequencyDomain);
    DhtHelper.checkFrequencyDomain(isFrequencyDomain);
  }

  /**
   * Converts this DHT to a discrete Fourier transform (DFT) and returns it as a two image pair. The
   * image is assumed to be in the frequency domain.
   *
   * @param real the buffer to use for the real component (can be null)
   * @param imaginary the buffer to use for the imaginary component (can be null)
   * @return [real, imaginary]
   * @throws IllegalArgumentException if not in the frequency domain
   * @see <A href=
   *      "https://en.wikipedia.org/wiki/Hartley_transform#Relation_to_Fourier_transform">https://en.wikipedia.org/
   *      wiki/Hartley_transform#Relation_to_Fourier_transform</a>
   */
  public FloatImage3D[] toDft(float[] real, float[] imaginary) {
    DhtHelper.checkFrequencyDomain(isFrequencyDomain);

    final float[] h1 = this.data;
    if (real == null || real.length != h1.length) {
      real = new float[h1.length];
    }
    if (imaginary == null || imaginary.length != h1.length) {
      imaginary = new float[h1.length];
    }

    for (int s = 0, nsMinusS = 0, i = 0; s < ns; s++, nsMinusS = ns - s) {
      for (int r = 0, nrMinusR = 0; r < nr; r++, nrMinusR = nr - r) {
        for (int c = 0, ncMinusC = 0; c < nc; c++, ncMinusC = nc - c, i++) {
          // This is a copy of the getComplexTransform operation in ij.process.FHT
          final int j = nsMinusS * nrByNc + nrMinusR * nc + ncMinusC;
          real[i] = (h1[i] + h1[j]) * 0.5f;
          imaginary[i] = (-h1[i] + h1[j]) * 0.5f;
        }
      }
    }

    return new FloatImage3D[] {new FloatImage3D(nc, nr, ns, nrByNc, real),
        new FloatImage3D(nc, nr, ns, nrByNc, imaginary)};
  }

  /**
   * Convert a discrete Fourier transform (DFT) to a DHT.
   *
   * @param real the real component
   * @param imaginary the imaginary component
   * @param tmp the tmp buffer to use for the result
   * @return the DHT
   * @throws IllegalArgumentException If there is a dimension mismatch
   */
  public static FloatDht3D fromDft(FloatImage3D real, FloatImage3D imaginary, float[] tmp) {
    if (real.ns != imaginary.ns || real.nr != imaginary.nr || real.nc != imaginary.nc) {
      throw new IllegalArgumentException("Dimension mismatch");
    }

    final float[] re = real.getData();
    final float[] im = imaginary.getData();
    if (tmp == null || tmp.length != re.length) {
      tmp = new float[re.length];
    }

    final int nc = real.nc;
    final int nr = real.nr;
    final int ns = real.ns;
    final int nrByNc = real.nrByNc;

    for (int s = 0, nsMinusS = 0, i = 0; s < ns; s++, nsMinusS = ns - s) {
      for (int r = 0, nrMinusR = 0; r < nr; r++, nrMinusR = nr - r) {
        for (int c = 0, ncMinusC = 0; c < nc; c++, ncMinusC = nc - c, i++) {
          final int j = nsMinusS * nrByNc + nrMinusR * nc + ncMinusC;
          // Reverse the toDFT() method
          // re = (a+b)/2
          // im = (-a+b)/2
          // b = re + im
          // a = 2*re - b
          tmp[j] = re[i] + im[i];
          tmp[i] = 2 * re[i] - tmp[j];
        }
      }
    }

    return new FloatDht3D(nc, nr, ns, nrByNc, tmp, true, new FloatDHT_3D(ns, nr, nc));
  }

  /**
   * Returns the absolute value (amplitude) of the Hartley transform. The image is assumed to be in
   * the frequency domain.
   *
   * @param tmp the tmp buffer to use for the result
   * @return the result
   * @throws IllegalArgumentException if not in the frequency domain
   */
  public FloatImage3D getAbsoluteValue(float[] tmp) {
    DhtHelper.checkFrequencyDomain(isFrequencyDomain);

    final float[] h1 = this.data;
    if (tmp == null || tmp.length != h1.length) {
      tmp = new float[h1.length];
    }

    for (int s = 0, nsMinusS = 0, i = 0; s < ns; s++, nsMinusS = ns - s) {
      for (int r = 0, nrMinusR = 0; r < nr; r++, nrMinusR = nr - r) {
        for (int c = 0, ncMinusC = 0; c < nc; c++, ncMinusC = nc - c, i++) {
          // This is a copy of the amplitude operation in ij.process.FHT
          final int j = nsMinusS * nrByNc + nrMinusR * nc + ncMinusC;
          tmp[i] = (float) Math.sqrt(h1[i] * h1[i] + h1[j] * h1[j]);
        }
      }
    }

    return new FloatImage3D(nc, nr, ns, nrByNc, tmp);
  }

  /**
   * Swap octants so the power spectrum origin is at the centre of the image.
   *
   * <pre>
   * 1 +++ &lt;=&gt; 7 ---
   * 2 -++ &lt;=&gt; 8 +--
   * 3 --+ &lt;=&gt; 5 ++-
   * 4 +-+ &lt;=&gt; 6 -+-
   * </pre>
   *
   * <p>Requires even dimensions.
   *
   * @throws IllegalArgumentException If not even dimensions
   * @see <a href=
   *      "https://en.m.wikipedia.org/wiki/Octant_(solid_geometry)">https://en.m.wikipedia.org/wiki/Octant_(solid_geometry)</a>
   */
  public void swapOctants() {
    swapOctants(this);
  }

  /**
   * Swap octants of the specified image stack so the power spectrum origin is at the centre of the
   * image.
   *
   * <pre>
   * 1 +++ &lt;=&gt; 7 ---
   * 2 -++ &lt;=&gt; 8 +--
   * 3 --+ &lt;=&gt; 5 ++-
   * 4 +-+ &lt;=&gt; 6 -+-
   * </pre>
   *
   * <p>Requires even dimensions.
   *
   * @param image the image
   * @throws IllegalArgumentException If not even dimensions
   * @see <a href=
   *      "https://en.m.wikipedia.org/wiki/Octant_(solid_geometry)">https://en.m.wikipedia.org/wiki/Octant_(solid_geometry)</a>
   */
  public static void swapOctants(FloatImage3D image) {
    final int ns = image.ns;
    final int nr = image.nr;
    final int nc = image.nc;

    if ((ns & 1) == 1 || (nr & 1) == 1 || (nc & 1) == 1) {
      throw new IllegalArgumentException("Require even dimensions");
    }

    final int nsOver2 = ns / 2;
    final int nrOver2 = nr / 2;
    final int ncOver2 = nc / 2;

    final float[] tmp = new float[nc];

    final int nrByNc = image.nrByNc;
    final float[] a = image.data;

    for (int s = 0; s < nsOver2; s++) {
      // Insert points
      final int ia = s * nrByNc;
      final int ib = (s + nsOver2) * nrByNc;

      //@formatter:off
      // We swap: 0 <=> ncOver2, 0 <=> ncOver2
      // 1 <=> 7
      swap(a, ia, a, ib, nc, ncOver2,       0,       0, nrOver2, ncOver2, nrOver2, tmp);
      // 2 <=> 8
      swap(a, ia, a, ib, nc,       0,       0, ncOver2, nrOver2, ncOver2, nrOver2, tmp);
      // 3 <=> 5
      swap(a, ia, a, ib, nc,       0, nrOver2, ncOver2,    0,    ncOver2, nrOver2, tmp);
      // 4 <=> 6
      swap(a, ia, a, ib, nc, ncOver2, nrOver2,       0,    0,    ncOver2, nrOver2, tmp);
      //@formatter:on
    }
  }

  /**
   * Swap octants of the specified image stack so the power spectrum origin is at the centre of the
   * image.
   *
   * <pre>
   * 1 +++ &lt;=&gt; 7 ---
   * 2 -++ &lt;=&gt; 8 +--
   * 3 --+ &lt;=&gt; 5 ++-
   * 4 +-+ &lt;=&gt; 6 -+-
   * </pre>
   *
   * <p>Requires even dimensions in a 32-bit float stack.
   *
   * @param stack the stack
   * @throws IllegalArgumentException If not a float stack with even dimensions
   * @see <a href=
   *      "https://en.m.wikipedia.org/wiki/Octant_(solid_geometry)">https://en.m.wikipedia.org/wiki/Octant_(solid_geometry)</a>
   */
  public static void swapOctants(ImageStack stack) {
    if (stack.getBitDepth() != 32) {
      throw new IllegalArgumentException("Require float stack");
    }
    final int ns = stack.getSize();
    final int nr = stack.getHeight();
    final int nc = stack.getWidth();
    if ((ns & 1) == 1 || (nr & 1) == 1 || (nc & 1) == 1) {
      throw new IllegalArgumentException("Require even dimensions");
    }

    final int nsOver2 = ns / 2;
    final int nrOver2 = nr / 2;
    final int ncOver2 = nc / 2;

    final float[] tmp = new float[nc];

    for (int s = 0; s < nsOver2; s++) {
      // slice index is 1-based
      final float[] a = (float[]) stack.getPixels(1 + s);
      final float[] b = (float[]) stack.getPixels(1 + s + nsOver2);
      //@formatter:off
      // We swap: 0 <=> ncOver2, 0 <=> ncOver2
      // 1 <=> 7
      Fht.swap(a, b, nc, ncOver2,       0,       0, nrOver2, ncOver2, nrOver2, tmp);
      // 2 <=> 8
      Fht.swap(a, b, nc,       0,       0, ncOver2, nrOver2, ncOver2, nrOver2, tmp);
      // 3 <=> 5
      Fht.swap(a, b, nc,       0, nrOver2, ncOver2,       0, ncOver2, nrOver2, tmp);
      // 4 <=> 6
      Fht.swap(a, b, nc, ncOver2, nrOver2,       0,       0, ncOver2, nrOver2, tmp);
      //@formatter:on
    }
  }

  /**
   * Swap the rectangle pixel values from a with b.
   *
   * <p>No bounds checks are performed so use with care!
   *
   * @param apixels the a pixels
   * @param ia the insert position for a
   * @param bpixels the b pixels (must match a.length)
   * @param ib the insert position for b
   * @param width the width of each set of XY pixels
   * @param ax the x origin from a
   * @param ay the y origin from a
   * @param bx the x origin from b
   * @param by the b origin from b
   * @param rw the width of the rectangle to swap
   * @param rh the height of the rectangle to swap
   * @param tmp the tmp buffer (must be at least width in length)
   */
  private static void swap(float[] apixels, int ia, float[] bpixels, int ib, int width, int ax,
      int ay, int bx, int by, int rw, int rh, float[] tmp) {
    for (int ayy = ay + rh, byy = by + rh - 1; ayy-- > ay; byy--) {
      final int ai = ia + ayy * width + ax;
      final int bi = ib + byy * width + bx;
      System.arraycopy(apixels, ai, tmp, 0, rw);
      System.arraycopy(bpixels, bi, apixels, ai, rw);
      System.arraycopy(tmp, 0, bpixels, bi, rw);
    }
  }
}
