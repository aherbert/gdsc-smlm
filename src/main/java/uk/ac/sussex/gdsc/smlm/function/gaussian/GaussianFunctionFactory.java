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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiAstigmatismErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiFixedErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiNbCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.MultiNbFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleAstigmatismErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFixedErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleFreeCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleNbCircularErfGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.erf.SingleNbFreeCircularErfGaussian2DFunction;

/**
 * Creates the appropriate Gaussian function.
 *
 * <p>Note that currently all functions support computing gradients for x/y position.
 */
public class GaussianFunctionFactory {
  /**
   * Compute gradients for background.
   */
  public static final int FIT_BACKGROUND = 0x00000001;
  /**
   * Compute gradients for rotation angle.
   */
  public static final int FIT_ANGLE = 0x00000002;
  /**
   * Compute gradients for x width.
   */
  public static final int FIT_X_WIDTH = 0x00000004;
  /**
   * Compute gradients for y width.
   */
  public static final int FIT_Y_WIDTH = 0x00000008;
  /**
   * Compute gradients for signal.
   */
  public static final int FIT_SIGNAL = 0x00000010;
  /**
   * Compute gradients for z position.
   */
  public static final int FIT_Z = 0x00000020;

  /**
   * An elliptical 2D Gaussian with gradients for background, signal, rotation angle, x/y position,
   * x/y width.
   */
  public static final int FIT_ELLIPTICAL =
      FIT_BACKGROUND | FIT_ANGLE | FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL;

  /** An elliptical 2D Gaussian with gradients for background, signal, x/y position, x/y width. */
  public static final int FIT_FREE_CIRCLE = FIT_BACKGROUND | FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL;

  /** A 2D Gaussian with gradients for background, signal, x/y position, width. */
  public static final int FIT_CIRCLE = FIT_BACKGROUND | FIT_X_WIDTH | FIT_SIGNAL;

  /** A 2D Gaussian with gradients for background, signal, x/y position. */
  public static final int FIT_FIXED = FIT_BACKGROUND | FIT_SIGNAL;
  /**
   * An elliptical 2D Gaussian with gradients for background, signal, x/y/z position. The z position
   * determines the x/y width using an astigmatism model.
   */
  public static final int FIT_ASTIGMATISM = FIT_BACKGROUND | FIT_Z | FIT_SIGNAL;

  // -=-=-=-=-=-=-=-=-=-=-=-=-
  // Flags for ERF Gaussian functions.
  // These are evaluated as a full integration over the pixel using the error function (erf)
  // -=-=-=-=-=-=-=-=-=-=-=-=-

  /**
   * Use a full integration over the pixel using the error function (erf). A rotation angle is not
   * supported.
   */
  public static final int FIT_ERF = 0x00000100;

  /**
   * An elliptical 2D Gaussian (full integration per pixel) with gradients for background, signal,
   * x/y position, x/y width.
   */
  public static final int FIT_ERF_FREE_CIRCLE = FIT_FREE_CIRCLE | FIT_ERF;

  /**
   * A 2D Gaussian (full integration per pixel) with gradients for background, signal, x/y position,
   * width.
   */
  public static final int FIT_ERF_CIRCLE = FIT_CIRCLE | FIT_ERF;

  /**
   * A 2D Gaussian (full integration per pixel) with gradients for background, signal, x/y position.
   */
  public static final int FIT_ERF_FIXED = FIT_FIXED | FIT_ERF;
  /**
   * An elliptical 2D Gaussian (full integration per pixel) with gradients for background, signal,
   * x/y/z position. The z position determines the x/y width using an astigmatism model.
   */
  public static final int FIT_ERF_ASTIGMATISM = FIT_BACKGROUND | FIT_Z | FIT_SIGNAL | FIT_ERF;

  // -=-=-=-=-=-=-=-=-=-=-=-=-
  // Flags for simple Gaussian functions.
  // These are evaluated using a single exponential at the centre of the pixel.
  // They support rotating the X/Y elliptical Gaussian (if X and Y are different).
  // -=-=-=-=-=-=-=-=-=-=-=-=-

  /**
   * Use a single exponential at the centre of the pixel. A rotation angle is supported.
   */
  public static final int FIT_SIMPLE = 0x00000200;

  /**
   * An elliptical 2D Gaussian (single exponential per pixel) with gradients for background, signal,
   * rotation angle, x/y position, x/y width.
   */
  public static final int FIT_SIMPLE_ELLIPTICAL = FIT_ELLIPTICAL | FIT_SIMPLE;

  /**
   * An elliptical 2D Gaussian (single exponential per pixel) with gradients for background, signal,
   * x/y position, x/y width.
   */
  public static final int FIT_SIMPLE_FREE_CIRCLE = FIT_FREE_CIRCLE | FIT_SIMPLE;

  /**
   * A 2D Gaussian (single exponential per pixel) with gradients for background, signal, x/y
   * position, width.
   */
  public static final int FIT_SIMPLE_CIRCLE = FIT_CIRCLE | FIT_SIMPLE;

  /**
   * A 2D Gaussian (single exponential per pixel) with gradients for background, signal, x/y
   * position.
   */
  public static final int FIT_SIMPLE_FIXED = FIT_FIXED | FIT_SIMPLE;

  // Extra support for functions without background

  /**
   * An elliptical 2D Gaussian (full integration per pixel) with gradients for signal, x/y position,
   * x/y width.
   */
  public static final int FIT_ERF_NB_FREE_CIRCLE = FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL | FIT_ERF;

  /** An 2D Gaussian (full integration per pixel) with gradients for signal, x/y position, width. */
  public static final int FIT_ERF_NB_CIRCLE = FIT_X_WIDTH | FIT_SIGNAL | FIT_ERF;

  /**
   * An elliptical 2D Gaussian (single exponential per pixel) with gradients for signal, rotation
   * angle, x/y position, x/y width.
   */
  public static final int FIT_SIMPLE_NB_ELLIPTICAL =
      FIT_ANGLE | FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL | FIT_SIMPLE;

  /**
   * An elliptical 2D Gaussian (single exponential per pixel) with gradients for signal, x/y
   * position, x/y width.
   */
  public static final int FIT_SIMPLE_NB_FREE_CIRCLE =
      FIT_X_WIDTH | FIT_Y_WIDTH | FIT_SIGNAL | FIT_SIMPLE;

  /**
   * A 2D Gaussian (single exponential per pixel) with gradients for signal, x/y position, width.
   */
  public static final int FIT_SIMPLE_NB_CIRCLE = FIT_X_WIDTH | FIT_SIGNAL | FIT_SIMPLE;

  /** A 2D Gaussian (single exponential per pixel) with gradients for signal, x/y position. */
  public static final int FIT_SIMPLE_NB_FIXED = FIT_SIGNAL | FIT_SIMPLE;

  /** A 2D Gaussian (single exponential per pixel) with gradients for background, x/y position. */
  public static final int FIT_SIMPLE_NS_FIXED = FIT_BACKGROUND | FIT_SIMPLE;

  /** A 2D Gaussian (single exponential per pixel) with gradients for x/y position. */
  public static final int FIT_SIMPLE_NS_NB_FIXED = FIT_SIMPLE;

  /** The flags. */
  private final int flags;

  /** The z model. */
  private final AstigmatismZModel zModel;

  /**
   * Instantiates a new gaussian function factory.
   *
   * <p>This is a simple instance that calls
   * {@link #create2D(int, int, int, int, AstigmatismZModel)} with a fixed set of flags and z-model.
   *
   * @param flags the flags
   * @param zModel the z model
   */
  public GaussianFunctionFactory(int flags, AstigmatismZModel zModel) {
    this.flags = flags;
    this.zModel = zModel;
  }

  /**
   * Create the correct 2D Gaussian function for the specified parameters.
   *
   * <p>Defaults to using the ERF Gaussian functions if the user has not requested a simple Gaussian
   * or angle fitting.
   *
   * @param numberOfPeaks The number of peaks (N)
   * @param maxx The maximum X-dimension
   * @param maxy The maximum Y-dimension
   * @return The function
   */
  public Gaussian2DFunction create2D(int numberOfPeaks, int maxx, int maxy) {
    return create2D(numberOfPeaks, maxx, maxy, flags, zModel);
  }

  /**
   * Create the correct 2D Gaussian function for the specified parameters.
   *
   * <p>Defaults to using the ERF Gaussian functions if the user has not requested a simple Gaussian
   * or angle fitting.
   *
   * @param numberOfPeaks The number of peaks (N)
   * @param maxx The maximum X-dimension
   * @param maxy The maximum Y-dimension
   * @param flags Enable all the parameters that should evaluate gradient
   * @param zModel the z model
   * @return The function
   */
  public static Gaussian2DFunction create2D(int numberOfPeaks, int maxx, int maxy, int flags,
      AstigmatismZModel zModel) {
    // Default to using the ERF functions if the user has not requested a simple Gaussian or angle
    // fitting
    if ((flags & (FIT_SIMPLE | FIT_ANGLE)) == 0) {
      if (numberOfPeaks == 1) {
        if ((flags & FIT_BACKGROUND) != 0) {
          // Independent X/Y width
          if ((flags & FIT_Y_WIDTH) != 0) {
            return new SingleFreeCircularErfGaussian2DFunction(maxx, maxy);
          }
          // Combined X/Y width
          if ((flags & FIT_X_WIDTH) != 0) {
            return new SingleCircularErfGaussian2DFunction(maxx, maxy);
          }
          // Z-depth function
          if ((flags & FIT_Z) != 0) {
            return new SingleAstigmatismErfGaussian2DFunction(maxx, maxy, zModel);
          }
          // Fixed width
          if ((flags & FIT_SIGNAL) != 0) {
            return new SingleFixedErfGaussian2DFunction(maxx, maxy);
          }
        }
        if ((flags & FIT_Y_WIDTH) != 0) {
          return new SingleNbFreeCircularErfGaussian2DFunction(maxx, maxy);
        }
        if ((flags & FIT_X_WIDTH) != 0) {
          return new SingleNbCircularErfGaussian2DFunction(maxx, maxy);
        }
      } else {
        if ((flags & FIT_BACKGROUND) != 0) {
          // Independent X/Y width
          if ((flags & FIT_Y_WIDTH) != 0) {
            return new MultiFreeCircularErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
          }
          // Combined X/Y width
          if ((flags & FIT_X_WIDTH) != 0) {
            return new MultiCircularErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
          }
          // Z-depth function
          if ((flags & FIT_Z) != 0) {
            return new MultiAstigmatismErfGaussian2DFunction(numberOfPeaks, maxx, maxy, zModel);
          }
          // Fixed width
          if ((flags & FIT_SIGNAL) != 0) {
            return new MultiFixedErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
          }
        }
        if ((flags & FIT_Y_WIDTH) != 0) {
          return new MultiNbFreeCircularErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
        }
        if ((flags & FIT_X_WIDTH) != 0) {
          return new MultiNbCircularErfGaussian2DFunction(numberOfPeaks, maxx, maxy);
        }
      }
    }

    // Legacy simple Gaussian functions
    if (numberOfPeaks == 1) {
      if ((flags & FIT_BACKGROUND) != 0) {
        if ((flags & FIT_ANGLE) != 0) {
          return new SingleEllipticalGaussian2DFunction(maxx, maxy);
        }
        if ((flags & FIT_Y_WIDTH) != 0) {
          return new SingleFreeCircularGaussian2DFunction(maxx, maxy);
        }
        if ((flags & FIT_X_WIDTH) != 0) {
          return new SingleCircularGaussian2DFunction(maxx, maxy);
        }

        // Fixed function
        if ((flags & FIT_SIGNAL) != 0) {
          return new SingleFixedGaussian2DFunction(maxx, maxy);
        }

        return new SingleNsFixedGaussian2DFunction(maxx, maxy);
      }

      if ((flags & FIT_ANGLE) != 0) {
        return new SingleNbEllipticalGaussian2DFunction(maxx, maxy);
      }
      if ((flags & FIT_Y_WIDTH) != 0) {
        return new SingleNbFreeCircularGaussian2DFunction(maxx, maxy);
      }
      if ((flags & FIT_X_WIDTH) != 0) {
        return new SingleNbCircularGaussian2DFunction(maxx, maxy);
      }

      // Fixed function
      if ((flags & FIT_SIGNAL) != 0) {
        return new SingleNbFixedGaussian2DFunction(maxx, maxy);
      }

      return new SingleNsNbFixedGaussian2DFunction(maxx, maxy);
    }

    if ((flags & FIT_BACKGROUND) != 0) {
      if ((flags & FIT_ANGLE) != 0) {
        return new EllipticalGaussian2DFunction(numberOfPeaks, maxx, maxy);
      }
      if ((flags & FIT_Y_WIDTH) != 0) {
        return new FreeCircularGaussian2DFunction(numberOfPeaks, maxx, maxy);
      }
      if ((flags & FIT_X_WIDTH) != 0) {
        return new CircularGaussian2DFunction(numberOfPeaks, maxx, maxy);
      }

      // Fixed function
      if ((flags & FIT_SIGNAL) != 0) {
        return new FixedGaussian2DFunction(numberOfPeaks, maxx, maxy);
      }

      return new NsFixedGaussian2DFunction(numberOfPeaks, maxx, maxy);
    }

    if ((flags & FIT_ANGLE) != 0) {
      return new NbEllipticalGaussian2DFunction(numberOfPeaks, maxx, maxy);
    }
    if ((flags & FIT_Y_WIDTH) != 0) {
      return new NbFreeCircularGaussian2DFunction(numberOfPeaks, maxx, maxy);
    }
    if ((flags & FIT_X_WIDTH) != 0) {
      return new NbCircularGaussian2DFunction(numberOfPeaks, maxx, maxy);
    }

    // Fixed function
    if ((flags & FIT_SIGNAL) != 0) {
      return new NbFixedGaussian2DFunction(numberOfPeaks, maxx, maxy);
    }

    return new NsNbFixedGaussian2DFunction(numberOfPeaks, maxx, maxy);
  }

  /**
   * Freeze the Gaussian function parameters to a set of parameters and function flags that will
   * evaluate faster. The widths of the input parameters may be modified. The shape will remain
   * unchanged even if the function flags are altered.
   *
   * @param flags The function flags (defines the type of function)
   * @param zModel the z model
   * @param a the parameters of the function (modified in place)
   * @return the flags for the new function
   * @throws IllegalArgumentException If the input array is null or not the correct size for a
   *         multiple of Gaussian peak parameters.
   */
  public static int freeze(int flags, AstigmatismZModel zModel, double[] a) {
    if (a == null) {
      throw new IllegalArgumentException("Parameter array is null");
    }
    final int numberOfPeaks = a.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
    // Check the array length is correct (including the background at position 0)
    if (a.length != 1 + numberOfPeaks * Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      throw new IllegalArgumentException("Incorrect array size for the number of peaks");
    }
    if (numberOfPeaks == 1) {
      return freezeSingle(flags, zModel, a);
    }
    // This supports zero peaks too
    return freezeMulti(numberOfPeaks, flags, zModel, a);
  }

  /**
   * Freeze the Gaussian function parameters to a set of parameters and function flags that will
   * evaluate faster. The widths of the input parameters may be modified. The shape will remain
   * unchanged even if the function flags are altered.
   *
   * @param numberOfPeaks The number of peaks (N)
   * @param flags The function flags (defines the type of function)
   * @param zModel the z model
   * @param a the parameters of the function (modified in place)
   * @return the flags for the new function
   */
  private static int freezeMulti(int numberOfPeaks, int flags, AstigmatismZModel zModel,
      double[] a) {
    // Default to using the ERF functions if the user has not requested a simple Gaussian or angle
    // fitting
    if ((flags & (FIT_SIMPLE | FIT_ANGLE)) == 0) {
      // Check for z-model
      if ((flags & FIT_Z) != 0 && zModel != null) {
        // Convert the parameters for a free-circle function
        for (int i = 0, j = 0; i < numberOfPeaks;
            i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
          final double tz = a[j + Gaussian2DFunction.Z_POSITION];
          a[j + Gaussian2DFunction.X_SD] = zModel.getSx(tz);
          a[j + Gaussian2DFunction.Y_SD] = zModel.getSy(tz);
        }
        return FIT_ERF_FREE_CIRCLE;
      }

      // Check the need for X/Y widths
      if ((flags & FIT_Y_WIDTH) == 0) {
        for (int i = 0, j = 0; i < numberOfPeaks;
            i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
          a[j + Gaussian2DFunction.Y_SD] = a[j + Gaussian2DFunction.X_SD];
        }
        return FIT_ERF_FIXED;
      }

      return FIT_ERF_FREE_CIRCLE;
    }

    // Legacy simple Gaussian functions

    // Check the need for X/Y widths
    if ((flags & FIT_Y_WIDTH) == 0) {
      for (int i = 0, j = 0; i < numberOfPeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
        a[j + Gaussian2DFunction.Y_SD] = a[j + Gaussian2DFunction.X_SD];
      }
      return FIT_SIMPLE_FIXED;
    }

    // Default to free circle (which can use the angle if present)
    return FIT_SIMPLE_FREE_CIRCLE;
  }

  /**
   * Freeze the Gaussian function parameters to a set of parameters and function flags that will
   * evaluate faster. The widths of the input parameters may be modified. The shape will remain
   * unchanged even if the function flags are altered.
   *
   * @param flags The function flags (defines the type of function)
   * @param zModel the z model
   * @param a the parameters of the function (modified in place)
   * @return the flags for the new function
   */
  private static int freezeSingle(int flags, AstigmatismZModel zModel, double[] a) {
    // Default to using the ERF functions if the user has not requested a simple Gaussian or angle
    // fitting
    if ((flags & (FIT_SIMPLE | FIT_ANGLE)) == 0) {
      // Check for z-model
      if ((flags & FIT_Z) != 0 && zModel != null) {
        // Convert the parameters for a free-circle function
        final double tz = a[Gaussian2DFunction.Z_POSITION];
        a[Gaussian2DFunction.X_SD] = zModel.getSx(tz);
        a[Gaussian2DFunction.Y_SD] = zModel.getSy(tz);
        return FIT_ERF_FREE_CIRCLE;
      }

      // Check the need for X/Y widths
      if ((flags & FIT_Y_WIDTH) == 0) {
        a[Gaussian2DFunction.Y_SD] = a[Gaussian2DFunction.X_SD];
        return FIT_ERF_FIXED;
      }

      return FIT_ERF_FREE_CIRCLE;
    }

    // Legacy simple Gaussian functions

    // Check the need for X/Y widths
    if ((flags & FIT_Y_WIDTH) == 0) {
      a[Gaussian2DFunction.Y_SD] = a[Gaussian2DFunction.X_SD];
      return FIT_SIMPLE_FIXED;
    }

    // Default to free circle (which can use the angle if present)
    return FIT_SIMPLE_FREE_CIRCLE;
  }
}
