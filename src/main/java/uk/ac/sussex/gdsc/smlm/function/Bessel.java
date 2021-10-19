/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

// All code has been ported from the Boost C++ implementation and is subject to
// the terms of the Boost Software Licence:

// Copyright (c) 2006 Xiaogang Zhang
// Copyright (c) 2017 John Maddock
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

package uk.ac.sussex.gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;

/**
 * Class for computing various Bessel functions
 *
 * <p>The implementation is based upon the Boost C++ implementation.
 *
 * @see <a href="https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/bessel.html">
 *      Boost C++ Bessel functions</a>
 */
public final class Bessel {
  /** 1 / sqrt(pi). */
  private static final double ONE_DIV_ROOT_PI = 5.641895835477562869480794515607725858e-01;
  /** sqrt(pi). */
  private static final double ROOT_PI = 1.772453850905516027298167483341145182;

  /** No public constructor. */
  private Bessel() {}

  // Implementation notes
  //
  // Polynomials in the Bessel I functions have been unrolled using
  // a second order Horner's rule.
  //
  // The functions are only called in Boost with a positive argument.
  // An explicit check for negative x is now performed and the functions
  // computed with abs(x) and the result returned appropriately for
  // odd/even functions

  /**
   * Compute the zero th order modified Bessel function of the first kind.
   *
   * @param x the x value
   * @return the modified Bessel function I0
   */
  public static double i0(double x) {
    // even function
    if (x <= 0) {
      if (x == 0) {
        return 1;
      }
      x = -x;
    }

    // boost/math/special_functions/detail/bessel_i0.hpp

    // Modified Bessel function of the first kind of order zero
    // we use the approximating forms derived in:
    // "Rational Approximations for the Modified Bessel Function of the First Kind - I0(x) for
    // Computations with Double Precision" by Pavel Holoborodko, see
    // http://www.advanpix.com/2015/11/11/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i0-computations-double-precision
    // The actual coefficients used are our own, and extend Pavel's work to precision's other than
    // double.

    if (x < 7.75) {
      // Bessel I0 over[10 ^ -16, 7.75]
      // Max error in interpolated form : 3.042e-18
      // Max Error found at double precision = Poly : 5.106609e-16 Cheb : 5.239199e-16
      final double a = x * x / 4;
      final double a2 = a * a;
      double p0;
      double p1;
      p0 = 9.07926920085624812e-25 * a2 + 2.63417742690109154e-20;
      p1 = 1.13943037744822825e-22 * a2 + 4.34709704153272287e-18;
      p0 = p0 * a2 + 6.27767773636292611e-16;
      p1 = p1 * a2 + 7.59389793369836367e-14;
      p0 = p0 * a2 + 7.59407002058973446e-12;
      p1 = p1 * a2 + 6.15118672704439289e-10;
      p0 = p0 * a2 + 3.93675991102510739e-08;
      p1 = p1 * a2 + 1.92901234513219920e-06;
      p0 = p0 * a2 + 6.94444444453352521e-05;
      p1 = p1 * a2 + 1.73611111111023792e-03;
      p0 = p0 * a2 + 2.77777777777782257e-02;
      p1 = p1 * a2 + 2.49999999999999909e-01;
      p0 = p0 * a2 + 1.00000000000000000e+00;
      final double p = p0 + p1 * a;
      return a * p + 1;
    } else if (x < 500) {
      // Max error in interpolated form : 1.685e-16
      // Max Error found at double precision = Poly : 2.575063e-16 Cheb : 2.247615e+00
      final double a = 1 / x;
      final double a2 = a * a;
      double p0;
      double p1;
      p0 = 2.17587543863819074e+15 * a2 + 2.02391097391687777e+15;
      p1 = -3.08675715295370878e+15 * a2 + -8.13426467865659318e+14;
      p0 = p0 * a2 + 2.24155239966958995e+14;
      p1 = p1 * a2 + -4.49034849696138065e+13;
      p0 = p0 * a2 + 6.76825737854096565e+12;
      p1 = p1 * a2 + -7.84261082124811106e+11;
      p0 = p0 * a2 + 7.08029243015109113e+10;
      p1 = p1 * a2 + -5.01883999713777929e+09;
      p0 = p0 * a2 + 2.80231938155267516e+08;
      p1 = p1 * a2 + -1.23157028595698731e+07;
      p0 = p0 * a2 + 4.24057674317867331e+05;
      p1 = p1 * a2 + -1.13366350697172355e+04;
      p0 = p0 * a2 + 2.33025711583514727e+02;
      p1 = p1 * a2 + -3.35052280231727022e+00;
      p0 = p0 * a2 + 1.30970574605856719e-01;
      p1 = p1 * a2 + 4.44207299493659561e-02;
      p0 = p0 * a2 + 2.92211225166047873e-02;
      p1 = p1 * a2 + 2.80506233928312623e-02;
      p0 = p0 * a2 + 4.98677850604961985e-02;
      p1 = p1 * a2 + 3.98942280401425088e-01;
      final double p = p0 * a + p1;
      return FastMath.exp(x) * p / Math.sqrt(x);
    } else {
      // Max error in interpolated form : 2.437e-18
      // Max Error found at double precision = Poly : 1.216719e-16
      final double a = 1 / x;
      final double a2 = a * a;
      double p0;
      double p1;
      p0 = 4.53371208762579442e-02 * a2 + 2.80506308916506102e-02;
      p1 = 2.92179096853915176e-02 * a2 + 4.98677850491434560e-02;
      p0 = p0 * a2 + 3.98942280401432905e-01;
      final double p = p0 + p1 * a;
      final double ex = FastMath.exp(x / 2);
      double result = ex * p / Math.sqrt(x);
      result *= ex;
      return result;
    }
  }

  /**
   * Compute the first order modified Bessel function of the first kind.
   *
   * @param x the x value
   * @return the modified Bessel function I1
   */
  public static double i1(double x) {
    // odd function
    double s;
    if (x <= 0) {
      if (x == 0) {
        return x;
      }
      x = -x;
      s = -1;
    } else {
      s = 1;
    }

    // boost/math/special_functions/detail/bessel_i1.hpp

    // Modified Bessel function of the first kind of order one
    // minimax rational approximations on intervals, see
    // Blair and Edwards, Chalk River Report AECL-4928, 1974
    if (x < 7.75) {
      // Bessel I0 over[10 ^ -16, 7.75]
      // Max error in interpolated form: 5.639e-17
      // Max Error found at double precision = Poly: 1.795559e-16

      final double a = x * x / 4;
      final double a2 = a * a;
      double p0;
      double p1;
      p0 = 1.332898928162290861e-23 * a2 + 3.410720494727771276e-19;
      p1 = 1.625212890947171108e-21 * a2 + 5.220157095351373194e-17;
      p0 = p0 * a2 + 6.904822652741917551e-15;
      p1 = p1 * a2 + 7.593969849687574339e-13;
      p0 = p0 * a2 + 6.834657311305621830e-11;
      p1 = p1 * a2 + 4.920949692800671435e-09;
      p0 = p0 * a2 + 2.755731926254790268e-07;
      p1 = p1 * a2 + 1.157407407354987232e-05;
      p0 = p0 * a2 + 3.472222222225921045e-04;
      p1 = p1 * a2 + 6.944444444444341983e-03;
      p0 = p0 * a2 + 8.333333333333333803e-02;
      final double p = p0 + p1 * a;
      final double q = 1 + (0.5 + p * a) * a;
      return s * x * q / 2;
    } else if (x < 500) {
      // Max error in interpolated form: 1.796e-16
      // Max Error found at double precision = Poly: 2.898731e-16

      final double a = 1 / x;
      final double a2 = a * a;
      double p0;
      double p1;
      p0 = -2.213318202179221945e+15 * a2 + -2.067285045778906105e+15;
      p1 = 3.146401654361325073e+15 * a2 + 8.325554073334618015e+14;
      p0 = p0 * a2 + -2.298849639457172489e+14;
      p1 = p1 * a2 + 4.614040809616582764e+13;
      p0 = p0 * a2 + -6.967602516005787001e+12;
      p1 = p1 * a2 + 8.087824484994859552e+11;
      p0 = p0 * a2 + -7.313784438967834057e+10;
      p1 = p1 * a2 + 5.192386898222206474e+09;
      p0 = p0 * a2 + -2.903390398236656519e+08;
      p1 = p1 * a2 + 1.277677779341446497e+07;
      p0 = p0 * a2 + -4.404655582443487334e+05;
      p1 = p1 * a2 + 1.178785865993440669e+04;
      p0 = p0 * a2 + -2.426181371595021021e+02;
      p1 = p1 * a2 + 3.458284470977172076e+00;
      p0 = p0 * a2 + -1.528189554374492735e-01;
      p1 = p1 * a2 + -5.719036414430205390e-02;
      p0 = p0 * a2 + -4.090895951581637791e-02;
      p1 = p1 * a2 + -4.675104253598537322e-02;
      p0 = p0 * a2 + -1.496033551613111533e-01;
      p1 = p1 * a2 + 3.989422804014406054e-01;
      final double p = p0 * a + p1;
      return s * FastMath.exp(x) * p / Math.sqrt(x);
    } else {
      // Max error in interpolated form: 1.320e-19
      // Max Error found at double precision = Poly: 7.065357e-17
      final double a = 1 / x;
      final double a2 = a * a;
      double p0;
      double p1;
      p0 = -5.843630344778927582e-02 * a2 + -4.675105322571775911e-02;
      p1 = -4.090421597376992892e-02 * a2 + -1.496033551467584157e-01;
      p0 = p0 * a2 + 3.989422804014314820e-01;
      final double p = p0 + p1 * a;
      final double ex = FastMath.exp(x / 2);
      double result = ex * p / Math.sqrt(x);
      result *= ex;
      return s * result;
    }
  }

  /**
   * Compute the zero th order Bessel function of the first kind.
   *
   * @param x the x value
   * @return the Bessel function J0
   */
  public static double j0(double x) {
    // even function
    if (x <= 0) {
      if (x == 0) {
        return 1;
      }
      x = -x;
    }

    // boost/math/special_functions/detail/bessel_j0.hpp

    // Bessel function of the first kind of order zero
    // x <= 8, minimax rational approximations on root-bracketing intervals
    // x > 8, Hankel asymptotic expansion in Hart, Computer Approximations, 1968

    // @formatter:off
    final double[] P1 = {
        -4.1298668500990866786e+11,
        2.7282507878605942706e+10,
        -6.2140700423540120665e+08,
        6.6302997904833794242e+06,
        -3.6629814655107086448e+04,
        1.0344222815443188943e+02,
        -1.2117036164593528341e-01
    };
    final double[] Q1 = {
         2.3883787996332290397e+12,
         2.6328198300859648632e+10,
         1.3985097372263433271e+08,
         4.5612696224219938200e+05,
         9.3614022392337710626e+02,
         1.0,
         0.0
    };
    final double[] P2 = {
         -1.8319397969392084011e+03,
         -1.2254078161378989535e+04,
         -7.2879702464464618998e+03,
         1.0341910641583726701e+04,
         1.1725046279757103576e+04,
         4.4176707025325087628e+03,
         7.4321196680624245801e+02,
         4.8591703355916499363e+01
    };
    final double[] Q2 = {
         -3.5783478026152301072e+05,
         2.4599102262586308984e+05,
         -8.4055062591169562211e+04,
         1.8680990008359188352e+04,
         -2.9458766545509337327e+03,
         3.3307310774649071172e+02,
         -2.5258076240801555057e+01,
         1.0
    };
    final double[] PC = {
         2.2779090197304684302e+04,
         4.1345386639580765797e+04,
         2.1170523380864944322e+04,
         3.4806486443249270347e+03,
         1.5376201909008354296e+02,
         8.8961548424210455236e-01
    };
    final double[] QC = {
         2.2779090197304684318e+04,
         4.1370412495510416640e+04,
         2.1215350561880115730e+04,
         3.5028735138235608207e+03,
         1.5711159858080893649e+02,
         1.0
    };
    final double[] PS = {
        -8.9226600200800094098e+01,
        -1.8591953644342993800e+02,
        -1.1183429920482737611e+02,
        -2.2300261666214198472e+01,
        -1.2441026745835638459e+00,
        -8.8033303048680751817e-03
    };
    final double[] QS = {
         5.7105024128512061905e+03,
         1.1951131543434613647e+04,
         7.2642780169211018836e+03,
         1.4887231232283756582e+03,
         9.0593769594993125859e+01,
         1.0
    };
    // @formatter:on
    final double x1 = 2.4048255576957727686e+00;
    final double x2 = 5.5200781102863106496e+00;
    final double x11 = 6.160e+02;
    final double x12 = -1.42444230422723137837e-03;
    final double x21 = 1.4130e+03;
    final double x22 = 5.46860286310649596604e-04;

    double factor;
    double r;
    double rc;
    double rs;

    if (x <= 4) {
      // x in (0, 4]
      final double y = x * x;
      r = evaluateRational7(P1, Q1, y);
      factor = (x + x1) * ((x - x11 / 256) - x12);
      return factor * r;
    } else if (x <= 8.0) {
      // x in (4, 8]
      final double y = 1 - (x * x) / 64;
      r = evaluateRational8(P2, Q2, y);
      factor = (x + x2) * ((x - x21 / 256) - x22);
      return factor * r;
    } else {
      // x in (8, \infty)
      final double y = 8 / x;
      final double y2 = y * y;
      rc = evaluateRational6(PC, QC, y2);
      rs = evaluateRational6(PS, QS, y2);
      factor = ONE_DIV_ROOT_PI / Math.sqrt(x);
      //
      // What follows is really just:
      //
      // double z = x - pi/4;
      // return factor * (rc * Math.cos(z) - y * rs * Math.sin(z));
      //
      // But using the addition formulae for Math.sin and Math.cos, plus
      // the special values for Math.sin/Math.cos of pi/4.
      //
      final double sx = Math.sin(x);
      final double cx = Math.cos(x);
      return factor * (rc * (cx + sx) - y * rs * (sx - cx));
    }
  }

  /**
   * Compute the first order Bessel function of the first kind.
   *
   * @param x the x value
   * @return the Bessel function J1
   */
  public static double j1(double x) {
    // odd function
    double s;
    if (x <= 0) {
      if (x == 0) {
        return x;
      }
      x = -x;
      s = -1;
    } else {
      s = 1;
    }

    // boost/math/special_functions/detail/bessel_j1.hpp

    // Bessel function of the first kind of order one
    // x <= 8, minimax rational approximations on root-bracketing intervals
    // x > 8, Hankel asymptotic expansion in Hart, Computer Approximations, 1968

    // @formatter:off
    final double[] P1 = {
        -1.4258509801366645672e+11,
        6.6781041261492395835e+09,
        -1.1548696764841276794e+08,
        9.8062904098958257677e+05,
        -4.4615792982775076130e+03,
        1.0650724020080236441e+01,
        -1.0767857011487300348e-02
    };
    final double[] Q1 = {
         4.1868604460820175290e+12,
         4.2091902282580133541e+10,
         2.0228375140097033958e+08,
         5.9117614494174794095e+05,
         1.0742272239517380498e+03,
         1.0,
         0.0
    };
    final double[] P2 = {
         -1.7527881995806511112e+16,
         1.6608531731299018674e+15,
         -3.6658018905416665164e+13,
         3.5580665670910619166e+11,
         -1.8113931269860667829e+09,
         5.0793266148011179143e+06,
         -7.5023342220781607561e+03,
         4.6179191852758252278e+00
    };
    final double[] Q2 = {
         1.7253905888447681194e+18,
         1.7128800897135812012e+16,
         8.4899346165481429307e+13,
         2.7622777286244082666e+11,
         6.4872502899596389593e+08,
         1.1267125065029138050e+06,
         1.3886978985861357615e+03,
         1.0
    };
    final double[] PC = {
        -4.4357578167941278571e+06,
        -9.9422465050776411957e+06,
        -6.6033732483649391093e+06,
        -1.5235293511811373833e+06,
        -1.0982405543459346727e+05,
        -1.6116166443246101165e+03,
        0.0
    };
    final double[] QC = {
        -4.4357578167941278568e+06,
        -9.9341243899345856590e+06,
        -6.5853394797230870728e+06,
        -1.5118095066341608816e+06,
        -1.0726385991103820119e+05,
        -1.4550094401904961825e+03,
        1.0
    };
    final double[] PS = {
         3.3220913409857223519e+04,
         8.5145160675335701966e+04,
         6.6178836581270835179e+04,
         1.8494262873223866797e+04,
         1.7063754290207680021e+03,
         3.5265133846636032186e+01,
         0.0
    };
    final double[] QS = {
         7.0871281941028743574e+05,
         1.8194580422439972989e+06,
         1.4194606696037208929e+06,
         4.0029443582266975117e+05,
         3.7890229745772202641e+04,
         8.6383677696049909675e+02,
         1.0
    };
    // @formatter:on
    final double x1 = 3.8317059702075123156e+00;
    final double x2 = 7.0155866698156187535e+00;
    final double x11 = 9.810e+02;
    final double x12 = -3.2527979248768438556e-04;
    final double x21 = 1.7960e+03;
    final double x22 = -3.8330184381246462950e-05;

    double factor;
    double r;
    double rc;
    double rs;

    if (x <= 4) {
      // x in (0, 4]
      final double y = x * x;
      r = evaluateRational7(P1, Q1, y);
      factor = x * (x + x1) * ((x - x11 / 256) - x12);
      return s * factor * r;
    } else if (x <= 8) {
      // x in (4, 8]
      final double y = x * x;
      r = evaluateRational8(P2, Q2, y);
      factor = x * (x + x2) * ((x - x21 / 256) - x22);
      return s * factor * r;
    } else {
      // x in (8, \infty)
      final double y = 8 / x;
      final double y2 = y * y;
      rc = evaluateRational7(PC, QC, y2);
      rs = evaluateRational7(PS, QS, y2);
      factor = 1 / (Math.sqrt(x) * ROOT_PI);
      //
      // What follows is really just:
      //
      // double z = x - 0.75f * Math.PI;
      // return s * factor * (rc * Math.cos(z) - y * rs * Math.sin(z));
      //
      // but using the Math.sin/Math.cos addition rules plus constants
      // for the values of Math.sin/Math.cos of 3PI/4 which then cancel
      // out with corresponding terms in "factor".
      //
      // Updated to use x in-place of x
      final double sx = Math.sin(x);
      final double cx = Math.cos(x);
      return s * factor * (rc * (sx - cx) + y * rs * (sx + cx));
    }
  }

  /**
   * Evaluate the rational number (numerator / denominator) when each is a polynomial evaluation of
   * length 6.
   *
   * <p>Adapted from boost/math/tools/detail/rational_horner3_6.hpp
   *
   * @param a First set of coefficients to compute numerator sum s1
   * @param b Second set of coefficients to compute denominator sum s2
   * @param x Point to evaluate
   * @return Polynomial rational s1 / s2
   */
  private static double evaluateRational6(double[] a, double[] b, double x) {
    double t0;
    double t1;
    double t2;
    double t3;

    // Note: This assumes x <= 1.
    // It is only called when x > 8 using (8 / x)^2.

    final double x2 = x * x;
    t0 = a[5] * x2 + a[3];
    t1 = a[4] * x2 + a[2];
    t2 = b[5] * x2 + b[3];
    t3 = b[4] * x2 + b[2];
    t0 *= x2;
    t1 *= x2;
    t2 *= x2;
    t3 *= x2;
    t0 += a[1];
    t1 += a[0];
    t2 += b[1];
    t3 += b[0];
    t0 *= x;
    t2 *= x;
    return (t0 + t1) / (t2 + t3);
  }

  /**
   * Evaluate the rational number (numerator / denominator) when each is a polynomial evaluation of
   * length 7.
   *
   * <p>Adapted from boost/math/tools/detail/rational_horner3_7.hpp
   *
   * @param a First set of coefficients to compute numerator sum s1
   * @param b Second set of coefficients to compute denominator sum s2
   * @param x Point to evaluate
   * @return Polynomial rational s1 / s2
   */
  private static double evaluateRational7(double[] a, double[] b, double x) {
    double t0;
    double t1;
    double t2;
    double t3;
    if (x <= 1) {
      final double x2 = x * x;
      t0 = a[6] * x2 + a[4];
      t1 = a[5] * x2 + a[3];
      t2 = b[6] * x2 + b[4];
      t3 = b[5] * x2 + b[3];
      t0 *= x2;
      t1 *= x2;
      t2 *= x2;
      t3 *= x2;
      t0 += a[2];
      t1 += a[1];
      t2 += b[2];
      t3 += b[1];
      t0 *= x2;
      t2 *= x2;
      t0 += a[0];
      t2 += b[0];
      t1 *= x;
      t3 *= x;
      return (t0 + t1) / (t2 + t3);
    }
    final double z = 1 / x;
    final double z2 = 1 / (x * x);
    t0 = a[0] * z2 + a[2];
    t1 = a[1] * z2 + a[3];
    t2 = b[0] * z2 + b[2];
    t3 = b[1] * z2 + b[3];
    t0 *= z2;
    t1 *= z2;
    t2 *= z2;
    t3 *= z2;
    t0 += a[4];
    t1 += a[5];
    t2 += b[4];
    t3 += b[5];
    t0 *= z2;
    t2 *= z2;
    t0 += a[6];
    t2 += b[6];
    t1 *= z;
    t3 *= z;
    return (t0 + t1) / (t2 + t3);
  }

  /**
   * Evaluate the rational number (numerator / denominator) when each is a polynomial evaluation of
   * length 8.
   *
   * <p>Adapted from boost/math/tools/detail/rational_horner3_8.hpp
   *
   * @param a First set of coefficients to compute numerator sum s1
   * @param b Second set of coefficients to compute denominator sum s2
   * @param x Point to evaluate
   * @return Polynomial rational s1 / s2
   */
  private static double evaluateRational8(double[] a, double[] b, double x) {
    double t0;
    double t1;
    double t2;
    double t3;
    if (x <= 1) {
      final double x2 = x * x;
      t0 = a[7] * x2 + a[5];
      t1 = a[6] * x2 + a[4];
      t2 = b[7] * x2 + b[5];
      t3 = b[6] * x2 + b[4];
      t0 *= x2;
      t1 *= x2;
      t2 *= x2;
      t3 *= x2;
      t0 += a[3];
      t1 += a[2];
      t2 += b[3];
      t3 += b[2];
      t0 *= x2;
      t1 *= x2;
      t2 *= x2;
      t3 *= x2;
      t0 += a[1];
      t1 += a[0];
      t2 += b[1];
      t3 += b[0];
      t0 *= x;
      t2 *= x;
      return (t0 + t1) / (t2 + t3);
    }
    final double z = 1 / x;
    final double z2 = 1 / (x * x);
    t0 = a[0] * z2 + a[2];
    t1 = a[1] * z2 + a[3];
    t2 = b[0] * z2 + b[2];
    t3 = b[1] * z2 + b[3];
    t0 *= z2;
    t1 *= z2;
    t2 *= z2;
    t3 *= z2;
    t0 += a[4];
    t1 += a[5];
    t2 += b[4];
    t3 += b[5];
    t0 *= z2;
    t1 *= z2;
    t2 *= z2;
    t3 *= z2;
    t0 += a[6];
    t1 += a[7];
    t2 += b[6];
    t3 += b[7];
    t0 *= z;
    t2 *= z;
    return (t0 + t1) / (t2 + t3);
  }
}
