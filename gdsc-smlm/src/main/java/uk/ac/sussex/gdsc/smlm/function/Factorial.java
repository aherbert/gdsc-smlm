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

import org.apache.commons.math3.special.Gamma;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

/**
 * Compute the factorial function ({@code n!}).
 */
public final class Factorial {
  /** All factorials that can be represented as a double (values up to 170!). */
  private static final double[] FACTORIALS = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880,
      3628800, 3.99168E7, 4.790016E8, 6.2270208E9, 8.71782912E10, 1.307674368E12, 2.0922789888E13,
      3.55687428096E14, 6.402373705728E15, 1.21645100408832E17, 2.43290200817664E18,
      5.109094217170944E19, 1.1240007277776077E21, 2.585201673888498E22, 6.204484017332394E23,
      1.5511210043330986E25, 4.0329146112660565E26, 1.0888869450418352E28, 3.0488834461171387E29,
      8.841761993739702E30, 2.6525285981219107E32, 8.222838654177922E33, 2.631308369336935E35,
      8.683317618811886E36, 2.9523279903960416E38, 1.0333147966386145E40, 3.7199332678990125E41,
      1.3763753091226346E43, 5.230226174666011E44, 2.0397882081197444E46, 8.159152832478977E47,
      3.345252661316381E49, 1.40500611775288E51, 6.041526306337383E52, 2.658271574788449E54,
      1.1962222086548019E56, 5.502622159812089E57, 2.5862324151116818E59, 1.2413915592536073E61,
      6.082818640342675E62, 3.0414093201713376E64, 1.5511187532873822E66, 8.065817517094388E67,
      4.2748832840600255E69, 2.308436973392414E71, 1.2696403353658276E73, 7.109985878048635E74,
      4.0526919504877214E76, 2.3505613312828785E78, 1.3868311854568984E80, 8.32098711274139E81,
      5.075802138772248E83, 3.146997326038794E85, 1.98260831540444E87, 1.2688693218588417E89,
      8.247650592082472E90, 5.443449390774431E92, 3.647111091818868E94, 2.4800355424368305E96,
      1.711224524281413E98, 1.1978571669969892E100, 8.504785885678623E101, 6.1234458376886085E103,
      4.4701154615126844E105, 3.307885441519386E107, 2.48091408113954E109, 1.8854947016660504E111,
      1.4518309202828587E113, 1.1324281178206297E115, 8.946182130782976E116, 7.156945704626381E118,
      5.797126020747368E120, 4.753643337012842E122, 3.945523969720659E124, 3.314240134565353E126,
      2.81710411438055E128, 2.4227095383672734E130, 2.107757298379528E132, 1.8548264225739844E134,
      1.650795516090846E136, 1.4857159644817615E138, 1.352001527678403E140, 1.2438414054641308E142,
      1.1567725070816416E144, 1.087366156656743E146, 1.032997848823906E148, 9.916779348709496E149,
      9.619275968248212E151, 9.426890448883248E153, 9.332621544394415E155, 9.332621544394415E157,
      9.42594775983836E159, 9.614466715035127E161, 9.90290071648618E163, 1.0299016745145628E166,
      1.081396758240291E168, 1.1462805637347084E170, 1.226520203196138E172, 1.324641819451829E174,
      1.4438595832024937E176, 1.588245541522743E178, 1.7629525510902446E180, 1.974506857221074E182,
      2.2311927486598138E184, 2.5435597334721877E186, 2.925093693493016E188, 3.393108684451898E190,
      3.969937160808721E192, 4.684525849754291E194, 5.574585761207606E196, 6.689502913449127E198,
      8.094298525273444E200, 9.875044200833601E202, 1.214630436702533E205, 1.506141741511141E207,
      1.882677176888926E209, 2.372173242880047E211, 3.0126600184576594E213, 3.856204823625804E215,
      4.974504222477287E217, 6.466855489220474E219, 8.47158069087882E221, 1.1182486511960043E224,
      1.4872707060906857E226, 1.9929427461615188E228, 2.6904727073180504E230, 3.659042881952549E232,
      5.012888748274992E234, 6.917786472619489E236, 9.615723196941089E238, 1.3462012475717526E241,
      1.898143759076171E243, 2.695364137888163E245, 3.854370717180073E247, 5.5502938327393044E249,
      8.047926057471992E251, 1.1749972043909107E254, 1.727245890454639E256, 2.5563239178728654E258,
      3.80892263763057E260, 5.713383956445855E262, 8.62720977423324E264, 1.3113358856834524E267,
      2.0063439050956823E269, 3.0897696138473508E271, 4.789142901463394E273, 7.471062926282894E275,
      1.1729568794264145E278, 1.853271869493735E280, 2.9467022724950384E282, 4.7147236359920616E284,
      7.590705053947219E286, 1.2296942187394494E289, 2.0044015765453026E291, 3.287218585534296E293,
      5.423910666131589E295, 9.003691705778438E297, 1.503616514864999E300, 2.5260757449731984E302,
      4.269068009004705E304, 7.257415615307999E306};
  /** Precomputed factor LANCZOS_G + 1.5. This has a value of 799 / 128. */
  private static final double LANCZOS_G_PLUS_1_5 = Gamma.LANCZOS_G + 1.5;
  /** sqrt(2 pi). */
  private static final double SQRT_TWO_PI = 2.506628274631000502415765284811;

  /** No instances. */
  private Factorial() {}

  /**
   * Get the factorial of n.
   *
   * <p>The maximum value representable as a double is {@code 170!}.
   *
   * @param n argument (must be positive)
   * @return n!
   * @throws ArrayIndexOutOfBoundsException if n is negative
   */
  public static double value(int n) {
    if (n < FACTORIALS.length) {
      return FACTORIALS[n];
    }
    return Double.POSITIVE_INFINITY;
  }

  /**
   * Get the factorial of n.
   *
   * <p>The maximum value representable as a double is {@code 170!}.
   *
   * @param n argument (must be positive)
   * @return n!
   * @throws ArrayIndexOutOfBoundsException if n is negative or above 170
   */
  static double uncheckedValue(int n) {
    return FACTORIALS[n];
  }

  /**
   * Get the factorial of n as a real number.
   *
   * <p>This is computed using {@code gamma(1+n)}.
   *
   * @param n argument (must be positive)
   * @return n!
   * @throws IllegalArgumentException if n is negative
   */
  public static double value(double n) {
    ValidationUtils.checkPositive(n);
    if (n < FACTORIALS.length) {
      if (Math.rint(n) == n) {
        return FACTORIALS[(int) n];
      }
      return gamma1p(n);
    }
    return Double.POSITIVE_INFINITY;
  }

  /**
   * Compute {@code Gamma(1+x)}.
   *
   * <p>Adapted from {@code org.apache.commons.math3.special.Gamma} gamma and gamma1p. Removed
   * support for NaN and negative x. Updated the code to use x or (x+1) where appropriate.
   *
   * <p>For {@code x <= 19}, the implementation is based on the double precision implementation in
   * the <em>NSWC Library of Mathematics Subroutines</em>; otherwise the implementation is based on
   * </p>
   *
   * <ul>
   *
   * <li><a href="https://mathworld.wolfram.com/LanczosApproximation.html"> Lanczos
   * Approximation</a>, equations (1) through (5).</li>
   *
   * <li><a href="http://www.numericana.com/answer/info/godfrey.htm">Paul Godfrey, A note on the
   * computation of the convergent Lanczos complex Gamma approximation</a></li>
   *
   * </ul>
   *
   * @param x argument
   * @return the value of {@code log(Gamma(1+x))}
   */
  private static double gamma1p(double x) {
    if (x <= 19) {
      // @formatter:off
      // From the recurrence relation
      // Gamma(x + 1) = x * (x - 1) * ... * (x - n + 1) * Gamma(x - n + 1),
      // and
      // Gamma(t + 1) = 1 / [1 + invGamma1pm1(t)],
      // Shift t into range of invGamma1pm1. This requires t <= 1.5.
      // @formatter:on
      double prod = 1;
      double t = x;
      while (t > 1.5) {
        prod *= t;
        t -= 1;
      }
      return prod / (1 + Gamma.invGamma1pm1(t));
    }

    final double x1p = x + 1;
    final double y = x + LANCZOS_G_PLUS_1_5;
    return SQRT_TWO_PI / x1p * Math.pow(y, x + 1.5) * Math.exp(-y) * Gamma.lanczos(x1p);
  }
}
