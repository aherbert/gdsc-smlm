/*
 * BesselFunctions.java
 *
 * Created on February, 16, 2009 
 * Christopher K. Allen
 */

/*
 * Copyright (c) 2003-2006 jMonkeyEngine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * Neither the name of 'jMonkeyEngine' nor the names of its contributors 
 *   may be used to endorse or promote products derived from this software 
 *   without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package xal.tools.math;


/**
 *<h1>Bessel Function Implementations</h1> 
 * <p>
 * Utility case for computing various Bessel functions that are not
 * included in the Java API.
 * </p>
 * <p>
 * Implementations for integer-order cylindrical Bessel functions of the
 * first and second kind (i.e., <i>J<sub>n</sub></i>(<i>x</i>) and 
 * <i>Y<sub>n</sub></i>(<i>x</i>), <i>n</i> = 1,2,3,... )were taken from
 * <a href="http://www.koders.com">www.koders.com</a> and the copyright
 * notice is included in the Java source file.  The implementation is based
 * upon that presented in 
 * <a href="http://www.nr.com"><i>Numerical Recipes</i></a> 
 * by W.H. Press, <i>et. al.</i>
 * </p>
 * <p>
 * Spherical Bessel functions <i>j<sub>n</sub></i>(<i>x</i>)
 * can be represented as cylindrical Bessel functions
 * with half-integer order.  Specifically, we have
 * <br/>
 * <br/>
 * &nbsp; <i>j<sub>n</sub></i>(<i>x</i>) = (&pi;/2<i>x</i>)<sup>1/2</sup><i>J<sub>n</i>+&frac12;</sub>(<i>x</i>).
 * <br/>
 * <br/>
 * However, since half-order cylindrical Bessel functions are not included in 
 * this class
 * (they are more difficult to implement), an implementation based upon the
 * above formula is not feasible.  Instead, 
 * spherical Bessel functions <i>j<sub>n</sub></i>(<i>x</i>)
 * for <i>n</i> = 0, 1, 2, 3, 4 are implemented using their trigonometric 
 * representations. 
 * </p>
 * <p>
 * <strong>NOTES</strong>: (CKA)
 * <br/>
 * &middot; There exists a recurrence relationship between Bessel functions of
 * different orders.  For example, if <i>B<sub>n</sub></i>(<i>x</i>) is any
 * cylindrical Bessel function of order <i>n</i>, we have the following:
 * <br/>
 * <br/>
 * &nbsp;  <i>B<sub>n</i>+1</sub>(<i>x</i>) = (2<i>n</i>/<i>x</i>)<i>B<sub>n</sub></i>(<i>x</i>) - 
 *                                            <i>B<sub>n</i>-1</sub>(<i>x</i>).
 * <br/>
 * <br/>
 * However, this formula is numerically unstable for Bessel functions of the first
 * kind <i>J<sub>n</sub></i>(<i>x</i>). Thus, it cannot be used to compute the 
 * higher order <i>J<sub>n</sub></i>(<i>x</i>) using recursion over lower orders.
 * <br/>
 * &middot; We can apply the above recurrence relation to spherical Bessel 
 * functions by expressing them in terms of cylindrical Bessel functions.  Letting
 * <i>b<sub>n</sub></i>(<i>x</i>) represent any spherical Bessel function we have
 * <br/>
 * <br/>
 * &nbsp; <i>b<sub>n</i>+1</sub>(<i>x</i>) = [(2<i>n</i>+1)/<i>x</i>]<i>b<sub>n</sub></i>(<i>x</i>)
 *                                         - <i>b<sub>n</i>-1</sub>(<i>x</i>).
 * <br/>
 * <br/>
 * Again, this formula is numerically unstable for the Bessel functions of the first 
 * kind <i>j<sub>n</sub></i>(<i>x</i>).  However, it can be used to determine the
 * representation of each Bessel function in terms of trigonometric functions.
 * </p>
 * <p>
 * <h2>References</h2>
 * [1] <a href="http://www.koders.com">www.koders.com</a>
 * <br/>
 * [2]<a href="http://www.nr.com/"><i>Numerical Recipes, The Art of Scientific
 *                                    Computing, Third Edition,</i>
 *                                    W.H. Press,
 *                                    S.A. Teukolsky, 
 *                                    W.T. Vetterling, 
 *                                    B.P. Flannery
 *                                    (Cambridge University Press, Cambridge, 2007).
 *                                    </a>
 * <br/>
 * </p>
 *                                    
 * @author  Christopher K. Allen
 * 
 * @see java.lang.Math
 * @see gov.sns.tools.math.ElementaryFunction
 * 
 */
public final class BesselFunction {
    
    /*
     *  Global Attributes
     */


    /** number of Units in the Last Place (ULPs) used for bracketing approximately equal values */
    public static final int     ULPS_BRACKET = 100;

    /** the value PI/2 */
    public static final double PI_BY_2 = Math.PI/2;

    /** smallest tolerance value - not used at the moment */
    public static final double EPS = 1000.0*Double.MIN_VALUE;
    
    /** Conditional value where polynomial expansions are employed */
    public static final double SMALL_ARG = 0.1;




    /**
     * Compute the zero<sup>th</sup> order Bessel function of the
     * first kind, <i>J</i><sub>0</sub>(<i>x</i>).
     * 
     * @param x a double value
     * 
     * @return    the Bessel function of order zero in the argument,
     *            <i>J</i><sub>0</sub>(<i>x</i>).
     */
    public static double J0(final double x) {
        double ax;

        if ( (ax = Math.abs(x)) < 8.0) {
            double y = x * x;
            double ans1 = 57568490574.0 + y * ( -13362590354.0 + y * (651619640.7
                    + y * ( -11214424.18 + y * (77392.33017 + y * ( -184.9052456)))));
            double ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                    + y * (59272.64853 + y * (267.8532712 + y * 1.0))));

            return ans1 / ans2;

        }

        double z = 8.0 / ax;
        double y = z * z;
        double xx = ax - 0.785398164;
        double ans1 = 1.0 + y * ( -0.1098628627e-2 + y * (0.2734510407e-4
                + y * ( -0.2073370639e-5 + y * 0.2093887211e-6)));
        double ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                +
                y * ( -0.6911147651e-5 + y * (0.7621095161e-6
                        - y * 0.934935152e-7)));

        return Math.sqrt(0.636619772 / ax) *
        (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);    
    }

    /**
     * Compute the first-order Bessel function of the
     * first kind, <i>J</i><sub>1</sub>(<i>x</i>).
     * 
     * @param x a double value
     * @return the Bessel function of order 1 of the argument,
     *         <i>J</i><sub>1</sub>(<i>x</i>).
     */
    static public double J1(final double x) {

        double ax;
        double y;
        double ans1, ans2;

        if ( (ax = Math.abs(x)) < 8.0) {
            y = x * x;
            ans1 = x * (72362614232.0 + y * ( -7895059235.0 + y * (242396853.1
                    + y * ( -2972611.439 + y * (15704.48260 + y * ( -30.16036606))))));
            ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
                    + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
            return ans1 / ans2;
        }

        double z = 8.0 / ax;
        double xx = ax - 2.356194491;
        y = z * z;

        ans1 = 1.0 + y * (0.183105e-2 + y * ( -0.3516396496e-4
                +
                y * (0.2457520174e-5 + y * ( -0.240337019e-6))));
        ans2 = 0.04687499995 + y * ( -0.2002690873e-3
                + y * (0.8449199096e-5 + y * ( -0.88228987e-6
                        + y * 0.105787412e-6)));
        double ans = Math.sqrt(0.636619772 / ax) *
        (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
        if (x < 0.0) ans = -ans;
        return ans;    
    }

    /**
     * <h2>Arbitrary Order Bessel Function of First Kind</h2>
     * <p>
     * Computes <strong>integer</strong> order Bessel function 
     * of the first kind, <i>J</i><sub><i>n</i></sub>(<i>x</i>).
     * </p>
     * <p>
     * This implementation relies upon evaluation of the zero
     * and first order Bessel functions 
     * <i>J</i><sub>0</sub>(<i>x</i>)
     * and <i>J</i><sub>1</sub>(<i>x</i>)..
     * </p>
     * 
     * @param n integer order
     * @param x a double value
     * 
     * @return  the Bessel function of order <i>n</i> of the argument,
     *          <i>J</i><sub><i>n</i></sub>(<i>x</i>).
     * 
     * @see BesselFunction#J0(double)
     * @see BesselFunction#J1(double)
     */
    static public double Jn(final int n, final double x) {
        int j, m;
        double ax, bj, bjm, bjp, sum, tox, ans;
        boolean jsum;

        double ACC = 40.0;
        double BIGNO = 1.0e+10;
        double BIGNI = 1.0e-10;

        if (n == 0)return J0(x);
        if (n == 1)return J1(x);

        ax = Math.abs(x);
        if (ax == 0.0)return 0.0;
        else
            if (ax > n) {
                tox = 2.0 / ax;
                bjm = J0(ax);
                bj = J1(ax);
                for (j = 1; j < n; j++) {
                    bjp = j * tox * bj - bjm;
                    bjm = bj;
                    bj = bjp;
                }
                ans = bj;
            } else {
                tox = 2.0 / ax;
                m = 2 * ( (n + (int) Math.sqrt(ACC * n)) / 2);
                jsum = false;
                bjp = ans = sum = 0.0;
                bj = 1.0;
                for (j = m; j > 0; j--) {
                    bjm = j * tox * bj - bjp;
                    bjp = bj;
                    bj = bjm;
                    if (Math.abs(bj) > BIGNO) {
                        bj *= BIGNI;
                        bjp *= BIGNI;
                        ans *= BIGNI;
                        sum *= BIGNI;
                    }
                    if (jsum) sum += bj;
                    jsum = !jsum;
                    if (j == n) ans = bjp;
                }
                sum = 2.0 * sum - bj;
                ans /= sum;
            }
        return x < 0.0 && n % 2 == 1 ? -ans : ans;
    }

    /**
     * <p>
     * Compute the zero<sup>th</sup> order Bessel function of the
     * second kind, <i>Y</i><sub>0</sub>(<i>x</i>).
     * </p>
     * <p>
     * This implementation relies upon evaluation of the zero
     * order Bessel function <i>Y</i><sub>0</sub>(<i>x</i>).
     * </p>
     * 
     * @param x a double value
     * 
     * @return the Bessel function of the second kind,
     *          of order 0 of the argument, <i>Y</i><sub>0</sub>(<i>x</i>).
     *          
     * @see BesselFunction#J0(double)
     */
    static public double Y0(final double x) {

        if (x < 8.0) {
            double y = x * x;

            double ans1 = -2957821389.0 + y * (7062834065.0 + y * ( -512359803.6
                    + y * (10879881.29 + y * ( -86327.92757 + y * 228.4622733))));
            double ans2 = 40076544269.0 + y * (745249964.8 + y * (7189466.438
                    + y * (47447.26470 + y * (226.1030244 + y * 1.0))));

            return (ans1 / ans2) + 0.636619772 * J0(x) * Math.log(x);
        } 

        double z = 8.0 / x;
        double y = z * z;
        double xx = x - 0.785398164;

        double ans1 = 1.0 + y * ( -0.1098628627e-2 + y * (0.2734510407e-4
                + y * ( -0.2073370639e-5 + y * 0.2093887211e-6)));
        double ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                +
                y * ( -0.6911147651e-5 + y * (0.7621095161e-6
                        + y * ( -0.934945152e-7))));
        return Math.sqrt(0.636619772 / x) *
        (Math.sin(xx) * ans1 + z * Math.cos(xx) * ans2);   
    }

    /**
     * <p>
     * Compute the first-order Bessel function of the
     * second kind, <i>Y</i><sub>1</sub>(<i>x</i>).
     * </p>
     * <p>
     * This implementation relies upon evaluation of the first
     * order Bessel function <i>Y</i><sub>1</sub>(<i>x</i>).
     * </p>
     * 
     * @param x a double value
     * 
     * @return the Bessel function of the second kind,
     *          of order 1 of the argument, <i>Y</i><sub>1</sub>(<i>x</i>).
     *  
     * @see BesselFunction#J1(double)
     */
    static public double Y1(final double x) {

        if (x < 8.0) {
            double y = x * x;
            double ans1 = x * ( -0.4900604943e13 + y * (0.1275274390e13
                    +
                    y * ( -0.5153438139e11 + y * (0.7349264551e9
                            + y * ( -0.4237922726e7 + y * 0.8511937935e4)))));
            double ans2 = 0.2499580570e14 + y * (0.4244419664e12
                    +
                    y * (0.3733650367e10 + y * (0.2245904002e8
                            + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));
            return (ans1 / ans2) + 0.636619772 * (J1(x) * Math.log(x) - 1.0 / x);
        } 

        double z = 8.0 / x;
        double y = z * z;
        double xx = x - 2.356194491;
        double ans1 = 1.0 + y * (0.183105e-2 + y * ( -0.3516396496e-4
                +
                y * (0.2457520174e-5 + y * ( -0.240337019e-6))));
        double ans2 = 0.04687499995 + y * ( -0.2002690873e-3
                +
                y * (0.8449199096e-5 + y * ( -0.88228987e-6
                        + y * 0.105787412e-6)));
        return Math.sqrt(0.636619772 / x) *
        (Math.sin(xx) * ans1 + z * Math.cos(xx) * ans2);    
    }

    /**
     * <h2>Arbitrary Order Bessel Function of Second Kind</h2>
     * <p>
     * Computes <strong>integer</strong> order Bessel function 
     * of the second kind, <i>Y</i><sub><i>n</i></sub>(<i>x</i>).
     * </p>
     * <p>
     * This implementation relies upon evaluation of the zero
     * and first order Bessel functions 
     * <i>Y</i><sub>0</sub>(<i>x</i>)
     * and <i>Y</i><sub>1</sub>(<i>x</i>).
     * </p>
     * 
     * @param n integer order
     * @param x a double value
     * 
     * @return  the Bessel function of the second kind,
     *          of order <i>n</i> of the argument, 
     *          <i>Y</i><sub><i>n</i></sub>(<i>x</i>).
     */
    static public double Yn(final int n, final double x) {
        double by, bym, byp, tox;

        if (n == 0)return Y0(x);
        if (n == 1)return Y1(x);

        tox = 2.0 / x;
        by = Y1(x);
        bym = Y0(x);
        for (int j = 1; j < n; j++) {
            byp = j * tox * by - bym;
            bym = by;
            by = byp;
        }
        return by;
    }



    /*
     * Spherical Bessel Functions
     */
    
    /**
     *  <h2>Implementation of the sinc Function</h2>
     *  <p>
     *  The sinc function arises frequency in engineering
     *  applications, especially communications and signal
     *  processing.  It is defined as follows:
     *  <br/>
     *  <br/>
     *  &nbsp; sinc(<i>x</i>) &equiv; sin(<i>x</i>)/<i>x</i>.
     *  <br/>
     *  <br/>
     *  The function is <strong>not</strong> singular at 
     *  <i>x</i> = 0, which may easily be verified with 
     *  L'hopital's rule.
     *  </p>
     *  <p>
     *  To avoid numerical instability, for small values of 
     *  <i>x</i> we Taylor expand sinc(<i>x</i>) 
     *  to sixth order about <i>x</i> = 0.  
     *  <br/>
     *  <br/>
     *  &nbsp;  sinc(x) &asymp; 1 - <i>x</i><sup>2</sup>/6 + 
     *                               <i>x</i><sup>4</sup>/120 - 
     *                               <i>x</i><sup>6</sup>/5040 + 
     *                               <i>O</i>(<i>x</i><sup>8</sup>).
     *  <br/>
     *  <br/>
     *  otherwise we return sin(<i>x</i>)/<i>x</i>.
     *  </p>
     *  
     *  <p>
     *  <strong>NOTE</strong>: (CKA)
     *  <br/>
     *  &middot; The sinc function is the zero<sup>th</sup> order
     *  spherical bessel function <i>j</i><sub>0</sub>.
     *  </p>
     * 
     * @param   x   any real number
     * 
     * @return  sinc(<var>x</var>) &equiv; sin(<var>x</var>)/<var>x</var>
     */
    public static double sinc(final double x) {
        
        if (Math.abs(x) < BesselFunction.SMALL_ARG) {    // avoid singularity at zero
            double      x_2 = x*x;
            double      x_4 = x_2*x_2;
            
            return 1.0 - x_2/6.0 + x_4/120.0 - x_2*x_4/5040.0;
        } else {
        
            return Math.sin(x)/x;
        }
    }
    
    /**
     * <h2>Zero-Order Spherical Bessel Function of the First Kind</h2>
     * 
     * <p>
     * This function is simply an alias for the <code>sinc</code>
     * function, which is the first spherical Bessel function,
     * <i>j</i><sub>0</sub>(<i>x</i>).
     * </p>
     * 
     * @param x     any real number (double value)
     *  
     * @return      spherical Bessel function of the first kind
     *              of order 0, <i>j</i><sub>0</sub>(<i>x</i>)
     * 
     * @see ElementaryFunction#sinc(double)
     * @see BesselFunction#sinc(double)
     */
    public static double j0(final double x)   {
        
        return sinc(x);
    }
    
    /**
     * <h2>First-Order Spherical Bessel Function of the First Kind</h2>
     * 
     * <p>
     * Direct computation of the first-order spherical Bessel function of the
     * first kind, <i>j</i><sub>1</sub>(<i>x</i>) in terms of trigonometric
     * functions.
     * <p>
     * To avoid numerical instability, for small values of 
     * <i>x</i> we Taylor expand <i>j</i><sub>1</sub>(<i>x</i>) 
     * to seventh order about <i>x</i> = 0.  
     * <br/>
     * <br/>
     * &nbsp;  <i>j</i><sub>1</sub>(<i>x</i>) &asymp; 
     *                         <i>x</i>/3 - 
     *                         <i>x</i><sup>3</sup>/30 + 
     *                         <i>x</i><sup>5</sup>/840 - 
     *                         <i>x</i><sup>7</sup>/45360 + 
     *                         <i>O</i>(<i>x</i><sup>9</sup>).
     * <br/>
     * <br/>
     * otherwise we return 
     * <br/>
     * <br/>
     * <i>j</i><sub>1</sub>(<i>x</i>) = sin(<i>x</i>)/<i>x</i><sup>2</sup> - 
     *                                  cos(<i>x</i>)/<i>x</i>.
     * </p>
     *  
     * @param   x   any real number (double value)
     * 
     * @return first-order spherical Bessel function of the first kind,
     *         <i>j</i><sub>1</sub>(<i>x</i>)  
     */
    public static double j1(final double x) {

        if (Math.abs(x) < BesselFunction.SMALL_ARG) {    // avoid singularity at zero
            
            // Numerically unstable at x=0, compute expansion
            double      x_2 = x*x;
            double      x_3 = x_2*x;
            double      x_5 = x_2*x_3;
            double      x_7 = x_2*x_5;

            
            return x/3.0 - x_3/30.0 + x_5/840.0 - x_7/45360.0;
            
        } else {
        
            // Numerically stable use exact expression
            double      x_2 = x*x;
            
            return Math.sin(x)/x_2 - Math.cos(x)/x;
        }
    }
    
    /**
     * <h2>Second-Order Spherical Bessel Function of the First Kind</h2>
     * 
     * <p>
     * Direct computation of the second-order spherical Bessel function of the
     * first kind, <i>j</i><sub>2</sub>(<i>x</i>) in terms of trigonometric
     * functions.
     * <p>
     * To avoid numerical instability, for small values of 
     * <i>x</i> we Taylor expand <i>j</i><sub>2</sub>(<i>x</i>) 
     * to eighth order about <i>x</i> = 0.  
     * <br/>
     * <br/>
     * &nbsp;  <i>j</i><sub>2</sub>(<i>x</i>) &asymp; 
     *                         <i>x</i><sup>2</sup>/15 - 
     *                         <i>x</i><sup>4</sup>/210 + 
     *                         <i>x</i><sup>6</sup>/7560 - 
     *                         <i>x</i><sup>8</sup>/498960 + 
     *                         <i>O</i>(<i>x</i><sup>10</sup>).
     * <br/>
     * <br/>
     * otherwise we return 
     * <br/>
     * <br/>
     * <i>j</i><sub>2</sub>(<i>x</i>) = (3/<i>x</i> - 1)sin(<i>x</i>)/<i>x</i> - 
     *                                  3cos(<i>x</i>)/<i>x</i><sup>2</sup>.
     * </p>
     *  
     * @param   x   any real number (double value)
     * 
     * @return second-order spherical Bessel function of the first kind,
     *         <i>j</i><sub>2</sub>(<i>x</i>)  
     */
    public static double j2(final double x) {

        if (Math.abs(x) < BesselFunction.SMALL_ARG) {    // avoid singularity at zero

            // Numerically unstable at x=0, compute expansion
            double      x_2 = x*x;
            double      x_4 = x_2*x_2;
            double      x_6 = x_2*x_4;
            double      x_8 = x_2*x_6;
            
            return x_2/15.0 - x_4/210.0 + x_6/7560.0 - x_8/498960.0;
            
        } else {
        
            // Numerically stable use exact expression
            double      x_2 = x*x;

            return  (3.0/x_2 - 1.0)*Math.sin(x)/x - 3.0*Math.cos(x)/x_2;
        }
    }
    
    /**
     * <h2>Third-Order Spherical Bessel Function of the First Kind</h2>
     * 
     * <p>
     * Direct computation of the third-order spherical Bessel function of the
     * first kind, <i>j</i><sub>3</sub>(<i>x</i>) in terms of trigonometric
     * functions.
     * <p>
     * To avoid numerical instability, for small values of 
     * <i>x</i> we Taylor expand <i>j</i><sub>3</sub>(<i>x</i>) 
     * to seventh order about <i>x</i> = 0.  
     * <br/>
     * <br/>
     * &nbsp;  <i>j</i><sub>3</sub>(<i>x</i>) &asymp; 
     *                         <i>x</i><sup>3</sup>/105 + 
     *                         <i>x</i><sup>5</sup>/1890 - 
     *                         <i>x</i><sup>7</sup>/83160 + 
     *                         <i>O</i>(<i>x</i><sup>9</sup>).
     * <br/>
     * <br/>
     * otherwise we return 
     * <br/>
     * <br/>
     * <i>j</i><sub>3</sub>(<i>x</i>) = 
     *            (15/<i>x</i><sup>3</sup> - 6/<i>x</i>)sin(<i>x</i>)/<i>x</i> - 
     *            (1 - 15/<i>x</i><sup>2</sup>)cos(<i>x</i>)/<i>x</i>.
     * </p>
     *  
     * @param   x   any real number (double value)
     * 
     * @return third-order spherical Bessel function of the first kind,
     *         <i>j</i><sub>3</sub>(<i>x</i>)  
     */
    public static double j3(final double x) {

        if (Math.abs(x) < BesselFunction.SMALL_ARG) {    // avoid singularity at zero
            
            // Numerically unstable at x=0, compute expansion
            double      x_2 = x*x;
            double      x_3 = x_2*x;
            double      x_5 = x_2*x_3;
            double      x_7 = x_2*x_5;

            return x_3/105.0 - x_5/1890.0 + x_7/83160.0;
            
        } else {

            // Numerically stable, use exact expression
            double      x_2 = x*x;
            double      x_3 = x_2*x;
            
            return (15.0/x_3 - 6.0/x)*Math.sin(x)/x + (1.0 - 15.0/x_2)*Math.cos(x)/x;
        }
    }
    
    /**
     * <h2>Fourth-Order Spherical Bessel Function of the First Kind</h2>
     * 
     * <p>
     * Direct computation of the fourth-order spherical Bessel function of the
     * first kind, <i>j</i><sub>4</sub>(<i>x</i>) in terms of trigonometric
     * functions.
     * <p>
     * To avoid numerical instability, for small values of 
     * <i>x</i> we Taylor expand <i>j</i><sub>4</sub>(<i>x</i>) 
     * to eighth order about <i>x</i> = 0.  
     * <br/>
     * <br/>
     * &nbsp;  <i>j</i><sub>2</sub>(<i>x</i>) &asymp; 
     *                         <i>x</i><sup>4</sup>/945 - 
     *                         <i>x</i><sup>6</sup>/20790 + 
     *                         <i>x</i><sup>8</sup>/1081080 + 
     *                         <i>O</i>(<i>x</i><sup>10</sup>).
     * <br/>
     * <br/>
     * otherwise we return 
     * <br/>
     * <br/>
     * <i>j</i><sub>4</sub>(<i>x</i>) = 
     *    (1 - 45/<i>x</i><sup>2</sup> + 105/<i>x</i><sup>4</sup>)sin(<i>x</i>)/<i>x</i> 
     *  + (10/<i>x</i> - 105/<i>x</i><sup>3</sup>)cos(<i>x</i>)/<i>x</i>.
     * </p>
     *  
     * @param   x   any real number (double value)
     * 
     * @return fourth-order spherical Bessel function of the first kind,
     *         <i>j</i><sub>4</sub>(<i>x</i>)  
     */
    public static double j4(final double x) {

        if (Math.abs(x) < BesselFunction.SMALL_ARG) {    // avoid singularity at zero
            
            // Numerically unstable at x=0, compute expansion
            double      x_2 = x*x;
            double      x_3 = x_2*x;
            double      x_4 = x_2*x_2;
            double      x_6 = x_2*x_4;
            double      x_8 = x_2*x_6;

            return x_4/945.0 - x_6/20790.0 + x_8/1081080.0;
            
        } else {
        
            // Numerically stable, use exact expression
            double      x_2 = x*x;
            double      x_3 = x_2*x;
            double      x_4 = x_2*x_2;

            return (1.0 - 45.0/x_2 + 105.0/x_4)*Math.sin(x)/x 
                 + (10.0/x - 105.0/x_3)*Math.cos(x)/x;
        }
    }
    
    
    
    
};
