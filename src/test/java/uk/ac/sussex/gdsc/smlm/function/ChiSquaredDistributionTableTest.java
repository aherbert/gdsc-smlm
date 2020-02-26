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

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.rng.UniformRandomProvider;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.GdscSmlmTestUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
public class ChiSquaredDistributionTableTest {
  // Taken from:
  // https://en.wikipedia.org/wiki/Chi-squared_distribution#Table_of_.CF.872_values_vs_p-values
  //@formatter:off
  double[] p = {0.95,0.90,0.80,0.70,0.50,0.30,0.20,0.10,0.05,0.01,0.001};
  double[][] chi2 = new double[][] {
        {0.004,0.02,0.06,0.15,0.46,1.07,1.64,2.71,3.84,6.64,10.83},
        {0.10,0.21,0.45,0.71,1.39,2.41,3.22,4.60,5.99,9.21,13.82},
        {0.35,0.58,1.01,1.42,2.37,3.66,4.64,6.25,7.82,11.34,16.27},
        {0.71,1.06,1.65,2.20,3.36,4.88,5.99,7.78,9.49,13.28,18.47},
        {1.14,1.61,2.34,3.00,4.35,6.06,7.29,9.24,11.07,15.09,20.52},
        {1.63,2.20,3.07,3.83,5.35,7.23,8.56,10.64,12.59,16.81,22.46},
        {2.17,2.83,3.82,4.67,6.35,8.38,9.80,12.02,14.07,18.48,24.32},
        {2.73,3.49,4.59,5.53,7.34,9.52,11.03,13.36,15.51,20.09,26.12},
        {3.32,4.17,5.38,6.39,8.34,10.66,12.24,14.68,16.92,21.67,27.88},
        {3.94,4.87,6.18,7.27,9.34,11.78,13.44,15.99,18.31,23.21,29.59},
  };
  //@formatter:on

  @SeededTest
  public void canComputeProbability() {
    for (final int df : new int[] {5, 10}) {
      double obs;
      double exp;
      double chi = 0;
      final ChiSquaredDistribution d = new ChiSquaredDistribution(null, df);

      obs = ChiSquaredDistributionTable.computePValue(chi, df);
      exp = d.cumulativeProbability(chi);
      Assertions.assertEquals(exp, obs, 1e-10);

      chi = 1;
      for (int i = 0; i < 10; i++, chi *= 2) {
        obs = ChiSquaredDistributionTable.computePValue(chi, df);
        exp = d.cumulativeProbability(chi);
        Assertions.assertEquals(exp, obs, 1e-10);

        obs = ChiSquaredDistributionTable.computeQValue(chi, df);
        exp = 1 - exp;
        Assertions.assertEquals(exp, obs, 1e-10);
      }
    }
  }

  @SeededTest
  public void canComputeChiSquared() {
    // We have to use the transpose of the table
    final DenseMatrix64F m = new DenseMatrix64F(chi2);
    CommonOps.transpose(m);
    final int max = m.numCols;
    final double[] et = m.data;
    for (int i = 0, j = 0; i < p.length; i++) {
      final ChiSquaredDistributionTable upperTable =
          ChiSquaredDistributionTable.createUpperTailed(p[i], max);
      // Use 1-p as the significance level to get the same critical values
      final ChiSquaredDistributionTable lowerTable =
          ChiSquaredDistributionTable.createLowerTailed(1 - p[i], max);
      for (int df = 1; df <= max; df++) {
        final double o = upperTable.getCrititalValue(df);
        final double e = et[j++];
        // logger.fine(FunctionUtils.getSupplier("p=%.3f,df=%d = %f", p[i], df, o);
        Assertions.assertEquals(e, o, 1e-2);

        // The test only stores 2 decimal places so use the computed value to set upper/lower
        final double upper = o * 1.01;
        final double lower = o * 0.99;

        Assertions.assertTrue(upperTable.reject(upper, df), "Upper did not reject higher");
        Assertions.assertFalse(upperTable.reject(o, df), "Upper did not reject actual value");
        Assertions.assertFalse(upperTable.reject(lower, df), "Upper did not accept lower");

        Assertions.assertTrue(lowerTable.reject(lower, df), "Lower did not reject lower");
        Assertions.assertFalse(lowerTable.reject(o, df), "Lower did not accept actual value");
        Assertions.assertFalse(lowerTable.reject(upper, df), "Lower did not accept higher");
      }
    }
  }

  @SeededTest
  public void canPerformChiSquaredTest(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());
    final ChiSquareTest test = new ChiSquareTest();
    for (final int n : new int[] {10, 50, 100}) {
      final double[] x = SimpleArrayUtils.newArray(n, 0.5, 1.0);
      final long[] l = new long[x.length];
      for (int i = 0; i < x.length; i++) {
        l[i] = GdscSmlmTestUtils.createPoissonSampler(rng, x[i]).sample();
      }
      final double chi2 = test.chiSquare(x, l);
      final double ep = test.chiSquareTest(x, l);
      final int df = x.length - 1;
      final double o = ChiSquaredDistributionTable.computeQValue(chi2, df);
      Assertions.assertEquals(ep, o, 1e-10);

      final ChiSquaredDistributionTable upperTable =
          ChiSquaredDistributionTable.createUpperTailed(o, df);

      final double upper = chi2 * 1.01;
      final double lower = chi2 * 0.99;

      Assertions.assertTrue(upperTable.reject(upper, df), "Upper did not reject higher");
      Assertions.assertFalse(upperTable.reject(o, df), "Upper did not reject actual value");
      Assertions.assertFalse(upperTable.reject(lower, df), "Upper did not accept lower");
    }
  }
}
