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

package uk.ac.sussex.gdsc.smlm.math3.distribution.fitting;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;
import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousUniformSampler;
import org.apache.commons.rng.sampling.distribution.SharedStateContinuousSampler;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.MixtureMultivariateGaussianDistribution;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.MixtureMultivariateGaussianDistribution.MultivariateGaussianDistribution;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class MultivariateGaussianMixtureExpectationMaximizationTest {
  /** Convergence checker to match the Commons Math 3 default on absolute value. */
  // @formatter:off
  private static final MultivariateGaussianMixtureExpectationMaximization.DoubleDoubleBiPredicate
      DEFAULT_CONVERGENCE_CHECKER = Predicates.doublesAreAbsolutelyClose(1e-5)::test;
  // @formatter:on

  @SeededTest
  void canComputeCovariance(RandomSeed seed) {
    final UniformRandomProvider rng = RngFactory.create(seed.get());

    final int rows = 20;
    final int cols = 3;
    final double[][] data = new double[rows][cols];
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        data[i][j] = rng.nextDouble();
      }
    }
    // Compute the same mean as the Covariance class
    final double[] means = getColumnMeans(data);

    final double[][] obs =
        MultivariateGaussianMixtureExpectationMaximization.covariance(means, data);
    final double[][] exp = new Covariance(data).getCovarianceMatrix().getData();

    Assertions.assertArrayEquals(exp, obs);
  }

  @Test
  void testCreateMultivariateGaussianDistributionThrows() {
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianDistribution.create(new double[2], new double[3][3]);
    });
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianDistribution.create(new double[2], new double[2][3]);
    });
    // A non-positive definite matrix
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      final double[] means = {3.4, 5.6};
      final double[][] covariances = {{1.1, 2.3}, {2.3, 2.0}};
      MultivariateGaussianDistribution.create(means, covariances);
    });
    // Failure of Eigen decomposition
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      final double[] means =
          {3.2394652282161466E-159, 0.9874091290326361, 0.3148579631286939, 0.3718812084490761};
      final double[][] covariances = {
          {4.7656408265957816E-164, 1.1864890699789289E-160, 3.486364577999923E-160,
              2.799059273090178E-160},
          {1.1864890699789289E-160, 0.24556306971562533, 0.03418519684554344, 0.04106282423640166},
          {3.486364577999923E-160, 0.03418519684554344, 0.005991172243435974, 0.007750670597351695},
          {2.799059273090178E-160, 0.04106282423640166, 0.007750670597351695, 0.03725445634104687}};
      MultivariateGaussianDistribution.create(means, covariances);
    });
  }

  @Test
  void canCreateMultivariateGaussianDistribution() {
    final double[][] data = {{1, 2}, {2.5, 1.5}, {3.5, 1.0}};
    final double[] means = getColumnMeans(data);
    final double[][] covariances = getCovariance(data);
    final MultivariateGaussianDistribution dist =
        MultivariateGaussianDistribution.create(means, covariances);
    Assertions.assertSame(means, dist.getMeans());
    Assertions.assertSame(covariances, dist.getCovariances());
    final double[] sd = dist.getStandardDeviations();
    Assertions.assertEquals(covariances.length, sd.length);
    for (int i = 0; i < sd.length; i++) {
      Assertions.assertEquals(Math.sqrt(covariances[i][i]), sd[i]);
    }
    // Test against Apache commons
    final MultivariateNormalDistribution expDist =
        new MultivariateNormalDistribution(means, covariances);
    for (final double[] x : data) {
      Assertions.assertEquals(expDist.density(x), dist.density(x));
    }
  }

  @Test
  void testCreateUnmixedMultivariateGaussianDistributionThrows() {
    Assertions.assertThrows(NullPointerException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createUnmixed(null);
    });
    // Not enough data
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createUnmixed(new double[0][]);
    });
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createUnmixed(new double[1][]);
    });
    // Not enough columns
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createUnmixed(new double[2][1]);
    });
  }

  @Test
  void canCreateUnmixedMultivariateGaussianDistribution() {
    final double[][] data = {{1, 2}, {2.5, 1.5}, {3.5, 1.0}};
    final double[] means = getColumnMeans(data);
    final double[][] covariances = getCovariance(data);
    final MultivariateGaussianDistribution exp =
        MultivariateGaussianDistribution.create(means, covariances);
    final MultivariateGaussianDistribution obs =
        MultivariateGaussianMixtureExpectationMaximization.createUnmixed(data);
    final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-6);
    TestAssertions.assertArrayTest(exp.getMeans(), obs.getMeans(), test);
    TestAssertions.assertArrayTest(exp.getCovariances(), obs.getCovariances(), test);
  }

  @Test
  void testCreateMixedMultivariateGaussianDistributionThrows() {
    final double[][] data = new double[2][2];
    final int[] component = {0, 1};
    Assertions.assertThrows(NullPointerException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(null, component);
    });
    Assertions.assertThrows(NullPointerException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(data, null);
    });
    // Not enough data
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(new double[0][2], component);
    });
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(new double[1][2], component);
    });
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(new double[3][2], component);
    });
    // Not enough columns
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(new double[2][1], component);
    });
    // Array length mismatch
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(data, new int[1]);
    });
    // Not enough components
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.createMixed(data, new int[2]);
    });
  }

  @SeededTest
  void canCreateMixedMultivariateGaussianDistribution(RandomSeed seed) {
    // Will be normalised
    final double[] weights = {1, 1};
    final double[][] means = new double[2][];
    final double[][][] covariances = new double[2][][];
    final double[][] data1 = {{1, 2}, {2.5, 1.5}, {3.5, 1.0}};
    final double[][] data2 = {{4, 2}, {3.5, -1.5}, {-3.5, 1.0}};
    means[0] = getColumnMeans(data1);
    covariances[0] = getCovariance(data1);
    means[1] = getColumnMeans(data2);
    covariances[1] = getCovariance(data2);

    // Create components. This does not have to be zero based.
    final LocalList<double[]> list = new LocalList<>();
    list.addAll(Arrays.asList(data1));
    list.addAll(Arrays.asList(data2));
    final double[][] data = list.toArray(new double[0][]);
    final int[] components = {-1, -1, -1, 3, 3, 3};

    // Randomise the data
    for (int n = 0; n < 3; n++) {
      final long start = n + seed.getAsLong();
      // This relies on the shuffle being the same
      RandomUtils.shuffle(data, RngFactory.create(start));
      RandomUtils.shuffle(components, RngFactory.create(start));

      final MixtureMultivariateGaussianDistribution dist =
          MultivariateGaussianMixtureExpectationMaximization.createMixed(data, components);

      Assertions.assertArrayEquals(new double[] {0.5, 0.5}, dist.getWeights());
      final MultivariateGaussianDistribution[] distributions = dist.getDistributions();
      Assertions.assertEquals(weights.length, distributions.length);
      final DoubleDoubleBiPredicate test = Predicates.doublesAreRelativelyClose(1e-8);
      for (int i = 0; i < means.length; i++) {
        TestAssertions.assertArrayTest(means[i], distributions[i].getMeans(), test);
        TestAssertions.assertArrayTest(covariances[i], distributions[i].getCovariances(), test);
      }
    }
  }

  @Test
  void testCreateMixtureMultivariateGaussianDistributionThrows() {
    // Will be normalised
    final double[] weights = {1, 3};
    final double[][] means = new double[2][];
    final double[][][] covariances = new double[2][][];
    final double[][] data = {{1, 2}, {2.5, 1.5}, {3.5, 1.0}};
    final double[][] data2 = {{4, 2}, {3.5, -1.5}, {-3.5, 1.0}};
    means[0] = getColumnMeans(data);
    covariances[0] = getCovariance(data);
    means[1] = getColumnMeans(data2);
    covariances[1] = getCovariance(data2);

    // Test with means and covariances

    // Non positive weights
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MixtureMultivariateGaussianDistribution.create(new double[] {-1, 1}, means, covariances);
    });
    // Weights sum is not finite
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MixtureMultivariateGaussianDistribution
          .create(new double[] {Double.MAX_VALUE, Double.MAX_VALUE}, means, covariances);
    });
    // Incorrect size mean
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      final double[][] means2 = means.clone();
      means2[0] = Arrays.copyOf(means2[0], means2[0].length + 1);
      MixtureMultivariateGaussianDistribution.create(weights, means2, covariances);
    });
    // Bad covariance matrix
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      final double[][][] covariances2 = covariances.clone();
      covariances2[0] = new double[3][2];
      MixtureMultivariateGaussianDistribution.create(weights, means, covariances2);
    });

    // Test with the weights and distributions.
    // Create a valid mixture to use the distributions.
    final MixtureMultivariateGaussianDistribution dist =
        MixtureMultivariateGaussianDistribution.create(weights, means, covariances);

    // Weights sum is not finite
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MixtureMultivariateGaussianDistribution
          .create(new double[] {Double.MAX_VALUE, Double.MAX_VALUE}, dist.getDistributions());
    });
    // Weights and distributions length mismatch
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MixtureMultivariateGaussianDistribution.create(new double[] {1}, dist.getDistributions());
    });
  }

  @Test
  void canCreateMixtureMultivariateGaussianDistribution() {
    // Will be normalised
    final double[] weights = {1, 3};
    final double[][] means = new double[2][];
    final double[][][] covariances = new double[2][][];
    final double[][] data = {{1, 2}, {2.5, 1.5}, {3.5, 1.0}};
    final double[][] data2 = {{4, 2}, {3.5, -1.5}, {-3.5, 1.0}};
    means[0] = getColumnMeans(data);
    covariances[0] = getCovariance(data);
    means[1] = getColumnMeans(data2);
    covariances[1] = getCovariance(data2);
    final MixtureMultivariateGaussianDistribution dist =
        MixtureMultivariateGaussianDistribution.create(weights, means, covariances);
    Assertions.assertArrayEquals(new double[] {0.25, 0.75}, dist.getWeights());
    final MultivariateGaussianDistribution[] distributions = dist.getDistributions();
    Assertions.assertEquals(weights.length, distributions.length);
    for (int i = 0; i < means.length; i++) {
      Assertions.assertArrayEquals(means[i], distributions[i].getMeans());
      Assertions.assertArrayEquals(covariances[i], distributions[i].getCovariances());
    }

    // Test against Apache commons
    final MixtureMultivariateNormalDistribution expDist =
        new MixtureMultivariateNormalDistribution(weights, means, covariances);
    for (final double[] x : data) {
      Assertions.assertEquals(expDist.density(x), dist.density(x), 1e-10);
    }

    // Test the package private create method normalises the weights
    Assertions.assertArrayEquals(new double[] {1, 3}, weights);
    final MixtureMultivariateGaussianDistribution dist2 =
        MixtureMultivariateGaussianDistribution.create(weights, distributions);
    // Stored by reference
    Assertions.assertArrayEquals(weights, dist2.getWeights());
    // Normalised in-place
    Assertions.assertArrayEquals(new double[] {0.25, 0.75}, weights);
  }

  @Test
  void testEstimateInitialMixtureThrows() {
    // Not enough data
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.estimate(new double[1][2], 2);
    });
    // Not enough components
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.estimate(new double[2][2], 1);
    });
    // Components > data length
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.estimate(new double[2][2], 3);
    });
    // Not enough columns
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      MultivariateGaussianMixtureExpectationMaximization.estimate(new double[2][1], 2);
    });
  }

  @SeededTest
  void canEstimateInitialMixture(RandomSeed seed) {
    // Test verses the Commons Math estimation
    final UniformRandomProvider rng = RngFactory.create(seed.get());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreClose(1e-5, 1e-16);
    // Number of components
    for (int n = 2; n <= 3; n++) {
      final double[] sampleWeights = createWeights(n, rng);
      final double[][] sampleMeans = create(n, 2, rng, -5, 5);
      final double[][] sampleStdDevs = create(n, 2, rng, 1, 10);
      final double[] sampleCorrelations = create(n, rng, -0.9, 0.9);
      final double[][] data =
          createData2d(1000, rng, sampleWeights, sampleMeans, sampleStdDevs, sampleCorrelations);

      final MixtureMultivariateGaussianDistribution model1 =
          MultivariateGaussianMixtureExpectationMaximization.estimate(data, n);
      final MixtureMultivariateNormalDistribution model2 =
          MultivariateNormalMixtureExpectationMaximization.estimate(data, n);
      final List<Pair<Double, MultivariateNormalDistribution>> comp = model2.getComponents();
      final double[] weights = model1.getWeights();
      final MultivariateGaussianDistribution[] distributions = model1.getDistributions();
      Assertions.assertEquals(n, comp.size());
      Assertions.assertEquals(n, weights.length);
      Assertions.assertEquals(n, distributions.length);
      for (int i = 0; i < n; i++) {
        // Must be binary equal for estimated model
        Assertions.assertEquals(comp.get(i).getFirst(), weights[i], "weight");
        final MultivariateNormalDistribution d = comp.get(i).getSecond();
        TestAssertions.assertArrayTest(d.getMeans(), distributions[i].getMeans(), test, "means");
        TestAssertions.assertArrayTest(d.getCovariances().getData(),
            distributions[i].getCovariances(), test, "covariances");
      }
    }
  }

  @SuppressWarnings("unused")
  @Test
  void testCreateMultivariateGaussianMixtureExpectationMaximizationThrows() {
    // No data
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new MultivariateGaussianMixtureExpectationMaximization(new double[0][2]);
    });
    // Data of dimension 1
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new MultivariateGaussianMixtureExpectationMaximization(new double[2][1]);
    });
    // Non rectangular matrix
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      new MultivariateGaussianMixtureExpectationMaximization(new double[][] {{0, 1, 2}, {3, 4}});
    });
  }

  @Test
  void testFitThrows() {
    // This does not matter for the initial checks
    final MixtureMultivariateGaussianDistribution initialMixture = null;
    final MultivariateGaussianMixtureExpectationMaximization fitter =
        new MultivariateGaussianMixtureExpectationMaximization(new double[2][3]);
    // Test initial model is null and the likelihood is zero
    Assertions.assertEquals(0, fitter.getLogLikelihood());
    Assertions.assertEquals(0, fitter.getIterations());
    Assertions.assertNull(fitter.getFittedModel());
    // Valid parameters
    final int maxIterations = 10;
    // Not positive iterations
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      fitter.fit(initialMixture, 0, DEFAULT_CONVERGENCE_CHECKER);
    });
    // Not convergence checker
    Assertions.assertThrows(NullPointerException.class, () -> {
      fitter.fit(initialMixture, maxIterations, null);
    });
    // Null mixture
    Assertions.assertThrows(NullPointerException.class, () -> {
      fitter.fit(initialMixture, maxIterations, DEFAULT_CONVERGENCE_CHECKER);
    });
    // Incorrect data dimensions. Create a 50-50 mixture of 2D Gaussians
    final MixtureMultivariateGaussianDistribution initialMixture2 =
        new MixtureMultivariateGaussianDistribution(new double[] {0.5, 0.5},
            new MultivariateGaussianDistribution[] {
                new MultivariateGaussianDistribution(new double[2],
                    new double[][] {{1, 0}, {0, 2}}),
                new MultivariateGaussianDistribution(new double[2],
                    new double[][] {{1, 0}, {0, 2}})});
    Assertions.assertThrows(IllegalArgumentException.class, () -> {
      fitter.fit(initialMixture2, maxIterations, DEFAULT_CONVERGENCE_CHECKER);
    });
  }

  @SeededTest
  void canFit(RandomSeed seed) {
    // Test verses the Commons Math estimation
    final UniformRandomProvider rng = RngFactory.create(seed.get());
    final DoubleDoubleBiPredicate test = Predicates.doublesAreClose(1e-5, 1e-16);
    final int sampleSize = 1000;
    // Number of components
    for (int n = 2; n <= 3; n++) {
      final double[] sampleWeights = createWeights(n, rng);
      final double[][] sampleMeans = create(n, 2, rng, -5, 5);
      final double[][] sampleStdDevs = create(n, 2, rng, 1, 10);
      final double[] sampleCorrelations = create(n, rng, -0.9, 0.9);
      final double[][] data = createData2d(sampleSize, rng, sampleWeights, sampleMeans,
          sampleStdDevs, sampleCorrelations);

      final MixtureMultivariateGaussianDistribution initialModel1 =
          MultivariateGaussianMixtureExpectationMaximization.estimate(data, n);
      final MultivariateGaussianMixtureExpectationMaximization fitter1 =
          new MultivariateGaussianMixtureExpectationMaximization(data);
      Assertions.assertTrue(fitter1.fit(initialModel1));

      final MultivariateNormalMixtureExpectationMaximization fitter2 =
          new MultivariateNormalMixtureExpectationMaximization(data);
      fitter2.fit(MultivariateNormalMixtureExpectationMaximization.estimate(data, n));

      final double ll1 = fitter1.getLogLikelihood() / sampleSize;
      Assertions.assertNotEquals(0, ll1);
      final double ll2 = fitter2.getLogLikelihood();
      TestAssertions.assertTest(ll2, ll1, test);

      final MixtureMultivariateGaussianDistribution model1 = fitter1.getFittedModel();
      Assertions.assertNotNull(model1);
      final MixtureMultivariateNormalDistribution model2 = fitter2.getFittedModel();

      // Check fitted models are the same
      final List<Pair<Double, MultivariateNormalDistribution>> comp = model2.getComponents();
      final double[] weights = model1.getWeights();
      final MultivariateGaussianDistribution[] distributions = model1.getDistributions();
      Assertions.assertEquals(n, comp.size());
      Assertions.assertEquals(n, weights.length);
      Assertions.assertEquals(n, distributions.length);
      for (int i = 0; i < n; i++) {
        TestAssertions.assertTest(comp.get(i).getFirst(), weights[i], test, "weight");
        final MultivariateNormalDistribution d = comp.get(i).getSecond();
        TestAssertions.assertArrayTest(d.getMeans(), distributions[i].getMeans(), test, "means");
        TestAssertions.assertArrayTest(d.getCovariances().getData(),
            distributions[i].getCovariances(), test, "covariances");
      }

      final int iterations = fitter1.getIterations();
      Assertions.assertNotEquals(0, iterations);
      // Test without convergence
      if (iterations > 2) {
        Assertions.assertFalse(fitter1.fit(initialModel1, 2, DEFAULT_CONVERGENCE_CHECKER));
      }
    }
  }

  /**
   * Gets the column means. This is done using the same method as the means in the Apache Commons
   * Math Covariance class.
   *
   * @param data the data
   * @return the column means
   */
  private static double[] getColumnMeans(double[][] data) {
    final Array2DRowRealMatrix m = new Array2DRowRealMatrix(data);
    final Mean mean = new Mean();
    return IntStream.range(0, data[0].length).mapToDouble(i -> mean.evaluate(m.getColumn(i)))
        .toArray();
  }

  /**
   * Gets the covariance.
   *
   * @param data the data
   * @return the covariance
   */
  private static double[][] getCovariance(double[][] data) {
    return new Covariance(data).getCovarianceMatrix().getData();
  }

  /**
   * Creates the weights.
   *
   * @param n the number of weights
   * @param rng the random generator
   * @return the weights
   */
  private static double[] createWeights(int n, UniformRandomProvider rng) {
    final double[] weights = new double[n];
    double sum = 0;
    for (int i = 0; i < n; i++) {
      weights[i] = rng.nextDouble();
      sum += weights[i];
    }
    if (sum == 0) {
      return createWeights(n, rng);
    }
    SimpleArrayUtils.multiply(weights, 1.0 / sum);
    return weights;
  }

  /**
   * Creates m n-dimension arrays with values in the given range.
   *
   * @param n the number to generate
   * @param n the array dimension
   * @param rng the random generator
   * @param lower the lower
   * @param upper the upper
   * @return the data
   */
  private static double[][] create(int m, int n, UniformRandomProvider rng, double lower,
      double upper) {
    final double[][] data = new double[m][];
    for (int i = 0; i < m; i++) {
      data[i] = create(n, rng, lower, upper);
    }
    return data;
  }

  /**
   * Creates a n-dimension array with values in the given range.
   *
   * @param n the array dimension
   * @param rng the random generator
   * @param lower the lower
   * @param upper the upper
   * @return the data
   */
  private static double[] create(int n, UniformRandomProvider rng, double lower, double upper) {
    final SharedStateContinuousSampler sampler = ContinuousUniformSampler.of(rng, lower, upper);
    final double[] data = new double[n];
    for (int i = 0; i < n; i++) {
      data[i] = sampler.sample();
    }
    return data;
  }

  /**
   * Creates the data from a mixture of n 2D Gaussian distributions. The length of the weights array
   * (and all other arrays) is the number of mixture components.
   *
   * @param sampleSize the sample size
   * @param rng the random generator
   * @param weights the weights for each component
   * @param means the means for the x and y dimensions
   * @param stdDevs the std devs for the x and y dimensions
   * @param correlations the correlations between the x and y dimensions
   * @return the double[][]
   */
  private static double[][] createData2d(int sampleSize, UniformRandomProvider rng,
      double[] weights, double[][] means, double[][] stdDevs, double[] correlations) {
    // Use Commons Math for sampling
    final ArrayList<Pair<Double, MultivariateNormalDistribution>> components = new ArrayList<>();
    for (int i = 0; i < weights.length; i++) {
      // Create covariance matrix
      final double sx = stdDevs[i][0];
      final double sy = stdDevs[i][1];
      final double sxsy = correlations[i] * sx * sy;
      final double[][] covar = new double[][] {{sx * sx, sxsy}, {sxsy, sy * sy}};
      components.add(new Pair<>(weights[i], new MultivariateNormalDistribution(means[i], covar)));
    }
    final MixtureMultivariateNormalDistribution dist =
        new MixtureMultivariateNormalDistribution(new RandomGeneratorAdapter(rng), components);
    return dist.sample(sampleSize);
  }
}
