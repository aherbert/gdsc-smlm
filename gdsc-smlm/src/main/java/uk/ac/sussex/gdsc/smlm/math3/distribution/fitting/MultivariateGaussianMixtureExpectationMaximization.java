/*
 * Licensed to the Apache Software Foundation (ASF) under one or more contributor license
 * agreements. See the NOTICE file distributed with this work for additional information regarding
 * copyright ownership. The ASF licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License. You may obtain a
 * copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */

package uk.ac.sussex.gdsc.smlm.math3.distribution.fitting;

import java.util.Arrays;
import java.util.function.ToDoubleFunction;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.SingularMatrixException;
import uk.ac.sussex.gdsc.core.data.VisibleForTesting;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.MixtureMultivariateGaussianDistribution.MultivariateGaussianDistribution;

/**
 * Expectation-Maximization algorithm for fitting the parameters of multivariate Gaussian mixture
 * model distributions.
 *
 * <p>This implementation is an adaption of the Apache Commons Math 3 class
 * {@link org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization
 * MultivariateNormalMixtureExpectationMaximization}. The core functionality has been reimplemented
 * to avoid generic matrix operations that provide a new instance for the result; all matrix
 * operations are done in place. The input and output classes have been reimplemented to avoid
 * unnecessary code used for the sampling API of the Commons Math distributions. All data is passed
 * by reference and precomputation is used where possible. All changes create a speed increase over
 * the original implementation of approximately 8 to 10-fold for mixtures of 2-4 components of 2-4
 * dimension Gaussians.
 *
 * <p>The convergence criteria is specified using a predicate allowing customisation of the
 * convergence, for example using a relative threshold.
 *
 * @see <a
 *      href="https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm">Expectation–maximization
 *      algorithm (Wikipedia)</a>
 * @since 1.0
 */
public class MultivariateGaussianMixtureExpectationMaximization {
  /** Default maximum number of iterations allowed per fitting process. */
  private static final int DEFAULT_MAX_ITERATIONS = 1000;
  /** Default convergence threshold for fitting. */
  private static final double DEFAULT_THRESHOLD = 1e-5;

  /** The data to fit. */
  private final double[][] data;
  /** The model fit against the data. */
  private MixtureMultivariateGaussianDistribution fittedModel;
  /** The log likelihood of the data given the fitted model. */
  private double logLikelihood;
  /** The iterations for the most recent fitted model. */
  private int iterations;

  /**
   * Multivariate Gaussian mixture distribution.
   */
  public static final class MixtureMultivariateGaussianDistribution {
    /** The normalized weight of each mixture component. */
    final double[] weights;
    /** The mixture components. */
    final MultivariateGaussianDistribution[] distributions;

    /**
     * Implementation of the multivariate Gaussian distribution.
     *
     * @see <a href="http://en.wikipedia.org/wiki/Multivariate_normal_distribution"> Multivariate
     *      normal distribution (Wikipedia)</a>
     */
    public static final class MultivariateGaussianDistribution {
      /** The means. */
      private final double[] means;
      /** The covariance matrix. */
      private final double[][] covarianceMatrix;
      /** The matrix inverse of the covariance matrix. */
      private final double[][] covarianceMatrixInverse;
      /** The density prefactor. */
      private final double densityPrefactor;

      /**
       * Creates a multivariate Gaussian distribution with the given mean vector and covariance
       * matrix. The data is stored by reference and not copied.
       *
       * @param means the means
       * @param covariances the covariances
       * @throws SingularMatrixException if the inverse cannot be performed on the provided
       *         covariance matrix (the matrix is singular).
       * @throws NonPositiveDefiniteMatrixException if the Eigen values of the matrix are not
       *         positive
       */
      MultivariateGaussianDistribution(double[] means, double[][] covariances) {
        this.means = means;
        this.covarianceMatrix = covariances;

        // Note: Simple testing with a version based on EJML show that there is no
        // difference to the speed of the Expectation-Maximization algorithm.
        // This is left using the Commons Math inverse and determinant so the
        // distribution should match the source implementation.

        // Covariance matrix eigen decomposition.
        final EigenDecomposition covMatDec;
        try {
          covMatDec = new EigenDecomposition(new Array2DRowRealMatrix(covariances));
        } catch (MaxCountExceededException | MathArithmeticException ex) {
          throw (SingularMatrixException) new SingularMatrixException().initCause(ex);
        }

        // Compute and store the inverse.
        covarianceMatrixInverse = covMatDec.getSolver().getInverse().getData();
        // Compute and store the determinant.
        final double determinant = covMatDec.getDeterminant();

        // Eigenvalues of the covariance matrix.
        final double[] covMatEigenvalues = covMatDec.getRealEigenvalues();

        for (int i = 0; i < covMatEigenvalues.length; i++) {
          if (covMatEigenvalues[i] < 0) {
            throw new NonPositiveDefiniteMatrixException(covMatEigenvalues[i], i, 0);
          }
        }

        densityPrefactor = Math.pow(2 * Math.PI, -0.5 * means.length) * Math.pow(determinant, -0.5);
      }

      /**
       * Creates a multivariate Gaussian distribution with the given mean vector and covariance
       * matrix. The data is stored by reference and not copied.
       *
       * @param means the means
       * @param covariances the covariances
       * @return the multivariate gaussian distribution
       * @throws IllegalArgumentException if the arrays length are inconsistent, or if the inverse
       *         cannot be performed on the provided covariance matrix (the matrix is singular).
       */
      public static MultivariateGaussianDistribution create(double[] means,
          double[][] covariances) {
        final int dim = means.length;
        ValidationUtils.checkArgument(dim == covariances.length,
            "Mean and covariance matrix size mismatch: %d != %d", dim, covariances.length);
        for (int i = 0; i < dim; i++) {
          ValidationUtils.checkArgument(dim == covariances[i].length,
              "Covariance matrix size is not square: %d != %d", dim, covariances[i].length);
        }
        return new MultivariateGaussianDistribution(means, covariances);
      }

      /**
       * Gets the mean vector. Returns a reference.
       *
       * @return the mean vector.
       */
      public double[] getMeans() {
        return means;
      }

      /**
       * Gets the covariance matrix. Returns a reference.
       *
       * @return the covariance matrix.
       */
      public double[][] getCovariances() {
        return covarianceMatrix;
      }

      /**
       * Gets the square root of each element on the diagonal of the covariance matrix. Returns a
       * new array.
       *
       * @return the standard deviations.
       */
      public double[] getStandardDeviations() {
        final double[][] s = covarianceMatrix;
        final int dim = s.length;
        final double[] std = new double[dim];
        for (int i = 0; i < dim; i++) {
          std[i] = Math.sqrt(s[i][i]);
        }
        return std;
      }

      /**
       * Returns the probability density function (PDF) of this distribution evaluated at the
       * specified point x.
       *
       * @param x the values for point x
       * @return PDF(x)
       * @throws ArrayIndexOutOfBoundsException if there is a dimension mismatch between the
       *         Gaussian distribution and the point
       */
      public double density(double[] x) {
        return densityPrefactor * getExponentTerm(x);
      }

      /**
       * Computes the term used in the exponent (see definition of the distribution).
       *
       * @param values values at which to compute density.
       * @return the multiplication factor of density calculations.
       */
      private double getExponentTerm(double[] values) {
        final int n = means.length;
        final double[] centered = new double[n];
        for (int i = 0; i < n; i++) {
          centered[i] = values[i] - means[i];
        }
        final double[] preMultiplied = new double[n];
        final double[][] data = covarianceMatrixInverse;
        for (int col = 0; col < n; ++col) {
          double sum = 0;
          for (int i = 0; i < n; ++i) {
            sum += data[i][col] * centered[i];
          }
          preMultiplied[col] = sum;
        }
        double sum = 0;
        for (int i = 0; i < preMultiplied.length; i++) {
          sum += preMultiplied[i] * centered[i];
        }
        return Math.exp(-0.5 * sum);
      }
    }

    /**
     * Creates a multivariate Gaussian mixture distribution.
     *
     * @param weights Weights of each component
     * @param distributions Mixture components
     */
    MixtureMultivariateGaussianDistribution(double[] weights,
        MultivariateGaussianDistribution[] distributions) {
      this.weights = weights;
      this.distributions = distributions;
    }

    /**
     * Creates a multivariate Gaussian mixture distribution. This normalises the input weights
     * in-place to sum to 1.
     *
     * @param weights weights of each component
     * @param distributions the distributions
     * @return the mixture multivariate gaussian distribution
     * @throws IllegalArgumentException if the arrays length are inconsistent, or if the sum of the
     *         weights is infinite
     */
    public static MixtureMultivariateGaussianDistribution create(double[] weights,
        MultivariateGaussianDistribution[] distributions) {
      double sum = 0;
      final int numComp = weights.length;
      ValidationUtils.checkArgument(numComp == distributions.length,
          "Weights and distributions size mismatch: %d != %d", numComp, distributions.length);
      for (int i = 0; i < numComp; i++) {
        sum += weights[i];
      }

      // Check for overflow.
      ValidationUtils.checkArgument(Double.isFinite(sum), "sum of weights is infinite");

      // Store normalised weights.
      for (int i = 0; i < numComp; i++) {
        weights[i] /= sum;
      }
      return new MixtureMultivariateGaussianDistribution(weights, distributions);
    }

    /**
     * Creates a multivariate Gaussian mixture distribution.
     *
     * @param weights weights of each component
     * @param means means for each component
     * @param covariances covariance matrix for each component
     * @return the mixture multivariate gaussian distribution
     * @throws IllegalArgumentException if the arrays length are inconsistent, if the sum of the
     *         weights is infinite, or if the inverse cannot be performed on the provided covariance
     *         matrix (the matrix is singular).
     */
    public static MixtureMultivariateGaussianDistribution create(double[] weights, double[][] means,
        double[][][] covariances) {
      // Validate inputs
      double sum = 0;
      final int numComp = weights.length;
      final int dim = means[0].length;
      final MultivariateGaussianDistribution[] distributions =
          new MultivariateGaussianDistribution[numComp];
      for (int i = 0; i < numComp; i++) {
        final double w = weights[i];
        ValidationUtils.checkPositive(w, "weight");
        sum += w;
        final int dim2 = means[i].length;
        ValidationUtils.checkArgument(dim == dim2, "Incorrect size mean on component %d: ", i, dim);
        distributions[i] = new MultivariateGaussianDistribution(means[i], covariances[i]);
      }

      // Check for overflow.
      ValidationUtils.checkArgument(Double.isFinite(sum), "sum of weights is infinite");

      // Store normalised weights.
      final double[] normWeights = new double[numComp];
      for (int i = 0; i < numComp; i++) {
        normWeights[i] = weights[i] / sum;
      }
      return new MixtureMultivariateGaussianDistribution(normWeights, distributions);
    }

    /**
     * Returns the probability density function (PDF) of this distribution evaluated at the
     * specified point x.
     *
     * @param x the values for point x
     * @return PDF(x)
     */
    public double density(double[] x) {
      double p = 0;
      for (int i = 0; i < weights.length; i++) {
        p += weights[i] * distributions[i].density(x);
      }
      return p;
    }

    /**
     * Gets the weights of each mixture component.
     *
     * @return the weights
     */
    public double[] getWeights() {
      return weights.clone();
    }

    /**
     * Gets the distributions of each mixture component.
     *
     * @return the distributions
     */
    public MultivariateGaussianDistribution[] getDistributions() {
      return distributions.clone();
    }
  }

  /**
   * Class used for sorting user-supplied data.
   */
  private static class ProjectedData {
    /** The data. */
    final double[] data;
    /** The projected value used for the sort. */
    final double value;

    /**
     * Create an instance.
     *
     * @param data the data
     * @param value the projected value to use for the sort
     */
    ProjectedData(final double[] data, double value) {
      this.data = data;
      this.value = value;
    }
  }

  /**
   * Class used for sorting user-supplied data.
   */
  private static class ClassifiedData {
    /** The data. */
    final double[] data;
    /** The classification. */
    final int value;

    /**
     * Create an instance.
     *
     * @param data the data
     * @param value the value to use for the sort
     */
    ClassifiedData(final double[] data, int value) {
      this.data = data;
      this.value = value;
    }
  }

  /**
   * Represents a predicate (boolean-valued function) of two primitive-valued arguments. This is a
   * {@code double},{@code double}-consuming primitive type specialisation of
   * {@link java.util.function.BiPredicate}.
   *
   * <p>This is a {@link FunctionalInterface} whose functional method is
   * {@link #test(double, double)}.
   *
   * @see java.util.function.BiPredicate
   */
  @FunctionalInterface
  public interface DoubleDoubleBiPredicate {
    /**
     * Evaluates this predicate on the given arguments.
     *
     * @param value1 the first input value
     * @param value2 the second input value
     * @return {@code true} if the input arguments match the predicate, otherwise {@code false}
     */
    boolean test(double value1, double value2);
  }

  /**
   * Creates an object to fit a multivariate Gaussian mixture model to data. The data is stored by
   * reference and not copied.
   *
   * @param data data to use in fitting procedure
   * @throws IllegalArgumentException if data has no rows, if rows of data have different numbers of
   *         columns, or if the number of columns in the data is less than 2
   */
  public MultivariateGaussianMixtureExpectationMaximization(double[][] data) {
    ValidationUtils.checkStrictlyPositive(data.length, "data length");
    final int size = data[0].length;
    ValidationUtils.checkArgument(size >= 2,
        "Multivariate Gaussian requires at least 2 data columns: %d", size);
    // Check for rectangular matrix
    for (int i = 1; i < data.length; i++) {
      final int size2 = data[i].length;
      ValidationUtils.checkArgument(size == size2, "Incorrect size data on row %d: ", i, size);
    }
    this.data = data;
  }

  /**
   * Fit a mixture model to the data supplied to the constructor.
   *
   * <p>The quality of the fit depends on the concavity of the data provided to the constructor and
   * the initial mixture provided to this function. If the data has many local optima, multiple runs
   * of the fitting function with different initial mixtures may be required to find the optimal
   * solution. If a SingularMatrixException is encountered, it is possible that another
   * initialization would work.
   *
   * @param initialMixture model containing initial values of weights and multivariate normals
   * @return true if converged within the threshold
   * @throws SingularMatrixException if any component's covariance matrix is singular during fitting
   * @throws IllegalArgumentException if maxIterations is less than one or if initialMixture mean
   *         vector and data number of columns are not equal
   * @throws SingularMatrixException if any component's covariance matrix is singular during
   *         fitting.
   * @throws NonPositiveDefiniteMatrixException if any component's covariance matrix is not positive
   *         definite during fitting.
   */
  public boolean fit(MixtureMultivariateGaussianDistribution initialMixture) {
    // Normalise the log-likelihood by the number of data points.
    final int n = data.length;
    return fit(initialMixture, DEFAULT_MAX_ITERATIONS,
        (prev, curr) -> Math.abs(prev - curr) / n < DEFAULT_THRESHOLD);
  }

  /**
   * Fit a mixture model to the data supplied to the constructor.
   *
   * <p>The quality of the fit depends on the concavity of the data provided to the constructor and
   * the initial mixture provided to this function. If the data has many local optima, multiple runs
   * of the fitting function with different initial mixtures may be required to find the optimal
   * solution. If a SingularMatrixException is encountered, it is possible that another
   * initialisation would work.
   *
   * @param initialMixture model containing initial values of weights and multivariate normals
   * @param maxIterations maximum iterations allowed for fit
   * @param convergencePredicate convergence predicated used to test the logLikelihoods between
   *        successive iterations
   * @return true if converged within the threshold
   * @throws IllegalArgumentException if maxIterations is less than one or if initialMixture mean
   *         vector and data number of columns are not equal
   * @throws SingularMatrixException if any component's covariance matrix is singular during
   *         fitting.
   * @throws NonPositiveDefiniteMatrixException if any component's covariance matrix is not positive
   *         definite during fitting.
   */
  public boolean fit(MixtureMultivariateGaussianDistribution initialMixture, int maxIterations,
      DoubleDoubleBiPredicate convergencePredicate) {
    ValidationUtils.checkStrictlyPositive(maxIterations, "maxIterations");
    ValidationUtils.checkNotNull(convergencePredicate, "convergencePredicate");
    ValidationUtils.checkNotNull(initialMixture, "initialMixture");

    final int n = data.length;
    final int k = initialMixture.weights.length;

    // Number of data columns.
    final int numCols = data[0].length;
    final int numMeanColumns = initialMixture.distributions[0].means.length;
    ValidationUtils.checkArgument(numCols == numMeanColumns,
        "Mixture model dimension mismatch with data columns: %d != %d", numCols, numMeanColumns);

    logLikelihood = -Double.MAX_VALUE;
    iterations = 0;

    // Initialize model to fit to initial mixture.
    fittedModel = initialMixture;

    while (iterations++ <= maxIterations) {
      final double previousLogLikelihood = logLikelihood;
      double sumLogLikelihood = 0;

      // Weight and distribution of each component
      final double[] weights = fittedModel.weights;
      final MultivariateGaussianDistribution[] mvns = fittedModel.distributions;

      // E-step: compute the data dependent parameters of the expectation function.
      // The percentage of row's total density between a row and a component
      final double[][] gamma = new double[n][k];

      // Sum of gamma for each component
      final double[] gammaSums = new double[k];

      // Sum of gamma times its row for each each component
      final double[][] gammaDataProdSums = new double[k][numCols];

      // Cache for the weight multiplied by the distribution density
      final double[] mvnDensity = new double[k];

      for (int i = 0; i < n; i++) {
        final double[] point = data[i];
        // Compute densities for each component and the row density
        double rowDensity = 0;
        for (int j = 0; j < k; j++) {
          final double d = weights[j] * mvns[j].density(point);
          mvnDensity[j] = d;
          rowDensity += d;
        }

        sumLogLikelihood += Math.log(rowDensity);

        for (int j = 0; j < k; j++) {
          gamma[i][j] = mvnDensity[j] / rowDensity;
          gammaSums[j] += gamma[i][j];

          for (int col = 0; col < numCols; col++) {
            gammaDataProdSums[j][col] += gamma[i][j] * point[col];
          }
        }
      }

      logLikelihood = sumLogLikelihood;

      // M-step: compute the new parameters based on the expectation function.
      final double[] newWeights = new double[k];
      final double[][] newMeans = new double[k][numCols];

      for (int j = 0; j < k; j++) {
        newWeights[j] = gammaSums[j] / n;
        for (int col = 0; col < numCols; col++) {
          newMeans[j][col] = gammaDataProdSums[j][col] / gammaSums[j];
        }
      }

      // Compute new covariance matrices.
      // These are symmetric so we compute the triangular half.
      final double[][][] newCovMats = new double[k][numCols][numCols];
      final double[] vec = new double[numCols];
      for (int i = 0; i < n; i++) {
        final double[] point = data[i];
        for (int j = 0; j < k; j++) {
          subtract(point, newMeans[j], vec);
          final double g = gamma[i][j];
          // covariance = vec * vecT
          // covariance[ii][jj] = vec[ii] * vec[jj] * gamma[i][j]
          final double[][] covar = newCovMats[j];
          for (int ii = 0; ii < numCols; ii++) {
            // pre-compute
            final double vig = vec[ii] * g;
            final double[] covari = covar[ii];
            for (int jj = 0; jj <= ii; jj++) {
              covari[jj] += vig * vec[jj];
            }
          }
        }
      }

      // Converting to arrays for use by fitted model
      final MultivariateGaussianDistribution[] distributions =
          new MultivariateGaussianDistribution[k];
      for (int j = 0; j < k; j++) {
        // Make symmetric and normalise by gamma sum
        final double norm = 1.0 / gammaSums[j];
        final double[][] covar = newCovMats[j];
        for (int ii = 0; ii < numCols; ii++) {
          // diagonal
          covar[ii][ii] *= norm;
          // elements
          for (int jj = 0; jj < ii; jj++) {
            final double tmp = covar[ii][jj] * norm;
            covar[ii][jj] = tmp;
            covar[jj][ii] = tmp;
          }
        }
        distributions[j] = new MultivariateGaussianDistribution(newMeans[j], covar);
      }

      // Update current model
      fittedModel = MixtureMultivariateGaussianDistribution.create(newWeights, distributions);

      // Check convergence
      if (convergencePredicate.test(previousLogLikelihood, logLikelihood)) {
        return true;
      }
    }

    // No convergence
    return false;
  }

  /**
   * Perform element-by-element subtraction into the provided array:
   * {@code result = minuend - subtrahend}.
   *
   * @param minuend the minuend
   * @param subtrahend the subtrahend
   * @param result the result
   */
  private static void subtract(double[] minuend, double[] subtrahend, double[] result) {
    for (int i = 0; i < result.length; i++) {
      result[i] = minuend[i] - subtrahend[i];
    }
  }

  /**
   * Helper method to create a multivariate Gaussian mixture model which can be used to initialize
   * {@link #fit(MixtureMultivariateGaussianDistribution)}.
   *
   * <p>This sums the values of each data point to create a metric for ranking. It produces an
   * equivalent rank order to a projection of the point onto the unit vector where all n-dimensions
   * of the vector are {@code 1.0 / Math.sqrt(n)} and then measuring the distance of the projected
   * point from the origin, i.e. the dot product with the described unit vector.
   *
   * @param data data to estimate distribution
   * @param numComponents number of components for estimated mixture
   * @return Multivariate Gaussian mixture model estimated from the data
   * @throws IllegalArgumentException if data has less than 2 rows, if {@code numComponents < 2}, or
   *         if {@code numComponents} is greater than the number of data rows.
   * @see #estimate(double[][], int, ToDoubleFunction)
   */
  public static MixtureMultivariateGaussianDistribution estimate(double[][] data,
      int numComponents) {
    return estimate(data, numComponents, MathUtils::sum);
  }

  /**
   * Helper method to create a multivariate Gaussian mixture model which can be used to initialize
   * {@link #fit(MixtureMultivariateGaussianDistribution)}.
   *
   * <p>The function is used to extract a value for ranking all points. These are then split
   * uniformly into the specified number of components. This is equivalent to projection of the
   * points onto a line and partitioning the points along the line uniformly; this cuts the points
   * into sets using hyper-planes defined by the vector of the line. A good ranking metric is a dot
   * product with a random unit vector uniformly sampled from the surface of an n-dimension
   * hypersphere.
   *
   * <pre>
   * double[] vec = ...;
   * ToDoubleFunction&lt;double[]&gt; rankingMetric =
   *   values -&gt; {
   *     double d = 0;
   *     for (int i = 0; i &lt; vec.length; i++) {
   *       d += vec[i] * values[i];
   *     }
   *     return d;
   *   };
   * </pre>
   *
   * <p>This method can be used with the data supplied to the instance constructor to try to
   * determine a good mixture model at which to start the fit, but it is not guaranteed to supply a
   * model which will find the optimal solution or even converge.
   *
   * @param data data to estimate distribution
   * @param numComponents number of components for estimated mixture
   * @param rankingMetric the function to generate the ranking metric
   * @return Multivariate Gaussian mixture model estimated from the data
   * @throws IllegalArgumentException if data has less than 2 rows, if {@code numComponents < 2}, or
   *         if {@code numComponents} is greater than the number of data rows.
   */
  public static MixtureMultivariateGaussianDistribution estimate(double[][] data, int numComponents,
      ToDoubleFunction<double[]> rankingMetric) {
    ValidationUtils.checkArgument(data.length >= 2,
        "Estimation requires at least 2 data points: %d", data.length);
    ValidationUtils.checkArgument(numComponents >= 2,
        "Multivariate Gaussian mixture requires at least 2 components: %d", numComponents);
    ValidationUtils.checkArgument(numComponents <= data.length,
        "Number of components %d greater than data length %d", numComponents, data.length);
    ValidationUtils.checkNotNull(rankingMetric, "rankingMetric");

    final int numRows = data.length;
    final int numCols = data[0].length;
    ValidationUtils.checkArgument(numCols >= 2,
        "Multivariate Gaussian requires at least 2 data columns: %d", numCols);

    // Sort the data
    final ProjectedData[] sortedData = new ProjectedData[numRows];
    for (int i = 0; i < numRows; i++) {
      sortedData[i] = new ProjectedData(data[i], rankingMetric.applyAsDouble(data[i]));
    }
    Arrays.sort(sortedData, (o1, o2) -> Double.compare(o1.value, o2.value));

    // components of mixture model to be created
    final MultivariateGaussianDistribution[] distributions =
        new MultivariateGaussianDistribution[numComponents];

    // create a component based on data in each bin
    for (int binIndex = 0; binIndex < numComponents; binIndex++) {
      // minimum index (inclusive) from sorted data for this bin
      final int minIndex = (binIndex * numRows) / numComponents;

      // maximum index (exclusive) from sorted data for this bin
      final int maxIndex = ((binIndex + 1) * numRows) / numComponents;

      // number of data records that will be in this bin
      final int numBinRows = maxIndex - minIndex;

      // data for this bin
      final double[][] binData = new double[numBinRows][];

      // mean of each column for the data in the this bin
      final double[] columnMeans = new double[numCols];

      // populate bin and create component
      for (int i = minIndex, iBin = 0; i < maxIndex; i++, iBin++) {
        final double[] values = sortedData[i].data;
        binData[iBin] = values;
        for (int j = 0; j < numCols; j++) {
          columnMeans[j] += values[j];
        }
      }

      SimpleArrayUtils.multiply(columnMeans, 1.0 / numBinRows);

      distributions[binIndex] =
          new MultivariateGaussianDistribution(columnMeans, covariance(columnMeans, binData));
    }

    // uniform weight for each bin
    final double[] weights = SimpleArrayUtils.newDoubleArray(numComponents, 1.0 / numComponents);

    return new MixtureMultivariateGaussianDistribution(weights, distributions);
  }

  /**
   * Helper method to create a multivariate Gaussian model from the multivariate data.
   *
   * <p>This method can be used to create an unmixed model that can be compared with any mixture
   * model output from the expectation maximisation fitting algorithm.
   *
   * @param data data for the distribution
   * @return Multivariate Gaussian model created from the data
   * @throws IllegalArgumentException if data has less than 2 rows.
   */
  public static MultivariateGaussianDistribution createUnmixed(double[][] data) {
    ValidationUtils.checkArgument(data.length >= 2,
        "Estimation requires at least 2 data points: %d", data.length);
    final int numRows = data.length;
    final int numCols = data[0].length;
    ValidationUtils.checkArgument(numCols >= 2,
        "Multivariate Gaussian requires at least 2 data columns: %d", numCols);

    // mean of each column for the data
    final double[] columnMeans = new double[numCols];

    for (final double[] values : data) {
      for (int j = 0; j < numCols; j++) {
        columnMeans[j] += values[j];
      }
    }

    SimpleArrayUtils.multiply(columnMeans, 1.0 / numRows);

    return new MultivariateGaussianDistribution(columnMeans, covariance(columnMeans, data));
  }

  /**
   * Helper method to create a multivariate Gaussian mixture model which can be used to initialize
   * {@link #fit(MixtureMultivariateGaussianDistribution)}.
   *
   * <p>This method can be used with the data supplied to the instance constructor to try to
   * determine a good mixture model at which to start the fit, but it is not guaranteed to supply a
   * model which will find the optimal solution or even converge.
   *
   * <p>The weights for each component will be uniform.
   *
   * @param data data for the distribution
   * @param component the component for each data point
   * @return Multivariate Gaussian mixture model estimated from the data
   * @throws IllegalArgumentException if data has less than 2 rows, if there is a size mismatch
   *         between the data and components length, or if the number of components is less than 2.
   */
  public static MixtureMultivariateGaussianDistribution createMixed(double[][] data,
      int[] component) {
    ValidationUtils.checkArgument(data.length >= 2,
        "Estimation requires at least 2 data points: %d", data.length);
    ValidationUtils.checkArgument(data.length == component.length,
        "Data and component size mismatch: %d != %d", data.length, component.length);
    final int numRows = data.length;
    final int numCols = data[0].length;
    ValidationUtils.checkArgument(numCols >= 2,
        "Multivariate Gaussian requires at least 2 data columns: %d", numCols);

    // Sort the data
    final ClassifiedData[] sortedData = new ClassifiedData[data.length];
    for (int i = 0; i < numRows; i++) {
      sortedData[i] = new ClassifiedData(data[i], component[i]);
    }
    Arrays.sort(sortedData, (o1, o2) -> Integer.compare(o1.value, o2.value));

    ValidationUtils.checkArgument(sortedData[0].value != sortedData[sortedData.length - 1].value,
        "Mixture model requires at least 2 data components");

    // components of mixture model to be created
    final LocalList<MultivariateGaussianDistribution> distributions = new LocalList<>();


    int from = 0;
    while (from < sortedData.length) {
      // Find the end
      int to = from + 1;
      final int comp = sortedData[from].value;
      while (to < sortedData.length && sortedData[to].value == comp) {
        to++;
      }

      // number of data records that will be in this component
      final int numCompRows = to - from;

      // data for this component
      final double[][] compData = new double[numCompRows][];

      // mean of each column for the data
      final double[] columnMeans = new double[numCols];

      // populate and create component
      int count = 0;
      for (int i = from; i < to; i++) {
        final double[] values = sortedData[i].data;
        compData[count++] = values;
        for (int j = 0; j < numCols; j++) {
          columnMeans[j] += values[j];
        }
      }

      SimpleArrayUtils.multiply(columnMeans, 1.0 / numCompRows);
      distributions.add(
          new MultivariateGaussianDistribution(columnMeans, covariance(columnMeans, compData)));

      from = to;
    }

    // uniform weight for each bin
    final int numComponents = distributions.size();
    final double[] weights = SimpleArrayUtils.newDoubleArray(numComponents, 1.0 / numComponents);

    return new MixtureMultivariateGaussianDistribution(weights,
        distributions.toArray(new MultivariateGaussianDistribution[0]));
  }

  /**
   * Compute the covariance of the data with columns representing covariates.
   *
   * @param means the column means
   * @param data the data
   * @return the covariance
   */
  @VisibleForTesting
  static double[][] covariance(double[] means, double[][] data) {
    final int cols = data[0].length;
    final double[][] covar = new double[cols][cols];
    for (int i = 0; i < cols; i++) {
      for (int j = 0; j < i; j++) {
        final double cov = covariance(means, data, i, j);
        covar[i][j] = cov;
        covar[j][i] = cov;
      }
      covar[i][i] = variance(means, data, i);
    }
    return covar;
  }

  /**
   * Compute the covariance between column 1 and 2.
   *
   * @param means the means of the columns
   * @param data the data
   * @param col1 the column 1
   * @param col2 the column 2
   * @return the covariance
   */
  private static double covariance(double[] means, double[][] data, int col1, int col2) {
    final int length = data.length;
    final double mean1 = means[col1];
    final double mean2 = means[col2];
    double result = 0;
    for (int row = 0; row < length; row++) {
      final double x = data[row][col1] - mean1;
      final double y = data[row][col2] - mean2;
      // Rolling mean of (X-meanX)*(Y-meanY).
      // This avoids overflow of a rolling sum then normalisation to the mean at the end.
      result += (x * y - result) / (row + 1);
    }
    // Bias corrected
    return result * (length / (length - 1.0));
  }

  /**
   * Compute the variance within the column.
   *
   * @param means the means of the columns
   * @param data the data
   * @param col the column
   * @return the variance
   */
  private static double variance(double[] means, double[][] data, int col) {
    final int length = data.length;
    final double mean = means[col];
    double sum = 0;
    double sumSq = 0;
    for (int row = 0; row < length; row++) {
      final double x = data[row][col] - mean;
      sum += x;
      sumSq += x * x;
    }
    // Bias corrected
    return (sumSq - (sum * sum / length)) / (length - 1.0);
  }

  /**
   * Gets the log likelihood of the data under the fitted model. This is the sum of the log
   * likelihood of each point.
   *
   * <p>Note: In contrast to the the Apache Commons Math class this is not divided by the number of
   * points (i.e. the mean log likelihood).
   *
   * @return log likelihood of data or zero of no data has been fit
   */
  public double getLogLikelihood() {
    return logLikelihood;
  }

  /**
   * Gets the number of iterations used for the last fit.
   *
   * @return number of iterations or zero of no data has been fit
   */
  public int getIterations() {
    return iterations;
  }

  /**
   * Gets the fitted model.
   *
   * @return fitted model or {@code null} if no fit has been performed yet.
   */
  public MixtureMultivariateGaussianDistribution getFittedModel() {
    return fittedModel;
  }
}
