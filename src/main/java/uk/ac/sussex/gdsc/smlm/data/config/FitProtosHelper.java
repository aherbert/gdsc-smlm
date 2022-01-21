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

package uk.ac.sussex.gdsc.smlm.data.config;

import uk.ac.sussex.gdsc.core.utils.NoiseEstimator.Method;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterType;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FilterSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitEngineSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolverSettings;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.LineSearchMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.NoiseEstimatorMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.RelativeParameter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.SearchMethod;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.FastMleSteppingFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;

/**
 * Contains helper functions for the FitProtos class.
 */
public final class FitProtosHelper {

  /** The default FitSolverSettings. */
  public static final FitSolverSettings defaultFitSolverSettings;

  static {
    final FitSolverSettings.Builder builder = FitSolverSettings.newBuilder();
    builder.setFixedPsf(false);
    builder.setDisableBackgroundFitting(false);
    builder.setDisableSignalFitting(false);

    builder.setFitSolver(FitSolver.LVM_LSE);
    builder.setFixedIterations(false);
    builder.setMaxIterations(20);
    builder.setRelativeThreshold(1e-6);
    builder.setAbsoluteThreshold(1e-16);
    // Disabled. These are rarely used
    //builder.setParameterRelativeThreshold(1e-3);
    //builder.setParameterAbsoluteThreshold(1e-6);

    builder.setLambda(10);

    builder.setSearchMethod(SearchMethod.POWELL_BOUNDED);
    builder.setGradientLineMinimisation(false);
    builder.setModelCamera(false);
    builder.setMaxFunctionEvaluations(2000);

    builder.setUseClamping(false);
    builder.setUseDynamicClamping(false);
    // Add defaults for a two-axis and theta Gaussian 2D function.
    // The units are photons and pixels.
    // B (3D DAOSTORM uses 100 which is high compared to the expected
    // background of a 'clean' image)
    builder.addClampValues(10);
    builder.addClampValues(1000); // I
    builder.addClampValues(1); // X
    builder.addClampValues(1); // Y
    // TODO: determine what this should be
    builder.addClampValues(10); // Z
    builder.addClampValues(3); // Sx
    builder.addClampValues(3); // Sy
    builder.addClampValues(Math.PI); // A

    builder.setLineSearchMethod(LineSearchMethod.PARTIAL_IGNORE);

    defaultFitSolverSettings = builder.build();
  }

  /** The default FilterSettings. */
  public static final FilterSettings defaultFilterSettings;

  static {
    final FilterSettings.Builder builder = FilterSettings.newBuilder();
    builder.setShiftFactor(1);
    builder.setSignalStrength(5);
    builder.setMinPhotons(30);
    builder.setPrecisionThreshold(40);
    builder.setPrecisionMethod(PrecisionMethod.MORTENSEN);
    builder.setMinWidthFactor(0.5);
    builder.setMaxWidthFactor(2);
    builder.setDisableSimpleFilter(false);
    builder.setSmartFilter(true);
    builder.setSmartFilterString("");
    defaultFilterSettings = builder.build();
  }

  /** The default FitSettings. */
  public static final FitSettings defaultFitSettings;

  static {
    final FitSettings.Builder builder = FitSettings.newBuilder();
    builder.setFitSolverSettings(defaultFitSolverSettings);
    builder.setFilterSettings(defaultFilterSettings);
    defaultFitSettings = builder.build();
  }

  /** The default FitEngineSettings. */
  public static final FitEngineSettings defaultFitEngineSettings;

  static {
    // Analysis* shows the best area-under-precision-recall curve (AUC) using a mean filter or
    // a Gaussian filter with ~1.2 SD smoothing. The Gaussian filter is more robust to width
    // mismatch but
    // the mean filter will be faster as it uses a smaller block size. The Gaussian filter has
    // higher
    // recall but lower precision as it identifies more spots due to the shape of the smoothing
    // filter.
    // The overall AUC is very similar.
    //
    // Note: Setting the parameter at a higher level allows the filter to work on out-of-focus spots
    // which
    // will have a wider PSF.
    //
    // *Analysis was performed on simulated data using a Image PSF with spots of 20-100 photons at a
    // depth of up to 1380nm (the PSF limit).

    final FitEngineSettings.Builder builder = FitEngineSettings.newBuilder();
    builder.setFitSettings(defaultFitSettings);

    builder.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES);

    final RelativeParameter.Builder rp = RelativeParameter.newBuilder();

    final DataFilterSettings.Builder dfs = builder.getDataFilterSettingsBuilder();
    dfs.setDataFilterType(DataFilterType.SINGLE);
    final DataFilter.Builder dfb = dfs.addDataFiltersBuilder();
    dfb.setDataFilterMethod(DataFilterMethod.MEAN);
    rp.setAbsolute(false);
    rp.setValue(1.2);
    dfb.addParameters(rp.build());

    rp.setAbsolute(false);
    rp.setValue(1);
    builder.setSearch(rp.build());
    builder.setBorder(rp.build());
    rp.setValue(3);
    builder.setFitting(rp.build());

    builder.setIncludeNeighbours(true);
    builder.setNeighbourHeightThreshold(0.3);
    builder.setResidualsThreshold(1);
    rp.setValue(0.5);
    rp.setAbsolute(true);
    builder.setDuplicateDistance(rp.build());

    // The new pass rate parameter should be more adaptable to different image sizes.
    // XXX Revisit this when more is known about how to set the pass rate.
    builder.setFailuresLimit(3);
    builder.setPassRate(0.5);

    defaultFitEngineSettings = builder.build();
  }

  /** No public constructor. */
  private FitProtosHelper() {}

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(FitSolver value) {
    switch (value) {
      case BACKTRACKING_FAST_MLE:
        return "Backtracking Fast MLE";
      case FAST_MLE:
        return "Fast MLE";
      case LVM_LSE:
        return "LVM LSE";
      case LVM_MLE:
        return "LVM MLE";
      case LVM_WLSE:
        return "LVM WLSE";
      case MLE:
        return "MLE";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownNameMessage(value));
    }
  }

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(SearchMethod value) {
    switch (value) {
      case BFGS:
        return "BFGS";
      case BOBYQA:
        return "BOBYQA";
      case CMAES:
        return "CMAES";
      case CONJUGATE_GRADIENT_FR:
        return "Conjugate Gradient FR";
      case CONJUGATE_GRADIENT_PR:
        return "Conjugate Gradient PR";
      case POWELL:
        return "Powell";
      case POWELL_ADAPTER:
        return "Powell (adpator)";
      case POWELL_BOUNDED:
        return "Powell (bounded)";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownNameMessage(value));
    }
  }

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(DataFilterType value) {
    switch (value) {
      case DIFFERENCE:
        return "Difference";
      case JURY:
        return "Jury";
      case SINGLE:
        return "Single";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownNameMessage(value));
    }
  }

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(DataFilterMethod value) {
    switch (value) {
      case BLOCK_MEAN:
        return "Block Mean";
      case CIRCULAR_MEAN:
        return "Circular Mean";
      case GAUSSIAN:
        return "Gaussian";
      case MEAN:
        return "Mean";
      case MEDIAN:
        return "Median";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownNameMessage(value));
    }
  }

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(NoiseEstimatorMethod value) {
    switch (value) {
      case ALL_PIXELS:
        return "All pixels";
      case LOWEST_PIXELS:
        return "Lowest pixels";
      case QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES:
        return "Quick residuals least mean of squares";
      case QUICK_RESIDUALS_LEAST_MEDIAN_OF_SQUARES:
        return "Quick residuals least median of squares";
      case QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES:
        return "Quick residuals least trimmed of squares";
      case RESIDUALS_LEAST_MEAN_OF_SQUARES:
        return "Residuals least mean of squares";
      case RESIDUALS_LEAST_MEDIAN_OF_SQUARES:
        return "Residuals least median of squares";
      case RESIDUALS_LEAST_TRIMMED_OF_SQUARES:
        return "Residuals least trimmed of squares";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownNameMessage(value));
    }
  }

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(PrecisionMethod value) {
    switch (value) {
      case PRECISION_METHOD_NA:
        return "NA";
      case MORTENSEN:
        return "Mortensen";
      case MORTENSEN_LOCAL_BACKGROUND:
        return "Mortensen (local background)";
      case POISSON_CRLB:
        return "Poisson CRLB";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownNameMessage(value));
    }
  }

  /**
   * Convert noise estimator method to the
   * {@link uk.ac.sussex.gdsc.core.utils.NoiseEstimator.Method }.
   *
   * @param method the method
   * @return the method
   */
  public static Method convertNoiseEstimatorMethod(NoiseEstimatorMethod method) {
    switch (method) {
      case ALL_PIXELS:
        return Method.ALL_PIXELS;
      case LOWEST_PIXELS:
        return Method.LOWEST_PIXELS;
      case QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES:
        return Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES;
      case QUICK_RESIDUALS_LEAST_MEDIAN_OF_SQUARES:
        return Method.QUICK_RESIDUALS_LEAST_MEDIAN_OF_SQUARES;
      case QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES:
        return Method.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES;
      case RESIDUALS_LEAST_MEAN_OF_SQUARES:
        return Method.RESIDUALS_LEAST_MEAN_OF_SQUARES;
      case RESIDUALS_LEAST_MEDIAN_OF_SQUARES:
        return Method.RESIDUALS_LEAST_MEDIAN_OF_SQUARES;
      case RESIDUALS_LEAST_TRIMMED_OF_SQUARES:
        return Method.RESIDUALS_LEAST_TRIMMED_OF_SQUARES;
      case UNRECOGNIZED:
        break;
      default:
        break;
    }
    throw new IllegalArgumentException(ProtosHelperUtils.unknownMethodMessage(method));
  }

  /**
   * Convert search method to the
   * {@link uk.ac.sussex.gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter.SearchMethod}.
   *
   * @param searchMethod the search method
   * @return the search method
   */
  public static MaximumLikelihoodFitter.SearchMethod
      convertSearchMethod(SearchMethod searchMethod) {
    switch (searchMethod) {
      case BOBYQA:
        return MaximumLikelihoodFitter.SearchMethod.BOBYQA;
      case CMAES:
        return MaximumLikelihoodFitter.SearchMethod.CMAES;
      case CONJUGATE_GRADIENT_FR:
        return MaximumLikelihoodFitter.SearchMethod.CONJUGATE_GRADIENT_FR;
      case CONJUGATE_GRADIENT_PR:
        return MaximumLikelihoodFitter.SearchMethod.CONJUGATE_GRADIENT_PR;
      case POWELL:
        return MaximumLikelihoodFitter.SearchMethod.POWELL;
      case POWELL_ADAPTER:
        return MaximumLikelihoodFitter.SearchMethod.POWELL_ADAPTER;
      case POWELL_BOUNDED:
        return MaximumLikelihoodFitter.SearchMethod.POWELL_BOUNDED;
      // BFGS is unsupported
      case BFGS:
      case UNRECOGNIZED:
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownMethodMessage(searchMethod));
    }
  }

  /**
   * Convert line search method to the
   * {@link uk.ac.sussex.gdsc.smlm.fitting.nonlinear.FastMleSteppingFunctionSolver.LineSearchMethod}.
   *
   * @param lineSearchMethod the line search method
   * @return the line search method
   */
  public static FastMleSteppingFunctionSolver.LineSearchMethod
      convertLineSearchMethod(LineSearchMethod lineSearchMethod) {
    switch (lineSearchMethod) {
      case IGNORE:
        return FastMleSteppingFunctionSolver.LineSearchMethod.IGNORE;
      case NONE:
        return FastMleSteppingFunctionSolver.LineSearchMethod.NONE;
      case PARTIAL_IGNORE:
        return FastMleSteppingFunctionSolver.LineSearchMethod.PARTIAL_IGNORE;
      case UNRECOGNIZED:
      default:
        throw new IllegalArgumentException(
            ProtosHelperUtils.unknownMethodMessage(lineSearchMethod));
    }
  }
}
