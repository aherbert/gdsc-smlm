package gdsc.smlm.data.config;

import gdsc.core.utils.NoiseEstimator.Method;
import gdsc.smlm.data.config.FitConfig.DataFilter;
import gdsc.smlm.data.config.FitConfig.DataFilterMethod;
import gdsc.smlm.data.config.FitConfig.DataFilterSettings;
import gdsc.smlm.data.config.FitConfig.DataFilterType;
import gdsc.smlm.data.config.FitConfig.FilterSettings;
import gdsc.smlm.data.config.FitConfig.FitEngineSettings;
import gdsc.smlm.data.config.FitConfig.FitSettings;
import gdsc.smlm.data.config.FitConfig.FitSolver;
import gdsc.smlm.data.config.FitConfig.FitSolverSettings;
import gdsc.smlm.data.config.FitConfig.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitConfig.RelativeParameter;
import gdsc.smlm.data.config.FitConfig.SearchMethod;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains helper functions for the FitConfig class.
 */
public class FitConfigHelper
{
	/** The default FitSolverSettings */
	public static final FitSolverSettings defaultFitSolverSettings;
	static
	{
		FitSolverSettings.Builder builder = FitSolverSettings.newBuilder();
		builder.setFixedPsf(false);
		builder.setDisableBackgroundFitting(false);
		builder.setDisableSignalFitting(false);
		
		builder.setFitSolver(FitSolver.LVM_LSE);
		builder.setFixedIterations(false);
		builder.setMaxIterations(20);
		builder.setRelativeThreshold(1e-6);
		builder.setAbsoluteThreshold(1e-16);
		builder.setParameterRelativeThreshold(1e-3);
		builder.setParameterAbsoluteThreshold(1e-6);
		
		builder.setLambda(10);
		
		builder.setSearchMethod(SearchMethod.POWELL_BOUNDED);
		builder.setGradientLineMinimisation(false);
		builder.setModelCamera(false);
		builder.setMaxFunctionEvaluations(2000);
		
		builder.setUseClamping(false);
		builder.setUseDynamicClamping(false);
		// Add defaults for a two-axis and theta Gaussian 2D function.
		// The units are photons and pixels.
		builder.addClampValue(10); // B (3D DAOSTORM uses 100 which is high compared to the expected background of a 'clean' image)
		builder.addClampValue(1000); // I
		builder.addClampValue(1); // X
		builder.addClampValue(1); // Y
		// TODO: determine what this should be
		builder.addClampValue(10); // Z
		builder.addClampValue(3); // Sx
		builder.addClampValue(3); // Sy
		builder.addClampValue(Math.PI); // A
		
		defaultFitSolverSettings = builder.build();
	}

	/** The default FilterSettings */
	public static final FilterSettings defaultFilterSettings;
	static
	{
		FilterSettings.Builder builder = FilterSettings.newBuilder();
		builder.setShiftFactor(1);
		builder.setSignalStrength(30);
		builder.setMinPhotons(30);
		builder.setPrecisionThreshold(40);
		builder.setPrecisionUsingBackground(false);
		builder.setMinWidthFactor(0.5);
		builder.setMaxWidthFactor(2);
		builder.setDisableSimpleFilter(false);
		builder.setSmartFilter(true);
		builder.setSmartFilterString("");
		defaultFilterSettings = builder.build();
	}

	/** The default FitSettings */
	public static final FitSettings defaultFitSettings;
	static
	{
		FitSettings.Builder builder = FitSettings.newBuilder();
		builder.setFitSolverSettings(defaultFitSolverSettings);
		builder.setFilterSettings(defaultFilterSettings);
		defaultFitSettings = builder.build();
	}

	/** The default FitEngineSettings */
	public static final FitEngineSettings defaultFitEngineSettings;
	static
	{
		// Analysis* shows the best area-under-precision-recall curve (AUC) using a mean filter or
		// a Gaussian filter with ~1.2 SD smoothing. The Gaussian filter is more robust to width mismatch but
		// the mean filter will be faster as it uses a smaller block size. The Gaussian filter has higher 
		// recall but lower precision as it identifies more spots due to the shape of the smoothing filter.
		// The overall AUC is very similar.
		//
		// Note: Setting the parameter at a higher level allows the filter to work on out-of-focus spots which
		// will have a wider PSF.
		//
		// *Analysis was performed on simulated data using a Image PSF with spots of 20-100 photons at a 
		// depth of up to 1380nm (the PSF limit).

		FitEngineSettings.Builder builder = FitEngineSettings.newBuilder();
		builder.setFitSettings(defaultFitSettings);
		
		builder.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES);

		
		RelativeParameter.Builder rp =  RelativeParameter.newBuilder();
		
		DataFilterSettings.Builder dfs = builder.getDataFilterSettingsBuilder();
		dfs.setDataFilterType(DataFilterType.SINGLE);
		DataFilter.Builder dfb =  dfs.addDataFilterBuilder();
		dfb.setDataFilterMethod(DataFilterMethod.MEAN);
		rp.setAbsolute(false);
		rp.setValue(1.2);
		dfb.addParameter(rp.build());
		
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

		builder.setFailuresLimit(3);
		
		defaultFitEngineSettings = builder.build();
	}

	public static String getName(FitSolver value)
	{
		switch (value)
		{
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
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}

	public static String getName(SearchMethod value)
	{
		switch (value)
		{
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
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}

	public static String getName(DataFilterType value)
	{
		switch (value)
		{
			case DIFFERENCE:
				return "Difference";
			case JURY:
				return "Jury";
			case SINGLE:
				return "Single";
			case UNRECOGNIZED:
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}
	public static String getName(DataFilterMethod value)
	{
		switch (value)
		{
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
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}

	public static String getName(NoiseEstimatorMethod value)
	{
		switch (value)
		{
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
				return "Unknown";
			default:
				throw new IllegalStateException("Unknown name: " + value);
		}
	}

	public static Method convertNoiseEstimatorMethod(NoiseEstimatorMethod method)
	{
		switch (method)
		{
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
		throw new IllegalStateException("Unknown method: " + method);
	}
	
	public static gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter.SearchMethod convertSearchMethod(SearchMethod searchMethod)
	{
		switch (searchMethod)
		{
			case BFGS:
				return MaximumLikelihoodFitter.SearchMethod.BFGS;
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
			case UNRECOGNIZED:
			default:
				throw new IllegalStateException("Unknown method: " + searchMethod);
		}
	}
}
