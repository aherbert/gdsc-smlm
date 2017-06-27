package gdsc.smlm.data.config;

import gdsc.core.utils.Maths;
import gdsc.smlm.data.config.FitConfig.DataFilter;
import gdsc.smlm.data.config.FitConfig.DataFilterType;
import gdsc.smlm.data.config.FitConfig.FilterSettings;
import gdsc.smlm.data.config.FitConfig.FitDistance;
import gdsc.smlm.data.config.FitConfig.FitEngineSettings;
import gdsc.smlm.data.config.FitConfig.FitSettings;
import gdsc.smlm.data.config.FitConfig.FitSolver;
import gdsc.smlm.data.config.FitConfig.FitSolverSettings;
import gdsc.smlm.data.config.FitConfig.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitConfig.SearchMethod;

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
		builder.setMinIterations(0);
		builder.setMaxIterations(20);
		builder.setRelativeThreshold(1e-6);
		builder.setAbsoluteThreshold(1e-16);
		builder.setParameterRelativeThreshold(1e-3);
		builder.setParameterAbsoluteThreshold(1e-6);
		
		builder.setLambda(10);
		
		builder.setSearchMethod(SearchMethod.POWELL_BOUNDED);
		builder.setGradientLineMinimisation(false);
		builder.setModelCamera(false);
		builder.setMaxFunctionEvaluations(100);
		
		builder.setUseClamping(false);
		builder.setUseDynamicClamping(false);
		// Add defaults for a two-axis and theta Gaussian 2D function.
		// The units are photons and pixels.
		builder.addClampValue(100); // B
		builder.addClampValue(1000); // I
		builder.addClampValue(1); // X
		builder.addClampValue(1); // Y
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
		builder.setPrecisionThreshold(Maths.pow2(40));
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
		// No calibration by default 
		//builder.setCalibration(CalibrationConfigHelper.defaultCalibration);
		builder.setPsf(PSFConfigHelper.defaultOneAxisGaussian2DPSF);
		builder.setFitSolverSettings(defaultFitSolverSettings);
		builder.setFilterSettings(defaultFilterSettings);
		defaultFitSettings = builder.build();
	}

	/** The default FitEngineSettings */
	public static final FitEngineSettings defaultFitEngineSettings;
	static
	{
		FitEngineSettings.Builder builder = FitEngineSettings.newBuilder();
		builder.setFitSettings(defaultFitSettings);
		
		builder.setNoiseMethod(NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES);

		FitDistance.Builder fd = FitDistance.newBuilder();
		DataFilter.Builder df = DataFilter.newBuilder();
		df.setDataFilterType(DataFilterType.SINGLE);
		fd.setAbsolute(false);
		fd.setValue(1.2);
		df.addDistance(fd.build());
		
		fd.setAbsolute(false);
		fd.setValue(1);
		builder.setSearch(fd.build());
		builder.setBorder(fd.build());
		fd.setValue(3);
		builder.setFitting(fd.build());
		
		builder.setIncludeNeighbours(true);
		builder.setNeighbourHeightThreshold(0.3);
		builder.setResidualsThreshold(1);
		fd.setValue(0.5);
		fd.setAbsolute(true);
		builder.setDuplicateDistance(fd.build());

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

}
