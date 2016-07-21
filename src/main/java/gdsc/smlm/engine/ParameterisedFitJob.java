package gdsc.smlm.engine;

import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.results.PeakResult;

import java.awt.Rectangle;
import java.util.List;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specifies a job for peak fitting.
 */
public class ParameterisedFitJob extends FitJob
{
	private FitParameters parameters;
	private List<PeakResult> peakResults;
	private int[] indices = new int[0];
	private FitResult[] fitResults = new FitResult[0];
	private FitResult[] fitResultsWithNeighbours = null;

	/**
	 * Constructor with data. Exceptions are thrown if invalid bounds or data are passed
	 * 
	 * @param id
	 * @param parameters
	 * @param slice
	 * @param data
	 * @param bounds
	 */
	public ParameterisedFitJob(int id, FitParameters parameters, int slice, float[] data, Rectangle bounds)
	{
		super(id, slice, data, bounds);
		this.parameters = parameters;
	}

	/**
	 * Constructor with data. Exceptions are thrown if invalid bounds or data are passed
	 * 
	 * @param parameters
	 * @param slice
	 * @param data
	 * @param bounds
	 */
	public ParameterisedFitJob(FitParameters parameters, int slice, float[] data, Rectangle bounds)
	{
		super(slice, slice, data, bounds);
		this.parameters = parameters;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.engine.FitJob#getFitParameters()
	 */
	public FitParameters getFitParameters()
	{
		return parameters;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.engine.FitJob#setResults(java.util.List)
	 */
	public void setResults(List<PeakResult> results)
	{
		this.peakResults = results;
	}

	/**
	 * @return The results
	 */
	public List<PeakResult> getResults()
	{
		return peakResults;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.engine.FitJob#setIndices(int[])
	 */
	@Override
	public void setIndices(int[] indices)
	{
		this.indices = indices;
		fitResults = new FitResult[indices.length];
	}

	@Override
	public void setFitResult(int n, FitResult fitResult)
	{
		if (n < fitResults.length)
			fitResults[n] = fitResult;
	}

	@Override
	public void setFitResultWithNeighbours(int n, FitResult fitResult)
	{
		// Dynamically initialise as we may not fit neighbours
		if (fitResultsWithNeighbours == null)
			fitResultsWithNeighbours = new FitResult[indices.length]; 
		if (n < fitResultsWithNeighbours.length)
			fitResultsWithNeighbours[n] = fitResult;
	}
	
	/**
	 * @return The indices of the data that were fitted
	 */
	public int[] getIndices()
	{
		return indices;
	}

	/**
	 * The fit result of the specified index in the array of fitted indices
	 * 
	 * @param n
	 * @return
	 */
	public FitResult getFitResult(int n)
	{
		return fitResults[n];
	}

	/**
	 * The fit result of the specified index in the array of fitted indices. This is only set if the fit was attempted
	 * with neighbours.
	 * 
	 * @param n
	 * @return
	 */
	public FitResult getFitResultWithNeighbours(int n)
	{
		if (fitResultsWithNeighbours != null && n < fitResultsWithNeighbours.length)
			return fitResultsWithNeighbours[n];
		return null;
	}
}
