/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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
package gdsc.smlm.engine;

import gdsc.smlm.fitting.FitResult;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.filter.MultiPathFitResult;

import java.awt.Rectangle;
import java.util.List;


/**
 * Specifies a job for peak fitting.
 */
public class FitJob
{
	public enum Status
	{
		//@formatter:off
		PENDING{ public String getName() { return "Pending"; }}, 
		IN_PROGRESS{ public String getName() { return "In-progress"; }}, 
		FINISHED{ public String getName() { return "Finished"; }};
		//@formatter:on

		@Override
		public String toString()
		{
			return getName();
		}

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();
	}

	private int id = 0;
	int slice;
	float[] data;
	Rectangle bounds;
	Status status = Status.PENDING;

	/**
	 * Constructor with data. Exceptions are thrown if invalid bounds or data are passed
	 * 
	 * @param id
	 * @param slice
	 * @param data
	 * @param bounds
	 */
	public FitJob(int id, int slice, float[] data, Rectangle bounds)
	{
		if (data == null)
			throw new NullPointerException("Data must not be null");
		if (bounds == null)
			throw new NullPointerException("Bounds must not be null");
		if (bounds.width < 1 || bounds.height < 1)
			throw new IllegalArgumentException("Bounds width and height must be positive");
		if (data.length < bounds.width * bounds.height)
			throw new IllegalArgumentException("Data must be at least equal to the width multiplied by the height");

		this.id = id;
		this.slice = slice;
		this.data = data;
		this.bounds = bounds;
	}

	/**
	 * Constructor with data. Exceptions are thrown if invalid bounds or data are passed
	 * 
	 * @param slice
	 * @param data
	 * @param bounds
	 */
	public FitJob(int slice, float[] data, Rectangle bounds)
	{
		this(slice, slice, data, bounds);
	}

	/**
	 * Constructor with data. Exceptions are thrown if invalid bounds or data are passed
	 * 
	 * @param id
	 * @param slice
	 * @param data
	 * @param width
	 * @param height
	 */
	public FitJob(int id, int slice, float[] data, int width, int height)
	{
		this(id, slice, data, new Rectangle(width, height));
	}

	/**
	 * Constructor with data. Exceptions are thrown if invalid bounds or data are passed
	 * 
	 * @param slice
	 * @param data
	 * @param width
	 * @param height
	 */
	public FitJob(int slice, float[] data, int width, int height)
	{
		this(slice, slice, data, new Rectangle(width, height));
	}

	/**
	 * Empty constructor with no data
	 */
	public FitJob()
	{
		this.id = 0;
		this.slice = 0;
		this.data = null;
		this.bounds = null;
	}

	/**
	 * @return the id
	 */
	public int getId()
	{
		return id;
	}

	/**
	 * @return the slice
	 */
	public int getSlice()
	{
		return slice;
	}

	/**
	 * @return the data
	 */
	public float[] getData()
	{
		return data;
	}

	/**
	 * @return the bounds
	 */
	public Rectangle getBounds()
	{
		return bounds;
	}

	/**
	 * Called to indicate that processing of the job is in progress
	 */
	public void start()
	{
		status = Status.IN_PROGRESS;
	}

	/**
	 * Called to indicate that processing of the job has finished
	 */
	public void finished()
	{
		status = Status.FINISHED;
	}

	/**
	 * @return The additional parameters for the job (or null)
	 */
	public FitParameters getFitParameters()
	{
		return null;
	}

	/**
	 * @return The job status
	 */
	public Status getStatus()
	{
		return status;
	}

	/**
	 * Set the indices of the the data that were fitted
	 * <p>
	 * This method is not implemented within this class but can be overridden if necessary to allow access to the
	 * results from the {@link FitWorker}.
	 * 
	 * @param indices
	 */
	public void setIndices(int[] indices)
	{

	}

	/**
	 * Set the fitting result for the specified fitted index
	 * <p>
	 * This method is not implemented within this class but can be overridden if necessary to allow access to the
	 * results from the {@link FitWorker}.
	 * 
	 * @param n
	 * @param fitResult
	 */
	public void setFitResult(int n, FitResult fitResult)
	{

	}

	/**
	 * Used to store the results of the fitting job.
	 * <p>
	 * This method is not implemented within this class but can be overridden if necessary to allow access to the
	 * results from the {@link FitWorker}.
	 * 
	 * @param results
	 */
	public void setResults(List<PeakResult> results)
	{
	}

	/**
	 * Set the multi-path fitting result for the specified fitted index
	 * <p>
	 * This method is not implemented within this class but can be overridden if necessary to allow access to the
	 * benchmarking results from the {@link FitWorker}.
	 * 
	 * @param n
	 * @param fitResult
	 */
	public void setMultiPathFitResult(int n, MultiPathFitResult fitResult)
	{

	}
}
