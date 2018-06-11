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
package gdsc.smlm.filters;

/**
 * Contains common functionality for weighted filters.
 */
public abstract class BaseWeightedFilter extends BaseFilter
{
	/** The weights. */
	protected float[] weights;

	protected int weightWidth, weightHeight;

	/**
	 * Sets the weights of the data. This should be called before filtering data samples.
	 *
	 * @param weights
	 *            the weights of the data (can be null)
	 * @param width
	 *            The width of the data
	 * @param height
	 *            The height of the data
	 */
	public void setWeights(final float[] weights, final int width, final int height)
	{
		this.weights = weights;
		this.weightWidth = width;
		this.weightHeight = height;
		newWeights();
	}

	/**
	 * Checks for weights.
	 *
	 * @return true, if successful
	 */
	public boolean hasWeights()
	{
		return weights != null;
	}

	/**
	 * Signal that new weight parameters have been set. Sub-classes can re-initialise for the new weights.
	 */
	protected abstract void newWeights();

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public BaseWeightedFilter clone()
	{
		return (BaseWeightedFilter) super.clone();
	}
}
