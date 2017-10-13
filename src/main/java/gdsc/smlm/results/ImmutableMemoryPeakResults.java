package gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.Collection;

import gdsc.core.data.DataException;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.results.predicates.PeakResultPredicate;

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
 * Wraps peak results in memory and prevents modification of the results size.
 * <p>
 * Any method that modifies the size of the results set will throw a data exception.
 */
public class ImmutableMemoryPeakResults extends MemoryPeakResults
{
	private boolean built = false;

	/**
	 * Instantiates a new immutable memory peak results with the original results store.
	 *
	 * @param results
	 *            the results
	 */
	public ImmutableMemoryPeakResults(MemoryPeakResults results)
	{
		this(results, false);
	}

	/**
	 * Instantiates a new immutable memory peak results with an optional copy of the results store.
	 *
	 * @param results
	 *            the results
	 * @param copy
	 *            Set to true to copy the original results store
	 */
	public ImmutableMemoryPeakResults(MemoryPeakResults results, boolean copy)
	{
		super((copy) ? results.results.copy() : results.results);
		copySettings(results);
		built = true;
	}

	@Override
	public void setSource(ImageSource source)
	{
		if (built)
			throw new DataException("This results set is immutable");
		super.setSource(source);
	}

	@Override
	public void setBounds(Rectangle bounds)
	{
		if (built)
			throw new DataException("This results set is immutable");
		super.setBounds(bounds);
	}

	@Override
	public void setCalibration(Calibration calibration)
	{
		if (built)
			throw new DataException("This results set is immutable");
		super.setCalibration(calibration);
	}

	@Override
	public void setPSF(PSF psf)
	{
		if (built)
			throw new DataException("This results set is immutable");
		super.setPSF(psf);
	}

	@Override
	public void setConfiguration(String configuration)
	{
		if (built)
			throw new DataException("This results set is immutable");
		super.setConfiguration(configuration);
	}

	@Override
	public void setName(String name)
	{
		if (built)
			throw new DataException("This results set is immutable");
		super.setName(name);
	}

	@Override
	public void add(PeakResult result)
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public void addAll(Collection<PeakResult> results)
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public void addAll(PeakResult[] results)
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public void addAll(PeakResultStore results)
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public void add(MemoryPeakResults results)
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public void removeNullResults()
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public boolean removeIf(PeakResultPredicate filter)
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public void begin()
	{
		throw new DataException("This results set is immutable");
	}
	
	@Override
	public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise, float[] params,
			float[] paramsStdDev)
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public void end()
	{
		throw new DataException("This results set is immutable");
	}

	@Override
	public boolean isActive()
	{
		throw new DataException("This results set is immutable");
	}
}
