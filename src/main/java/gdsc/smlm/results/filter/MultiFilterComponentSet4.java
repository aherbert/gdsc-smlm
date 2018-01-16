package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains a set of components of the multi filter.
 */
public class MultiFilterComponentSet4 extends MultiFilterComponentSet
{
	private MultiFilterComponent component0;
	private MultiFilterComponent component1;
	private MultiFilterComponent component2;
	private MultiFilterComponent component3;

	public MultiFilterComponentSet4(MultiFilterComponent[] components)
	{
		this.component0 = components[0];
		this.component1 = components[1];
		this.component2 = components[2];
		this.component3 = components[3];
	}

	@Override
	public int getValidationFlags()
	{
		return component0.getType() | component1.getType() | component2.getType() | component3.getType();
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		//@formatter:off
		if (component0.fail(peak)) return component0.getType();
		if (component1.fail(peak)) return component1.getType();
		if (component2.fail(peak)) return component2.getType();
		if (component3.fail(peak)) return component3.getType();
		//@formatter:on
		return 0;
	}

	@Override
	void replace0(MultiFilterComponent c)
	{
		component0 = c;
	}
}