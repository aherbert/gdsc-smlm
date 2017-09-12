package gdsc.smlm.results.predicates;

import gdsc.smlm.results.PeakResult;

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
 * Test a result using the id
 */
public class IdPeakResultPredicate implements PeakResultPredicate
{
	/** The id. */
	private final int id;

	/**
	 * Instantiates a new id peak result predicate.
	 *
	 * @param id
	 *            the id
	 */
	public IdPeakResultPredicate(int id)
	{
		this.id = id;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.PeakResultPredicate#test(gdsc.smlm.results.PeakResult)
	 */
	public boolean test(PeakResult t)
	{
		return t.getId() == id;
	}
}