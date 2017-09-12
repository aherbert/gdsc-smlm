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
 * Test a result using the frame
 */
public class FramePeakResultPredicate implements PeakResultPredicate
{
	/** The frame. */
	private final int frame;

	/**
	 * Instantiates a new frame peak result predicate.
	 *
	 * @param frame
	 *            the frame
	 */
	public FramePeakResultPredicate(int frame)
	{
		this.frame = frame;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.PeakResultPredicate#test(gdsc.smlm.results.PeakResult)
	 */
	public boolean test(PeakResult t)
	{
		return t.getFrame() == frame;
	}
}