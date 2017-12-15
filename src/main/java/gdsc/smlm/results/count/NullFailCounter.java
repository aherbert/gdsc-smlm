package gdsc.smlm.results.count;

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
 * A fail counter that always returns OK. Failures are ignored.
 */
public class NullFailCounter implements FailCounter
{
	/** An instance */
	public static final NullFailCounter INSTANCE = new NullFailCounter();

	public String getDescription()
	{
		return "ignoreAllFailures";
	}
	
	public void pass()
	{
	}

	public void pass(int n)
	{
	}

	public void fail()
	{
	}

	public void fail(int n)
	{
	}

	public boolean isOK()
	{
		return true;
	}

	public FailCounter newCounter()
	{
		return this; // This doesn't matter
	}

	public void reset()
	{
	}
}
