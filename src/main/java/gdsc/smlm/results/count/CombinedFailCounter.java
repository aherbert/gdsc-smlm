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
package gdsc.smlm.results.count;

/**
 * Combine the result of two fail counters
 */
public abstract class CombinedFailCounter extends BaseFailCounter
{
	/** The first fail counter . */
	protected final FailCounter c1;
	/** The second fail counter . */
	protected final FailCounter c2;

	/**
	 * Instantiates a new combined fail counter.
	 *
	 * @param c1
	 *            the first counter
	 * @param c2
	 *            the second counter
	 */
	public CombinedFailCounter(FailCounter c1, FailCounter c2)
	{
		if (c1 == null || c2 == null)
			throw new NullPointerException();
		this.c1 = c1;
		this.c2 = c2;
	}

	@Override
	protected String generateDescription()
	{
		StringBuilder sb = new StringBuilder();
		add(sb, c1);
		sb.append(" ").append(getOperator()).append(" ");
		add(sb, c2);
		return sb.toString();
	}

	private static void add(StringBuilder sb, FailCounter c)
	{
		if (c instanceof CombinedFailCounter)
		{
			sb.append("(");
			sb.append(c.getDescription());
			sb.append(")");
		}
		else
		{
			sb.append(c.getDescription());
		}
	}

	/**
	 * Get the string representation of the operator used to combine the two fail counters. This is used in the filter
	 * name.
	 *
	 * @return The operator
	 */
	protected abstract String getOperator();

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#pass()
	 */
	@Override
	public void pass()
	{
		c1.pass();
		c2.pass();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#pass(int)
	 */
	@Override
	public void pass(int n)
	{
		c1.pass(n);
		c2.pass(n);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#fail()
	 */
	@Override
	public void fail()
	{
		c1.fail();
		c2.fail();
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#fail(int)
	 */
	@Override
	public void fail(int n)
	{
		c1.fail(n);
		c2.fail(n);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#reset()
	 */
	@Override
	public void reset()
	{
		c1.reset();
		c2.reset();
	}
}
