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
 * Class to encapsulate counting
 */
public class Counter
{
	private int count;

	/**
	 * Instantiates a new counter.
	 */
	public Counter()
	{
		this.count = 0;
	}

	/**
	 * Instantiates a new counter.
	 *
	 * @param count
	 *            the count
	 */
	public Counter(int count)
	{
		this.count = count;
	}

	/**
	 * Reset the count
	 */
	public void reset()
	{
		count = 0;
	}

	/**
	 * Increment the count.
	 */
	public void increment()
	{
		count++;
	}

	/**
	 * Increment the count by the value.
	 *
	 * @param value
	 *            the value
	 */
	public void increment(int value)
	{
		count += value;
	}

	/**
	 * Increment the count.
	 *
	 * @return the new count
	 */
	public int incrementAndGet()
	{
		return ++count;
	}

	/**
	 * Increment the count by the value.
	 *
	 * @param value
	 *            the value
	 * @return the new count
	 */
	public int incrementAndGet(int value)
	{
		return count += value;
	}

	/**
	 * Increment the count.
	 *
	 * @return the old count
	 */
	public int getAndIncrement()
	{
		return count++;
	}

	/**
	 * Increment the count by the value.
	 *
	 * @param value
	 *            the value
	 * @return the old count
	 */
	public int getAndIncrement(int value)
	{
		final int old = count;
		count += value;
		return old;
	}

	/**
	 * Decrement the count.
	 */
	public void decrement()
	{
		count--;
	}

	/**
	 * Decrement the count by the value.
	 *
	 * @param value
	 *            the value
	 */
	public void decrement(int value)
	{
		count -= value;
	}

	/**
	 * Decrement the count.
	 *
	 * @return the new count
	 */
	public int decrementAndGet()
	{
		return --count;
	}

	/**
	 * Decrement the count by the value.
	 *
	 * @param value
	 *            the value
	 * @return the new count
	 */
	public int decrementAndGet(int value)
	{
		return count -= value;
	}

	/**
	 * Decrement the count.
	 *
	 * @return the old count
	 */
	public int getAndDecrement()
	{
		return count--;
	}

	/**
	 * Decrement the count by the value.
	 *
	 * @param value
	 *            the value
	 * @return the old count
	 */
	public int getAndDecrement(int value)
	{
		final int old = count;
		count -= value;
		return old;
	}

	/**
	 * Gets the count.
	 *
	 * @return the count
	 */
	public int getCount()
	{
		return count;
	}
}
