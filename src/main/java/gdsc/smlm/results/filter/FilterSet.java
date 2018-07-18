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
package gdsc.smlm.results.filter;

import java.util.Collections;
import java.util.List;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Specify a named set of filters
 */
public class FilterSet
{
	@XStreamAsAttribute
	private String name;
	private final List<Filter> filters;

	@XStreamOmitField
	private Filter weakest;
	@XStreamOmitField
	private boolean initialisedWeakest;
	@XStreamOmitField
	private int allSameType = 0;

	/**
	 * Instantiates a new filter set.
	 *
	 * @param name
	 *            the name
	 * @param filters
	 *            the filters
	 */
	public FilterSet(String name, List<Filter> filters)
	{
		this.name = name;
		this.filters = filters;
	}

	/**
	 * Instantiates a new filter set.
	 *
	 * @param filters
	 *            the filters
	 */
	public FilterSet(List<Filter> filters)
	{
		this.filters = filters;
	}

	/**
	 * Get the size of the set
	 *
	 * @return the size
	 */
	public int size()
	{
		return filters.size();
	}

	/**
	 * The name of the set. If empty return the name of the first filter in the set.
	 *
	 * @return the name
	 */
	public String getName()
	{
		if (name == null || name.length() == 0)
			if (filters != null && !filters.isEmpty())
				name = filters.get(0).getName();
		return name;
	}

	/**
	 * Return the name of the value that will be returned from the first Filter getNumericalValueName() method
	 *
	 * @return the valueName
	 */
	public String getValueName()
	{
		if (filters != null && !filters.isEmpty())
			return filters.get(0).getNumericalValueName();
		return "";
	}

	/**
	 * @return the filters
	 */
	public List<Filter> getFilters()
	{
		return filters;
	}

	/**
	 * @return An XML representation of this object
	 */
	public String toXML()
	{
		return XStreamWrapper.toXML(this);
	}

	/**
	 * Create the filter set from the XML representation.
	 *
	 * @param xml
	 *            the xml
	 * @return the filter set
	 */
	public static FilterSet fromXML(String xml)
	{
		try
		{
			return (FilterSet) XStreamWrapper.fromXML(xml);
		}
		catch (final ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		return null;
	}

	/**
	 * Sort the filters
	 */
	public void sort()
	{
		Collections.sort(filters);
	}

	/**
	 * Check all filters in the set are the same type. If so find the parameters that are the weakest and generate a new
	 * filter.
	 *
	 * @return The weakest filter
	 */
	public Filter createWeakestFilter()
	{
		if (!initialisedWeakest)
		{
			weakest = createWeakest();
			initialisedWeakest = true;
		}
		return weakest;
	}

	private Filter createWeakest()
	{
		if (!allSameType())
			return null;

		// Initialise the parameters
		final Filter f1 = filters.get(0);
		final double[] parameters = f1.getParameters();

		// Find the weakest
		for (final Filter f : filters)
			f.weakestParameters(parameters);

		return f1.create(parameters);
	}

	/**
	 * @return True if all the filters are the same type
	 */
	public boolean allSameType()
	{
		if (allSameType == 0)
			allSameType = checkAllSameType();
		return (allSameType == 1);
	}

	/**
	 * @return 1 if all the filters are the same type, -1 otherwise
	 */
	private int checkAllSameType()
	{
		if (size() == 0)
			return -1;

		// Check for the same type
		final String type = filters.get(0).getType();
		for (final Filter f : filters)
			// Use the != since the Strings should be immutable
			//if (f.getType() != type)
			if (!f.getType().equals(type))
				return -1;
		return 1;
	}
}
