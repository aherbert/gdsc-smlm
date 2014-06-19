package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import java.util.Collections;
import java.util.List;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;

/**
 * Specify a named set of filters
 */
public class FilterSet
{
	@XStreamAsAttribute
	private String name;
	private List<Filter> filters;

	public FilterSet(String name, List<Filter> filters)
	{
		this.name = name;
		this.filters = filters;
	}

	public FilterSet(List<Filter> filters)
	{
		this.filters = filters;
	}

	public int size()
	{
		return filters.size();
	}

	/**
	 * The name of the set. If empty return the name of the first filter in the set.
	 * @return the name
	 */
	public String getName()
	{
		if (name == null || name.length() == 0)
		{
			if (filters != null && !filters.isEmpty())
				name = filters.get(0).getName();
		}
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
	 * Create the filter set from the XML representation
	 * 
	 * @param xml
	 * @return the filter set
	 */
	public static FilterSet fromXML(String xml)
	{
		try
		{
			return (FilterSet) XStreamWrapper.fromXML(xml);
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		return null;
	}

	public void sort()
	{
		Collections.sort(filters);
	}
}