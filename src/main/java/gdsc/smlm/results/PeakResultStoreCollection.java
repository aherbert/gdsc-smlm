package gdsc.smlm.results;

import java.util.Collection;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Marker interface for PeakResultStore objects that use a collection.
 */
public interface PeakResultStoreCollection
{
	/**
	 * Gets the collection.
	 *
	 * @return the collection
	 */
	public Collection<PeakResult> getCollection();
	
	/**
	 * Gets the collection by reference.
	 *
	 * @return the collection
	 */
	Collection<PeakResult> getCollectionReference();
}
