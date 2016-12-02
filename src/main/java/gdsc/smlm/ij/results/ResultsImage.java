package gdsc.smlm.ij.results;

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

public enum ResultsImage
{
	//@formatter:off
	NONE{ public String getName() { return "None"; }}, 
	LOCALISATIONS{ public String getName() { return "Localisations"; }}, 
	SIGNAL_INTENSITY{ public String getName() { return "Signal intensity"; }},
	FRAME_NUMBER{ public String getName() { return "Frame number"; }},
	PSF{ public String getName() { return "PSF"; }}, 
	LOCALISATIONS_PRECISION{ public String getName() { return "Localisations (width=precision)"; }}, 
	SIGNAL_PRECISION{ public String getName() { return "Signal (width=precision)"; }}, 
	LOCALISATIONS_AV_PRECISION{ public String getName() { return "Localisations (width=av.precision)"; }}, 
	SIGNAL_AV_PRECISION{ public String getName() { return "Signal (width=av.precision)"; }}, 
	ERROR{ public String getName() { return "Fit error"; }};
	//@formatter:on

	@Override
	public String toString()
	{
		return getName();
	}

	/**
	 * Gets the name.
	 *
	 * @return the name
	 */
	abstract public String getName();
}
