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
	NONE("None"), LOCALISATIONS("Localisations"), SIGNAL_INTENSITY("Signal intenisty"), FRAME_NUMBER("Frame number"), PSF(
			"PSF"), LOCALISATIONS_PRECISION("Localisations (width=precision)"), SIGNAL_PRECISION(
			"Signal (width=precision)"), LOCALISATIONS_AV_PRECISION("Localisations (width=av.precision)"), SIGNAL_AV_PRECISION(
			"Signal (width=av.precision)"), ERROR("Fit error");

	private String name;

	private ResultsImage(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}
