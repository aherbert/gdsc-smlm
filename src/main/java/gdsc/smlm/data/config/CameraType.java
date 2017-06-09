package gdsc.smlm.data.config;

import gdsc.smlm.data.NamedObject;

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
 * Define a camera type.
 */
public enum CameraType implements NamedObject
{
	/**
	 * Charge Coupled Device (CCD).
	 */
	CCD
	{
		public String getName()
		{
			return "CCD";
		}
	},

	/**
	 * Electron Multiplying Charge Coupled Device (EM CCD).
	 */
	EM_CCD
	{
		public String getName()
		{
			return "EM-CCD";
		}
	},

	/** Scientific Complementary Metal-Oxide-Semiconductor (sCMOS). */
	SCMOS
	{
		public String getName()
		{
			return "sCMOS";
		}
	},;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.NamedObject#getShortName()
	 */
	public String getShortName()
	{
		return getName();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Enum#toString()
	 */
	public String toString()
	{
		return getName();
	}
}