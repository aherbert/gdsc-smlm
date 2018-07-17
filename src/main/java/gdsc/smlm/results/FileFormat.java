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
package gdsc.smlm.results;

/**
 * Specify the file format when reading results from file
 */
public enum FileFormat
{
	//@formatter:off
	/** SMLM Text */
    SMLM_TEXT{ 
    	@Override public String getName() { return "SMLM Text"; } 
        @Override public boolean isSMLM(){return true;}},
    /** SMLM Binary */
    SMLM_BINARY{ 
        @Override public String getName() { return "SMLM Binary"; } 
        @Override public boolean isSMLM(){return true;}},
    /** RapidSTORM */
    RAPID_STORM{ @Override public String getName() { return "RapidSTORM"; }},
    /** NSTORM */
    NSTORM{ @Override public String getName() { return "NSTORM"; }},
    /** SMLM Table */
    SMLM_TABLE{ @Override public String getName() { return "SMLM Table"; }},
    /** MALK */
    MALK{ @Override public String getName() { return "MALK"; }},
    /** TSF Binary */
    TSF_BINARY{ @Override public String getName() { return "TSF Binary"; }},
    /** Unknown */
    UNKNOWN{ @Override public String getName() { return "Unknown"; }};
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

	/**
	 * Checks if is a GDSC SMLM format.
	 *
	 * @return true, if is a GDSC SMLM format
	 */
	public boolean isSMLM()
	{
		return false;
	}
}
