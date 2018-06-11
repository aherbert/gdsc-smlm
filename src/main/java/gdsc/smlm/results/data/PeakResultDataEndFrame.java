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
package gdsc.smlm.results.data;

import gdsc.smlm.results.PeakResult;

/**
 * Gets the end frame from a result.
 */
public class PeakResultDataEndFrame extends PeakResultDataInteger
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultData#getValue(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public Integer getValue(PeakResult result)
	{
		return result.getEndFrame();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResultData#getValueName()
	 */
	@Override
	public String getValueName()
	{
		return "End Frame";
	}
}
