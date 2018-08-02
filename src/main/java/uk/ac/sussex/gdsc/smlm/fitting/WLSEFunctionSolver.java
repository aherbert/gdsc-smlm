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
package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Defines methods to fit a function with coefficients (a) using weighted least-squares estimation.
 */
public interface WLSEFunctionSolver extends FunctionSolver
{
    /**
     * Gets the chi squared.
     *
     * @return the chi squared
     */
    public double getChiSquared();

    /**
     * Gets the probability Q that a value of chi-square as poor as the value should occur by chance.
     * <p>
     * A low value indicates greater statistical significance, i.e. greater confidence that the observed deviation from
     * the null hypothesis is significant, with the null model being that the fit is good. The confidence in rejecting
     * the null hypothesis is 100 * (1 - q) percent.
     *
     * @return the q-value
     */
    public double getQ();
}
