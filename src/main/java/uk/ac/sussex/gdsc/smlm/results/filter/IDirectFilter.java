/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Support direct filtering of PreprocessedPeakResult objects.
 *
 * <p>The decision to support for filtering as both a DirectFilter and Filter concurrently is left
 * to the implementing class. It is not a requirement.
 */
public interface IDirectFilter {
  /**
   * Gets the flags indicating all the fields that are used during validation. These flags may be
   * returned by the filter {@link #validate(PreprocessedPeakResult)} method if the result fails
   * validation.
   *
   * @return the validation flags
   */
  int getValidationFlags();

  /**
   * Called before the accept method is called for PreprocessedPeakResult.
   *
   * <p>This should be called once to initialise the filter before processing a batch of results.
   *
   * @see #validate(PreprocessedPeakResult)
   */
  void setup();

  /**
   * Called before the accept method is called for PreprocessedPeakResult. The flags can control the
   * type of filtering requested. Filters are asked to respect the flags defined in this class.
   *
   * <p>This should be called once to initialise the filter before processing a batch of results.
   *
   * @param flags Flags used to control the filter
   * @see #validate(PreprocessedPeakResult)
   */
  void setup(final int flags);

  /**
   * Called before the accept method is called for PreprocessedPeakResult. The filter data can
   * control the type of filtering requested.
   *
   * <p>This should be called once to initialise the filter before processing a batch of results.
   *
   * @param flags Flags used to control the filter
   * @param filterSetupData Data used to control the filter
   * @see #validate(PreprocessedPeakResult)
   */
  void setup(final int flags, final FilterSetupData... filterSetupData);

  /**
   * Gets the flags required to reinitialise the current filter state using {@link #setup(int)} or
   * {@link #setup(int, FilterSetupData...)}.
   *
   * @return the flags
   * @throws IllegalStateException If setup has not been called and the flags cannot be created
   */
  int getFilterSetupFlags();

  /**
   * Gets the filter setup data required to reinitialise the current filter state using
   * {@link #setup(int, FilterSetupData...)}.
   *
   * @return the filter setup data (can be null)
   * @throws IllegalStateException If setup has not been called and the data cannot be created
   */
  FilterSetupData[] getFilterSetupData();

  /**
   * Filter the peak result.
   *
   * <p>Calls {@link #validate(PreprocessedPeakResult)} and stores the result. This can be obtained
   * using {@link #getResult()}.
   *
   * @param peak The peak result
   * @return true if the peak should be accepted
   */
  boolean accept(final PreprocessedPeakResult peak);

  /**
   * Filter the peak result.
   *
   * @param peak The peak result
   * @return zero if the peak should be accepted, otherwise set to flags indicating the field that
   *         failed validation.
   */
  int validate(final PreprocessedPeakResult peak);

  /**
   * Return the type of filter. This should be a DirectFilter.
   *
   * @return Should return DirectFilter
   */
  FilterType getFilterType();

  /**
   * Return the result flag generated during the last call to
   * {@link #accept(PreprocessedPeakResult)}.
   *
   * @return the validation result from the last call to {@link #accept(PreprocessedPeakResult)}
   */
  int getResult();

  /**
   * Copy this filter.
   *
   * @return the copy
   */
  IDirectFilter copy();
}
