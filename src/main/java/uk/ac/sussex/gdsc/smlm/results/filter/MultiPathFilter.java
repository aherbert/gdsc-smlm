/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.match.FractionClassificationResult;
import uk.ac.sussex.gdsc.core.match.FractionalAssignment;
import uk.ac.sussex.gdsc.core.match.RankedScoreCalculator;
import uk.ac.sussex.gdsc.smlm.results.count.ConsecutiveFailCounter;
import uk.ac.sussex.gdsc.smlm.results.count.FailCounter;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult.FitResult;

/**
 * Filter a multi-path set of peak results into accepted/rejected.
 */
public class MultiPathFilter {

  private static final FailCounter defaultFailCounter = ConsecutiveFailCounter.create(0);

  /** The direct filter to apply to the results. */
  final IDirectFilter filter;

  /**
   * The minimal direct filter to apply to the results.
   *
   * <p>This is applied if the result fails the primary filter. It is used to indicate that the
   * result achieves a minimum set of criteria.
   */
  final IDirectFilter minFilter;

  /**
   * The residuals threshold to consider the residuals Quadrant Analysis (QA) score of a single for
   * doublet fitting. The score should range from 0 to 1. A score equal or above 1 will ignore
   * doublet fitting.
   */
  @XStreamAsAttribute
  public final double residualsThreshold;

  /** The validation results used when filtering a result. */
  @XStreamOmitField
  private int[] validationResults;
  @XStreamOmitField
  private boolean failExisting;
  @XStreamOmitField
  private boolean failNew;

  /**
   * Stores the results that were accepted when filtering a multi-path result. Also stores the fit
   * result that was used to select the results.
   */
  public static class SelectedResult {

    /** The results. */
    public final PreprocessedPeakResult[] results;

    /** The fit result. */
    public final FitResult fitResult;

    /**
     * Instantiates a new selected result.
     *
     * @param results the results
     * @param fitResult the fit result
     */
    public SelectedResult(PreprocessedPeakResult[] results, FitResult fitResult) {
      this.results = results;
      this.fitResult = fitResult;
    }
  }

  /**
   * Allows storage of results that have been selected during multi-path filtering.
   */
  public interface SelectedResultStore {
    /**
     * Add a selected result to the store.
     *
     * @param selectedResult the selected result
     */
    void add(SelectedResult selectedResult);

    /**
     * Checks if is fit. Any candidate that has already been fit will not be stored.
     *
     * @param candidateId the candidate id
     * @return True if the candidate has been fit
     */
    boolean isFit(int candidateId);

    /**
     * Checks if is valid.
     *
     * <p>Return true if this candidate should definitely be filtered.
     *
     * @param candidateId the candidate id
     * @return true, if is valid
     */
    boolean isValid(int candidateId);

    /**
     * A result that passed the primary filter.
     *
     * @param result the result
     */
    void pass(PreprocessedPeakResult result);

    /**
     * A result that passed the minimal filter.
     *
     * @param result the result
     */
    void passMin(PreprocessedPeakResult result);
  }

  /**
   * Allow tracking of candidates that have been fit.
   */
  private static class SimpleSelectedResultStore implements SelectedResultStore {
    /** The is fit. */
    boolean[] isFit;

    /** The is valid. */
    boolean[] isValid;

    /**
     * Instantiates a new simple selected result store.
     */
    SimpleSelectedResultStore() {
      isFit = new boolean[0];
      isValid = new boolean[0];
    }

    /**
     * Instantiates a new simple selected result store.
     *
     * @param totalCandidates the total candidates
     */
    SimpleSelectedResultStore(int totalCandidates) {
      isFit = new boolean[totalCandidates];
      isValid = new boolean[totalCandidates];
    }

    @Override
    public void add(SelectedResult selectedResult) {
      final PreprocessedPeakResult[] results = selectedResult.results;
      if (results == null) {
        return;
      }
      for (int i = 0; i < results.length; i++) {
        if (results[i].isNewResult()) {
          isFit[results[i].getCandidateId()] = true;
        }
      }
    }

    @Override
    public boolean isFit(int candidateId) {
      return isFit[candidateId];
    }

    @Override
    public boolean isValid(int candidateId) {
      return isValid[candidateId];
    }

    @Override
    public void pass(PreprocessedPeakResult result) {
      //// If we are debugging
      // if (result.getFrame() == 46 && result.getCandidateId() == 70)
      // {
      // System.out.println("Follow this...");
      // }

      // This an existing result or candidate. Mark as valid so candidates will be processed
      isValid[result.getCandidateId()] = true;
    }

    @Override
    public void passMin(PreprocessedPeakResult result) {
      // Passing the minimal filter does not mean it is valid. This would be used to store
      // a fit estimate during processing for this candidate.
    }

    /**
     * Resize.
     *
     * @param totalCandidates the total candidates
     */
    public void resize(int totalCandidates) {
      if (isFit.length < totalCandidates) {
        isFit = new boolean[totalCandidates];
        isValid = new boolean[totalCandidates];
      } else {
        for (int i = 0; i < totalCandidates; i++) {
          isFit[i] = false;
          isValid[i] = false;
        }
      }
    }
  }

  /**
   * Used to return default behaviour for acceptAny/acceptAll.
   */
  private static class NullSelectedResultStore implements SelectedResultStore {
    /** The null selected result store. */
    static final NullSelectedResultStore INSTANCE = new NullSelectedResultStore();

    @Override
    public void add(SelectedResult selectedResult) {
      // Do nothing
    }

    @Override
    public boolean isFit(int candidateId) {
      // Make sure non-candidate fits are ignored.
      return true;
    }

    @Override
    public boolean isValid(int candidateId) {
      return false;
    }

    @Override
    public void pass(PreprocessedPeakResult result) {
      // Do nothing
    }

    @Override
    public void passMin(PreprocessedPeakResult result) {
      // Do nothing
    }
  }

  /**
   * Allows signalling of results that have been selected during multi-path filter scoring.
   */
  @FunctionalInterface
  public interface FractionScoreStore {
    /**
     * Add the unique Id of a result that was selected.
     *
     * @param uniqueId the unique id
     */
    void add(int uniqueId);
  }

  /**
   * Used to return default behaviour.
   */
  private static class NullFractionScoreStore implements FractionScoreStore {
    /** The null fraction result store. */
    static final NullFractionScoreStore INSTANCE = new NullFractionScoreStore();

    @Override
    public void add(int uniqueId) {
      // Do nothing
    }
  }

  private class FilterSetupState {
    /** Store the initial state of the filter flags. */
    int setupFlags;

    FilterSetupState(int flags) {
      this.setupFlags = flags;
    }

    void restoreState() {
      filter.setup(setupFlags);
    }
  }

  private class FilterSetupDataState extends FilterSetupState {
    /** Store the initial state of the filter setup data. */
    FilterSetupData[] setupData;

    FilterSetupDataState(int flags, FilterSetupData[] setupData) {
      super(flags);
      this.setupData = setupData;
    }

    @Override
    void restoreState() {
      filter.setup(setupFlags, setupData);
    }
  }

  private FilterSetupState filterSetupState;

  /**
   * Save the initial setup state of the main filter.
   */
  private void saveState() {
    final int setupFlags = filter.getFilterSetupFlags();
    final FilterSetupData[] setupData = filter.getFilterSetupData();
    filterSetupState = (setupData == null) ? new FilterSetupState(setupFlags)
        : new FilterSetupDataState(setupFlags, setupData);
  }

  /**
   * Sets up the main filter.
   *
   * @param flags the filter flags
   */
  private void setupFilter(int flags) {
    filter.setup(flags);
  }

  /**
   * Restore the initial setup state of the main filter.
   */
  private void restoreFilterState() {
    Objects.requireNonNull(filterSetupState,
        "Unknown initial filter state. setup(...) must be called before using the filter methods.");
    filterSetupState.restoreState();
  }

  /**
   * Create a new MultiPathFilter.
   *
   * @param filter the direct filter for filtering the results
   * @param residualsThreshold The residuals threshold to consider a single fit for doublet fitting
   */
  public MultiPathFilter(IDirectFilter filter, double residualsThreshold) {
    this(filter, null, residualsThreshold);
  }

  /**
   * Create a new MultiPathFilter.
   *
   * @param filter the direct filter for filtering the results
   * @param minFilter the minimal direct filter for filtering the results
   * @param residualsThreshold The residuals threshold to consider a single fit for doublet fitting
   */
  public MultiPathFilter(IDirectFilter filter, IDirectFilter minFilter, double residualsThreshold) {
    this.filter = filter;
    this.minFilter = minFilter;
    this.residualsThreshold = residualsThreshold;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected MultiPathFilter(MultiPathFilter source) {
    this.filter = copyFilter(source.filter);
    this.minFilter = copyFilter(source.minFilter);
    this.residualsThreshold = source.residualsThreshold;
  }

  /**
   * Return a deep copy of this object with a copy of the configured filters.
   *
   * @return the copy
   */
  public MultiPathFilter copy() {
    return new MultiPathFilter(this);
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    final MultiPathFilter other = (MultiPathFilter) obj;
    if (this.residualsThreshold != other.residualsThreshold) {
      return false;
    }
    if ((this.filter == null) ? (other.filter != null) : !this.filter.equals(other.filter)) {
      return false;
    }
    if (this.minFilter == null) {
      return other.minFilter == null;
    }
    return this.minFilter.equals(other.minFilter);
  }

  @Override
  public int hashCode() {
    int hash = Double.hashCode(residualsThreshold);
    if (filter != null) {
      hash = hash * 31 + filter.hashCode();
    }
    if (minFilter != null) {
      hash = hash * 31 + minFilter.hashCode();
    }
    return hash;
  }

  private static IDirectFilter copyFilter(IDirectFilter filter) {
    return (filter == null) ? null : filter.copy();
  }

  /**
   * Gets the filter.
   *
   * @return the filter
   */
  public IDirectFilter getFilter() {
    return copyFilter(filter);
  }

  /**
   * Gets the minimal filter.
   *
   * @return the minimal filter
   */
  public IDirectFilter getMinimalFilter() {
    return copyFilter(minFilter);
  }

  /**
   * Called before the accept method is called for PreprocessedPeakResult. This calls the setup()
   * method in the DirectFilter. It also saves the initial state of the filter so that it can be
   * restored when the filter is updated, for example disabling
   *
   * <p>This should be called once to initialise the filter before processing a batch of results.
   *
   * @see #accept(PreprocessedPeakResult)
   */
  public void setup() {
    filter.setup();
    if (minFilter != null) {
      minFilter.setup();
    }
    saveState();
  }

  /**
   * Called before the accept method is called for PreprocessedPeakResult. The flags can control the
   * type of filtering requested. Filters are asked to respect the flags defined in this class. This
   * calls the setup(int) method in the DirectFilter.
   *
   * <p>This should be called once to initialise the filter before processing a batch of results.
   *
   * @param flags Flags used to control the filter
   * @see #accept(PreprocessedPeakResult)
   */
  public void setup(final int flags) {
    filter.setup(flags);
    if (minFilter != null) {
      minFilter.setup(flags);
    }
    saveState();
  }

  /**
   * Called before the accept method is called for PreprocessedPeakResult. The flags can control the
   * type of filtering requested. Filters are asked to respect the flags defined in this class. This
   * calls the setup(int) method in the DirectFilter.
   *
   * <p>This should be called once to initialise the filter before processing a batch of results.
   *
   * @param flags Flags used to control the filter
   * @param data the data used to control the filter
   * @see #accept(PreprocessedPeakResult)
   */
  public void setup(final int flags, FilterSetupData... data) {
    filter.setup(flags, data);
    if (minFilter != null) {
      minFilter.setup(flags, data);
    }
    saveState();
  }

  /**
   * Filter the peak result. This calls the accept() method in the DirectFilter.
   *
   * @param peak The peak result
   * @return true if the peak should be accepted, otherwise false to reject.
   */
  public boolean accept(final PreprocessedPeakResult peak) {
    return filter.accept(peak);
  }

  /**
   * Filter a multi-path set of peak results into a set that are accepted.
   *
   * <p>Any existing or new results must pass the {@link #accept(PreprocessedPeakResult)} method.
   * Any other results are assumed to be candidates that were fitted but will not be validated
   * unless required.
   *
   * <p>Note that new results may not be for the candidate identified by the MultiPathFitResult.
   * This can happen when multi-fitting has fit another candidate that previously did not have a
   * result. The SelectedResultStore is used to determine if that result has been fit already. If
   * not it is added to the output list.
   *
   * <p>The SelectedResultStore will be passed any result that passes the configured filters.
   *
   * @param multiPathResult the multi path result
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @return The new peak results that are accepted (and any valid candidates if found); or null
   */
  public final PreprocessedPeakResult[] accept(final MultiPathFitResult multiPathResult,
      boolean validateCandidates, SelectedResultStore store) {
    return accept(multiPathResult, validateCandidates, store, false);
  }

  /**
   * Filter a multi-path set of peak results into a set that are accepted.
   *
   * <p>Any existing or new results must pass the {@link #accept(PreprocessedPeakResult)} method.
   * Any other results are assumed to be candidates that were fitted but will not be validated
   * unless required.
   *
   * <p>Note that new results may not be for the candidate identified by the MultiPathFitResult.
   * This can happen when multi-fitting has fit another candidate that previously did not have a
   * result. The SelectedResultStore is used to determine if that result has been fit already. If
   * not it is added to the output list.
   *
   * <p>The SelectedResultStore will be passed any result that passes the configured filters.
   *
   * @param multiPathResult the multi path result
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @param precomputed True if this is a subset with pre-computed validation results
   * @return The new peak results that are accepted (and any valid candidates if found); or null
   */
  private final PreprocessedPeakResult[] accept(final MultiPathFitResult multiPathResult,
      boolean validateCandidates, SelectedResultStore store, boolean precomputed) {
    final int candidateId = multiPathResult.getCandidateId();

    // Ensure we don't have to check the store in acceptAll/acceptAny
    if (store == null) {
      store = NullSelectedResultStore.INSTANCE;
    }

    // The aim is to obtain a new result for the current candidate Id.
    // acceptAll/acceptAny will return all new results, even if they do not match the candidate.
    // So we check the candidate Id and return when we have a new result for the candidate.
    // We accept the doublet fit over the single fit if we are performing doublet fitting.
    // If nothing matches then pick the result with the most new results, or use the default
    // order we processed the fits.

    boolean doDoublet = false;

    // Filter multi-fit

    // Accept all and then check if we can perform a doublet fit
    // final PreprocessedPeakResult[] multiResults = acceptAll(candidateId,
    // multiPathResult.getMultiFitResult(),
    // validateCandidates, store);
    // if (multiResults == null)
    // {
    // // The fit was not accepted. However it may have been rejected for being too wide
    // // and is suitable for a doublet fit.
    // doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(),
    // false);
    // }
    // else
    // {
    // doDoublet = (residualsThreshold < 1 && multiPathResult.getMultiQAScore() >
    // residualsThreshold);
    // }

    // Accept any and then check if we can perform a doublet fit
    final PreprocessedPeakResult[] multiResults = acceptAny(candidateId,
        multiPathResult.getMultiFitResult(), validateCandidates, store, precomputed);
    doDoublet =
        isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(), false);

    final PreprocessedPeakResult[] multiDoubletResults;
    if (doDoublet) {
      multiDoubletResults =
          acceptAnyDoublet(multiPathResult, validateCandidates, store, candidateId);
      // Check we have a new result for the candidate
      if (multiDoubletResults != null && contains(multiDoubletResults, candidateId)) {
        return multiDoubletResults;
      }
    } else {
      multiDoubletResults = null;
    }

    // Check if the multi result is to the correct candidate.
    if (multiResults != null && contains(multiResults, candidateId)) {
      return multiResults;
    }

    // We reached here with:
    // a multi fit that failed or matched a different candidate
    // a doublet multi fit that failed or matched a different candidate

    // Filter single-fit
    final PreprocessedPeakResult[] singleResults = acceptAll(candidateId,
        multiPathResult.getSingleFitResult(), validateCandidates, store, precomputed);
    if (singleResults == null) {
      // The fit was not accepted. However it may have been rejected for being too wide
      // and is suitable for a doublet fit.
      doDoublet =
          isSuitableForDoubletFit(multiPathResult, multiPathResult.getSingleFitResult(), true);
    } else {
      // The single fit is OK.
      doDoublet =
          (residualsThreshold < 1 && multiPathResult.getSingleQaScore() > residualsThreshold);
    }

    // We reached here with:
    // a multi fit that failed or matched a different candidate
    // a doublet multi fit that failed or matched a different candidate
    // a single fit that is eligible for doublet fitting, it may be null (if it passed without width
    // filtering)

    final PreprocessedPeakResult[] singleDoubletResults;
    if (doDoublet) {
      // We must validate the spot without shift filtering. Doublets may drift further than single
      // spot candidates.
      setupFilter(FilterValidationOption.NO_SHIFT);
      singleDoubletResults = acceptAny(candidateId, multiPathResult.getDoubletFitResult(),
          validateCandidates, store, precomputed);
      restoreFilterState();
      // Check we have a new result for the candidate
      if (singleDoubletResults != null && contains(singleDoubletResults, candidateId)) {
        return singleDoubletResults;
      }
    } else {
      singleDoubletResults = null;
    }

    // Check if the single result is to the correct candidate.
    if (singleResults != null && contains(singleResults, candidateId)) {
      return singleResults;
    }

    // We reached here with:
    // a multi fit that failed or matched a different candidate
    // a multi doublet fit that failed or matched a different candidate
    // a single fit that failed or matched a different candidate
    // a doublet fit that failed or matched a different candidate
    return rank(multiResults, multiDoubletResults, singleResults, singleDoubletResults);
  }

  /**
   * Filter the peak result. This calls the accept() method in the minimal DirectFilter.
   *
   * @param peak The peak result
   * @return true if the peak should be accepted, otherwise false to reject.
   */
  private boolean minAccept(final PreprocessedPeakResult peak) {
    return minFilter.accept(peak);
  }

  /**
   * Allows results to be ranked.
   */
  private static class ResultRank {
    /** The results. */
    final PreprocessedPeakResult[] results;

    /** The default rank (when the count of new result is the same). */
    final int rank;

    /** The count of new results. */
    final int count;

    /**
     * Instantiates a new result rank.
     *
     * @param results the results
     * @param rank the rank
     */
    public ResultRank(PreprocessedPeakResult[] results, int rank) {
      this.results = results;
      this.rank = rank;
      if (results == null) {
        // Negative so null results are ranked below not-null results with no new results
        count = -1;
      } else {
        count = countNewResult(results);
      }
    }

    static int compare(ResultRank o1, ResultRank o2) {
      if (o2.count < o1.count) {
        return -1;
      }
      if (o2.count > o1.count) {
        return 1;
      }
      return Integer.compare(o1.rank, o2.rank);
    }
  }

  /**
   * Rank the results. It is assumed that each result is either null or has results that do not
   * match the current candidate Id. In this case we will return the result that has the highest
   * number of new results. In the event of a tie we order as multi then single, with doublet
   * results ahead of non-doublet results. Results are discounted if null.
   *
   * @param multiResults the multi results
   * @param multiDoubletResults the multi doublet results
   * @param singleResults the single results
   * @param singleDoubletResults the doublet results
   * @return the preprocessed peak result[]
   */
  private static @Nullable PreprocessedPeakResult[] rank(PreprocessedPeakResult[] multiResults,
      PreprocessedPeakResult[] multiDoubletResults, PreprocessedPeakResult[] singleResults,
      PreprocessedPeakResult[] singleDoubletResults) {
    if (multiResults == null && multiDoubletResults == null && singleResults == null
        && singleDoubletResults == null) {
      return null;
    }
    final ResultRank[] rank = new ResultRank[4];
    // rank[0] = new ResultRank(multiResults, 1);
    // rank[1] = new ResultRank(multiDoubletResults, 2);
    // rank[2] = new ResultRank(singleResults, 3);
    // rank[3] = new ResultRank(singleDoubletResults, 4);
    rank[0] = new ResultRank(multiDoubletResults, 1);
    rank[1] = new ResultRank(multiResults, 2);
    rank[2] = new ResultRank(singleDoubletResults, 3);
    rank[3] = new ResultRank(singleResults, 4);
    Arrays.sort(rank, ResultRank::compare);
    return rank[0].results;
  }

  /**
   * Filter a multi-path set of peak results into a set that are accepted.
   *
   * <p>Any existing or new results must pass the {@link #accept(PreprocessedPeakResult)} method.
   * Any other results are assumed to be candidates that were fitted but will not be validated
   * unless required.
   *
   * <p>Note that new results may not be for the candidate identified by the MultiPathFitResult.
   * This can happen when multi-fitting has fit another candidate that previously did not have a
   * result. The SelectedResultStore is used to determine if that result has been fit already. If
   * not it is added to the output list.
   *
   * <p>The method returns the the same results as
   * {@link #accept(MultiPathFitResult, boolean, SelectedResultStore)} but includes the FitResult
   * that the data originated from.
   *
   * <p>The SelectedResultStore will be passed any result that passes the configured filters. It
   * will not be passed the returned SelectedResult as the results will be duplicates of those
   * passed to the store individually. They may also contain validated candidates. The returned
   * results must thus be filtered for new results (e.g. not existing or candidate results).
   *
   * @param multiPathResult the multi path result
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @return The results that are accepted; or null
   */
  public final @Nullable SelectedResult select(final MultiPathFitResult multiPathResult,
      boolean validateCandidates, SelectedResultStore store) {
    final int candidateId = multiPathResult.getCandidateId();

    // Ensure we don't have to check the store in acceptAll/acceptAny
    if (store == null) {
      store = NullSelectedResultStore.INSTANCE;
    }

    // The aim is to obtain a new result for the current candidate Id.
    // acceptAll/acceptAny will return all new results, even if they do not match the candidate.
    // So we check the candidate Id and return when we have a new result for the candidate.
    // We accept the doublet fit over the single fit if we are performing doublet fitting.
    // If nothing matches then pick the result with the most new results, or use the default
    // order we processed the fits.

    boolean doDoublet = false;

    // Filter multi-fit
    // Accept all and then check if we can perform a doublet fit
    // final PreprocessedPeakResult[] multiResults = acceptAll(candidateId,
    // multiPathResult.getMultiFitResult(),
    // validateCandidates, store);
    // if (multiResults == null)
    // {
    // // The fit was not accepted. However it may have been rejected for being too wide
    // // and is suitable for a doublet fit.
    // doDoublet = isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(),
    // false);
    // }
    // else
    // {
    // doDoublet = (residualsThreshold < 1 && multiPathResult.getMultiQAScore() >
    // residualsThreshold);
    // }

    // Accept any and then check if we can perform a doublet fit
    final PreprocessedPeakResult[] multiResults =
        acceptAny(candidateId, multiPathResult.getMultiFitResult(), validateCandidates, store);
    doDoublet =
        isSuitableForDoubletFit(multiPathResult, multiPathResult.getMultiFitResult(), false);

    final PreprocessedPeakResult[] multiDoubletResults;
    if (doDoublet) {
      multiDoubletResults =
          acceptAnyDoublet(multiPathResult, validateCandidates, store, candidateId);
      if (multiDoubletResults != null) {
        // Check we have a new result for the candidate
        if (contains(multiDoubletResults, candidateId)) {
          return new SelectedResult(multiDoubletResults,
              multiPathResult.getMultiDoubletFitResult());
        }
      }
    } else {
      multiDoubletResults = null;
    }

    // Check if the multi result is to the correct candidate.
    if (multiResults != null && contains(multiResults, candidateId)) {
      return new SelectedResult(multiResults, multiPathResult.getMultiFitResult());
    }

    // We reached here with:
    // a multi fit that failed or matched a different candidate
    // a doublet multi fit that failed or matched a different candidate

    doDoublet = false;

    // Filter single-fit
    final PreprocessedPeakResult[] singleResults =
        acceptAll(candidateId, multiPathResult.getSingleFitResult(), validateCandidates, store);
    if (singleResults == null) {
      // The fit was not accepted. However it may have been rejected for being too wide
      // and is suitable for a doublet fit.
      doDoublet =
          isSuitableForDoubletFit(multiPathResult, multiPathResult.getSingleFitResult(), true);
    } else {
      // The single fit is OK.
      doDoublet =
          (residualsThreshold < 1 && multiPathResult.getSingleQaScore() > residualsThreshold);
    }

    // We reached here with:
    // a multi fit that failed or matched a different candidate
    // a doublet multi fit that failed or matched a different candidate
    // a single fit that is eligible for doublet fitting, it may be null (if it passed without width
    // filtering)

    final PreprocessedPeakResult[] singleDoubletResults;
    if (doDoublet) {
      // We must validate the spot without shift filtering. Doublets may drift further than single
      // spot candidates.
      setupFilter(FilterValidationOption.NO_SHIFT);
      singleDoubletResults =
          acceptAny(candidateId, multiPathResult.getDoubletFitResult(), validateCandidates, store);
      restoreFilterState();
      if (singleDoubletResults != null) {
        // Check we have a new result for the candidate
        if (contains(singleDoubletResults, candidateId)) {
          return new SelectedResult(singleDoubletResults, multiPathResult.getDoubletFitResult());
        }
      }
    } else {
      singleDoubletResults = null;
    }

    // Check if the single result is to the correct candidate.
    if (singleResults != null && contains(singleResults, candidateId)) {
      return new SelectedResult(singleResults, multiPathResult.getSingleFitResult());
    }

    // We reached here with:
    // a multi fit that failed or matched a different candidate
    // a multi doublet fit that failed or matched a different candidate
    // a single fit that failed or matched a different candidate
    // a doublet fit that failed or matched a different candidate
    final PreprocessedPeakResult[] result =
        rank(multiResults, multiDoubletResults, singleResults, singleDoubletResults);
    if (result == null) {
      return null;
    }
    //@formatter:off
    if (result == multiResults) {
      return new SelectedResult(multiResults,        multiPathResult.getMultiFitResult());
    }
    if (result == multiDoubletResults) {
      return new SelectedResult(multiDoubletResults, multiPathResult.getMultiDoubletFitResult());
    }
    if (result == singleResults) {
      return new SelectedResult(singleResults,       multiPathResult.getSingleFitResult());
    }
    return new SelectedResult(singleDoubletResults,  multiPathResult.getDoubletFitResult());
    //@formatter:on
  }

  /**
   * Select a set of peak results.
   *
   * <p>The results are processed in order. Results are only processed if either the fail counter
   * {@link FailCounter#isOk() isOK} or the store determines the result
   * {@link SelectedResultStore#isValid(int) isValid}. Results are selected using the configured
   * filters and may be excluded using the coordinate store. Selected results are passed to the
   * selected results store. The fail counter is sent a pass/fail signal for each multi-path result.
   *
   * <p>A selected result will be stored for each MultiPathFitResult that is assessed, even if the
   * fitting failed. In this case the list of accepted results will be null.
   *
   * <p>The SelectedResultStore can be used to track results that pass validation. If this is null
   * then the default behaviour is to track candidates that were successfully fit and that pass
   * validation. These will be marked as valid and will be processed even if the fail counter is not
   * OK. This allows filtering of candidate known to be good that may have been fit as a multi-fit
   * result for a previous candidate.
   *
   * <p>The coordinate store can be used to check for duplicates.
   *
   * @param multiPathResults the multi path results
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param setup Set to true to run the {@link #setup()} method
   * @param store the store (can be used to track results that pass validation)
   * @param coordinateStore the coordinate store (can be null)
   */
  public void select(final IMultiPathFitResults multiPathResults, FailCounter failCounter,
      boolean setup, SelectedResultStore store, CoordinateStore coordinateStore) {
    if (setup) {
      setup();
    }

    failCounter = replaceIfNull(failCounter);
    if (store == null) {
      store = new SimpleSelectedResultStore(multiPathResults.getTotalCandidates());
    }
    coordinateStore = NullCoordinateStore.replaceIfNull(coordinateStore);

    // // Debugging the results that are scored
    // java.io.OutputStreamWriter out = null;
    // if (debugFilename != null && multiPathResults.getFrame() == 46)
    // {
    // try
    // {
    // out = new java.io.OutputStreamWriter(new java.io.FileOutputStream(debugFilename), "UTF-8");
    // }
    // catch (Exception e)
    // {
    // }
    // }

    // Note - this could be made iterative. Any pass through the data may store estimates
    // using the SelectedResultStore and used to determine if

    failCounter.reset();
    final int total = multiPathResults.getNumberOfResults();
    // while (multiPathResults.begin())
    // {
    coordinateStore.clear();
    for (int c = 0; c < total; c++) {
      final MultiPathFitResult multiPathResult = multiPathResults.getResult(c);

      // if (multiPathResult == null)
      // {
      // multiPathResults.complete(c);
      // // Ignore this but do not count it as a failure
      // continue;
      // }

      final boolean evaluateFit = failCounter.isOk();
      if (evaluateFit || store.isValid(multiPathResult.getCandidateId())) {
        // if (out != null)
        // {
        // try
        // {
        // out.write(String.format("[%d] %d : %d %b %b\n", multiPathResults.getFrame(),
        // multiPathResult.candidateId, failCount, store.isValid(multiPathResult.candidateId),
        // isNewResult(accept(multiPathResult, true, null))));
        //
        // // todo - Write out the full set of initial and fitted parameters for the results...
        // }
        // catch (Exception e)
        // {
        // try
        // {
        // out.close();
        // }
        // catch (Exception ee)
        // {
        // }
        // finally
        // {
        // out = null;
        // }
        // }
        // }

        // Assess the result if we are below the fail limit or have an estimate
        final SelectedResult selectedResult = select(multiPathResult, true, store);
        boolean newResult = false;
        if (selectedResult != null) {
          int size = 0;
          final PreprocessedPeakResult[] result = selectedResult.results;
          final int[] ok = new int[result.length];
          for (int i = 0; i < ok.length; i++) {
            if (result[i].isNewResult()) {
              newResult = true;
              // if (out != null)
              // {
              // try
              // {
              // out.write(String.format("[%d] %d : %.2f %.2f\n", multiPathResults.getFrame(),
              // multiPathResult.candidateId, result.results[i].getX(),
              // result.results[i].getY()));
              // }
              // catch (Exception e)
              // {
              // try
              // {
              // out.close();
              // }
              // catch (Exception ee)
              // {
              // }
              // finally
              // {
              // out = null;
              // }
              // }
              // }

              // TODO - Check for duplicates
              if (result[i].isNotDuplicate() || !coordinateStore.contains(result[i].getX(),
                  result[i].getY(), result[i].getZ())) {
                coordinateStore.addToQueue(result[i].getX(), result[i].getY(), result[i].getZ());
                ok[size++] = i;
              }
              // else
              // {
              // double[] tmp = coordinateStore.find(result[i].getX(), result[i].getY());
              // System.out.printf("Duplicate [%d] %.2f,%.2f == %.2f,%.2f\n",
              // multiPathResult.candidateId,
              // result[i].getX(), result[i].getY(), tmp[0], tmp[1]);
              // }
            }
          }

          if (size != 0) {
            coordinateStore.flush();

            // This has valid results so add to the output subset only those that are new
            if (size == ok.length) {
              store.add(selectedResult);
            } else {
              final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[size];
              for (int i = 0; i < size; i++) {
                filtered[i] = result[ok[i]];
              }
              store.add(new SelectedResult(filtered, selectedResult.fitResult));
            }
          } else {
            // Add the selected result but with no new results (due to duplicates)
            store.add(new SelectedResult(null, selectedResult.fitResult));
          }
        } else {
          // This failed. Just return the single result
          store.add(new SelectedResult(null, multiPathResult.getSingleFitResult()));
        }
        if (evaluateFit) {
          if (newResult) {
            // More results were accepted so reset the fail count
            failCounter.pass();
          } else {
            // Nothing was accepted, increment fail count
            failCounter.fail();
          }
        }
      } else {
        // This was rejected, increment fail count
        failCounter.fail();
      }

      multiPathResults.complete(c);
    }
    // multiPathResults.end();
    // }

    // if (out != null)
    // {
    // try
    // {
    // out.write(String.format("[%d] %s\n", multiPathResults.getFrame(), "# selected = ?"));
    // }
    // catch (Exception e)
    // {
    // try
    // {
    // out.close();
    // }
    // catch (Exception ee)
    // {
    // }
    // finally
    // {
    // out = null;
    // }
    // }
    // // Close here to only do the first frame
    // if (out != null)
    // {
    // try
    // {
    // out.close();
    // }
    // catch (Exception ee)
    // {
    // }
    // }
    // }
  }

  private boolean isSuitableForDoubletFit(MultiPathFitResult multiPathResult, FitResult fitResult,
      boolean singleQa) {
    // Check there is a fit result
    if (fitResult == null || fitResult.status != 0 || fitResult.getResults() == null) {
      return false;
    }

    // Check if the residuals score is below the configured threshold
    if (residualsThreshold >= 1) {
      return false;
    }

    // Check the other results are OK. Candidates are allowed to fail. New and existing results must
    // pass.
    for (int i = 1; i < validationResults.length; i++) {
      if ((fitResult.getResults()[i].isNewResult() || fitResult.getResults()[i].isExistingResult())
          && validationResults[i] != 0) {
        return false;
      }
    }

    if (validationResults[0] == 0) {
      // The peak was valid so check the residuals
      return ((singleQa) ? multiPathResult.getSingleQaScore()
          : multiPathResult.getMultiQaScore()) > residualsThreshold;
    }

    // Check if it failed due to width
    if (!DirectFilter.anySet(validationResults[0],
        FilterValidationFlag.X_SD_FACTOR | FilterValidationFlag.X_SD_FACTOR)) {
      return false;
    }

    // Get the first spot
    final PreprocessedPeakResult firstResult = fitResult.getResults()[0];

    // Check the width is reasonable given the size of the fitted region.
    //@formatter:off
    if (  firstResult.getXSdFactor() < 1 || // Not a wide spot
        firstResult.getXSd() > multiPathResult.getWidth() || // width covers more than the region
        firstResult.getYSdFactor() < 1 || // Not a wide spot
        firstResult.getYSd() > multiPathResult.getHeight()
      )
     {
      return false;
    //@formatter:on
    }

    // Check the quadrant analysis on the fit residuals
    if (((singleQa) ? multiPathResult.getSingleQaScore()
        : multiPathResult.getMultiQaScore()) < residualsThreshold) {
      return false;
    }

    // We must validate the spot without width filtering. Do not change the min filter.
    setupFilter(FilterValidationOption.NO_WIDTH);

    try {
      if (!filter.accept(firstResult)) {
        // This is still a bad single result, without width filtering
        return false;
      }
    } finally {
      // reset
      restoreFilterState();
    }

    return true;
  }

  /**
   * Replace the fail counter with a default if null. The default allows no failures.
   *
   * @param failCounter the fail counter
   * @return the fail counter
   */
  public static FailCounter replaceIfNull(FailCounter failCounter) {
    return (failCounter == null) ? defaultFailCounter.newCounter() : failCounter;
  }

  /**
   * Check all new and all existing results are valid. Returns the new results.
   *
   * <p>New results and validated candidates that fail the primary filter can be filtered using the
   * minimal filter and sent to the store. The store can be used to determine if a fit for a
   * different candidate has been performed already.
   *
   * @param candidateId the candidate id
   * @param fitResult the results
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @return The new results that pass the filter
   */
  public @Nullable PreprocessedPeakResult[] acceptAll(int candidateId, final FitResult fitResult,
      boolean validateCandidates, SelectedResultStore store) {
    return acceptAll(candidateId, fitResult, validateCandidates, store, false);
  }

  /**
   * Check all new and all existing results are valid. Returns the new results.
   *
   * <p>New results and validated candidates that fail the primary filter can be filtered using the
   * minimal filter and sent to the store. The store can be used to determine if a fit for a
   * different candidate has been performed already.
   *
   * @param candidateId the candidate id
   * @param fitResult the results
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @param precomputed True if this is a subset with pre-computed validation results
   * @return The new results that pass the filter
   */
  public @Nullable PreprocessedPeakResult[] acceptAll(int candidateId, final FitResult fitResult,
      boolean validateCandidates, SelectedResultStore store, boolean precomputed) {
    if (fitResult == null || fitResult.getResults() == null) {
      return null;
    }
    final PreprocessedPeakResult[] results = fitResult.getResults();

    getValidationResults(precomputed, results);

    // All new and existing results should be valid
    int count = 0;
    final int[] ok = new int[results.length];

    // Support for testing using the minimal filter.
    // Note: We do not check the store is not null. This is private method
    // and we send in a null store if necessary.
    final boolean minimalFilter = minFilter != null;

    failExisting = false;
    failNew = false;

    for (int i = 0; i < results.length; i++) {
      if (results[i].isNewResult()) {
        if (results[i].getCandidateId() != candidateId) {
          // This is new result for a different candidate.
          // If a fit has already been accepted (or we don't know)
          // then it should be ignored.
          if (store.isFit(results[i].getCandidateId())) {
            continue;
          }
        }

        if (validationResults[i] == 0) {
          ok[count++] = i;
        } else {
          failNew = true;
          if (minimalFilter) {
            if (minAccept(results[i])) {
              store.passMin(results[i]);
            }
          }
        }
      } else if (results[i].isExistingResult()) {
        if (validationResults[i] != 0) {
          failExisting = true;
        }
      } else if (validateCandidates) {
        // Optionally candidates must pass
        if (validationResults[i] == 0) {
          ok[count++] = i;
        } else if (minimalFilter) {
          if (minAccept(results[i])) {
            store.passMin(results[i]);
          }
        }
      }
    }

    // All new results must pass
    // All existing results must pass
    if (count == 0 || failNew || failExisting) {
      return null;
    }

    // Return the new results
    final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[count];
    for (int i = 0; i < count; i++) {
      filtered[i] = results[ok[i]];
      store.pass(filtered[i]);
    }

    return filtered;
  }

  private void getValidationResults(boolean precomputed, final PreprocessedPeakResult[] results) {
    validationResults = new int[results.length];
    for (int i = 0; i < results.length; i++) {
      if (precomputed && results[i].getValidationResult() != 0) {
        // System.out.println("Ignored");
        validationResults[i] = results[i].getValidationResult();
        continue;
      }
      validationResults[i] = filter.validate(results[i]);
    }
  }

  /**
   * Check any new and all existing results are valid. Returns the new results.
   *
   * <p>New results and validated candidates that fail the primary filter can be filtered using the
   * minimal filter and sent to the store. The store can be used to determine if a fit for a
   * different candidate has been performed already.
   *
   * @param candidateId the candidate id
   * @param fitResult the results
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @return The new results that pass the filter
   */
  public @Nullable PreprocessedPeakResult[] acceptAny(int candidateId, final FitResult fitResult,
      boolean validateCandidates, SelectedResultStore store) {
    return acceptAny(candidateId, fitResult, validateCandidates, store, false);
  }

  /**
   * Check any new and all existing results are valid. Returns the new results.
   *
   * <p>New results and validated candidates that fail the primary filter can be filtered using the
   * minimal filter and sent to the store. The store can be used to determine if a fit for a
   * different candidate has been performed already.
   *
   * @param candidateId the candidate id
   * @param fitResult the results
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @param precomputed True if this is a subset with pre-computed validation results
   * @return The new results that pass the filter
   */
  public @Nullable PreprocessedPeakResult[] acceptAny(int candidateId, final FitResult fitResult,
      boolean validateCandidates, SelectedResultStore store, boolean precomputed) {
    if (fitResult == null || fitResult.getResults() == null) {
      return null;
    }
    final PreprocessedPeakResult[] results = fitResult.getResults();

    // Validate the results
    getValidationResults(precomputed, results);

    return acceptAnyInternal(candidateId, fitResult, validateCandidates, store);
  }

  /**
   * Check any new and all existing results are valid. Returns the new results.
   *
   * <p>New results and validated candidates that fail the primary filter can be filtered using the
   * minimal filter and sent to the store. The store can be used to determine if a fit for a
   * different candidate has been performed already.
   *
   * @param candidateId the candidate id
   * @param fitResult the results
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @return The new results that pass the filter
   */
  private @Nullable PreprocessedPeakResult[] acceptAnyInternal(int candidateId,
      final FitResult fitResult, boolean validateCandidates, SelectedResultStore store) {
    if (fitResult == null || fitResult.getResults() == null) {
      return null;
    }
    final PreprocessedPeakResult[] results = fitResult.getResults();

    // Results are already in the validationResults array

    // Any new and all existing results should be valid
    int count = 0;
    final int[] ok = new int[results.length];

    // Support for testing using the minimal filter
    // Note: We do not check the store is not null. This is private method
    // and we send in a null store if necessary.
    final boolean minimalFilter = minFilter != null;

    failExisting = false;
    failNew = false;

    for (int i = 0; i < results.length; i++) {
      if (results[i].isNewResult()) {
        if (results[i].getCandidateId() != candidateId) {
          // This is new result for a different candidate.
          // If a fit has already been accepted (or we don't know)
          // then it should be ignored.
          if (store.isFit(results[i].getCandidateId())) {
            continue;
          }
        }

        // Any new result that pass are OK
        if (validationResults[i] == 0) {
          ok[count++] = i;
        } else {
          failNew = true;
          if (minimalFilter) {
            if (minAccept(results[i])) {
              store.passMin(results[i]);
            }
          }
        }
      } else if (results[i].isExistingResult()) {
        // All existing results must pass
        if (validationResults[i] != 0) {
          failExisting = true;
        }
      } else if (validateCandidates) {
        // Optionally candidates must pass
        if (validationResults[i] == 0) {
          ok[count++] = i;
        } else if (minimalFilter) {
          if (minAccept(results[i])) {
            store.passMin(results[i]);
          }
        }
      }
    }

    // All existing results must pass
    if (count == 0 || failExisting) {
      return null;
    }

    // Return the new results
    final PreprocessedPeakResult[] filtered = new PreprocessedPeakResult[count];
    for (int i = 0; i < count; i++) {
      filtered[i] = results[ok[i]];
      store.pass(filtered[i]);
    }

    return filtered;
  }

  /**
   * Check any new and all existing results within the multi-doublet fit results are valid. Returns
   * the new results. Coordinate shift filter is disabled for the doublet results.
   *
   * <p>New results and validated candidates that fail the primary filter can be filtered using the
   * minimal filter and sent to the store. The store can be used to determine if a fit for a
   * different candidate has been performed already.
   *
   * @param multiPathResult the multi path result
   * @param validateCandidates Set to true to validate the candidates
   * @param store the store
   * @param candidateId the candidate id
   * @return The new results that pass the filter
   */
  private @Nullable PreprocessedPeakResult[] acceptAnyDoublet(
      final MultiPathFitResult multiPathResult, boolean validateCandidates,
      SelectedResultStore store, final int candidateId) {
    final FitResult multiDoubletFitResult = multiPathResult.getMultiDoubletFitResult();
    if (multiDoubletFitResult == null || multiDoubletFitResult.getResults() == null) {
      return null;
    }

    // Doublets may drift further than single spot candidates.
    // We must validate the doublet spot without shift filtering.
    // Note: Only disable shift for the doublet results.
    // doublets = len(multi-doublet) - len(multi) + 1

    final PreprocessedPeakResult[] results = multiDoubletFitResult.getResults();
    final int nDoublets =
        results.length - multiPathResult.getMultiFitResult().getResults().length + 1;

    setupFilter(FilterValidationOption.NO_SHIFT);

    validationResults = new int[results.length];
    for (int i = 0; i < nDoublets; i++) {
      validationResults[i] = filter.validate(results[i]);
    }

    restoreFilterState();

    for (int i = nDoublets; i < results.length; i++) {
      validationResults[i] = filter.validate(results[i]);
    }

    return acceptAnyInternal(candidateId, multiDoubletFitResult, validateCandidates, store);
  }

  /**
   * Checks if there is the given candidate in the results.
   *
   * @param results the results
   * @param candidateId the candidate id
   * @return true, if there is the given candidate
   */
  private static boolean contains(final PreprocessedPeakResult[] results, final int candidateId) {
    for (int i = 0; i < results.length; i++) {
      if (results[i].getCandidateId() == candidateId) {
        return true;
      }
    }
    return false;
  }

  /**
   * Counts the number of new results in the results.
   *
   * @param results the results
   * @return The count
   */
  private static int countNewResult(final PreprocessedPeakResult[] results) {
    int count = 0;
    if (results != null) {
      for (int i = 0; i < results.length; i++) {
        if (results[i].isNewResult()) {
          count++;
        }
      }
    }
    return count;
  }

  /**
   * Check if the results contain a new result.
   *
   * @param results the results
   * @return True if a new result
   */
  private static boolean isNewResult(final PreprocessedPeakResult[] results) {
    if (results != null) {
      for (int i = 0; i < results.length; i++) {
        if (results[i].isNewResult()) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Filter a set of multi-path results into a set of results.
   *
   * <p>The coordinate store is used to check for duplicates.
   *
   * @param results the results
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param subset True if a subset (the candidate Id will be used to determine the number of failed
   *        fits before the current candidate)
   * @param coordinateStore the coordinate store (can be null)
   * @return the filtered results
   */
  public final PreprocessedPeakResult[] filter(final MultiPathFitResults[] results,
      FailCounter failCounter, boolean subset, CoordinateStore coordinateStore) {
    setup();
    final SimpleSelectedResultStore store = new SimpleSelectedResultStore();

    failCounter = replaceIfNull(failCounter);
    coordinateStore = NullCoordinateStore.replaceIfNull(coordinateStore);

    final ArrayList<PreprocessedPeakResult> list = new ArrayList<>(results.length);
    for (int k = 0; k < results.length; k++) {
      final MultiPathFitResults multiPathResults = results[k];

      // Debugging the results that are scored
      // java.io.OutputStreamWriter out = null;
      // if (debugFilename != null && multiPathResults.getFrame() == 46)
      // {
      // try
      // {
      // out = new java.io.OutputStreamWriter(new java.io.FileOutputStream(debugFilename), "UTF-8");
      // }
      // catch (Exception e)
      // {
      // }
      // }

      failCounter.reset();
      int lastId = -1;
      final int length = multiPathResults.getNumberOfResults();
      store.resize(multiPathResults.getTotalCandidates());
      coordinateStore.clear();
      for (int c = 0; c < length; c++) {
        final MultiPathFitResult multiPathResult = multiPathResults.getResult(c);

        // Include the number of failures before this result from the larger set
        if (subset) {
          incrementFailures(failCounter, lastId, multiPathResult);
          lastId = multiPathResult.getCandidateId();
        }

        final boolean evaluateFit = failCounter.isOk();
        if (evaluateFit || store.isValid(multiPathResult.getCandidateId())) {
          // if (out != null)
          // {
          // try
          // {
          // out.write(String.format("[%d] %d : %d %b %b\n", multiPathResults.frame,
          // multiPathResult.candidateId, failCount, store.isValid(multiPathResult.candidateId),
          // isNewResult(accept(multiPathResult, true, null))));
          // }
          // catch (Exception e)
          // {
          // try
          // {
          // out.close();
          // }
          // catch (Exception ee)
          // {
          // }
          // finally
          // {
          // out = null;
          // }
          // }
          // }

          // Evaluate the result.
          // This allows storing more estimates in the store even if we are past the failures limit.
          final PreprocessedPeakResult[] result = accept(multiPathResult, true, store);
          boolean newResult = false;
          if (result != null) {
            for (int i = 0; i < result.length; i++) {
              if (result[i].isNewResult()) {
                newResult = true;
                // if (out != null)
                // {
                // try
                // {
                // out.write(String.format("[%d] %d : %.2f %.2f\n", multiPathResults.frame,
                // multiPathResult.candidateId, result[i].getX(), result[i].getY()));
                // }
                // catch (Exception e)
                // {
                // try
                // {
                // out.close();
                // }
                // catch (Exception ee)
                // {
                // }
                // finally
                // {
                // out = null;
                // }
                // }
                // }

                // TODO - Check for duplicates
                if (result[i].isNotDuplicate() || !coordinateStore.contains(result[i].getX(),
                    result[i].getY(), result[i].getZ())) {
                  coordinateStore.addToQueue(result[i].getX(), result[i].getY(), result[i].getZ());
                  // if (store.isFit[result[i].getCandidateId()] &&
                  // result[i].getCandidateId() != multiPathResult.candidateId)
                  // System.out.printf("Fitted candidate %d [%d] %f,%f ([%d])\n",
                  // multiPathResults.frame, multiPathResult.candidateId, result[i].getX(),
                  // result[i].getY(), result[i].getCandidateId());

                  list.add(result[i]);

                  // This is a new fitted result
                  store.isFit[result[i].getCandidateId()] = true;
                }
                // else
                // {
                // double[] tmp = cstore.find(result[i].getX(), result[i].getY());
                // System.out.printf("Duplicate %d [%d] %f,%f == %f,%f (%b [%d])\n",
                // multiPathResults.frame, multiPathResult.candidateId, result[i].getX(),
                // result[i].getY(), tmp[0], tmp[1], store.isFit[result[i].getCandidateId()],
                // result[i].getCandidateId());
                // // So we can see them
                // //list.add(result[i]);
                // }
              }
            }
            coordinateStore.flush();
          }
          if (evaluateFit) {
            if (newResult) {
              // More results were accepted so reset the fail count
              failCounter.pass();
            } else {
              // Nothing was accepted, increment fail count
              failCounter.fail();
            }
          }
        } else {
          // This was rejected, increment fail count
          failCounter.fail();
        }

        multiPathResults.complete(c);
      }

      // if (out != null)
      // {
      // try
      // {
      // out.write(String.format("[%d] %s\n", multiPathResults.getFrame(), "# selected = ?"));
      // }
      // catch (Exception e)
      // {
      // try
      // {
      // out.close();
      // }
      // catch (Exception ee)
      // {
      // }
      // finally
      // {
      // out = null;
      // }
      // }
      // // Close here to only do the first frame
      // if (out != null)
      // {
      // try
      // {
      // out.close();
      // }
      // catch (Exception ee)
      // {
      // }
      // }
      // }
    }
    return list.toArray(new PreprocessedPeakResult[list.size()]);
  }

  /**
   * Create a subset of multi-path results, i.e. all those that pass the filter.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>If the subset flag is set to true the candidate Id will be used to determine the number of
   * failed fits before the current candidate, assuming candidates start at zero and increment.
   *
   * @param multiPathResults the multi path results
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param subset True if a subset (the candidate Id will be used to determine the number of failed
   *        fits before the current candidate)
   * @return the filtered results
   */
  public MultiPathFitResult[] filter(final IMultiPathFitResults multiPathResults,
      final FailCounter failCounter, boolean subset) {
    return filter(multiPathResults, replaceIfNull(failCounter), true, subset);
  }

  /**
   * Create a subset of multi-path results, i.e. all those that pass the filter.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>If the subset flag is set to true the candidate Id will be used to determine the number of
   * failed fits before the current candidate, assuming candidates start at zero and increment.
   *
   * @param multiPathResults the multi path results
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param setup Set to true to run the {@link #setup()} method
   * @param subset True if a subset (the candidate Id will be used to determine the number of failed
   *        fits before the current candidate)
   * @return the filtered results
   */
  private @Nullable MultiPathFitResult[] filter(final IMultiPathFitResults multiPathResults,
      final FailCounter failCounter, boolean setup, boolean subset) {
    if (setup) {
      setup();
    }

    failCounter.reset();
    int lastId = -1;
    int size = 0;
    final MultiPathFitResult[] newMultiPathResults =
        new MultiPathFitResult[multiPathResults.getNumberOfResults()];
    final SimpleSelectedResultStore store =
        new SimpleSelectedResultStore(multiPathResults.getTotalCandidates());
    // if (multiPathResults.getFrame() == 12)
    // System.out.println("Debug");
    for (int c = 0; c < newMultiPathResults.length; c++) {
      final MultiPathFitResult multiPathResult = multiPathResults.getResult(c);

      // Include the number of failures before this result from the larger set
      if (subset) {
        incrementFailures(failCounter, lastId, multiPathResult);
        lastId = multiPathResult.getCandidateId();
      }

      final boolean evaluateFit = failCounter.isOk();
      if (evaluateFit || store.isValid(multiPathResult.getCandidateId())) {
        // Evaluate the result.
        // This allows storing more estimates in the store even if we are past the failures limit.
        final PreprocessedPeakResult[] result = accept(multiPathResult, false, store);

        // Note: Even if the actual result failed, the candidate may have passed and so
        // the entire multi-path result should be retained.

        // Also note that depending on the filter, different results can be selected and pushed
        // through
        // the store to set them valid. So we must push everything through the store to ensure
        // nothing
        // is removed that could be used.
        checkIsValid(multiPathResult.getSingleFitResult(), store);
        checkIsValid(multiPathResult.getMultiFitResult(), store);
        setupFilter(FilterValidationOption.NO_SHIFT);
        checkIsValid(multiPathResult.getDoubletFitResult(), store);

        // Fix to only disable shift filtering for the doublet results...
        final FitResult multiDoubletFitResult = multiPathResult.getMultiDoubletFitResult();
        if (multiDoubletFitResult != null && multiDoubletFitResult.getResults() != null) {
          // Note: Only disable shift for the doublet results.
          // doublets = len(multi-doublet) - len(multi) + 1
          final PreprocessedPeakResult[] results = multiDoubletFitResult.getResults();
          final int nDoublets =
              results.length - multiPathResult.getMultiFitResult().getResults().length + 1;
          checkIsValid(results, store, 0, nDoublets);
          restoreFilterState();
          checkIsValid(results, store, nDoublets, results.length);
        } else {
          restoreFilterState();
        }

        // This has valid results so add to the output subset
        newMultiPathResults[size++] = multiPathResult;

        if (evaluateFit) {
          if (isNewResult(result)) {
            // More results were accepted so reset the fail count
            failCounter.pass();
          } else {
            // Nothing was accepted, increment fail count
            failCounter.fail();
          }
        }
      } else {
        // This was rejected, increment fail count
        failCounter.fail();
      }

      multiPathResults.complete(c);
    }

    if (size != 0) {
      return Arrays.copyOf(newMultiPathResults, size);
    }

    return null;
  }

  /**
   * Increment the failures assuming that all the candidates between the current id and the last id
   * failed.
   *
   * @param failCounter the fail counter
   * @param lastId the last id
   * @param multiPathResult the multi path result
   */
  private static void incrementFailures(FailCounter failCounter, int lastId,
      final MultiPathFitResult multiPathResult) {
    final int n = multiPathResult.getCandidateId() - (lastId + 1);
    if (n > 0) {
      failCounter.fail(n);
    }
  }

  /**
   * Create a subset of multi-path results, i.e. all those that pass the filter.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>If the subset flag is set to true the candidate Id will be used to determine the number of
   * failed fits before the current candidate, assuming candidates start at zero and increment.
   *
   * <p>All results are validated with the filter and the result set in the PreprocessedPeakResult.
   * This can be reset using {@link #resetValidationFlag(MultiPathFitResults[])}. This result is
   * used when scoring a subset allowing results to be ignored from duplicate validation.
   *
   * @param results a set of results to analyse
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param subset True if a subset (the candidate Id will be used to determine the number of failed
   *        fits before the current candidate)
   * @return the filtered results
   */
  public MultiPathFitResults[] filterSubset(final MultiPathFitResults[] results,
      FailCounter failCounter, boolean subset) {
    final MultiPathFitResults[] newResults = new MultiPathFitResults[results.length];
    int size = 0;

    setup();
    failCounter = replaceIfNull(failCounter);
    for (int i = 0; i < results.length; i++) {
      final MultiPathFitResult[] newMultiPathResults =
          filter(results[i], failCounter, false, subset);
      if (newMultiPathResults != null) {
        newResults[size++] = new MultiPathFitResults(results[i].getFrame(), newMultiPathResults,
            results[i].getTotalCandidates(), results[i].getNumberOfActualResults());
      }
    }

    return Arrays.copyOf(newResults, size);
  }

  private void checkIsValid(FitResult fitResult, SimpleSelectedResultStore store) {
    if (fitResult == null || fitResult.getResults() == null) {
      return;
    }
    final PreprocessedPeakResult[] results = fitResult.getResults();

    for (int i = 0; i < results.length; i++) {
      // Validate everything
      final int r = filter.validate(results[i]);
      results[i].setValidationResult(r);

      // Mark as valid in the store
      if (r == 0) {
        store.isValid[results[i].getCandidateId()] = true;
      }
    }
  }

  private void checkIsValid(PreprocessedPeakResult[] results, SimpleSelectedResultStore store,
      int lower, int upper) {
    for (int i = lower; i < upper; i++) {
      // Validate everything
      final int r = filter.validate(results[i]);
      results[i].setValidationResult(r);

      // Mark as valid in the store
      if (r == 0) {
        store.isValid[results[i].getCandidateId()] = true;
      }
    }
  }

  /**
   * Score a set of multi-path results. Filter each multi-path result. Any output results that are
   * new results are assumed to be positives and their assignments used to score the results per
   * frame.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>Note: The fractional scores are totalled as well as the integer tp/fp scores. These are
   * returned in the positives and negatives fields of the result.
   *
   * @param results a set of results to analyse
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param n The number of actual results
   * @return the score
   */
  public FractionClassificationResult fractionScore(final MultiPathFitResults[] results,
      final FailCounter failCounter, final int n) {
    return computeFractionScore(results, replaceIfNull(failCounter), n, false, null, null, null);
  }

  /**
   * Score a set of multi-path results. Filter each multi-path result. Any output results that are
   * new results are assumed to be positives and their assignments used to score the results per
   * frame.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>Note: The fractional scores are totalled as well as the integer tp/fp scores. These are
   * returned in the positives and negatives fields of the result.
   *
   * <p>The coordinate store is used to check for duplicates.
   *
   * @param results a set of results to analyse
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param n The number of actual results
   * @param assignments the assignments
   * @param scoreStore the score store
   * @param coordinateStore the coordinate store (can be null)
   * @return the score
   */
  public FractionClassificationResult fractionScore(final MultiPathFitResults[] results,
      final FailCounter failCounter, final int n, List<FractionalAssignment[]> assignments,
      FractionScoreStore scoreStore, CoordinateStore coordinateStore) {
    return computeFractionScore(results, replaceIfNull(failCounter), n, false, assignments,
        scoreStore, coordinateStore);
  }

  /**
   * Score a subset of multi-path results. The subset should be created with
   * {@link #filterSubset(MultiPathFitResults[], FailCounter, boolean)}.
   *
   * <p>Filter each multi-path result. Any output results that are new results are assumed to be
   * positives and their assignments used to score the results per frame.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>Note: The fractional scores are totalled as well as the integer tp/fp scores. These are
   * returned in the positives and negatives fields of the result.
   *
   * @param results a set of results to analyse
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param n The number of actual results
   * @return the score
   */
  public FractionClassificationResult fractionScoreSubset(final MultiPathFitResults[] results,
      final FailCounter failCounter, final int n) {
    return computeFractionScore(results, replaceIfNull(failCounter), n, true, null, null, null);
  }

  /**
   * Score a subset of multi-path results. The subset can be created with
   * {@link #filterSubset(MultiPathFitResults[], FailCounter, boolean)}.
   *
   * <p>Filter each multi-path result. Any output results that are new results are assumed to be
   * positives and their assignments used to score the results per frame.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>Note: The fractional scores are totalled as well as the integer tp/fp scores. These are
   * returned in the positives and negatives fields of the result.
   *
   * <p>The coordinate store is used to check for duplicates.
   *
   * @param results a set of results to analyse
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param n The number of actual results
   * @param assignments the assignments
   * @param scoreStore the score store
   * @param coordinateStore the coordinate store (can be null)
   * @return the score
   */
  public FractionClassificationResult fractionScoreSubset(final MultiPathFitResults[] results,
      final FailCounter failCounter, final int n, List<FractionalAssignment[]> assignments,
      FractionScoreStore scoreStore, CoordinateStore coordinateStore) {
    return computeFractionScore(results, replaceIfNull(failCounter), n, true, assignments,
        scoreStore, coordinateStore);
  }

  /**
   * Score a set of multi-path results.
   *
   * <p>If the subset flag is set to true the candidate Id will be used to determine the number of
   * failed fits before the current candidate, assuming candidates start at zero and increment. The
   * validationResult property of the PreprocessedPeakResult will be used; if non-zero then the peak
   * will not be validated and will be counted as a fail.
   *
   * <p>Filter each multi-path result. Any output results that are new results are assumed to be
   * positives and their assignments used to score the results per frame.
   *
   * <p>The results are processed in order. Results are only processed if the fail counter
   * {@link FailCounter#isOk() isOK}.
   *
   * <p>Note: The fractional scores are totalled as well as the integer tp/fp scores. These are
   * returned in the positives and negatives fields of the result.
   *
   * <p>The coordinate store is used to check for duplicates.
   *
   * @param results a set of results to analyse
   * @param failCounter the counter to track the failures to allow per frame before all peaks are
   *        rejected
   * @param n The number of actual results
   * @param subset True if a subset (the candidate Id will be used to determine the number of failed
   *        fits before the current candidate)
   * @param allAssignments the assignments
   * @param scoreStore the score store
   * @param coordinateStore the coordinate store (can be null)
   * @return the score
   */
  private FractionClassificationResult computeFractionScore(final MultiPathFitResults[] results,
      final FailCounter failCounter, final int n, final boolean subset,
      List<FractionalAssignment[]> allAssignments, FractionScoreStore scoreStore,
      CoordinateStore coordinateStore) {
    final double[] score = new double[4];
    final ArrayList<FractionalAssignment> assignments = new ArrayList<>();

    final SimpleSelectedResultStore store = new SimpleSelectedResultStore();
    if (scoreStore == null) {
      scoreStore = NullFractionScoreStore.INSTANCE;
    }
    coordinateStore = NullCoordinateStore.replaceIfNull(coordinateStore);
    final boolean save = allAssignments != null;

    setup();
    for (int k = 0; k < results.length; k++) {
      final MultiPathFitResults multiPathResults = results[k];

      // Reset fail count for new frames
      failCounter.reset();
      int lastId = -1;
      final int length = multiPathResults.getMultiPathFitResults().length;
      int predicted = 0;
      store.resize(multiPathResults.getTotalCandidates());
      coordinateStore.clear();
      for (int c = 0; c < length; c++) {
        final MultiPathFitResult multiPathResult = multiPathResults.getMultiPathFitResults()[c];

        // Include the number of failures before this result from the larger set
        if (subset) {
          incrementFailures(failCounter, lastId, multiPathResult);
          lastId = multiPathResult.getCandidateId();
        }

        final boolean evaluateFit = failCounter.isOk();
        if (evaluateFit || store.isValid(multiPathResult.getCandidateId())) {
          // Assess the result if we are below the fail limit or have an estimate
          final PreprocessedPeakResult[] result = accept(multiPathResult, true, store, subset);
          boolean newResult = false;
          if (result != null) {
            // For all the results that were returned, check if any are classified results
            // and store the classifications
            for (int i = 0; i < result.length; i++) {
              if (result[i].isNewResult()) {
                newResult = true;

                if (result[i].ignore()) {
                  // Q. should this be passed to the scoreStore?
                } else if (result[i].isNotDuplicate() || !coordinateStore.contains(result[i].getX(),
                    result[i].getY(), result[i].getZ())) {
                  coordinateStore.addToQueue(result[i].getX(), result[i].getY(), result[i].getZ());
                  scoreStore.add(result[i].getUniqueId());
                  final FractionalAssignment[] a = result[i].getAssignments(predicted++);
                  if (a != null && a.length > 0) {
                    // assignments.addAll(Arrays.asList(a));
                    assignments.addAll(new DummyCollection(a));
                  }

                  // This is a new fitted result
                  store.isFit[result[i].getCandidateId()] = true;
                }
              }
            }
            coordinateStore.flush();
          }
          if (evaluateFit) {
            if (newResult) {
              // More results were accepted so reset the fail count
              failCounter.pass();
            } else {
              // Nothing was accepted, increment fail count
              failCounter.fail();
            }
          }
        } else {
          // This was rejected, increment fail count
          failCounter.fail();
        }
      }

      final FractionalAssignment[] tmp =
          score(assignments, score, predicted, save, multiPathResults.getNumberOfActualResults());
      if (allAssignments != null) {
        allAssignments.add(tmp);
      }

      // if (out != null)
      // {
      // try
      // {
      // out.write(String.format("[%d] %s\n", multiPathResults.frame, Arrays.toString(score)));
      // }
      // catch (Exception e)
      // {
      // try
      // {
      // out.close();
      // }
      // catch (Exception ee)
      // {
      // }
      // finally
      // {
      // out = null;
      // }
      // }
      // // Close here to only do the first frame
      // if (out != null)
      // {
      // try
      // {
      // out.close();
      // }
      // catch (Exception ee)
      // {
      // }
      // }
      // }
    }

    // if (out != null)
    // {
    // try
    // {
    // out.close();
    // }
    // catch (Exception ee)
    // {
    // }
    // }

    // Note: We are using the integer positives and negatives fields to actually store integer TP
    // and FP
    return new FractionClassificationResult(score[0], score[1], 0, n - score[0], (int) score[2],
        (int) score[3]);
  }

  /**
   * Score the assignments (TP/FP) and then clear the list.
   *
   * @param assignments The assignments
   * @param score Scores array to accumulate TP/FP scores
   * @param predicted The number of predictions
   * @param save Set to true to save the scored assignments
   * @param actual The number of actual results in the frame
   * @return the fractional assignments
   */
  private static @Nullable FractionalAssignment[] score(
      final ArrayList<FractionalAssignment> assignments, final double[] score, final int predicted,
      boolean save, int actual) {
    if (assignments.isEmpty()) {
      return null;
    }
    final FractionalAssignment[] tmp =
        assignments.toArray(new FractionalAssignment[assignments.size()]);
    final RankedScoreCalculator calc = RankedScoreCalculator.create(tmp, actual, predicted);
    final double[] result = calc.score(predicted, false, save);
    score[0] += result[0];
    score[1] += result[1];
    score[2] += result[2];
    score[3] += result[3];
    assignments.clear();
    return calc.getScoredAssignments();
  }

  /**
   * Reset validation flag.
   *
   * @param results the results
   */
  public static void resetValidationFlag(final MultiPathFitResults[] results) {
    if (results == null) {
      return;
    }
    for (int i = 0; i < results.length; i++) {
      resetValidationFlag(results[i]);
    }
  }

  /**
   * Reset validation flag.
   *
   * @param multiPathFitResults the multi path fit results
   */
  public static void resetValidationFlag(IMultiPathFitResults multiPathFitResults) {
    if (multiPathFitResults == null) {
      return;
    }
    final int size = multiPathFitResults.getNumberOfResults();
    for (int c = 0; c < size; c++) {
      resetValidationFlag(multiPathFitResults.getResult(c));
    }
  }

  /**
   * Reset validation flag.
   *
   * @param multiPathResult the multi path result
   */
  public static void resetValidationFlag(MultiPathFitResult multiPathResult) {
    if (multiPathResult == null) {
      return;
    }
    resetValidationFlag(multiPathResult.getSingleFitResult());
    resetValidationFlag(multiPathResult.getMultiFitResult());
    resetValidationFlag(multiPathResult.getDoubletFitResult());
    resetValidationFlag(multiPathResult.getMultiDoubletFitResult());
  }

  /**
   * Reset validation flag.
   *
   * @param fitResult the fit result
   */
  public static void resetValidationFlag(FitResult fitResult) {
    if (fitResult == null || fitResult.getResults() == null) {
      return;
    }
    final PreprocessedPeakResult[] results = fitResult.getResults();
    for (int i = 0; i < results.length; i++) {
      results[i].setValidationResult(0);
    }
  }

  /**
   * Create a dummy collection that implements toArray() without cloning for the addAll() method in
   * ArrayList.
   */
  private static class DummyCollection extends AbstractList<FractionalAssignment> {

    /** The assignments. */
    final FractionalAssignment[] assignments;

    /**
     * Instantiates a new dummy collection.
     *
     * @param assignments the assignments
     */
    DummyCollection(final FractionalAssignment[] assignments) {
      this.assignments = assignments;
    }

    @Override
    public int size() {
      return assignments.length;
    }

    @Override
    public Object[] toArray() {
      // Return by reference
      return assignments;
    }

    @Override
    public FractionalAssignment get(int index) {
      return assignments[index];
    }
  }

  /**
   * Convert to an XML representation.
   *
   * @return An XML representation of this object
   */
  public String toXml() {
    return FilterXStreamUtils.toXml(this);
  }

  /**
   * Create the filter from the XML representation.
   *
   * @param xml the xml
   * @return the filter (or null)
   */
  public static @Nullable MultiPathFilter fromXml(String xml) {
    try {
      return (MultiPathFilter) FilterXStreamUtils.fromXml(xml);
    } catch (final ClassCastException ex) {
      // Ignore
    }
    return null;
  }
}
