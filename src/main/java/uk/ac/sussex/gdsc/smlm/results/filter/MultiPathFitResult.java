/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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
 * Specifies a the result of fitting a position using different fitting methods.
 *
 * <p>The multi-path result can be evaluated by the MultiPathFilter to determine which result from
 * the different paths should be accepted.
 *
 * <p>This class is used for benchmarking the fitting path options in the PeakFit algorithm.
 */
public class MultiPathFitResult {
  /**
   * The frame containing the result.
   */
  private int frame;

  /**
   * The width of the fit region.
   */
  private int width;

  /**
   * The height of the fit region.
   */
  private int height;

  /**
   * The candidate Id of this result (i.e. the candidate used to identify this position for
   * fitting).
   */
  private int candidateId;

  /**
   * The score from residuals analysis on the residuals of the multi fit. This can be used to choose
   * if the doublet fit should be considered.
   */
  private double multiQaScore = -1;

  /**
   * The score from residuals analysis on the residuals of the single fit. This can be used to
   * choose if the doublet fit should be considered.
   */
  private double singleQaScore = -1;

  /**
   * The results from the multi-fit. It is expected that one result will be true for isNewResult()
   * and zero or more could be true for isExistingResult().
   */
  private FitResult multiFitResult;

  /**
   * The results from the doublet-fit on the multi-fit residuals. It is expected that this should be
   * one or two results that are true for isNewResult().
   */
  private FitResult multiDoubletFitResult;

  /**
   * The results from the single-fit. It is expected that this should be one result that is true for
   * isNewResult().
   */
  private FitResult singleFitResult;

  /**
   * The results from the doublet-fit. It is expected that this should be one or two results that
   * are true for isNewResult().
   */
  private FitResult doubletFitResult;

  /**
   * The fit result.
   */
  public static class FitResult {
    /**
     * Fitting status of the fit. Zero for OK.
     */
    public final int status;

    /**
     * The results from the fit.
     */
    private PreprocessedPeakResult[] results;

    /**
     * Allows storing any data associated with the fit result.
     */
    public final Object data;

    /**
     * Instantiates a new fit result.
     *
     * @param status the status
     */
    public FitResult(int status) {
      this(status, null);
    }

    /**
     * Instantiates a new fit result.
     *
     * @param status the status
     * @param data the data
     */
    public FitResult(int status, Object data) {
      this.status = status;
      this.data = data;
    }

    /**
     * Instantiates a new fit result.
     *
     * @param source the source
     */
    protected FitResult(FitResult source) {
      this.status = source.status;
      this.results = source.results;
      this.data = source.data;
    }

    /**
     * Create a shallow copy.
     *
     * @return the copy
     */
    public FitResult copy() {
      return new FitResult(this);
    }

    /**
     * Gets the status.
     *
     * @return the status
     */
    public int getStatus() {
      return status;
    }

    /**
     * Gets the results results from the fit. It is expected that one or more results will be true
     * for {@link PreprocessedPeakResult#isNewResult()} and zero or more could be true for
     * {@link PreprocessedPeakResult#isExistingResult()}.
     *
     * @return the results
     */
    public PreprocessedPeakResult[] getResults() {
      return results;
    }

    /**
     * Sets the results from the fit. It is expected that one or more results will be true for
     * {@link PreprocessedPeakResult#isNewResult()} and zero or more could be true for
     * {@link PreprocessedPeakResult#isExistingResult()}.
     *
     * @param results the new results
     */
    public void setResults(PreprocessedPeakResult[] results) {
      this.results = results;
    }

    /**
     * Gets the data.
     *
     * @return the data
     */
    public Object getData() {
      return data;
    }
  }

  /**
   * Instantiates a new multi path fit result.
   */
  public MultiPathFitResult() {
    // Do nothing
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   * @param copyResults the copy results
   */
  protected MultiPathFitResult(MultiPathFitResult source, boolean copyResults) {
    candidateId = source.candidateId;
    frame = source.frame;
    width = source.width;
    height = source.height;
    candidateId = source.candidateId;
    multiQaScore = source.multiQaScore;
    singleQaScore = source.singleQaScore;
    if (copyResults) {
      multiFitResult = copyResult(source.multiFitResult);
      multiDoubletFitResult = copyResult(source.multiDoubletFitResult);
      singleFitResult = copyResult(source.singleFitResult);
      doubletFitResult = copyResult(source.doubletFitResult);
    } else {
      multiFitResult = source.multiFitResult;
      multiDoubletFitResult = source.multiDoubletFitResult;
      singleFitResult = source.singleFitResult;
      doubletFitResult = source.doubletFitResult;
    }
  }

  /**
   * Copy the result if not null.
   *
   * @param result the result
   * @return the copy
   */
  private static FitResult copyResult(FitResult result) {
    return (result == null) ? null : result.copy();
  }

  /**
   * Copy the class level field values into a new object. Ignores the fail count fields.
   *
   * @param copyResults Set to true to do a copy of the FitResult objects. Their array objects will
   *        not be copied.
   * @return A copy
   */
  public MultiPathFitResult copy(boolean copyResults) {
    return new MultiPathFitResult(this, copyResults);
  }

  /**
   * Gets the frame containing the result.
   *
   * @return the frame
   */
  public int getFrame() {
    return frame;
  }

  /**
   * Sets the frame containing the result.
   *
   * @param frame the new frame
   */
  public void setFrame(int frame) {
    this.frame = frame;
  }

  /**
   * Gets the width of the fit rectangle.
   *
   * @return the width
   */
  public int getWidth() {
    return width;
  }

  /**
   * Sets the width of the fit rectangle.
   *
   * @param width the new width
   */
  public void setWidth(int width) {
    this.width = width;
  }

  /**
   * Gets the height of the fit rectangle.
   *
   * @return the height
   */
  public int getHeight() {
    return height;
  }

  /**
   * Sets the height of the fit rectangle.
   *
   * @param height the new height
   */
  public void setHeight(int height) {
    this.height = height;
  }

  /**
   * Gets the candidate Id of this result (i.e. the candidate used to identify this position for
   * fitting).
   *
   * @return the candidate id
   */
  public int getCandidateId() {
    return candidateId;
  }

  /**
   * Sets the candidate Id of this result (i.e. the candidate used to identify this position for
   * fitting).
   *
   * @param candidateId the new candidate id
   */
  public void setCandidateId(int candidateId) {
    this.candidateId = candidateId;
  }

  /**
   * Gets the multi fit result.
   *
   * @return the multi fit result
   */
  public FitResult getMultiFitResult() {
    return multiFitResult;
  }

  /**
   * Sets the multi fit result.
   *
   * @param multiFitResult the new multi fit result
   */
  protected void setMultiFitResult(FitResult multiFitResult) {
    this.multiFitResult = multiFitResult;
  }

  /**
   * Gets the multi quadrant analysis (QA) score.
   *
   * @return the multi QA score
   */
  public double getMultiQaScore() {
    return multiQaScore;
  }

  /**
   * Sets the multi quadrant analysis (QA) score.
   *
   * @param multiQaScore the new multi QA score
   */
  protected void setMultiQaScore(double multiQaScore) {
    this.multiQaScore = multiQaScore;
  }

  /**
   * Gets the multi doublet fit result.
   *
   * @return the multi doublet fit result
   */
  public FitResult getMultiDoubletFitResult() {
    return multiDoubletFitResult;
  }

  /**
   * Sets the multi doublet fit result.
   *
   * @param multiDoubletFitResult the new multi doublet fit result
   */
  protected void setMultiDoubletFitResult(FitResult multiDoubletFitResult) {
    this.multiDoubletFitResult = multiDoubletFitResult;
  }

  /**
   * Gets the single fit result.
   *
   * @return the single fit result
   */
  public FitResult getSingleFitResult() {
    return singleFitResult;
  }

  /**
   * Sets the single fit result.
   *
   * @param singleFitResult the new single fit result
   */
  protected void setSingleFitResult(FitResult singleFitResult) {
    this.singleFitResult = singleFitResult;
  }

  /**
   * Gets the single quadrant analysis (QA) score.
   *
   * @return the single QA score
   */
  public double getSingleQaScore() {
    return singleQaScore;
  }

  /**
   * Sets the single quadrant analysis (QA) score.
   *
   * @param singleQaScore the new single QA score
   */
  protected void setSingleQaScore(double singleQaScore) {
    this.singleQaScore = singleQaScore;
  }

  /**
   * Gets the doublet fit result.
   *
   * @return the doublet fit result
   */
  public FitResult getDoubletFitResult() {
    return doubletFitResult;
  }

  /**
   * Sets the doublet fit result.
   *
   * @param doubletFitResult the new doublet fit result
   */
  protected void setDoubletFitResult(FitResult doubletFitResult) {
    this.doubletFitResult = doubletFitResult;
  }
}
