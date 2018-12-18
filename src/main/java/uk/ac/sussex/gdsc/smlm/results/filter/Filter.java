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

package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.match.ClassificationResult;
import uk.ac.sussex.gdsc.core.match.FractionClassificationResult;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.ga.Chromosome;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import org.apache.commons.math3.util.FastMath;

import java.util.List;

/**
 * Filter a set of peak results into accepted/rejected.
 */
public abstract class Filter implements Comparable<Filter>, Chromosome<FilterScore>, Cloneable {
  private static final int FP = 0;
  private static final int FN = 1;
  private static final int TP = 2;
  private static final int TN = 3;

  @XStreamOmitField
  private String name;
  @XStreamOmitField
  private String type;
  @XStreamOmitField
  private FilterScore fitness;
  @XStreamOmitField
  private int hash;

  /**
   * Generate the name of the filter using the filter settings (defaults to the first parameter).
   *
   * @return The name of the filter
   */
  protected String generateName() {
    return getParameterName(0) + " " + getParameterValue(0);
  }

  /**
   * Generate the type of the filter using the filter settings (default to the class name with
   * 'Filter' removed).
   *
   * @return The type of the filter
   */
  protected String generateType() {
    return this.getClass().getSimpleName().replaceAll("Filter", "");
  }

  /**
   * Filter the results.
   *
   * @param results the results
   * @return the filtered results
   */
  public MemoryPeakResults filter(MemoryPeakResults results) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.copySettings(results);
    setup(results);
    results.forEach((PeakResultProcedure) peak -> {
      if (accept(peak)) {
        newResults.add(peak);
      }
    });
    end();
    return newResults;
  }

  /**
   * Filter the results.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * @param results the results
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @return the filtered results
   */
  public MemoryPeakResults filter(MemoryPeakResults results, final int failures) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    final FrameCounter counter = new FrameCounter();
    newResults.copySettings(results);
    setup(results);
    results.forEach((PeakResultProcedure) peak -> {
      counter.advanceAndReset(peak.getFrame());

      // Reject all peaks if we have exceeded the fail count
      final boolean isPositive;
      if (counter.getCount() > failures) {
        isPositive = false;
      } else {
        // Otherwise assess the peak
        isPositive = accept(peak);
      }

      if (isPositive) {
        counter.reset();
        newResults.add(peak);
      } else {
        counter.increment();
      }
    });
    end();
    return newResults;
  }

  /**
   * Filter the results.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * <p>Note that this method is to be used to score a set of results that may have been extracted
   * from a larger set since the number of consecutive failures before each peak are expected to be
   * stored in the origY property. Set this to zero and the results should be identical to
   * {@link #filter(MemoryPeakResults, int)}
   *
   * @param results the results
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @return the filtered results
   */
  public MemoryPeakResults filter2(MemoryPeakResults results, final int failures) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    final FrameCounter counter = new FrameCounter();
    newResults.copySettings(results);
    setup(results);
    results.forEach((PeakResultProcedure) peak -> {
      counter.advanceAndReset(peak.getFrame());

      counter.increment(peak.getOrigY());

      // Reject all peaks if we have exceeded the fail count
      final boolean isPositive;
      if (counter.getCount() > failures) {
        isPositive = false;
      } else {
        // Otherwise assess the peak
        isPositive = accept(peak);
      }

      if (isPositive) {
        counter.reset();
        newResults.add(peak);
      } else {
        counter.increment();
      }
    });
    end();
    return newResults;
  }

  /**
   * Filter the results.
   *
   * <p>Input PeakResults must be allocated a score for true positive, false positive, true negative
   * and false negative (accessed via the object property get methods). The filter is run and
   * results that pass accumulate scores for true positive and false positive, otherwise the scores
   * are accumulated for true negative and false negative. The simplest scoring scheme is to mark
   * valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
   *
   * <p>The number of failures before each peak is stored in the origX property of the PeakResult.
   *
   * @param results the results
   * @param score If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
   * @return the filtered results
   */
  public MemoryPeakResults filterSubset(MemoryPeakResults results, double[] score) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    final FrameCounter counter = new FrameCounter();
    newResults.copySettings(results);
    setup(results);
    final double[] s = new double[4];
    final Counter p = new Counter();
    results.forEach((PeakResultProcedure) peak -> {
      counter.advanceAndReset(peak.getFrame());

      // Reject all peaks if we have exceeded the fail count
      final boolean isPositive = accept(peak);

      if (isPositive) {
        peak.setOrigX(counter.getCount());
        counter.reset();
        newResults.add(peak);
      } else {
        counter.increment();
      }

      if (isPositive) {
        p.increment();
        s[TP] += peak.getTruePositiveScore();
        s[FP] += peak.getFalsePositiveScore();
      } else {
        s[FN] += peak.getFalseNegativeScore();
        s[TN] += peak.getTrueNegativeScore();
      }
    });
    end();

    if (score != null && score.length > 5) {
      score[0] = s[TP];
      score[1] = s[FP];
      score[2] = s[TN];
      score[3] = s[FN];
      score[4] = p.getCount();
      score[5] = results.size() - p.getCount();
    }

    return newResults;
  }

  /**
   * Filter the results.
   *
   * <p>Input PeakResults must be allocated a score for true positive, false positive, true negative
   * and false negative (accessed via the object property get methods). The filter is run and
   * results that pass accumulate scores for true positive and false positive, otherwise the scores
   * are accumulated for true negative and false negative. The simplest scoring scheme is to mark
   * valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * <p>The number of failures before each peak is stored in the origX property of the PeakResult.
   *
   * @param results the results
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @param score If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
   * @return the filtered results
   */
  public MemoryPeakResults filterSubset(MemoryPeakResults results, final int failures,
      double[] score) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    final FrameCounter counter = new FrameCounter();
    newResults.copySettings(results);
    setup(results);
    final double[] s = new double[4];
    results.forEach((PeakResultProcedure) peak -> {
      counter.advanceAndReset(peak.getFrame());

      // Reject all peaks if we have exceeded the fail count
      final boolean isPositive;
      if (counter.getCount() > failures) {
        isPositive = false;
      } else {
        // Otherwise assess the peak
        isPositive = accept(peak);
      }

      if (isPositive) {
        peak.setOrigX(counter.getCount());
        counter.reset();
        newResults.add(peak);
      } else {
        counter.increment();
      }

      if (isPositive) {
        s[TP] += peak.getTruePositiveScore();
        s[FP] += peak.getFalsePositiveScore();
      } else {
        s[FN] += peak.getFalseNegativeScore();
        s[TN] += peak.getTrueNegativeScore();
      }
    });
    end();

    if (score != null && score.length > 5) {
      score[0] = s[TP];
      score[1] = s[FP];
      score[2] = s[TN];
      score[3] = s[FN];
      score[4] = newResults.size();
      score[5] = results.size() - newResults.size();
    }

    return newResults;
  }

  /**
   * Filter the results.
   *
   * <p>Input PeakResults must be allocated a score for true positive, false positive, true negative
   * and false negative (accessed via the object property get methods). The filter is run and
   * results that pass accumulate scores for true positive and false positive, otherwise the scores
   * are accumulated for true negative and false negative. The simplest scoring scheme is to mark
   * valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * <p>Note that this method is to be used to score a set of results that may have been extracted
   * from a larger set since the number of consecutive failures before each peak are expected to be
   * stored in the origY property. Set this to zero and the results should be identical to
   * {@link #filterSubset(MemoryPeakResults, double[])}.
   *
   * <p>The number of failures before each peak is stored in the origX property of the PeakResult.
   *
   * @param results the results
   * @param score If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
   * @return the filtered results
   */
  public MemoryPeakResults filterSubset2(MemoryPeakResults results, double[] score) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    final FrameCounter counter = new FrameCounter();
    newResults.copySettings(results);
    setup(results);
    final double[] s = new double[4];
    final Counter p = new Counter();
    results.forEach((PeakResultProcedure) peak -> {
      counter.advanceAndReset(peak.getFrame());

      counter.increment(peak.getOrigY());

      // Reject all peaks if we have exceeded the fail count
      final boolean isPositive = accept(peak);

      if (isPositive) {
        peak.setOrigX(counter.getCount());
        counter.reset();
        newResults.add(peak);
      } else {
        counter.increment();
      }

      if (isPositive) {
        p.increment();
        s[TP] += peak.getTruePositiveScore();
        s[FP] += peak.getFalsePositiveScore();
      } else {
        s[FN] += peak.getFalseNegativeScore();
        s[TN] += peak.getTrueNegativeScore();
      }
    });
    end();

    if (score != null && score.length > 5) {
      score[0] = s[TP];
      score[1] = s[FP];
      score[2] = s[TN];
      score[3] = s[FN];
      score[4] = p.getCount();
      score[5] = results.size() - p.getCount();
    }

    return newResults;
  }

  /**
   * Filter the results.
   *
   * <p>Input PeakResults must be allocated a score for true positive, false positive, true negative
   * and false negative (accessed via the object property get methods). The filter is run and
   * results that pass accumulate scores for true positive and false positive, otherwise the scores
   * are accumulated for true negative and false negative. The simplest scoring scheme is to mark
   * valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * <p>Note that this method is to be used to score a set of results that may have been extracted
   * from a larger set since the number of consecutive failures before each peak are expected to be
   * stored in the origY property. Set this to zero and the results should be identical to
   * {@link #filterSubset(MemoryPeakResults, int, double[])}.
   *
   * <p>The number of failures before each peak is stored in the origX property of the PeakResult.
   *
   * @param results the results
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @param score If not null will be populated with the fraction score [ tp, fp, tn, fn, p, n ]
   * @return the filtered results
   */
  public MemoryPeakResults filterSubset2(MemoryPeakResults results, final int failures,
      double[] score) {
    final MemoryPeakResults newResults = new MemoryPeakResults();
    final FrameCounter counter = new FrameCounter();
    newResults.copySettings(results);
    setup(results);
    final double[] s = new double[4];
    results.forEach((PeakResultProcedure) peak -> {
      counter.advanceAndReset(peak.getFrame());

      counter.increment(peak.getOrigY());

      // Reject all peaks if we have exceeded the fail count
      final boolean isPositive;
      if (counter.getCount() > failures) {
        isPositive = false;
      } else {
        // Otherwise assess the peak
        isPositive = accept(peak);
      }

      if (isPositive) {
        peak.setOrigX(counter.getCount());
        counter.reset();
        newResults.add(peak);
      } else {
        counter.increment();
      }

      if (isPositive) {
        s[TP] += peak.getTruePositiveScore();
        s[FP] += peak.getFalsePositiveScore();
      } else {
        s[FN] += peak.getFalseNegativeScore();
        s[TN] += peak.getTrueNegativeScore();
      }
    });
    end();

    if (score != null && score.length > 5) {
      score[0] = s[TP];
      score[1] = s[FP];
      score[2] = s[TN];
      score[3] = s[FN];
      score[4] = newResults.size();
      score[5] = results.size() - newResults.size();
    }

    return newResults;
  }

  /**
   * Filter the results and return the performance score. Allows benchmarking the filter by marking
   * the results as true or false.
   *
   * <p>Any input PeakResult with an original value that is not zero will be treated as a true
   * result, all other results are false. The filter is run and the results are marked as true
   * positive, false negative and false positive.
   *
   * @param resultsList a list of results to analyse
   * @return the score
   */
  public ClassificationResult score(List<MemoryPeakResults> resultsList) {
    final int[] s = new int[4];
    for (final MemoryPeakResults peakResults : resultsList) {
      setup(peakResults);
      peakResults.forEach(new PeakResultProcedure() {
        @Override
        public void execute(PeakResult peak) {
          final boolean isTrue = peak.getOrigValue() != 0;
          final boolean isPositive = accept(peak);
          if (isTrue) {
            if (isPositive) {
              s[TP]++; // true positive
            } else {
              s[FN]++; // false negative
            }
          } else if (isPositive) {
            s[FP]++; // false positive
          } else {
            s[TN]++; // true negative
          }
        }
      });
      end();
    }
    return new ClassificationResult(s[TP], s[FP], s[TN], s[FN]);
  }

  /**
   * Filter the results and return the performance score. Allows benchmarking the filter by marking
   * the results as true or false.
   *
   * <p>Any input PeakResult with an original value that is not zero will be treated as a true
   * result, all other results are false. The filter is run and the results are marked as true
   * positive, false negative and false positive.
   *
   * @param resultsList a list of results to analyse
   * @param tn The initial true negatives (used when the results have been pre-filtered)
   * @param fn The initial false negatives (used when the results have been pre-filtered)
   * @return the classification result
   */
  public ClassificationResult score(List<MemoryPeakResults> resultsList, int tn, int fn) {
    final int[] s = new int[4];
    s[TN] = tn;
    s[FN] = fn;
    for (final MemoryPeakResults peakResults : resultsList) {
      setup(peakResults);
      peakResults.forEach((PeakResultProcedure) peak -> {
        final boolean isTrue = peak.getOrigValue() != 0;
        final boolean isPositive = accept(peak);
        if (isTrue) {
          if (isPositive) {
            s[TP]++; // true positive
          } else {
            s[FN]++; // false negative
          }
        } else if (isPositive) {
          s[FP]++; // false positive
        } else {
          s[TN]++; // true negative
        }
      });
      end();
    }
    return new ClassificationResult(s[TP], s[FP], s[TN], s[FN]);
  }

  /**
   * Filter the results and return the performance score. Allows benchmarking the filter by marking
   * the results as true or false.
   *
   * <p>Any input PeakResult with an original value that is not zero will be treated as a true
   * result, all other results are false. The filter is run and the results are marked as true
   * positive, false negative and false positive.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * @param resultsList a list of results to analyse
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @return the score
   */
  public ClassificationResult score(List<MemoryPeakResults> resultsList, final int failures) {
    final int[] s = new int[4];
    for (final MemoryPeakResults peakResults : resultsList) {
      setup(peakResults);

      final FrameCounter counter = new FrameCounter();
      peakResults.forEach((PeakResultProcedure) peak -> {
        counter.advanceAndReset(peak.getFrame());

        final boolean isTrue = peak.getOrigValue() != 0;

        // Reject all peaks if we have exceeded the fail count
        final boolean isPositive;
        if (counter.getCount() > failures) {
          isPositive = false;
        } else {
          // Otherwise assess the peak
          isPositive = accept(peak);
        }

        if (isPositive) {
          counter.reset();
        } else {
          counter.increment();
        }

        if (isTrue) {
          if (isPositive) {
            s[TP]++; // true positive
          } else {
            s[FN]++; // false negative
          }
        } else if (isPositive) {
          s[FP]++; // false positive
        } else {
          s[TN]++; // true negative
        }
      });
      end();
    }
    return new ClassificationResult(s[TP], s[FP], s[TN], s[FN]);
  }

  /**
   * Filter the results and return the performance score. Allows benchmarking the filter by marking
   * the results as true or false.
   *
   * <p>Any input PeakResult with an original value that is not zero will be treated as a true
   * result, all other results are false. The filter is run and the results are marked as true
   * positive, false negative and false positive.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * <p>Note that this method is to be used to score a subset that was generated using
   * {@link #filterSubset(MemoryPeakResults, int, double[])} since the number of consecutive
   * failures before each peak are expected to be stored in the origX property.
   *
   * @param resultsList a list of results to analyse
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @param tn The initial true negatives (used when the results have been pre-filtered)
   * @param fn The initial false negatives (used when the results have been pre-filtered)
   * @return the score
   */
  public ClassificationResult scoreSubset(List<MemoryPeakResults> resultsList, final int failures,
      int tn, int fn) {
    final int[] s = new int[4];
    s[TN] = tn;
    s[FN] = fn;
    for (final MemoryPeakResults peakResults : resultsList) {
      setup(peakResults);

      final FrameCounter counter = new FrameCounter();
      peakResults.forEach((PeakResultProcedure) peak -> {
        counter.advanceAndReset(peak.getFrame());

        final boolean isTrue = peak.getOrigValue() != 0;

        counter.increment(peak.getOrigX());

        // Reject all peaks if we have exceeded the fail count
        final boolean isPositive;
        if (counter.getCount() > failures) {
          isPositive = false;
        } else {
          // Otherwise assess the peak
          isPositive = accept(peak);
        }

        if (isPositive) {
          counter.reset();
        } else {
          counter.increment();
        }

        if (isTrue) {
          if (isPositive) {
            s[TP]++; // true positive
          } else {
            s[FN]++; // false negative
          }
        } else if (isPositive) {
          s[FP]++; // false positive
        } else {
          s[TN]++; // true negative
        }
      });
      end();
    }
    return new ClassificationResult(s[TP], s[FP], s[TN], s[FN]);
  }

  /**
   * Filter the results and return the performance score. Allows benchmarking the filter by marking
   * the results as true or false.
   *
   * <p>Input PeakResults must be allocated a score for true positive, false positive, true negative
   * and false negative (accessed via the object property get methods). The filter is run and
   * results that pass accumulate scores for true positive and false positive, otherwise the scores
   * are accumulated for true negative and false negative. The simplest scoring scheme is to mark
   * valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * @param resultsList a list of results to analyse
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @return the score
   */
  public FractionClassificationResult fractionScore(List<MemoryPeakResults> resultsList,
      final int failures) {
    final double[] s = new double[4];
    final Counter p = new Counter();
    int negatives = 0;
    for (final MemoryPeakResults peakResults : resultsList) {
      setup(peakResults);

      final FrameCounter counter = new FrameCounter();
      peakResults.forEach((PeakResultProcedure) peak -> {
        counter.advanceAndReset(peak.getFrame());

        // Reject all peaks if we have exceeded the fail count
        final boolean isPositive;
        if (counter.getCount() > failures) {
          isPositive = false;
        } else {
          // Otherwise assess the peak
          isPositive = accept(peak);
        }

        if (isPositive) {
          counter.reset();
        } else {
          counter.increment();
        }

        if (isPositive) {
          p.increment();
          s[TP] += peak.getTruePositiveScore();
          s[FP] += peak.getFalsePositiveScore();
        } else {
          s[FN] += peak.getFalseNegativeScore();
          s[TN] += peak.getTrueNegativeScore();
        }
      });
      negatives += peakResults.size();
      end();
    }
    negatives -= p.getCount();
    return new FractionClassificationResult(s[TP], s[FP], s[TN], s[FN], p.getCount(), negatives);
  }

  /**
   * Filter the results and return the performance score. Allows benchmarking the filter by marking
   * the results as true or false.
   *
   * <p>Input PeakResults must be allocated a score for true positive, false positive, true negative
   * and false negative (accessed via the object property get methods). The filter is run and
   * results that pass accumulate scores for true positive and false positive, otherwise the scores
   * are accumulated for true negative and false negative. The simplest scoring scheme is to mark
   * valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * <p>Note that this method is to be used to score a set of results that may have been extracted
   * from a larger set since the number of consecutive failures before each peak are expected to be
   * stored in the origY property. Set this to zero and the results should be identical to
   * {@link #fractionScore(List, int)}.
   *
   * @param resultsList a list of results to analyse
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @return the score
   */
  public FractionClassificationResult fractionScore2(List<MemoryPeakResults> resultsList,
      final int failures) {
    final double[] s = new double[4];
    final Counter p = new Counter();
    int negatives = 0;
    for (final MemoryPeakResults peakResults : resultsList) {
      setup(peakResults);

      final FrameCounter counter = new FrameCounter();
      peakResults.forEach((PeakResultProcedure) peak -> {
        counter.advanceAndReset(peak.getFrame());

        counter.increment(peak.getOrigY());

        // Reject all peaks if we have exceeded the fail count
        final boolean isPositive;
        if (counter.getCount() > failures) {
          isPositive = false;
        } else {
          // Otherwise assess the peak
          isPositive = accept(peak);
        }

        if (isPositive) {
          counter.reset();
        } else {
          counter.increment();
        }

        if (isPositive) {
          p.increment();
          s[TP] += peak.getTruePositiveScore();
          s[FP] += peak.getFalsePositiveScore();
        } else {
          s[FN] += peak.getFalseNegativeScore();
          s[TN] += peak.getTrueNegativeScore();
        }
      });
      negatives += peakResults.size();
      end();
    }
    negatives -= p.getCount();
    return new FractionClassificationResult(s[TP], s[FP], s[TN], s[FN], p.getCount(), negatives);
  }

  /**
   * Filter the results and return the performance score. Allows benchmarking the filter by marking
   * the results as true or false.
   *
   * <p>Input PeakResults must be allocated a score for true positive, false positive, true negative
   * and false negative (accessed via the object property get methods). The filter is run and
   * results that pass accumulate scores for true positive and false positive, otherwise the scores
   * are accumulated for true negative and false negative. The simplest scoring scheme is to mark
   * valid results as tp=fn=1 and fp=tn=0 and invalid results the opposite.
   *
   * <p>The number of consecutive rejections are counted per frame. When the configured number of
   * failures is reached all remaining results for the frame are rejected. This assumes the results
   * are ordered by the frame.
   *
   * <p>Note that this method is to be used to score a subset that was generated using
   * {@link #filterSubset(MemoryPeakResults, int, double[])} since the number of consecutive
   * failures before each peak are expected to be stored in the origX property.
   *
   * @param resultsList a list of results to analyse
   * @param failures the number of failures to allow per frame before all peaks are rejected
   * @param tn The initial true negatives (used when the results have been pre-filtered)
   * @param fn The initial false negatives (used when the results have been pre-filtered)
   * @param initialNegatives The initial negatives (used when the results have been pre-filtered)
   * @return the score
   */
  public FractionClassificationResult fractionScoreSubset(List<MemoryPeakResults> resultsList,
      final int failures, double tn, double fn, int initialNegatives) {
    final double[] s = new double[4];
    s[TN] = tn;
    s[FN] = fn;
    final Counter p = new Counter();
    int negatives = initialNegatives;

    for (final MemoryPeakResults peakResults : resultsList) {
      setup(peakResults);

      final FrameCounter counter = new FrameCounter();
      peakResults.forEach((PeakResultProcedure) peak -> {
        counter.advanceAndReset(peak.getFrame());

        counter.increment(peak.getOrigX());

        // Reject all peaks if we have exceeded the fail count
        final boolean isPositive;
        if (counter.getCount() > failures) {
          isPositive = false;
        } else {
          // Otherwise assess the peak
          isPositive = accept(peak);
        }

        if (isPositive) {
          counter.reset();
        } else {
          counter.increment();
        }

        if (isPositive) {
          p.increment();
          s[TP] += peak.getTruePositiveScore();
          s[FP] += peak.getFalsePositiveScore();
        } else {
          s[FN] += peak.getFalseNegativeScore();
          s[TN] += peak.getTrueNegativeScore();
        }
      });
      negatives += peakResults.size();
      end();
    }
    negatives -= p.getCount();
    return new FractionClassificationResult(s[TP], s[FP], s[TN], s[FN], p.getCount(), negatives);
  }

  /**
   * Called before the accept method is called for each peak in the results. Allows pre-processing
   * of the results.
   *
   * @param peakResults the new up
   */
  public abstract void setup(MemoryPeakResults peakResults);

  /**
   * Called for each peak in the results that are filtered.
   *
   * @param peak the peak
   * @return true if the peak should be accepted, otherwise false to reject.
   */
  public abstract boolean accept(PeakResult peak);

  /**
   * Called after the accept method has been called for each peak in the results. Allows memory
   * clean-up of the results.
   */
  public void end() {
    // Do nothing
  }

  /**
   * The numerical value of the filter (defaults to the first parameter).
   *
   * @return The numerical value of the filter. Used for plotting value against performance score.
   */
  public double getNumericalValue() {
    return getParameterValue(0);
  }

  /**
   * The name of the numerical value of the filter (defaults to the first parameter).
   *
   * @return The name of the numerical value of the filter. Used for plotting value against
   *         performance score.
   */
  public String getNumericalValueName() {
    return getParameterName(0);
  }

  /**
   * Gets the name.
   *
   * @return the name (including any parameter values)
   */
  public String getName() {
    if (name == null) {
      name = generateName();
    }
    return name;
  }

  /**
   * Gets the type.
   *
   * @return the type (excluding any parameter values)
   */
  public String getType() {
    if (type == null) {
      type = generateType();
    }
    return type;
  }

  /**
   * Gets the description.
   *
   * @return Describes the functionality of the filter
   */
  public abstract String getDescription();

  /**
   * Set to true if the filter requires parameter deviations within the PeakResult data. The default
   * just uses the result parameters.
   *
   * @return true, if the filter requires parameter deviations
   */
  public boolean requiresParameterDeviations() {
    return false;
  }

  /**
   * To XML.
   *
   * @return An XML representation of this object
   */
  public String toXML() {
    return FilterXStreamUtils.toXML(this);
  }

  /**
   * Create the filter from the XML representation.
   *
   * @param xml the xml
   * @return the filter
   */
  public static Filter fromXML(String xml) {
    try {
      final Filter f = (Filter) FilterXStreamUtils.fromXML(xml);
      f.initialiseState();
      return f;
    } catch (final ClassCastException ex) {
      // ex.printStackTrace();
    }
    return null;
  }

  /**
   * Run after the filter is deserialised using XStream or cloned. Overrride this method if the
   * filter has state that requires resetting.
   */
  protected void initialiseState() {
    // Do nothing
  }

  /**
   * Compare to the other filter, count the number of weakest parameters. If negative then this
   * filter has more weak parameters. If positive then this filter has less weak parameters. If the
   * same or the number of parameters do not match then return 0. If the other filter is null return
   * -1.
   *
   * @param other The other filter
   * @return the count difference
   */
  public int weakest(Filter other) {
    // Null to end of list
    if (other == null) {
      return -1;
    }

    // Use all the parameters
    int size = getNumberOfParameters();
    if (size == other.getNumberOfParameters()) {
      // Extract the parameters
      final double[] p1 = getParameters();
      final double[] p2 = other.getParameters();

      // Find the weakest

      final double[] weakest = p1.clone();
      other.weakestParameters(weakest);
      // Count the number of weakest
      int count = 0;
      while (size-- > 0) {
        if (p1[size] != p2[size]) {
          if (p1[size] == weakest[size]) {
            --count;
          } else {
            ++count;
          }
        }
      }
      return count;
    }

    return 0;
  }

  /**
   * Compare to the other filter, count the number of weakest parameters. If negative then this
   * filter has more weak parameters. If positive then this filter has less weak parameters.
   *
   * <p>This method does not check for null or if the other filter has a different number of
   * parameters.
   *
   * @param other The other filter
   * @return the count difference
   */
  public int weakestUnsafe(Filter other) {
    // Use all the parameters
    int size = getNumberOfParameters();

    // Extract the parameters
    final double[] p1 = getParameters();
    final double[] p2 = other.getParameters();

    // Find the weakest

    final double[] weakest = p1.clone();
    other.weakestParameters(weakest);
    // Count the number of weakest
    int count = 0;
    while (size-- > 0) {
      if (p1[size] != p2[size]) {
        if (p1[size] == weakest[size]) {
          --count;
        } else {
          ++count;
        }
      }
    }
    return count;
  }

  /**
   * Gets the number of parameters.
   *
   * @return The number of parameters for the filter
   */
  public abstract int getNumberOfParameters();

  /**
   * Check index.
   *
   * @param index the index
   */
  protected void checkIndex(final int index) {
    if (index < 0 || index >= getNumberOfParameters()) {
      throw new IndexOutOfBoundsException("Index must be >= 0 and < " + getNumberOfParameters());
    }
  }

  /**
   * Get the parameter value.
   *
   * @param index the index
   * @return The value of the specified parameter
   */
  public double getParameterValue(int index) {
    checkIndex(index);
    return getParameterValueInternal(index);
  }

  /**
   * Get the parameter value. The index should always be between 0 and
   * {@link #getNumberOfParameters()}
   *
   * @param index the index
   * @return The value of the specified parameter
   */
  protected abstract double getParameterValueInternal(int index);

  /**
   * Gets the parameters as an array.
   *
   * @return the parameters
   */
  public double[] getParameters() {
    final int n = getNumberOfParameters();
    final double[] p = new double[n];
    for (int i = 0; i < n; i++) {
      p[i] = getParameterValueInternal(i);
    }
    return p;
  }

  /**
   * Get the recommended minimum amount by which to increment the parameter.
   *
   * @param index the index
   * @return The increment value of the specified parameter
   */
  public abstract double getParameterIncrement(int index);

  /**
   * Return a value to use to disable the parameter.
   *
   * <p>Override this method if zero does not disable the parameter.
   *
   * @param index the index
   * @return The disabled value of the specified parameter
   */
  public double getDisabledParameterValue(int index) {
    checkIndex(index);
    return 0;
  }

  /**
   * Gets the parameter name.
   *
   * @param index the index
   * @return The name of the specified parameter
   */
  public String getParameterName(int index) {
    return getParameterType(index).toString();
  }

  /**
   * Gets the parameter type.
   *
   * @param index the index
   * @return The type of the specified parameter
   */
  public abstract ParameterType getParameterType(int index);

  /**
   * Create a new filter by adjusting the specified parameter.
   *
   * <p>A positive delta will adjust the parameter to be larger. A negative delta will adjust the
   * parameter to be smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
   *
   * <p>Filters can adjust the parameter by a different amount, e.g. by the delta multiplied by a
   * range expected to change the filter performance. This may be relevant in the case where the
   * value is presently zero since no relative change is possible.
   *
   * @param index The parameter index
   * @param delta The amount to adjust the parameter
   * @return The new filter
   */
  public abstract Filter adjustParameter(int index, double delta);

  /**
   * Adjust the specified parameter value.
   *
   * <p>A positive delta will adjust the parameter to be larger. A negative delta will adjust the
   * parameter to be smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
   *
   * @param value the value
   * @param delta the delta
   * @param defaultRange The default range to apply the delta to in the case where the value is zero
   *        and no relative adjustment is possible.
   * @return the double
   */
  protected double updateParameter(double value, double delta, double defaultRange) {
    if (value != 0) {
      return (value + value * delta);
    }
    return (value + defaultRange * delta);
  }

  /**
   * Adjust the specified parameter value.
   *
   * <p>A positive delta will adjust the parameter to be larger. A negative delta will adjust the
   * parameter to be smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
   *
   * @param value the value
   * @param delta the delta
   * @param defaultRange The default range to apply the delta to in the case where the value is zero
   *        and no relative adjustment is possible.
   * @return the float
   */
  protected float updateParameter(float value, double delta, double defaultRange) {
    if (value != 0) {
      return (float) (value + value * delta);
    }
    return (float) (value + defaultRange * delta);
  }

  /**
   * Adjust the specified parameter value.
   *
   * <p>A positive delta will adjust the parameter to be larger. A negative delta will adjust the
   * parameter to be smaller. The adjustment is relative to the parameter value, e.g. 0.1 is 10%.
   * The adjustment is rounded up to the next valid integer to ensure a new parameter value is
   * created.
   *
   * @param value the value
   * @param delta the delta
   * @param defaultRange The default range to apply the delta to in the case where the value is zero
   *        and no relative adjustment is possible.
   * @return the int
   */
  protected int updateParameter(int value, double delta, int defaultRange) {
    final int update;
    if (value != 0) {
      update = (int) Math.ceil(value * Math.abs(delta));
    } else {
      update = (int) Math.ceil(defaultRange * Math.abs(delta));
    }
    if (delta < 0) {
      return value - update;
    }
    return value + update;
  }

  /**
   * Sets the min.
   *
   * @param parameters the parameters
   * @param index the index
   * @param value the value
   */
  protected void setMin(double[] parameters, int index, double value) {
    if (parameters[index] > value) {
      parameters[index] = value;
    }
  }

  /**
   * Sets the max.
   *
   * @param parameters the parameters
   * @param index the index
   * @param value the value
   */
  protected void setMax(double[] parameters, int index, double value) {
    if (parameters[index] < value) {
      parameters[index] = value;
    }
  }

  /**
   * Create a new filter with the specified parameters.
   *
   * @param parameters the parameters
   * @return A new filter
   */
  public abstract Filter create(double... parameters);

  /**
   * Creates a new filter with only the specified parameters enabled.
   *
   * @param enable the enabled flags
   * @return the filter
   */
  public Filter create(boolean[] enable) {
    if (enable == null || enable.length != getNumberOfParameters()) {
      throw new IllegalArgumentException(
          "Enable array must match the number of parameters: " + getNumberOfParameters());
    }
    final double[] p = new double[enable.length];
    for (int i = 0; i < p.length; i++) {
      p[i] = (enable[i]) ? getParameterValueInternal(i) : getDisabledParameterValue(i);
    }
    return create(p);
  }

  /**
   * Update the input array if the Filter's parameters are weaker. This method can be used to find
   * the weakest parameters across a set of filters of the same type. The weakest filter can then be
   * used to create a subset of pre-filtered results to use for testing the filter set.
   *
   * @param parameters The parameters
   */
  public abstract void weakestParameters(double[] parameters);

  /**
   * Compare the two values and return a sort result for the minimum of the two.
   *
   * @param value1 the value 1
   * @param value2 the value 2
   * @return the result (-1 is value1 is lower, 0 is equal, 1 is value2 is lower)
   */
  public static int compareMin(double value1, double value2) {
    if (value1 < value2) {
      return -1;
    }
    if (value1 > value2) {
      return 1;
    }
    return 0;
  }

  /**
   * Compare the two values and return a sort result for the maximum of the two.
   *
   * @param value1 the value 1
   * @param value2 the value 2
   * @return the result (-1 is value1 is higher, 0 is equal, 1 is value2 is higher)
   */
  public static int compareMax(double value1, double value2) {
    if (value1 < value2) {
      return 1;
    }
    if (value1 > value2) {
      return -1;
    }
    return 0;
  }

  /**
   * Some filters requires all the data in a subset for scoring analysis. Others can create a subset
   * using the fail count parameter for a smaller subset that will evaluate faster. This method
   * returns true if the subset can be created using the fail count parameter that will be used to
   * score the subset.
   *
   * @return True if the {@link #filterSubset(MemoryPeakResults, int, double[])} is valid
   */
  public boolean subsetWithFailCount() {
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public int length() {
    // Assume all the parameters are included in the Chromosome
    return getNumberOfParameters();
  }

  /** {@inheritDoc} */
  @Override
  public double[] sequence() {
    // Assume all the parameters are included in the Chromosome
    return getParameters();
  }

  /** {@inheritDoc} */
  @Override
  public Chromosome<FilterScore> newChromosome(double[] sequence) {
    return create(sequence);
  }

  /** {@inheritDoc} */
  @Override
  public double[] lowerLimit() {
    // Set zero as the lower limit
    return new double[length()];
  }

  @Override
  public @Nullable double[] upperLimit() {
    // Default implementation has no upper limit
    return null;
  }

  /** {@inheritDoc} */
  @Override
  public void setFitness(FilterScore fitness) {
    this.fitness = fitness;
  }

  /** {@inheritDoc} */
  @Override
  public FilterScore getFitness() {
    return fitness;
  }

  /**
   * Return the Manhattan (city-block) distance between two chromosomes. This measure is intended to
   * return if the sequences are the same (zero distance) or not). It is not intended for use in
   * distance analysis.
   */
  @Override
  public double distance(Chromosome<FilterScore> other) {
    // NOTE: If the distance is required for a certain type of analysis then this could be done
    // using injection of an interface for calculating the distance.

    final int n = FastMath.min(length(), other.length());
    final double[] s1 = sequence();
    final double[] s2 = other.sequence();
    double distance = 0;
    for (int i = 0; i < n; i++) {
      distance += Math.abs(s1[i] - s2[i]);
    }
    return distance;
  }

  /** {@inheritDoc} */
  @Override
  public boolean equalTo(Chromosome<FilterScore> other) {
    if (length() != other.length()) {
      return false;
    }
    final int n = length();
    final double[] s1 = sequence();
    final double[] s2 = other.sequence();
    for (int i = 0; i < n; i++) {
      if (s1[i] != s2[i]) {
        return false;
      }
    }
    return true;
  }

  /**
   * Get the indices of the parameters that are included in the Chromosome interface. This can be
   * used to look up the name of the parameter using {@link #getParameterName(int)}.
   *
   * @return The indices of the parameters that are included in the Chromosome interface
   */
  public int[] getChromosomeParameters() {
    // Assume all the parameters are included in the Chromosome
    return SimpleArrayUtils.natural(getNumberOfParameters());
  }

  /**
   * Return the value or Float.POSITIVE_INFINITY if value is not positive.
   *
   * @param value the value
   * @return The limit
   */
  public static float getUpperLimit(double value) {
    if (value > 0) {
      return (float) value;
    }
    return Float.POSITIVE_INFINITY;
  }

  /**
   * Return the value squared or Float.POSITIVE_INFINITY if value is not positive.
   *
   * @param value the value
   * @return The squared limit
   */
  public static float getUpperSquaredLimit(double value) {
    if (value > 0) {
      return (float) (value * value);
    }
    return Float.POSITIVE_INFINITY;
  }

  /**
   * Return the value or Double.POSITIVE_INFINITY if value is not positive.
   *
   * @param value the value
   * @return The limit
   */
  public static double getDUpperLimit(double value) {
    if (value > 0) {
      return value;
    }
    return Double.POSITIVE_INFINITY;
  }

  /**
   * Return the value squared or Double.POSITIVE_INFINITY if value is not positive.
   *
   * @param value the value
   * @return The squared limit
   */
  public static double getDUpperSquaredLimit(double value) {
    if (value > 0) {
      return value * value;
    }
    return Double.POSITIVE_INFINITY;
  }

  /**
   * Checks if is finite strictly positive.
   *
   * @param value the value
   * @return true, if is finite strictly positive
   */
  public static boolean isFiniteStrictlyPositive(double value) {
    return value > 0 && value != Double.POSITIVE_INFINITY;
  }

  /**
   * Checks if is finite strictly negative.
   *
   * @param value the value
   * @return true, if is finite strictly negative
   */
  public static boolean isFiniteStrictlyNegative(double value) {
    return value < 0 && value != Double.NEGATIVE_INFINITY;
  }

  /**
   * Gets the filter type.
   *
   * @return The filter type
   */
  public FilterType getFilterType() {
    return FilterType.STANDARD;
  }

  @Override
  public Filter clone() {
    try {
      final Filter f = (Filter) super.clone();
      f.initialiseState();
      return f;
    } catch (final CloneNotSupportedException ex) {
      return null;
    }
  }

  /** {@inheritDoc} */
  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    // Use getClass() and this re-applies to sub-classes
    if (getClass() != obj.getClass()) {
      return false;
    }
    final Filter other = (Filter) obj;
    final int size = getNumberOfParameters();
    if (size != other.getNumberOfParameters()) {
      return false;
    }
    // Check the types are the same before a parameter comparison
    if (!this.getType().equals(other.getType())) {
      return false;
    }
    for (int i = 0; i < size; i++) {
      final double d1 = getParameterValueInternal(i);
      final double d2 = other.getParameterValueInternal(i);
      if (d1 != d2) {
        return false;
      }
    }
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public int hashCode() {
    // Use the cached hash.
    // This assumes that the filter parameter values cannot be modified.
    int result = hash;
    if (result == 0) {
      // Use all the fields used in equals
      result = getType().hashCode();
      final int size = getNumberOfParameters();
      for (int i = 0; i < size; i++) {
        result = result * 31 + Double.hashCode(getParameterValueInternal(i));
      }
      hash = result;
    }
    return result;
  }

  /** {@inheritDoc} */
  @Override
  public int compareTo(Filter obj) {
    // Anything is less than null
    if (obj == null) {
      return -1;
    }

    // Compare using the numerical value first
    final double v1 = getNumericalValue();
    final double v2 = obj.getNumericalValue();
    if (v1 < v2) {
      return -1;
    }
    if (v1 > v2) {
      return 1;
    }

    // Use all the parameters.
    final int size = getNumberOfParameters();
    if (size == obj.getNumberOfParameters()) {
      // Only do this if the same number of parameters
      for (int i = 0; i < size; i++) {
        final double d1 = getParameterValueInternal(i);
        final double d2 = obj.getParameterValueInternal(i);
        if (d1 < d2) {
          return -1;
        }
        if (d1 > d2) {
          return 1;
        }
      }
      return 0;
    }

    // Compare using the number of parameters
    return Integer.compare(size, obj.getNumberOfParameters());
  }
}
