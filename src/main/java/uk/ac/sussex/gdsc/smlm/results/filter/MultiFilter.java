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

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;
import java.util.Arrays;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift, precision and
 * z-depth.
 */
public class MultiFilter extends DirectFilter implements IMultiFilter {

  /** The signal. */
  @XStreamAsAttribute
  final double signal;

  /** The snr. */
  @XStreamAsAttribute
  final float snr;

  /** The min width. */
  @XStreamAsAttribute
  final double minWidth;

  /** The max width. */
  @XStreamAsAttribute
  final double maxWidth;

  /** The shift. */
  @XStreamAsAttribute
  final double shift;

  /** The eshift. */
  @XStreamAsAttribute
  final double eshift;

  /** The precision. */
  @XStreamAsAttribute
  final double precision;

  /** The min Z. */
  @XStreamAsAttribute
  final float minZ;

  /** The max Z. */
  @XStreamAsAttribute
  final float maxZ;

  /** The signal threshold. */
  @XStreamOmitField
  float signalThreshold;

  /** The lower sigma threshold. */
  @XStreamOmitField
  float lowerSigmaThreshold;

  /** The upper sigma threshold. */
  @XStreamOmitField
  float upperSigmaThreshold;

  /** The offsetx. */
  @XStreamOmitField
  float offsetx;

  /** The offsety. */
  @XStreamOmitField
  float offsety;

  /** The eoffset. */
  @XStreamOmitField
  float eoffset;

  /** The variance. */
  @XStreamOmitField
  double variance;

  /** The calculator. */
  @XStreamOmitField
  Gaussian2DPeakResultCalculator calculator;

  /** The width enabled. */
  @XStreamOmitField
  boolean widthEnabled;

  /** The flags. */
  @XStreamOmitField
  int flags;

  /** The z enabled. */
  @XStreamOmitField
  boolean zEnabled;

  /** The components. */
  @XStreamOmitField
  MultiFilterComponentSet components;

  /** The components for no width and shift. */
  @XStreamOmitField
  MultiFilterComponentSet componentsNoWidthShift;
  /** The components for width and shift. */
  @XStreamOmitField
  MultiFilterComponentSet componentsWidthShift;
  /** The components for no width and no shift. */
  @XStreamOmitField
  MultiFilterComponentSet componentsNoWidthNoShift;
  /** The components for width and no shift. */
  @XStreamOmitField
  MultiFilterComponentSet componentsWidthNoShift;
  /**
   * The components for width and no shift copied into a larger array so that the shift component
   * can be set at position 0.
   */
  @XStreamOmitField
  MultiFilterComponentSet componentsShift0;
  @XStreamOmitField
  private int filterSetupFlags;
  @XStreamOmitField
  private FilterSetupData[] filterSetupData;

  /**
   * Instantiates a new multi filter.
   *
   * @param signal the signal
   * @param snr the snr
   * @param minWidth the min width
   * @param maxWidth the max width
   * @param shift the shift
   * @param eshift the eshift
   * @param precision the precision
   * @param minZ the min Z
   * @param maxZ the max Z
   */
  public MultiFilter(double signal, float snr, double minWidth, double maxWidth, double shift,
      double eshift, double precision, float minZ, float maxZ) {
    this.signal = Math.max(0, signal);
    this.snr = Math.max(0, snr);
    // Only swap if max width is enabled
    if (maxWidth != 0 && maxWidth < minWidth) {
      final double f = maxWidth;
      maxWidth = minWidth;
      minWidth = f;
    }
    this.minWidth = Math.max(0, minWidth);
    this.maxWidth = Math.max(0, maxWidth);
    this.shift = Math.max(0, shift);
    this.eshift = Math.max(0, eshift);
    this.precision = Math.max(0, precision);
    this.minZ = minZ;
    this.maxZ = maxZ;
  }

  @Override
  protected String generateName() {
    return String.format(
        "Multi: Signal=%.1f, SNR=%.1f, Width=%.2f-%.2f, Shift=%.2f, EShift=%.2f, "
            + "Precision=%.1f, Width=%.2f-%.2f",
        signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ);
  }

  @Override
  public void setup(MemoryPeakResults peakResults) {
    setupCalculator(peakResults);

    signalThreshold = (float) (signal);

    // Set the width limit
    lowerSigmaThreshold = 0;
    upperSigmaThreshold = Float.POSITIVE_INFINITY;
    // Set the shift limit. The calculator can support both 1/2 axis widths
    // when extracting the Standard Deviation from the parameters.
    final double[] s = PsfHelper.getGaussian2DWxWy(peakResults.getPsf());
    final double s0 =
        (s[0] == s[1]) ? s[0] : Gaussian2DPeakResultHelper.getStandardDeviation(s[0], s[1]);
    lowerSigmaThreshold = (float) (s0 * minWidth);
    upperSigmaThreshold = Filter.getUpperLimit(s0 * maxWidth);
    offsetx = getUpperLimit(s[0] * shift);
    offsety = getUpperLimit(s[1] * shift);
    // Convert to squared distance
    eoffset = Filter.getUpperSquaredLimit(s0 * eshift);

    // Configure the precision limit
    variance = Filter.getDUpperSquaredLimit(precision);

    zEnabled = isZEnabled();
  }

  @Override
  public void setup() {
    filterSetupFlags = 0;
    this.filterSetupData = null;
    setupComponents(true, true, 0);
  }

  @Override
  public void setup(int flags) {
    filterSetupFlags = flags;
    this.filterSetupData = null;
    setupComponents(!areSet(flags, FilterValidationOption.NO_WIDTH),
        !areSet(flags, FilterValidationOption.NO_SHIFT),
        // Pass through the flags that are recognised
        flags & (FilterValidationOption.XY_WIDTH | FilterValidationOption.NO_Z));
  }

  @Override
  public void setup(int flags, FilterSetupData... filterSetupData) {
    setup(flags);
    for (int i = filterSetupData.length; i-- > 0;) {
      if (filterSetupData[i] instanceof ShiftFilterSetupData) {
        this.filterSetupData = getFilterSetupData(filterSetupData[i]);
        final double newShift = ((ShiftFilterSetupData) filterSetupData[i]).shift;
        if (newShift > 0) {
          componentsShift0.replace0(new MultiFilterShiftComponent(newShift));
          components = componentsShift0;
        } else {
          components = componentsWidthNoShift;
        }
        return;
      }
    }
  }

  private boolean isZEnabled() {
    return (minZ != 0 || maxZ != 0) && minZ <= maxZ;
  }

  /**
   * Sets up the calculator.
   *
   * @param peakResults the results
   */
  protected void setupCalculator(MemoryPeakResults peakResults) {
    calculator = Gaussian2DPeakResultHelper.create(peakResults.getPsf(),
        peakResults.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);
  }

  private void setupComponents(final boolean widthEnabled, final boolean shiftEnabled,
      final int flags) {
    // Note: The filter caches the combinations that are likely to be turned on/off:
    // width filtering and shift filtering
    // Other filters related to the PSF (XY widths, z-depth) are assumed to be constant.

    if (componentsWidthShift == null || this.flags != flags) {
      // Store this in case the filter is setup with different flags
      this.flags = flags;

      // Create the components we require
      final MultiFilterComponent[] components1 = new MultiFilterComponent[7];
      int s1 = 0;
      Class<? extends MultiFilterComponent> widthComponentClass = null;
      Class<? extends MultiFilterComponent> shiftComponentClass = null;

      // Current order of filter power obtained from BenchmarkFilterAnalysis:
      // SNR, Max Width, Precision, Shift, Min width
      if (isFiniteStrictlyPositive(snr)) {
        components1[s1++] = new MultiFilterSnrComponent(snr);
      }
      if ((maxWidth > 1 && maxWidth != Double.POSITIVE_INFINITY)
          || (minWidth > 0 && minWidth < 1)) {
        // Handle the width being 1/2 axis variable.
        if (areSet(flags, FilterValidationOption.XY_WIDTH)) {
          components1[s1++] = new MultiFilterXyWidthComponent(minWidth, maxWidth);
        } else {
          components1[s1++] = new MultiFilterWidthComponent(minWidth, maxWidth);
        }
        widthComponentClass = components1[s1 - 1].getClass();
      }
      if (isFiniteStrictlyPositive(precision)) {
        components1[s1++] = createPrecisionComponent();
      }
      if (isFiniteStrictlyPositive(shift)) {
        components1[s1++] = new MultiFilterShiftComponent(shift);
        shiftComponentClass = components1[s1 - 1].getClass();
      }
      if (isFiniteStrictlyPositive(signal)) {
        components1[s1++] = new MultiFilterSignalComponent(signal);
      }
      if (isFiniteStrictlyPositive(eshift)) {
        components1[s1++] = new MultiFilterEShiftComponent(eshift);
      }
      if (isZEnabled() && !areSet(flags, FilterValidationOption.NO_Z)) {
        components1[s1++] = new MultiFilterZComponent(minZ, maxZ);
      }

      final MultiFilterComponent[] components2 =
          MultiFilter.remove(components1, s1, widthComponentClass);
      final MultiFilterComponent[] components3 =
          MultiFilter.remove(components1, s1, shiftComponentClass);

      final MultiFilterComponent[] components4 =
          MultiFilter.remove(components2, components2.length, shiftComponentClass);

      componentsWidthShift = MultiFilterComponentSetUtils.create(components1, s1);
      componentsNoWidthShift = MultiFilterComponentSetUtils.create(components2, components2.length);
      componentsWidthNoShift = MultiFilterComponentSetUtils.create(components3, components3.length);
      componentsNoWidthNoShift =
          MultiFilterComponentSetUtils.create(components4, components4.length);

      final MultiFilterComponent[] data = new MultiFilterComponent[components3.length + 1];
      System.arraycopy(components3, 0, data, 1, components3.length);
      componentsShift0 = MultiFilterComponentSetUtils.create(data, data.length);
    }

    if (widthEnabled) {
      components = (shiftEnabled) ? componentsWidthShift : componentsNoWidthShift;
    } else {
      components = (shiftEnabled) ? componentsNoWidthShift : componentsNoWidthNoShift;
    }
  }

  /**
   * Creates the precision component.
   *
   * @return the precision component
   */
  protected MultiFilterComponent createPrecisionComponent() {
    return new MultiFilterVarianceComponent(precision);
  }

  /**
   * Removes the component from the input array using the class.
   *
   * @param in the input
   * @param size the size
   * @param clazz the clazz
   * @return the new array
   */
  static MultiFilterComponent[] remove(MultiFilterComponent[] in, int size,
      Class<? extends MultiFilterComponent> clazz) {
    if (clazz == null) {
      return Arrays.copyOf(in, size);
    }
    final MultiFilterComponent[] out = new MultiFilterComponent[size];
    int length = 0;
    for (int i = 0; i < size; i++) {
      if (in[i].getClass() != clazz) {
        out[length++] = in[i];
      }
    }
    return Arrays.copyOf(out, length);
  }

  @Override
  public int getFilterSetupFlags() {
    return filterSetupFlags;
  }

  @Override
  public FilterSetupData[] getFilterSetupData() {
    return filterSetupData;
  }

  @Override
  public boolean accept(PeakResult peak) {
    // Current order of filter power obtained from BenchmarkFilterAnalysis:
    // SNR, Max Width, Precision, Shift, Min width

    if (peak.getSnr() < this.snr) {
      return false;
    }

    // Width
    final float sd = calculator.getStandardDeviation(peak.getParameters());
    if (sd > upperSigmaThreshold || sd < lowerSigmaThreshold) {
      return false;
    }

    // Precision
    if (getVariance(peak) > variance) {
      return false;
    }

    // Shift
    if (Math.abs(peak.getXShift()) > offsetx || Math.abs(peak.getYShift()) > offsety) {
      return false;
    }

    if (peak.getIntensity() < signalThreshold) {
      return false;
    }

    if (zEnabled && (peak.getZPosition() < minZ || peak.getZPosition() > maxZ)) {
      return false;
    }

    // Euclidian shift
    final float dx = peak.getXPosition();
    final float dy = peak.getYPosition();
    return (dx * dx + dy * dy <= eoffset);
  }

  /**
   * Gets the variance.
   *
   * @param peak the peak
   * @return the variance
   */
  protected double getVariance(PeakResult peak) {
    return calculator.getLseVariance(peak.getParameters(), peak.getNoise());
  }

  @Override
  public int getValidationFlags() {
    if (components == null) {
      setup();
    }
    return components.getValidationFlags();
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    return components.validate(peak);
  }

  @Override
  public double getNumericalValue() {
    // This is not the first parameter so override
    return snr;
  }

  @Override
  public String getNumericalValueName() {
    // This is not the first parameter so override
    return ParameterType.SNR.toString();
  }

  @Override
  public String getDescription() {
    return "Filter results using multiple thresholds: Signal, SNR, width, shift, "
        + "Euclidian shift, precision and Z-depth";
  }

  @Override
  public int getNumberOfParameters() {
    return 9;
  }

  @Override
  protected double getParameterValueInternal(int index) {
    switch (index) {
      case 0:
        return signal;
      case 1:
        return snr;
      case 2:
        return minWidth;
      case 3:
        return maxWidth;
      case 4:
        return shift;
      case 5:
        return eshift;
      case 6:
        return precision;
      case 7:
        return minZ;
      default:
        return maxZ;
    }
  }

  @Override
  public double[] getParameters() {
    return new double[] {signal, snr, minWidth, maxWidth, shift, eshift, precision, minZ, maxZ};
  }

  @Override
  public double getParameterIncrement(int index) {
    checkIndex(index);
    switch (index) {
      case 0:
        return SignalFilter.DEFAULT_INCREMENT;
      case 1:
        return SnrFilter.DEFAULT_INCREMENT;
      case 2:
        return WidthFilter2.DEFAULT_MIN_INCREMENT;
      case 3:
        return WidthFilter.DEFAULT_INCREMENT;
      case 4:
        return ShiftFilter.DEFAULT_INCREMENT;
      case 5:
        return EShiftFilter.DEFAULT_INCREMENT;
      case 6:
        return PrecisionFilter.DEFAULT_INCREMENT;
      case 7:
        return ZCoordinateFilter.DEFAULT_INCREMENT;
      default:
        return ZCoordinateFilter.DEFAULT_INCREMENT;
    }
  }

  @Override
  public ParameterType getParameterType(int index) {
    checkIndex(index);
    switch (index) {
      case 0:
        return ParameterType.SIGNAL;
      case 1:
        return ParameterType.SNR;
      case 2:
        return ParameterType.MIN_WIDTH;
      case 3:
        return ParameterType.MAX_WIDTH;
      case 4:
        return ParameterType.SHIFT;
      case 5:
        return ParameterType.ESHIFT;
      case 6:
        return getPrecisionParamaterType();
      case 7:
        return ParameterType.MIN_Z;
      default:
        return ParameterType.MAX_Z;
    }
  }

  /**
   * Gets the precision paramater type.
   *
   * @return the precision paramater type
   */
  protected ParameterType getPrecisionParamaterType() {
    return ParameterType.PRECISION;
  }

  /** The default range. */
  static final double[] defaultRange = {SignalFilter.DEFAULT_RANGE, SnrFilter.DEFAULT_RANGE,
      WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE,
      EShiftFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE, ZCoordinateFilter.DEFAULT_RANGE,
      ZCoordinateFilter.DEFAULT_RANGE};

  @Override
  public Filter adjustParameter(int index, double delta) {
    checkIndex(index);
    final double[] params = getParameters();
    params[index] = updateParameter(params[index], delta, defaultRange[index]);
    return create(params);
  }

  @Override
  public Filter create(double... parameters) {
    return new MultiFilter(parameters[0], (float) parameters[1], parameters[2], parameters[3],
        parameters[4], parameters[5], parameters[6], (float) parameters[7], (float) parameters[8]);
  }

  @Override
  public void weakestParameters(double[] parameters) {
    setMin(parameters, 0, signal);
    setMin(parameters, 1, snr);
    setMin(parameters, 2, minWidth);
    setMax(parameters, 3, maxWidth);
    setMax(parameters, 4, shift);
    setMax(parameters, 5, eshift);
    setMax(parameters, 6, precision);
    setMin(parameters, 7, minZ);
    setMax(parameters, 8, maxZ);
  }

  @Override
  public int lowerBoundOrientation(int index) {
    return (index < 3 || index == 7) ? -1 : 1;
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
  public int weakest(MultiFilter other) {
    if (other == null) {
      return -1;
    }

    // Count the number of weakest
    //@formatter:off
    return compareMin(signal, other.signal)
        + compareMin(snr, other.snr)
        + compareMin(minWidth, other.minWidth)
        + compareMax(maxWidth, other.maxWidth)
        + compareMax(shift, other.shift)
        + compareMax(eshift, other.eshift)
        + compareMax(precision, other.precision)
        + compareMin(minZ, other.minZ)
        + compareMax(maxZ, other.maxZ);
    //@formatter:on
  }

  @Override
  public double[] upperLimit() {
    return new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
        WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT, ShiftFilter.UPPER_LIMIT,
        EShiftFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT, ZCoordinateFilter.UPPER_LIMIT,
        ZCoordinateFilter.UPPER_LIMIT};
  }

  @Override
  public double[] mutationStepRange() {
    return defaultRange;
  }

  @Override
  public double getSignal() {
    return signal;
  }

  @Override
  public double getSnr() {
    return snr;
  }

  @Override
  public double getMinWidth() {
    return minWidth;
  }

  @Override
  public double getMaxWidth() {
    return maxWidth;
  }

  @Override
  public double getShift() {
    return shift;
  }

  @Override
  public double getEShift() {
    return eshift;
  }

  @Override
  public double getPrecision() {
    return precision;
  }

  @Override
  public PrecisionType getPrecisionType() {
    return PrecisionType.ESTIMATE;
  }

  @Override
  public double getMinZ() {
    return minZ;
  }

  @Override
  public double getMaxZ() {
    return maxZ;
  }

  @Override
  protected void initialiseState() {
    // This is run after a clone() occurs.
    // Q. Can the setup state be maintained?

    // Replace any object that is manipulated by the instance
    if (componentsShift0 != null) {
      final boolean update = components == componentsShift0;
      componentsShift0 = componentsShift0.copy();
      if (update) {
        components = componentsShift0;
      }
    }
  }
}
