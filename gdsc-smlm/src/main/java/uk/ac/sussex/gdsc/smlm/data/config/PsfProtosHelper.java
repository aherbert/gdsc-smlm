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

package uk.ac.sussex.gdsc.smlm.data.config;

import java.util.List;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.AstigmatismModel;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFParameterUnit;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;

/**
 * Contains helper functions for the PSFProtos class.
 */
public final class PsfProtosHelper {
  /** The default one-axis Gaussian 2D PSF. */
  public static class DefaultOneAxisGaussian2dPsf {
    /** Default settings instance. */
    public static final PSF INSTANCE;

    static {
      final PSFParameter.Builder paramBuilder = PSFParameter.newBuilder();
      final PSF.Builder builder = PSF.newBuilder();

      builder.setPsfType(PSFType.ONE_AXIS_GAUSSIAN_2D);
      paramBuilder.setName("S");
      paramBuilder.setValue(1);
      paramBuilder.setUnit(PSFParameterUnit.DISTANCE);
      builder.addParameters(paramBuilder.build());
      INSTANCE = builder.build();
    }
  }

  /** The default two-axis Gaussian 2D PSF. */
  public static class DefaultTwoAxisGaussian2dPsf {
    /** Default settings instance. */
    public static final PSF INSTANCE;

    static {
      final PSFParameter.Builder paramBuilder = PSFParameter.newBuilder();
      final PSF.Builder builder = PSF.newBuilder();

      builder.setPsfType(PSFType.TWO_AXIS_GAUSSIAN_2D);
      paramBuilder.setName("Sx");
      paramBuilder.setValue(1);
      builder.addParameters(paramBuilder.build());
      paramBuilder.setName("Sy");
      builder.addParameters(paramBuilder.build());
      INSTANCE = builder.build();
    }
  }

  /** The default two-axis and theta Gaussian 2D PSF. */
  public static class DefaultTwoAxisAndThetaGaussian2dPsf {
    /** Default settings instance. */
    public static final PSF INSTANCE;

    static {
      final PSFParameter.Builder paramBuilder = PSFParameter.newBuilder();
      final PSF.Builder builder = PSF.newBuilder();

      builder.setPsfType(PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
      paramBuilder.setName("Sx");
      paramBuilder.setValue(1);
      builder.addParameters(paramBuilder.build());
      paramBuilder.setName("Sy");
      builder.addParameters(paramBuilder.build());
      paramBuilder.setName("Angle");
      paramBuilder.setUnit(PSFParameterUnit.ANGLE);
      paramBuilder.setValue(0);
      builder.addParameters(paramBuilder.build());
      INSTANCE = builder.build();
    }
  }

  /** No public constructor. */
  private PsfProtosHelper() {}

  /**
   * Gets the name.
   *
   * @param value the value
   * @return the name
   */
  public static String getName(PSFType value) {
    switch (value) {
      case ASTIGMATIC_GAUSSIAN_2D:
        return "Astigmatic Gaussian 2D";
      case CUSTOM:
        return "Custom";
      case ONE_AXIS_GAUSSIAN_2D:
        return "Circular Gaussian 2D";
      case PSF_TYPE_NA:
        return "NA";
      case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
        return "Rotating Elliptical Gaussian 2D";
      case TWO_AXIS_GAUSSIAN_2D:
        return "Elliptical Gaussian 2D";
      case UNRECOGNIZED:
        return ProtosHelperUtils.UNKNOWN;
      default:
        throw new IllegalArgumentException(ProtosHelperUtils.unknownNameMessage(value));
    }
  }

  /**
   * Gets the default PSF.
   *
   * @param value the value
   * @return the default PSF
   */
  public static PSF getDefaultPsf(PSFType value) {
    switch (value) {
      case ONE_AXIS_GAUSSIAN_2D:
        return DefaultOneAxisGaussian2dPsf.INSTANCE;
      case TWO_AXIS_GAUSSIAN_2D:
      case ASTIGMATIC_GAUSSIAN_2D:
        return DefaultTwoAxisGaussian2dPsf.INSTANCE;
      case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
        return DefaultTwoAxisAndThetaGaussian2dPsf.INSTANCE;

      case CUSTOM:
        return PSF.getDefaultInstance();
      case PSF_TYPE_NA:
        return PSF.getDefaultInstance();
      case UNRECOGNIZED:
      default:
        throw new IllegalArgumentException("No default PSF for type: " + value);
    }
  }

  /**
   * Convert the model to the given units.
   *
   * @param model the model
   * @param zDistanceUnit the desired input z distance unit
   * @param widthDistanceUnit the desired output width distance unit
   * @return the astigmatism model
   * @throws ConversionException if the units cannot be converted
   */
  public static AstigmatismModel convert(AstigmatismModel model, DistanceUnit zDistanceUnit,
      DistanceUnit widthDistanceUnit) {
    if (model.getZDistanceUnitValue() == zDistanceUnit.getNumber()
        && model.getSDistanceUnitValue() == widthDistanceUnit.getNumber()) {
      return model;
    }

    final AstigmatismModel.Builder builder = model.toBuilder();
    final TypeConverter<DistanceUnit> zc = UnitConverterUtils
        .createConverter(model.getZDistanceUnit(), zDistanceUnit, model.getNmPerPixel());
    final TypeConverter<DistanceUnit> sc = UnitConverterUtils
        .createConverter(model.getSDistanceUnit(), widthDistanceUnit, model.getNmPerPixel());
    builder.setZDistanceUnitValue(zDistanceUnit.getNumber());
    builder.setSDistanceUnitValue(widthDistanceUnit.getNumber());

    // Convert the input units
    builder.setGamma(zc.convert(model.getGamma()));
    builder.setD(zc.convert(model.getD()));
    builder.setAx(zc.convertBack(model.getAx()));
    builder.setAy(zc.convertBack(model.getAy()));
    builder.setBx(zc.convertBack(zc.convertBack(model.getBx())));
    builder.setBy(zc.convertBack(zc.convertBack(model.getBy())));

    // Convert the output units
    builder.setS0X(sc.convert(model.getS0X()));
    builder.setS0Y(sc.convert(model.getS0Y()));

    return builder.build();
  }

  /**
   * Creates the astigmatic Gaussian 2D PSF from an astigmatism model.
   *
   * <p>Note that the nm/pixel or the units are not contained within the PSF parameters.
   *
   * @param model the model
   * @param zDistanceUnit the desired input z distance unit
   * @param widthDistanceUnit the desired output width distance unit
   * @return the psf
   * @throws ConversionException if the units cannot be converted
   */
  public static PSF createPsf(AstigmatismModel model, DistanceUnit zDistanceUnit,
      DistanceUnit widthDistanceUnit) {
    model = convert(model, zDistanceUnit, widthDistanceUnit);

    final PSF.Builder psf = PSF.newBuilder();
    final PSFParameter.Builder param = PSFParameter.newBuilder();
    psf.setPsfTypeValue(PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE);

    // Add the widths first so the PSF parameters can be converted/used as a 2-axis Gaussian
    addParameter(psf, param, "s0x", model.getS0X(), PSFParameterUnit.DISTANCE);
    addParameter(psf, param, "s0y", model.getS0Y(), PSFParameterUnit.DISTANCE);
    addParameter(psf, param, "gamma", model.getGamma(), PSFParameterUnit.DISTANCE);
    addParameter(psf, param, "d", model.getD(), PSFParameterUnit.DISTANCE);
    // Just use NA for now. The unit is actually 1/distance_unit^2
    addParameter(psf, param, "Ax", model.getAx(), PSFParameterUnit.PSF_PARAMETER_UNIT_NA);
    addParameter(psf, param, "Bx", model.getBx(), PSFParameterUnit.PSF_PARAMETER_UNIT_NA);
    addParameter(psf, param, "Ay", model.getAy(), PSFParameterUnit.PSF_PARAMETER_UNIT_NA);
    addParameter(psf, param, "By", model.getBy(), PSFParameterUnit.PSF_PARAMETER_UNIT_NA);

    return psf.build();
  }

  private static void addParameter(PSF.Builder psf, PSFParameter.Builder param, String name,
      double value, PSFParameterUnit unit) {
    param.setName(name);
    param.setValue(value);
    param.setUnit(unit);
    psf.addParameters(param);
  }

  /**
   * Creates the astigmatism model from the astigmatic Gaussian 2D PSF.
   *
   * <p>Note that the nm/pixel and the units are not contained within the PSF parameters so these
   * must be input as parameters.
   *
   * @param psf the psf
   * @param zDistanceUnit the input z distance unit
   * @param widthDistanceUnit the output width distance unit
   * @param nmPerPixel the nm per pixel
   * @return the model
   * @throws ConfigurationException if the model cannot be created
   */
  public static AstigmatismModel createModel(PSF psf, DistanceUnit zDistanceUnit,
      DistanceUnit widthDistanceUnit, double nmPerPixel) {
    if (psf.getPsfTypeValue() != PSFType.ASTIGMATIC_GAUSSIAN_2D_VALUE) {
      throw new ConfigurationException("Not a " + getName(PSFType.ASTIGMATIC_GAUSSIAN_2D));
    }
    final List<PSFParameter> list = psf.getParametersList();
    final String[] names = {"s0x", "s0y", "gamma", "d", "Ax", "Bx", "Ay", "By"};
    if (list.size() != names.length) {
      throw new ConfigurationException("Invalid number of parameters");
    }
    for (int i = 0; i < names.length; i++) {
      if (!names[i].equals(list.get(i).getName())) {
        throw new ConfigurationException(
            "Invalid parameter name: " + list.get(i).getName() + " != " + names[i]);
      }
    }

    final AstigmatismModel.Builder model = AstigmatismModel.newBuilder();
    model.setSDistanceUnit(widthDistanceUnit);
    model.setZDistanceUnit(zDistanceUnit);
    model.setNmPerPixel(nmPerPixel);

    model.setS0X(list.get(0).getValue());
    model.setS0Y(list.get(1).getValue());
    model.setGamma(list.get(2).getValue());
    model.setD(list.get(3).getValue());
    model.setAx(list.get(4).getValue());
    model.setBx(list.get(5).getValue());
    model.setAy(list.get(6).getValue());
    model.setBy(list.get(7).getValue());

    return model.build();
  }
}
