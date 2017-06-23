// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: unit.proto

package gdsc.smlm.data.config;

public final class UnitConfig {
  private UnitConfig() {}
  public static void registerAllExtensions(
      com.google.protobuf.ExtensionRegistryLite registry) {
  }

  public static void registerAllExtensions(
      com.google.protobuf.ExtensionRegistry registry) {
    registerAllExtensions(
        (com.google.protobuf.ExtensionRegistryLite) registry);
  }
  /**
   * <pre>
   * Unit for measuring distance
   * </pre>
   *
   * Protobuf enum {@code gdsc.smlm.data.config.DistanceUnit}
   */
  public enum DistanceUnit
      implements com.google.protobuf.ProtocolMessageEnum {
    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>DISTANCE_UNIT_NA = 0;</code>
     */
    DISTANCE_UNIT_NA(0),
    /**
     * <pre>
     * Camera pixel units
     * </pre>
     *
     * <code>PIXEL = 1;</code>
     */
    PIXEL(1),
    /**
     * <pre>
     * Micrometer units
     * </pre>
     *
     * <code>UM = 2;</code>
     */
    UM(2),
    /**
     * <pre>
     * Nanometer units
     * </pre>
     *
     * <code>NM = 3;</code>
     */
    NM(3),
    UNRECOGNIZED(-1),
    ;

    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>DISTANCE_UNIT_NA = 0;</code>
     */
    public static final int DISTANCE_UNIT_NA_VALUE = 0;
    /**
     * <pre>
     * Camera pixel units
     * </pre>
     *
     * <code>PIXEL = 1;</code>
     */
    public static final int PIXEL_VALUE = 1;
    /**
     * <pre>
     * Micrometer units
     * </pre>
     *
     * <code>UM = 2;</code>
     */
    public static final int UM_VALUE = 2;
    /**
     * <pre>
     * Nanometer units
     * </pre>
     *
     * <code>NM = 3;</code>
     */
    public static final int NM_VALUE = 3;


    public final int getNumber() {
      if (this == UNRECOGNIZED) {
        throw new java.lang.IllegalArgumentException(
            "Can't get the number of an unknown enum value.");
      }
      return value;
    }

    /**
     * @deprecated Use {@link #forNumber(int)} instead.
     */
    @java.lang.Deprecated
    public static DistanceUnit valueOf(int value) {
      return forNumber(value);
    }

    public static DistanceUnit forNumber(int value) {
      switch (value) {
        case 0: return DISTANCE_UNIT_NA;
        case 1: return PIXEL;
        case 2: return UM;
        case 3: return NM;
        default: return null;
      }
    }

    public static com.google.protobuf.Internal.EnumLiteMap<DistanceUnit>
        internalGetValueMap() {
      return internalValueMap;
    }
    private static final com.google.protobuf.Internal.EnumLiteMap<
        DistanceUnit> internalValueMap =
          new com.google.protobuf.Internal.EnumLiteMap<DistanceUnit>() {
            public DistanceUnit findValueByNumber(int number) {
              return DistanceUnit.forNumber(number);
            }
          };

    public final com.google.protobuf.Descriptors.EnumValueDescriptor
        getValueDescriptor() {
      return getDescriptor().getValues().get(ordinal());
    }
    public final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptorForType() {
      return getDescriptor();
    }
    public static final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptor() {
      return gdsc.smlm.data.config.UnitConfig.getDescriptor().getEnumTypes().get(0);
    }

    private static final DistanceUnit[] VALUES = values();

    public static DistanceUnit valueOf(
        com.google.protobuf.Descriptors.EnumValueDescriptor desc) {
      if (desc.getType() != getDescriptor()) {
        throw new java.lang.IllegalArgumentException(
          "EnumValueDescriptor is not for this type.");
      }
      if (desc.getIndex() == -1) {
        return UNRECOGNIZED;
      }
      return VALUES[desc.getIndex()];
    }

    private final int value;

    private DistanceUnit(int value) {
      this.value = value;
    }

    // @@protoc_insertion_point(enum_scope:gdsc.smlm.data.config.DistanceUnit)
  }

  /**
   * <pre>
   * Unit for measuring intensity
   * </pre>
   *
   * Protobuf enum {@code gdsc.smlm.data.config.IntensityUnit}
   */
  public enum IntensityUnit
      implements com.google.protobuf.ProtocolMessageEnum {
    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>INTENSITY_UNIT_NA = 0;</code>
     */
    INTENSITY_UNIT_NA(0),
    /**
     * <pre>
     * Photon units
     * </pre>
     *
     * <code>PHOTON = 1;</code>
     */
    PHOTON(1),
    /**
     * <pre>
     * Camera count units
     * </pre>
     *
     * <code>COUNT = 2;</code>
     */
    COUNT(2),
    UNRECOGNIZED(-1),
    ;

    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>INTENSITY_UNIT_NA = 0;</code>
     */
    public static final int INTENSITY_UNIT_NA_VALUE = 0;
    /**
     * <pre>
     * Photon units
     * </pre>
     *
     * <code>PHOTON = 1;</code>
     */
    public static final int PHOTON_VALUE = 1;
    /**
     * <pre>
     * Camera count units
     * </pre>
     *
     * <code>COUNT = 2;</code>
     */
    public static final int COUNT_VALUE = 2;


    public final int getNumber() {
      if (this == UNRECOGNIZED) {
        throw new java.lang.IllegalArgumentException(
            "Can't get the number of an unknown enum value.");
      }
      return value;
    }

    /**
     * @deprecated Use {@link #forNumber(int)} instead.
     */
    @java.lang.Deprecated
    public static IntensityUnit valueOf(int value) {
      return forNumber(value);
    }

    public static IntensityUnit forNumber(int value) {
      switch (value) {
        case 0: return INTENSITY_UNIT_NA;
        case 1: return PHOTON;
        case 2: return COUNT;
        default: return null;
      }
    }

    public static com.google.protobuf.Internal.EnumLiteMap<IntensityUnit>
        internalGetValueMap() {
      return internalValueMap;
    }
    private static final com.google.protobuf.Internal.EnumLiteMap<
        IntensityUnit> internalValueMap =
          new com.google.protobuf.Internal.EnumLiteMap<IntensityUnit>() {
            public IntensityUnit findValueByNumber(int number) {
              return IntensityUnit.forNumber(number);
            }
          };

    public final com.google.protobuf.Descriptors.EnumValueDescriptor
        getValueDescriptor() {
      return getDescriptor().getValues().get(ordinal());
    }
    public final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptorForType() {
      return getDescriptor();
    }
    public static final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptor() {
      return gdsc.smlm.data.config.UnitConfig.getDescriptor().getEnumTypes().get(1);
    }

    private static final IntensityUnit[] VALUES = values();

    public static IntensityUnit valueOf(
        com.google.protobuf.Descriptors.EnumValueDescriptor desc) {
      if (desc.getType() != getDescriptor()) {
        throw new java.lang.IllegalArgumentException(
          "EnumValueDescriptor is not for this type.");
      }
      if (desc.getIndex() == -1) {
        return UNRECOGNIZED;
      }
      return VALUES[desc.getIndex()];
    }

    private final int value;

    private IntensityUnit(int value) {
      this.value = value;
    }

    // @@protoc_insertion_point(enum_scope:gdsc.smlm.data.config.IntensityUnit)
  }

  /**
   * <pre>
   * Unit for measuring angle
   * </pre>
   *
   * Protobuf enum {@code gdsc.smlm.data.config.AngleUnit}
   */
  public enum AngleUnit
      implements com.google.protobuf.ProtocolMessageEnum {
    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>ANGLE_UNIT_NA = 0;</code>
     */
    ANGLE_UNIT_NA(0),
    /**
     * <pre>
     * Radian units
     * </pre>
     *
     * <code>RADIAN = 1;</code>
     */
    RADIAN(1),
    /**
     * <pre>
     * Degree units
     * </pre>
     *
     * <code>DEGREE = 2;</code>
     */
    DEGREE(2),
    UNRECOGNIZED(-1),
    ;

    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>ANGLE_UNIT_NA = 0;</code>
     */
    public static final int ANGLE_UNIT_NA_VALUE = 0;
    /**
     * <pre>
     * Radian units
     * </pre>
     *
     * <code>RADIAN = 1;</code>
     */
    public static final int RADIAN_VALUE = 1;
    /**
     * <pre>
     * Degree units
     * </pre>
     *
     * <code>DEGREE = 2;</code>
     */
    public static final int DEGREE_VALUE = 2;


    public final int getNumber() {
      if (this == UNRECOGNIZED) {
        throw new java.lang.IllegalArgumentException(
            "Can't get the number of an unknown enum value.");
      }
      return value;
    }

    /**
     * @deprecated Use {@link #forNumber(int)} instead.
     */
    @java.lang.Deprecated
    public static AngleUnit valueOf(int value) {
      return forNumber(value);
    }

    public static AngleUnit forNumber(int value) {
      switch (value) {
        case 0: return ANGLE_UNIT_NA;
        case 1: return RADIAN;
        case 2: return DEGREE;
        default: return null;
      }
    }

    public static com.google.protobuf.Internal.EnumLiteMap<AngleUnit>
        internalGetValueMap() {
      return internalValueMap;
    }
    private static final com.google.protobuf.Internal.EnumLiteMap<
        AngleUnit> internalValueMap =
          new com.google.protobuf.Internal.EnumLiteMap<AngleUnit>() {
            public AngleUnit findValueByNumber(int number) {
              return AngleUnit.forNumber(number);
            }
          };

    public final com.google.protobuf.Descriptors.EnumValueDescriptor
        getValueDescriptor() {
      return getDescriptor().getValues().get(ordinal());
    }
    public final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptorForType() {
      return getDescriptor();
    }
    public static final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptor() {
      return gdsc.smlm.data.config.UnitConfig.getDescriptor().getEnumTypes().get(2);
    }

    private static final AngleUnit[] VALUES = values();

    public static AngleUnit valueOf(
        com.google.protobuf.Descriptors.EnumValueDescriptor desc) {
      if (desc.getType() != getDescriptor()) {
        throw new java.lang.IllegalArgumentException(
          "EnumValueDescriptor is not for this type.");
      }
      if (desc.getIndex() == -1) {
        return UNRECOGNIZED;
      }
      return VALUES[desc.getIndex()];
    }

    private final int value;

    private AngleUnit(int value) {
      this.value = value;
    }

    // @@protoc_insertion_point(enum_scope:gdsc.smlm.data.config.AngleUnit)
  }

  /**
   * <pre>
   * Unit for measuring time
   * </pre>
   *
   * Protobuf enum {@code gdsc.smlm.data.config.TimeUnit}
   */
  public enum TimeUnit
      implements com.google.protobuf.ProtocolMessageEnum {
    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>TIME_UNIT_NA = 0;</code>
     */
    TIME_UNIT_NA(0),
    /**
     * <pre>
     * Frame units
     * </pre>
     *
     * <code>FRAME = 1;</code>
     */
    FRAME(1),
    /**
     * <pre>
     * Second units
     * </pre>
     *
     * <code>SECOND = 2;</code>
     */
    SECOND(2),
    /**
     * <pre>
     * Millisecond units
     * </pre>
     *
     * <code>MILLISECOND = 3;</code>
     */
    MILLISECOND(3),
    UNRECOGNIZED(-1),
    ;

    /**
     * <pre>
     * Not available
     * </pre>
     *
     * <code>TIME_UNIT_NA = 0;</code>
     */
    public static final int TIME_UNIT_NA_VALUE = 0;
    /**
     * <pre>
     * Frame units
     * </pre>
     *
     * <code>FRAME = 1;</code>
     */
    public static final int FRAME_VALUE = 1;
    /**
     * <pre>
     * Second units
     * </pre>
     *
     * <code>SECOND = 2;</code>
     */
    public static final int SECOND_VALUE = 2;
    /**
     * <pre>
     * Millisecond units
     * </pre>
     *
     * <code>MILLISECOND = 3;</code>
     */
    public static final int MILLISECOND_VALUE = 3;


    public final int getNumber() {
      if (this == UNRECOGNIZED) {
        throw new java.lang.IllegalArgumentException(
            "Can't get the number of an unknown enum value.");
      }
      return value;
    }

    /**
     * @deprecated Use {@link #forNumber(int)} instead.
     */
    @java.lang.Deprecated
    public static TimeUnit valueOf(int value) {
      return forNumber(value);
    }

    public static TimeUnit forNumber(int value) {
      switch (value) {
        case 0: return TIME_UNIT_NA;
        case 1: return FRAME;
        case 2: return SECOND;
        case 3: return MILLISECOND;
        default: return null;
      }
    }

    public static com.google.protobuf.Internal.EnumLiteMap<TimeUnit>
        internalGetValueMap() {
      return internalValueMap;
    }
    private static final com.google.protobuf.Internal.EnumLiteMap<
        TimeUnit> internalValueMap =
          new com.google.protobuf.Internal.EnumLiteMap<TimeUnit>() {
            public TimeUnit findValueByNumber(int number) {
              return TimeUnit.forNumber(number);
            }
          };

    public final com.google.protobuf.Descriptors.EnumValueDescriptor
        getValueDescriptor() {
      return getDescriptor().getValues().get(ordinal());
    }
    public final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptorForType() {
      return getDescriptor();
    }
    public static final com.google.protobuf.Descriptors.EnumDescriptor
        getDescriptor() {
      return gdsc.smlm.data.config.UnitConfig.getDescriptor().getEnumTypes().get(3);
    }

    private static final TimeUnit[] VALUES = values();

    public static TimeUnit valueOf(
        com.google.protobuf.Descriptors.EnumValueDescriptor desc) {
      if (desc.getType() != getDescriptor()) {
        throw new java.lang.IllegalArgumentException(
          "EnumValueDescriptor is not for this type.");
      }
      if (desc.getIndex() == -1) {
        return UNRECOGNIZED;
      }
      return VALUES[desc.getIndex()];
    }

    private final int value;

    private TimeUnit(int value) {
      this.value = value;
    }

    // @@protoc_insertion_point(enum_scope:gdsc.smlm.data.config.TimeUnit)
  }


  public static com.google.protobuf.Descriptors.FileDescriptor
      getDescriptor() {
    return descriptor;
  }
  private static  com.google.protobuf.Descriptors.FileDescriptor
      descriptor;
  static {
    java.lang.String[] descriptorData = {
      "\n\nunit.proto\022\025gdsc.smlm.data.config*?\n\014D" +
      "istanceUnit\022\024\n\020DISTANCE_UNIT_NA\020\000\022\t\n\005PIX" +
      "EL\020\001\022\006\n\002UM\020\002\022\006\n\002NM\020\003*=\n\rIntensityUnit\022\025\n" +
      "\021INTENSITY_UNIT_NA\020\000\022\n\n\006PHOTON\020\001\022\t\n\005COUN" +
      "T\020\002*6\n\tAngleUnit\022\021\n\rANGLE_UNIT_NA\020\000\022\n\n\006R" +
      "ADIAN\020\001\022\n\n\006DEGREE\020\002*D\n\010TimeUnit\022\020\n\014TIME_" +
      "UNIT_NA\020\000\022\t\n\005FRAME\020\001\022\n\n\006SECOND\020\002\022\017\n\013MILL" +
      "ISECOND\020\003B\014B\nUnitConfigb\006proto3"
    };
    com.google.protobuf.Descriptors.FileDescriptor.InternalDescriptorAssigner assigner =
        new com.google.protobuf.Descriptors.FileDescriptor.    InternalDescriptorAssigner() {
          public com.google.protobuf.ExtensionRegistry assignDescriptors(
              com.google.protobuf.Descriptors.FileDescriptor root) {
            descriptor = root;
            return null;
          }
        };
    com.google.protobuf.Descriptors.FileDescriptor
      .internalBuildGeneratedFileFrom(descriptorData,
        new com.google.protobuf.Descriptors.FileDescriptor[] {
        }, assigner);
  }

  // @@protoc_insertion_point(outer_class_scope)
}
