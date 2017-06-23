@rem Build into the GDSC SMLM source tree
set DIR=..\src\main\java
@rem mkdir -p %DIR%
@rem protoc --java_out=%DIR% TSFProto.proto
@rem perl suppress.pl %DIR%\gdsc\smlm\tsf\TaggedSpotFile.java
protoc --java_out=%DIR% unit.proto psf.proto calibration.proto
perl suppress.pl %DIR%\gdsc\smlm\data\config\PSFConfig.java --unchecked --unused
perl suppress.pl %DIR%\gdsc\smlm\data\config\CalibrationConfig.java --unchecked --unused --deprecation

