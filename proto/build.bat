@rem Build into the GDSC SMLM source tree
@echo off
set DIR=..\src\main\java
IF EXIST %DIR%\NUL GOTO HAVEDIR
    mkdir %DIR%
:HAVEDIR

@rem protoc --java_out=%DIR% TSFProto.proto
@rem perl suppress.pl %DIR%\gdsc\smlm\tsf\TaggedSpotFile.java --unchecked --unused

protoc --java_out=%DIR% unit_config.proto psf_config.proto calibration_config.proto results_config.proto test_config.proto
perl suppress.pl %DIR%\gdsc\smlm\data\config\PSFConfig.java --unchecked --unused
perl suppress.pl %DIR%\gdsc\smlm\data\config\CalibrationConfig.java --unchecked --unused --deprecation
perl suppress.pl %DIR%\gdsc\smlm\data\config\ResultsConfig.java --unchecked --unused --deprecation
perl suppress.pl %DIR%\gdsc\smlm\data\config\TestConfig.java --unchecked --unused
