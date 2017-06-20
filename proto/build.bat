@rem Build into the GDSC SMLM source tree
set DIR=..\src\main\java
@rem mkdir -p %DIR%
protoc --java_out=%DIR% TSFProto.proto
perl suppress.pl %DIR%\gdsc\smlm\tsf\TaggedSpotFile.java
@rem protoc --java_out=%DIR% GDSC_SMLM.proto
@rem perl suppress.pl %DIR%\gdsc\smlm\data\config\SMLMSettings.java
