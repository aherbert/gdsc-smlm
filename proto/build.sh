#!/bin/sh
# Build into the GDSC SMLM source tree
DIR=../src/main/java
mkdir -p $DIR
#protoc --java_out=$DIR TSFProto.proto
#perl suppress.pl $DIR/gdsc/smlm/tsf/TaggedSpotFile.java
protoc --java_out=$DIR unit.proto psf.proto calibration.proto
perl suppress.pl $DIR/gdsc/smlm/data/config/PSFConfig.java --unchecked --unused
perl suppress.pl $DIR/gdsc/smlm/data/config/CalibrationConfig.java --unchecked --unused --deprecation
