#!/bin/sh
# Build into the GDSC SMLM source tree
DIR=../src/main/java
mkdir -p $DIR
#protoc --java_out=$DIR TSFProto.proto
#./suppress.pl $DIR/gdsc/smlm/tsf/TaggedSpotFile.java
protoc --java_out=$DIR GDSC_SMLM.proto
./suppress.pl $DIR/gdsc/smlm/data/config/SMLMSettings.java
