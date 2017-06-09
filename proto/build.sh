#!/bin/sh
# Build into the GDSC SMLM source tree
DIR=../src/main/java
mkdir -p $DIR
protoc --java_out=$DIR TSFProto.proto
./suppress.pl $DIR/gdsc/smlm/tsf/TaggedSpotFile.java
