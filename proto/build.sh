#!/bin/sh
# Build into the GDSC SMLM source tree
DIR=../src/main/java
mkdir -p $DIR

protoc --java_out=$DIR tsf.proto
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/tsf/TSFProtos.java --unchecked --unused --javadoc --static

protoc --java_out=$DIR unit.proto psf.proto calibration.proto results.proto test.proto fit.proto template.proto gui.proto molecule.proto fisher.proto
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/CalibrationProtos.java --unchecked --unused --deprecation --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/FisherProtos.java --unchecked --unused --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/FitProtos.java --unchecked --unused --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/GUIProtos.java --unchecked --unused --deprecation --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/MoleculeProtos.java --unchecked --unused --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/PSFProtos.java --unchecked --unused --deprecation --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/ResultsProtos.java --unchecked --unused --deprecation --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/TemplateProtos.java --unchecked --unused --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/TestProtos.java --unchecked --unused --javadoc --static
perl suppress.pl $DIR/uk/ac/sussex/gdsc/smlm/data/config/UnitProtos.java --unused --javadoc
