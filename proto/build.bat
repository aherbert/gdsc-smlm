@rem Build into the GDSC SMLM source tree
@echo off
set DIR=..\src\main\java
IF EXIST %DIR%\NUL GOTO HAVEDIR
    mkdir %DIR%
:HAVEDIR

@rem protoc --java_out=%DIR% tsf.proto
@rem perl suppress.pl %DIR%\gdsc\smlm\tsf\TSFProtos.java --unchecked --unused

protoc --java_out=%DIR% unit.proto psf.proto calibration.proto results.proto test.proto fit.proto template.proto gui.proto molecule.proto fisher.proto
perl suppress.pl %DIR%\gdsc\smlm\data\config\PSFProtos.java --unchecked --unused --deprecation
perl suppress.pl %DIR%\gdsc\smlm\data\config\CalibrationProtos.java --unchecked --unused --deprecation
perl suppress.pl %DIR%\gdsc\smlm\data\config\ResultsProtos.java --unchecked --unused --deprecation
perl suppress.pl %DIR%\gdsc\smlm\data\config\TestProtos.java --unchecked --unused
perl suppress.pl %DIR%\gdsc\smlm\data\config\FitProtos.java --unchecked --unused
perl suppress.pl %DIR%\gdsc\smlm\data\config\TemplateProtos.java --unchecked --unused
perl suppress.pl %DIR%\gdsc\smlm\data\config\GUIProtos.java --unchecked --unused --deprecation
perl suppress.pl %DIR%\gdsc\smlm\data\config\MoleculeProtos.java --unchecked --unused
perl suppress.pl %DIR%\gdsc\smlm\data\config\FisherProtos.java --unchecked --unused
