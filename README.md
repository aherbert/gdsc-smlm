GDSC Single Molecule Light Microscopy (SMLM) ImageJ Plugins
===========================================================

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://github.com/aherbert/gdsc-smlm/actions/workflows/build.yml/badge.svg)](https://github.com/aherbert/gdsc-smlm/actions/workflows/build.yml)
[![Coverage Status](https://codecov.io/gh/aherbert/gdsc-smlm/branch/master/graph/badge.svg)](https://app.codecov.io/gh/aherbert/gdsc-smlm)
[![Documentation Status](https://readthedocs.org/projects/gdsc-smlm/badge/?version=latest)](https://gdsc-smlm.readthedocs.io/en/latest/?badge=latest)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/uk.ac.sussex.gdsc/gdsc-smlm/badge.svg)](https://maven-badges.herokuapp.com/maven-central/uk.ac.sussex.gdsc/gdsc-smlm/)
[![Javadocs](https://javadoc.io/badge2/uk.ac.sussex.gdsc/gdsc-smlm/javadoc.svg)](https://javadoc.io/doc/uk.ac.sussex.gdsc/gdsc-smlm)

The GDSC Single Molecule Light Microscopy (SMLM) plugins provide various tools
for single molecule localisation analysis. This includes PALM, STORM and other
single molecule microscopy methods.

The plugins provide tools to:

- Fit an image, or series of images, using a 2D Guassian Point Spread Function
(PSF)
- Save results to a table, a file, an image and/or to memory
- Trace or cluster localisations through time to identify molecules
- Correct drift in long time course images
- Estimate fluorophore dark-time and blinking rate
- Create localisation density images
- Create custom PSFs from calibration images
- Create simulated data using a Gaussian or custom PSF with configurable
molecule populations and diffusion
- Calibrate the gain and read noise of your microscope camera
- Estimate noise in an image
- Estimate resolution using Fourier Ring Correlation

Results can be loaded from file for analysis using the following formats:

- SMLM Text/Binary
- RapidSTORM
- Nikon DSTORM

The code is split into modules:

- gdsc-smlm: Contains utilities for SMLM data analysis
- gdsc-smlm-ij: Contains SMLM plugins for ImageJ


Citation
--------

Etheridge TJ, Carr AM and Herbert AD.  
**GDSC SMLM: Single-molecule localisation microscopy software for ImageJ** \[version 1; peer review: 1 approved\].  
*Wellcome Open Res* 2022, **7:241**  
[doi.org/10.12688/wellcomeopenres.18327.1](https://doi.org/10.12688/wellcomeopenres.18327.1)


Documentation
-------------

See the latest documentation on [ReadTheDocs](https://gdsc-smlm.readthedocs.io).


Install
-------

The SMLM plugins are distributed using an ImageJ2/Fiji update site.

To install the plugins using Fiji (an ImageJ distribution) just follow the
instructions [How_to_follow_a_3rd_party_update_site](http://fiji.sc/How_to_follow_a_3rd_party_update_site)
and add the GDSC SMLM update site. All the plugins will appear under the 'Plugins > GDSC SMLM' menu.


Installation from source
------------------------

The source code is accessed using git and built using Maven.

The code depends on the gdsc-test, gdsc-ij-parent and gdsc-core artifacts so
you will have to install these to your local Maven repository before building:

1. Clone the required repositories

        git clone https://github.com/aherbert/gdsc-test.git
        git clone https://github.com/aherbert/gdsc-ij-parent.git
        git clone https://github.com/aherbert/gdsc-core.git
        git clone https://github.com/aherbert/gdsc-smlm.git

1. Build the code and install using Maven

        cd ../gdsc-test
        mvn install
        cd ../gdsc-ij-parent
        mvn install
        cd ../gdsc-core
        mvn install
        cd ../gdsc-smlm
        mvn package

   This will produce a gdsc-smlm-ij_-[VERSION].jar file in the target directory
   of the gdsc-smlm-ij module.  Dependencies can be copied into the
   target/dependencies directory using:

        mvn dependency:copy-dependencies

1. Installation into a Fiji/ImageJ2 install can be performed using the scijava
maven goal to populate the application:

        cd gdsc-smlm-ij
        mvn scijava:populate-app -Dscijava.app.directory=/usr/local/fiji
        cd ..

   where `/usr/local/fiji` is the root directory of the ImageJ install.

1. Manual installation must copy the gdsc-smlm-ij_* jar into the plugins
directory of ImageJ.

   Copy the dependencies into the plugins directory (or onto the Java
   classpath). Note that the Maven package routine puts all dependencies into
   the target/dependencies directory even if they are not required by the SMLM code
   (it does not check what functions are actually used by the code). The libraries
   required are:

        gdsc-core
        gdsc-core-ij
        gdsc-smlm
        commons-math3
        commons-lang3
        commons-rng-client-api
        commons-rng-core
        commons-rng-simple
        commons-rng-sampling
        ejml
        fastutil-core
        guava
        JTransforms
        JLargeArrays
        mfl-core
        mxparser
        protobuf-java
        protobuf-java-util
        quickhull3d
        xstream
        xpp3_min

   This excludes the 3D_Viewer functionality. View the dependencies using:

        mvn dependency:tree

1. The plugins will now appear under the 'Plugins > GDSC SMLM' menu in ImageJ.


Running from source
-------------------

Maven can be used to run ImageJ using a profile defined in the gdsc-ij-parent POM.
This profile compiles all classes and then executes ImageJ with the appropriate Java
classpath.

To run the GDSC SMLM plugins requires that they are detected by ImageJ. The plugins
use the ImageJ v1 plugin architecture which uses a search in the `plugins`
directory for classes with an underscore in the name. To avoid having to package the
classes to a jar file a bootstrap GDSC plugin loader is provided. To run
the GDSC SMLM plugins requires creating a link from the compilation target folder to
a `plugins` directory:


        cd gdsc-smlm-ij
        ln -s target/classes plugins
        mvn -P run-imagej

On start-up ImageJ will search in the `plugins` directory and find the bootstrap
plugin. This can then be run from the ImageJ Plugins menu to load the plugins
defined in the project's ImageJ plugins.config file.

Note: This file is normally detected by ImageJ when loading plugin jar files to
identify the available plugins. The bootstrap plugin has been written to duplicate
this functionality by reading the configuration and populating the ImageJ menu. It
allows ImageJ to be run with the latest classes using only the compile goal and
avoids having to package the classes to a jar. An alternative is to package the class
files and create a link from the jar location to a `plugins` directory:

        cd gdsc-smlm-ij
        ln -s target plugins
        mvn compile jar:jar
        mvn -P run-imagej

Note that when using this approach ImageJ may detect duplicate plugins if the sources
jar file has previously been created using the `package` goal. This can be avoided
by cleaning the target directory using `mvn clean`.

Using the jar approach will replicate the ImageJ `Plugins > GDSC SMLM` menu
structure.


Modifying the source
--------------------

The gdsc-smlm code was developed using the [Eclipse IDE](https://eclipse.org/).

Details of how to open the source code with Eclipse can be found in the eclipse
folder.


Legal
-----

See [LICENSE](LICENSE.txt)


# About #

###### Owner(s) ######
Alex Herbert

###### Institution ######
[Genome Damage and Stability Centre, University of Sussex](http://www.sussex.ac.uk/gdsc/)

###### URL ######
[GDSC SMLM ImageJ plugins](https://gdsc-smlm.readthedocs.io/en/latest/)
