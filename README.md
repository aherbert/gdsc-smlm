GDSC Single Molecule Light Microscopy (SMLM) ImageJ Plugins
===========================================================

The GDSC Single Molecule Light Microscopy (SMLM) plugins provide various tools
for single molecule localisation analysis. This includes PALM, STORM and other
single molecule microscopy methods.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/aherbert/gdsc-smlm.svg?branch=master)](https://travis-ci.com/aherbert/gdsc-smlm)
[![Coverage Status](https://coveralls.io/repos/github/aherbert/gdsc-smlm/badge.svg?branch=master)](https://coveralls.io/github/aherbert/gdsc-smlm?branch=master)
[![Documentation Status](https://readthedocs.org/projects/gdsc-smlm/badge/?version=latest)](https://gdsc-smlm.readthedocs.io/en/latest/?badge=latest)


Read the [SMLM User Manual](https://gdsc-smlm.readthedocs.io/en/latest) for full details.

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


Install
-------

The SMLM plugins are distributed using an ImageJ2/Fiji update site.

To install the plugins using Fiji (an ImageJ distribution) just follow the
instructions [How_to_follow_a_3rd_party_update_site](http://fiji.sc/How_to_follow_a_3rd_party_update_site)
and add the GDSC SMLM update site. All the plugins will appear under the 'Plugins > GDSC SMLM' menu.


Installation from source
------------------------

The source code is accessed using git and built using Maven.

The code depends on the gdsc-analytics, gdsc-test and gdsc-core artifacts so
you will have to install these to your local Maven repository before building:

1. Clone the required repositories

        git clone https://github.com/aherbert/gdsc-analytics.git
        git clone https://github.com/aherbert/gdsc-test.git
        git clone https://github.com/aherbert/gdsc-core.git
        git clone https://github.com/aherbert/gdsc-smlm.git

2. Build the code and install using Maven

        cd gdsc-analytics
        mvn install
        cd ../gdsc-test
        mvn install
        cd ../gdsc-core
        mvn install
        cd ../gdsc-smlm
        mvn package

	This will produce a gdsc_smlm-[VERSION].jar file in the target directory.
	All dependencies are copied into the target/dependencies directory.

3. Copy the gdsc_smlm* jar into the plugins directory of ImageJ.

4. Copy the dependencies into the plugins directory (or onto the Java
classpath). Note that the Maven package routine puts all dependencies into
the target/dependencies directory even if they are not required by the SMLM code
(it does not check what functions are actually used by the code). The libraries
required are:

        gdsc-analytics
        gdsc-core
        commons-rng-client-api
        commons-rng-core
        commons-rng-simple
        commons-rng-sampling
        xstream
        xpp3_min
        ejml
        JTransforms
        JLargeArrays
        protobuf-java
        protobuf-java-util

Libraries required if not using the Fiji distribution of ImageJ:

        commons-math3
        commons-lang3
        trove4j
        guava

5. To enable use of the 3D viewer also add the Java3D libraries and the native libraries for your
platform. This is unnecessary if using Fiji as the correct libraries are present. The libraries
required are:

        3D_Viewer
        j3dcore
        j3dutils
        jogl-all
        jogl-all-*-natives-[platform]
        gluegen-rt
        gluegen-rt-main
        gluegen-rt-*-natives-[platform]
        vecmath-*-scijava-2.jar
        VIB-lib


6. The plugins will now appear under the 'Plugins > GDSC SMLM' menu in ImageJ.


Running from source
-------------------

Maven can be used to run ImageJ using a profile defined in the gdsc-ij-parent POM:

        mvn -P run-imagej

This profile compiles all classes and then executes ImageJ with the appropriate Java classpath. A
default plugin can then be run from the ImageJ Plugins menu to load the plugins defined in the
project's ImageJ plugins.config file. 

Note: This file is normally detected by ImageJ when loading plugin jar files to identify the
available plugins. The default plugin has been written to duplicate this functionality by reading
the configuration and populating the ImageJ menu.


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
