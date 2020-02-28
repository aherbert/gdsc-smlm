.. index:: ! Installation

Installation
============

The plugin is designed to run within ``ImageJ``. This can be the original ``ImageJ`` (version 1) or any modified version of ``ImageJ`` such as ``ImageJ2`` or ``Fiji``.


.. index:: ! Install using ImageJ2/Fiji

Install using ImageJ2/Fiji
--------------------------

The SMLM plugins are distributed using an ``ImageJ2``/``Fiji`` update site. This allows the plugins to be easily installed and kept up-to-date. To install the plugins using ``Fiji`` just follow the instructions here:

* https://imagej.net/Following_an_update_site

Add the ``GDSC SMLM2`` update site. ``Fiji`` will automatically check for a new version during start-up and install it if desired.

All the plugins will appear under the ``Plugins > GDSC SMLM`` menu.


.. index:: ! Install using ImageJ version 1

Install using ImageJ version 1
------------------------------

The plugin is designed to run within ``ImageJ``. You can obtain the latest version of ``ImageJ`` from:

* https://imagej.nih.gov/ij/download.html

To get the plugins you can download the latest Jar files from the update site and put them in your ``ImageJ`` plugins folder. The jars can be found here:

* https://sites.imagej.net/GDSC-SMLM2/

You will also need to install the additional Apache Commons Math 3 library and EJML. These are already included in ``Fiji`` so are not on the update site. You can get the files here:

* `Apache Commons Math 3 <https://repo.maven.apache.org/maven2/org/apache/commons/commons-math3/3.6.1/>`_
* `EJML v0.24 <https://sourceforge.net/projects/ejml/files/v0.24/>`_

Place all of the following ``jar`` files into the ``ImageJ`` plugins directory:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Jar File
     - Description

   * - commons-math3
     - Contains math routines.

   * - commons-rng-*
     - Contains random number routines.

   * - gdsc-analytics
     - Contains the GDSC analytics library.

   * - gdsc-core
     - Contains the GDSC core library.

   * - gdsc-smlm
     - Contains the GDSC SMLM plugins.

   * - ejml
     - Efficient Java Matrix Library for linear algebra. Requires version 0.24. This is the version used by ``Fiji``.

   * - JTransforms
     - Library for multi-threaded Fourier transforms.

   * - JLargeArrays
     - Support library for JTransforms.

   * - protobuf-java-*
     - Google protocol buffers library for readin/writing data.

   * - xstream
     - A library for reading/writing XML.

The plugins will be visible under the ``Plugins > GDSC SMLM`` menu.

An ``ImageJ`` toolset can be installed using the ``GDSC SMLM > Toolset > Install SMLM Toolset`` plugin.
When selected the toolset adds a set of buttons on the ``ImageJ`` toolbar for commonly used plugins.
More details can be found in section :numref:`{number}: {name} <Toolset_Plugins:Toolset Plugins>`.
