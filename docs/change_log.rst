Change Log
==========

.. contents::
   :local:


Version 1.1
-----------

Minor release of GDSC SMLM.

* Fix the :numref:`{name} <analysis_plugins:Trace Diffusion>` plugin to use descending diffusion
  coefficients for multi-population models. This corrects selection of the model using the minimum
  relative difference between successive coefficients.
* Added the :numref:`{name} <analysis_plugins:Residence Time Analysis>` plugin to fit residence
  times of stationary (bound) molecules.
* Fix the :numref:`{name} <results_plugins:Save Localisations>` plugin error when the PSF is
  undefined.
* Update the :numref:`{name} <tools_plugins:TIFF Series Viewer>` plugin with a ``Validate`` output
  mode.
* Update the ``View a camera model`` action in the
  :numref:`{name} <calibration_plugins:Camera Model Manager>` plugin. Allow display of histograms
  of the model pixel data; and outlier pixels.
* Fix the :numref:`{name} <fitting_plugins:Fit Configuration>` plugin to allow templates with a
  sCMOS camera.
* Added :numref:`{name} <analysis_plugins:Trace Molecules (Multi)>` to support multiple dataset
  analysis.


Version 1.0.2
-------------

Patch release of GDSC SMLM.

* Update to GDSC Core 2.0.2


Version 1.0.1
-------------

Patch release of GDSC SMLM.

* Correct application of templates with PSF settings.


Version 1.0
-----------

First working version of GDSC SMLM.

Requires Java 8.
