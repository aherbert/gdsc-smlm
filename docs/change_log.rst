Change Log
==========

.. contents::
   :local:


Version 1.2
-----------

Minor release of GDSC SMLM.

* Update to GDSC Core 2.2
* Update the :numref:`{name} <analysis_plugins:Residence Time Analysis>` plugin to allow
  filtering of the counts histogram to remove long residence times.
* Update the tracing algorithm used by the :numref:`{name} <analysis_plugins:Trace Molecules>`
  plugin. This now performs tracing of localisations to existing tracks using the end time of the
  track only (i.e. its last known position).

  This corrects a bug in tracing using time thresholds above 1 where jump distances between
  adjacent localisations could be above the distance threshold if a localisation was matched to
  a prior position of the track. E.g. Given a track with localisations A and B, a new
  localisation C will only be compared to position B; the previous algorithm allowed comparison
  to A and B.

  The use of the exclusion distance has been updated to apply a stricter exclusion using all
  possible matches across the entire time distance. Previously this was applied using alternative
  matches in the same frame only.
* Update the interactive results table to allow a source image to be attached to the results.
  The results can be added to the source image as an overlay.
* Added the :numref:`{name} <results_plugins:Trace Viewer>` plugin to display a table of traced
  results.
* Update from Commons Math to Commons Statistics for inference testing and univariate
  descriptive statistics.
* Change use of EJML API from 0.24 to 0.41 (change mandated by scijava-pom for ImageJ
  compatibility). This is binary incompatible change to the public API.
* Change use of java3d API from 4.x to 5.x (change mandated by scijava-pom for ImageJ
  compatibility). This is binary incompatible change to the public API.
* Add ``Use stack`` option to :ref:`tools_plugins:Overlay Image`.


Version 1.1
-----------

Minor release of GDSC SMLM.

* Update to GDSC Core 2.1
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
