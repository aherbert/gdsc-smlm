.. index:: ! ImageJ Plugins Overview

ImageJ Plugins Overview
=======================

The various ``ImageJ`` plugins allow the processing of single-molecule light microscopy images and analysis of the results. This includes:

*   Finding and fitting spots on the image
*   Reconstructing an image from the list of localisations
*   Analysis of the blinking rate of fluorophores
*   Tracing of fluorophore molecules through time

The plugins used to analyse a set of localisations only require that the localisations be loaded into memory. The localisations do not have to be computed by the SMLM fitting plugins and can be generated by another program. For example it is possible to read the following file formats for analysis:

*   rapidSTORM
*   Nikon NSTORM
*   MicroManager Tagged Spot File (`TSF <https://micro-manager.org/wiki/Tagged_Spot_File_(tsf)_format>`_)

The plugins have been divided into the following sub-sets:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Set
     - Description

   * - :ref:`Fitting <fitting_plugins:Fitting Plugins>`
     - For identification of localisations on an image.

   * - :ref:`Results <results_plugins:Results Plugins>`
     - Allow loading and saving results in different formats. Results can be filtered to subsets and compared to a reference set, e.g. for benchmarking.

   * - :ref:`Analysis <analysis_plugins:Analysis Plugins>`
     - Perform analysis on localisations for example blinking rate estimation, molecule tracing and Fourier image resolution.

   * - :ref:`PC PALM <pc_palm_plugins:PC PALM Plugins>`
     - Plugins for Pair Correlation (PC) analysis.

   * - :ref:`Model <model_plugins:Model Plugins>`
     - Simulate single-molecule images.

   * - :ref:`Calibration <calibration_plugins:Calibration Plugins>`
     - Estimate PSF widths and allow calibration of the imaging camera noise and gain.

   * - :ref:`Tools <tools_plugins:Tools Plugins>`
     - Utility plugins for image manipulation.

   * - :ref:`Toolset <toolset_plugins:Toolset Plugins>`
     - For install of the SMLM Toolset and configuration of the SMLM Tools window.


.. index:: ! Usage Tracking

Usage Tracking
--------------

.. note::

    **No tracking** is performed by the GDSC plugins.

    This section details the tracking code that was used in a previous version. It is maintained to describe what and how anonymous user tracking was performed. Tracking was disabled in May 2022.

    Tracking had to be enabled via a dialog in a similar way to a user preferences pop-up dialog in a web browser. A key feature is that tracking was disabled by default and had to be enabled. If your usage was tracked in a previous version then this was only because you provided consent (and after being provided with similar information to that shown below).

    Tracking provided us with usage numbers to add to research grant applications for work that included the GDSC software. It was a useful facility but non-essential; we are grateful to those who participated.

Tracking allows a provider to understand its user base. In the context of ImageJ plugins it would help to improve the plugins. Example questions that can be answered are:

* Which plugins are popular?
* How often people update to the latest version?
* How many times people use different plugins in one ``ImageJ`` analysis session?
* What plugins are frequently used together, and in what order?
* What computer platform and software do we need to support?

To understand how the GDSC SMLM plugins were being used around the world we previously added some code to track usage.

**Tracking was disabled by default. You had to opt-in and could opt-out at any time.**

Tracking provided us with usage numbers to add to research grant applications for work that included the GDSC software. It also showed that some plugins we do not regularly use in GDSC research work had an uptake of users in the community. Thus we will continue to publish all plugins we develop in the hope that they are useful.

To track usage we used Google Analytics, a web analytics service provided by Google, Inc. ("Google"). The information about use of the plugins was transmitted to and stored by Google on servers in the United States.


.. index:: Usage Data

Usage Data
~~~~~~~~~~

Usage data was collected as if performing anonymous tracking of web page views by a web browser. The code identified each plugin as a page. The ImageJ instance acted as the web browser client. Each time a GDSC plugin was executed this event was recorded as if a user viewed a web page.

No personal information was transmitted to Google. There is no information that uniquely identifies a person. No data currently open in ``ImageJ`` was sent. Only data about the GDSC plugin names were sent. No data about any other actions within ``ImageJ`` were sent.

Data was only sent when a GDSC plugin was run. The following data was sent:

* Name of the plugin
* Plugin version
* ``ImageJ`` version
* Java version
* Operating system (e.g. Windows, Linux, Mac OS)
* Screen resolution\ :sup:`1`

:sup:`1` This can be used to design dialogs that fit on the screen.


.. index:: Tracking Identifiers

Tracking Identifiers
~~~~~~~~~~~~~~~~~~~~

All the usage information could have been collected without tracking the same individual. Each ``ImageJ`` session would count as a new individual using the software. To allow the distinction of new or repeat use a random identifier for the individual was generated. This was a 128-bit random `UUID <https://en.wikipedia.org/wiki/Universally_unique_identifier>`_ with a very low chance of being repeated on another ``ImageJ`` instance. This is stored in the ``ImageJ`` preferences file under the key ``.gdsc.ga.clientId``. It can only be read by programs with permission to read the ``ImageJ`` preferences file. This is usually in the user home directory and so would only be read by programs run by the user.

The tracking identifier was not created unless the user opted-in to tracking. The default was no tracking.


.. index:: Performance

Performance
~~~~~~~~~~~

Note that usage tracking should not slow down ``ImageJ``. The tracking was performed in the background only when the computer was doing nothing else. If it was always too busy then no tracking data was sent.

If not connected to the internet then the tracker identified that messages could not be sent and shutdown.


.. index:: User Preferences

User Preferences
~~~~~~~~~~~~~~~~

To notify users of analytics a dialog was shown when the user runs a GDSC plugin and their preferences were not known. This presented the user with the choice to opt-in.

Preferences were saved in the ``ImageJ`` preferences file. This is written when ``ImageJ`` closes and stores user preferences between ``ImageJ`` sessions. The following settings were stored:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Key
     - Description

   * - gdsc.ga.lastVersion
     - The version number of the most recently run GDSC analytics code. If a new version is released that does not match this stored version then the user preferences dialog would be shown again to ensure the preferences are correct.

   * - gdsc.ga.clientId
     - The random UUID for the user. This allows repeat sessions to be distinguished from new users.

   * - gdsc.ga.state
     - A flag indicating the user preference:

       * -1: Opt-out - no tracking is performed
       * 0: Unknown - show the user preferences dialog
       * 1: Opt-in - tracking is active

The user could change their options at any time by running the ``SMLM Usage Tracker`` plugin.


.. index:: Google Analytics Details

Google Analytics Details
~~~~~~~~~~~~~~~~~~~~~~~~

The code used the Analytics Measurement Protocol. This allows any web connected application to record simple usage data. The GDSC plugins fully complied with the Google protocol policy. In brief this means that no data should be sent to Google that allows the user to be personally identified, and the user can opt-out at any time. Breaking these rules results in deletion of both data and account by Google.

Data was sent using a secure HTTP connection. All information sent was non-personal however Google may have stored an IP address associated with the data.

Data used Google's Universal Analytics. This service will stop processing data in `July 2023 <https://support.google.com/analytics/answer/11583528>`_. It is currently unknown when the data will be made inaccessible.
