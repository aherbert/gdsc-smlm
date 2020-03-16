PC PALM Plugins
===============

The following plugins are used to analyse the auto-correlation of a set of localisations using Pair Correlation (PC) analysis. This can provide information on whether the localisations are randomly distributed or clustered ([Sengupta *et al* , 2011], [Sengupta, *et al*, 2013], [Puchnar, *et al*, 2013]).

The analysis uses a set of localisations assumed to represent single molecules or single fluorescent bursts. The data is compared to itself and a curve is computed as a function of the distance from the centre of the localisation. A flat curve is indicative that the distribution of localisations is no different from a random distribution. Shaped curves can be fit using models that apply to various distributions of localisations (random; fluctuations; or an emulsion).

The plugins are described in the following sections using the order presented on the
``Plugins > GDSC SMLM > PC PALM``
menu.


.. note::

    The PC-PALM plugins are under development and should be considered experimental.

    The plugins are subject to change and so no documentation has yet been produced to describe them. They are included in this package to assess interest from the community.

    If you wish to use the plugins please contact `Alex Herbert <a.herbert@sussex.ac.uk>`_ for more information.


PC-PALM Molecules
-----------------

Prepare results held in memory for analysis using pair correlation methods. This plugin either simulates results or filter results from a results set to a set of coordinates with time and photon signal information. The localisations are drawn on a binary image to allow regions of the data to be selected for analysis.


PC-PALM Analysis
----------------

Perform pair-correlation analysis in the frequency domain as per the paper by [Sengupta *et al* , 2011], [Sengupta, *et al*, 2013] to produce a *g(r)* correlation curve.


PC-PALM Spatial Analysis
------------------------

Perform pair-correlation spatial analysis as per the paper by [Puchnar, *et al*, 2013]. This methods plots the molecule density around each localisation as a function of distance from the localisation.


PC-PALM Save Results
--------------------

Saves all the PC-PALM results held in memory to a results folder.


PC-PALM Load Results
--------------------

Load all the PC-PALM results from a results folder to memory.


PC-PALM Fitting
---------------

Combines multiple correlation curves calculated by PC-PALM Analysis into an average curve and fits the curve using various models.


PC-PALM Clusters
----------------

Clusters localisations using a distance threshold and produces a histogram of cluster size. This can be fit using a zero-truncated negative Binomial distribution (with parameters *n*, *p*) to calculate the size of the clusters (*n*) and the probability of seeing a fluorophore (*p*).
