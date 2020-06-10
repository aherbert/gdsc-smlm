PC PALM Plugins
===============

The following plugins are used to analyse the auto-correlation of a set of localisations using Pair Correlation (PC) analysis. This can provide information on whether the localisations are randomly distributed or clustered ([Sengupta *et al* , 2011], [Sengupta, *et al*, 2013], [Puchnar, *et al*, 2013]).

The analysis uses a set of localisations assumed to represent single molecules or single fluorescent bursts. The data is compared to itself using auto-correlation and a curve is computed as a function of the distance from the centre of the localisation. A flat curve is indicative that the distribution of localisations is no different from a random distribution. Shaped curves can be fit using models that apply to various distributions of localisations (random; fluctuations; or an emulsion).

The plugins are described in the following sections using the order presented on the ``Plugins > GDSC SMLM > PC PALM`` menu.


.. note::

    The PC-PALM plugins are under development and should be considered experimental.

    The plugins are subject to change and so no documentation has yet been produced to describe them. They are included in this package to assess interest from the community.

    If you wish to use the plugins please contact `Alex Herbert <a.herbert@sussex.ac.uk>`_ for more information.


PC-PALM Molecules
-----------------

Prepare results held in memory for analysis using pair correlation methods. This plugin either simulates results or filters localisations from a results set to a set of coordinates with time and photon signal information. The molecules are drawn on an image to allow regions of the data to be selected for analysis.

The pair correlation analysis assumes that the molecules are not moving. The input data should be localisation positions from a fixed sample. The use of fiducial markers to correct drift during image acquisition is recommended by Sengupta, *et al* (2013).

The output of the ``PC-PALM Molecules`` plugin is a set of molecule positions stored in memory. Each position should represent a single continuous on time of a fluorophore over one or more frames. Blinking events resulting in intervening dark frames should not be joined and must remain as separate molecules in the output molecule positions. This is because the model correlation functions used during analysis assume a blinking component for the same molecules.

The plugin has the following modes to prepare the molecule positions:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Run mode
     - Description

   * - PC-PALM
     - Follows the PC-PALM protocol for estimating localisation precision and then tracing molecules.

       See :numref:`{number}: {name} <pc_palm_plugins:Run Mode: PC-PALM>`.

   * - Manual tracing
     - Trace molecules using a distance and time threshold. Any spot that occurred within time threshold and distance threshold of a previous spot is grouped into the same trace as that previous spot.

       The tracing uses nearest neighbour assignment and more recent frames have priority (i.e. localisations with a small time gap are joined in preference to a larger time gap and smaller distance). More control over the tracing can be performed using an external tracing analysis and the molecules loaded direct using the ``In-memory`` mode.

   * - In-memory results
     - Select molecules previously loaded into memory.

   * - Simulation
     - Simulates localisations using different molecule distributions.

       See :numref:`{number}: {name} <pc_palm_plugins:Run Mode: Simulation>`.

Once the molecules have been identified the plugin will construct a super-resolution image of the data. The image is binary with a value of 1 for any pixel that has one or more molecules, otherwise the pixels are 0. The image has an ROI drawn on the image. This represents a rectangular region for analysis by the ``PC-PALM Analysis`` plugin. The image is used to select regions of data for analysis. Ideally some level of non-random distribution (clustering) should be visible in the image.

Optionally a high resolution image of the data can be constructed. This will attempt to draw each molecule on a single pixel by defining the pixel pitch as the reciprocal of the minimum distance between any two molecules. In the event of colocalisation this may be a very small distance and the reciprocal which defines the output pixel pitch will tend towards infinity. Thus the plugin provides the option to limit the output pixel pitch. Note that this image allows debugging the molecule distribution that is used in the ``PC-PALM Analysis`` plugin. Ideally it should be possible to render all molecules to separate pixels on a high resolution image to maximise the information available during pair correlation. This may not be possible if colocalisation of blinks from the same molecules is present. However if the minimum distance between molecules is high then the conversion of localisations to molecules may have grouped together separate molecules.

Parameters
~~~~~~~~~~

The following parameters are available:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Input
     - The input localisations.

   * - Use ROI
     - Map the ROI from the currently selected image to the input localisations and crop the selected region. This options is only shown if the current image has an area ROI.

       This option can be used to dynamically crop results from a dataset using a ROI drawn on a super-resolution render of the data.

   * - Run mode
     - The mode used to map the localisations into molecules representing distinct blinks of fluorophores.

   * - Image size
     - The size (in pixels) of the output super-resolution image of the final molecules.

   * - ROI size
     - The size of the ROI to create on the output super resolution image.

   * - Show high res image
     - Set to **true** to show a high resolution image of the final molecules.

   * - nm per pixel limit
     - Set the minimum pixel pitch (in nm) for the high resolution image. A setting of 0 will attempt to create the largest image possible.

   * - Clear results
     - Set to **true** to remove any PC-PALM analysis results from memory. Use this option to clear old results when starting a new analysis of a different dataset.


Run Mode: PC-PALM
~~~~~~~~~~~~~~~~~

This mode follows the PC-PALM protocol of Sengupta, *et al* (2013), steps 13 to 18. The localisation precision of each localisation is used to build a histogram of precision. The stored precision associated with the localisation is used if available or it is computed using the Mortensen formula ([Mortensen *et al*, 2010]). The histogram of the precision is fit using a skewed Gaussian function to determine the average positional uncertainty. The histogram and fitted function will be displayed.

The localisations are then traced using a distance of 2.5 times the average positional uncertainty with a time gap of successive frames. Joined localisations correspond to a single molecule that is active over multiple frames. The centroid of each molecule is computed using the intensity weighted coordinates of the localisations. The localisation precision of the molecule is computed using the weighted distance from the centroid and the weighted precision of each localisation (see formula 7b from Sengupta, *et al* (2013)).

The precision of each molecule is used to build a histogram of precision. Optionally this can include the localisations that are singles, i.e. they only occur in one frame and are not traced into a molecule. The precision histogram is again fit with a skewed Gaussian to determine the average positional uncertainty of each molecule and the results displayed. Note that the distribution of the precision of singles may be very different from the distribution of the precision of multi-frame molecules. This will be evident on the displayed histogram as a bimodal function. Thus including the singles may not create a good fit of the histogram using a skewed Gaussian.

Any molecule with a positional uncertainty above 3 times the average positional uncertainty is discarded to create the final molecules dataset. If the singles were not included in the previous stage to create the molecule precision histogram then any singles can optionally be included in the final dataset if they have a precision lower than the threshold.

The following parameters are available:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Histogram bins
     - The number of bins to use for the histogram. Use zero for auto.

   * - Singles mode
     - Specify how to handle single localisations that cannot be traced into molecules:

       - ``Ignore``: Remove from the data.
       - ``Include in molecules histogram``: Include them as molecules and allow their precision to contribute to the molecules precision histogram.
       - ``Include in final filtering``: Add to the final output dataset if below the precision threshold set using the average positional uncertainty of the traced molecules.

   * - Simplex fit
     - Set to **true** to perform a simplex fit of the skewed Gaussian. The default is a least square optimisation using numerical gradients.

   * - Show histograms
     - Set to **true** to show the histograms.

   * - Binary image
     - Set to **true** to display the super resolution image of the molecules as a binary image. If **false** then the image is a histogram where the value of each pixel is the molecule count at that pixel. This will not effect the later analysis and is used for visualisation purposes of the molecule density.

   * - Blinking rate
     - Set the blinking rate. This only effects the protein density that is reported to the ``ImageJ`` log window. The protein density is the molecule density divided by the blinking rate. Blinking rate is of interest during later PC-PALM analysis.


Run Mode: Simulation
~~~~~~~~~~~~~~~~~~~~

This mode allows simulation of data using different spatial distributions. Note that the simulation was created to verify that the models used during PC-PALM analysis correctly fit the data. Thus the options are based around clusters of loosely associated molecules. These clusters have an average size (number of members) and cover a circular region that should not overlap other regions. This is the data that is fit by the emulsion model of PC-PALM. Parameters have been added to simulate fluorophore blinking of each molecule in the cluster.

The simulation creates positions randomly within a defined 2D region. The positions may be cluster centres or molecules. If cluster centres then molecules are created for each cluster. Each molecule may blink multiple times resulting in 0 or more localisations per molecule. The number of blinks is the number of localisations per molecule. The localisations are created with a specified positional uncertainty to simulate the fitting precision of a typical super-resolution experiment.

The following cluster simulations are available:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Cluster simulation
     - Description

   * - None
     - Molecule positions are sampled uniformly from the 2D region. Each molecule position generates zero or more localisations due to blinking.

   * - Circles
     - Cluster positions are sampled uniformly from the 2D region. Each cluster contains zero or more molecules within a circle of a specified radius. The cluster circles may overlap.

   * - Non-overlapping circles
     - Create a mask using randomly distributed non-overlapping circles of a specified radius. Each circle has only 1 cluster of molecules. This simulation may not create the desired number of clusters due to space restrictions. If this occurs a message is logged to the ``ImageJ`` log window.

       *Note: This is the distribution modelled by the emulsion model during PC-PALM analysis.*

   * - Circles Mask
     - Create a mask using non-overlapping circles of a specified radius. The region is filled with circles. Sample cluster positions from any circle, there may be more than 1 cluster per circle.

When molecules are simulated into localisations (i.e. blinking) the plugin can optionally compute data on the cluster sizes and the intra-molecule distances. A histogram of the distances is computed and summary statistics recorded to the ``ImageJ`` log window. The plugin also computes the mean distance from a cluster member to the cluster centroid and records this in the ``ImageJ`` log window. These distances should be analysed in conjunction to the simulation settings and also to the model produced by subsequent PC-PALM analysis. If the intra-molecule distances are computed an option is provided to run the paricle linkage algorithm to perform clustering. The clustering distance is the 99\ :sup:`th` percentile from the actual intra-molecule distances. During clustering the join distances between the same molecule (intra-molecule) and between different molecules (inter-molecule) are collected and these are displayed in a cumulative histogram. If the clusters are not dense relative to the localisation precision then there should be a larger frequency of intra-molecule links. As the clusters reduce in size different molecules will begin to be joined and the frequency of inter-molecule links will increase.

The following parameters are available:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Molecules
     - The number of molecules to simulate. When using a ``Binomial`` distribution this is the number of clusters.

   * - Simulation size
     - The size of the region (in |micro|\ m).

   * - Blinking rate
     - The average number of blinks per molecule. When using a ``Binomial`` distribution this is the number of molecules per cluster.

   * - Blinking distribution
     - The distribution of the blinks per molecule.

       - ``Poisson``: Use a Poisson distribution.
       - ``Geometric``: Use a geometric distribution.
       - ``None``: Use a fixed number of blinks.
       - ``Binomial``: Use a binomial distribution. The ``Blinking rate`` parameter is used as the number of trials and the p-value of a blink occurring is collected via a dialog.

   * - Average precision
     - Define the standard deviation (in nm) of the random Gaussian added to each molecule position when generating localisations to simulate localisation uncertainty.

   * - Show histograms
     - Set to **true** to display a histogram of the intra-molecule distances and the number of blinks per molecule.

   * - Distance analysis
     - Set to **true** to perform clustering and distance analysis on the final localisations. Requires that ``Show histograms`` is **true**.

   * - Cluster simulation
     - Specify the cluster simulation.

   * - Cluster number
     - Specify the number of molecules per cluster. This is called the cluster number in the PC-PALM analysis.

   * - Cluster variation
     - Specify the standard deviation of the cluster number to allow variation in cluster size.

   * - Cluster radius
     - Specify the cluster radius (in nm).

   * - Show cluster mask
     - Set to **true** to show a mask of the region where a molecule may occur. The actual molecule positions are shown on the mask image. Note: This is different from the output binary image from ``PC-PALM Molecules`` that shows the final molecule dataset, i.e. each blink of the simulated molecule. This option shows the actual coordinate of the molecule without blinking and can be used to inspect the number of molecules in each cluster.


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
