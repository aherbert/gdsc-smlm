.. index:: calibration plugins

Calibration Plugins
===================

The following plugins contain functionality to estimate the width of the point spread function (PSF) for an imaging set-up and analyse the noise and gain of the imaging camera.

The plugins are described in the following sections using the order presented on the
``Plugins > GDSC SMLM > Calibration``
menu.

.. index:: psf calculator

PSF Calculator
--------------

A simple plugin that estimates your Gaussian approximation to the PSF using the microscope imaging parameters and the wavelength of light (:numref:`Figure %s <fig_psf_calculator_dialog>`).

.. _fig_psf_calculator_dialog:
.. figure:: images/psf_calculator_dialog.png
    :align: center
    :figwidth: 80%

    PSF Calculator dialog

The calculator uses the following formula:

.. math::

    \mathit{Airy\:Width}=\frac{\lambda}{2\pi\mathit{NA}}

:math:`\mathit{Airy\:Width}` is the width of the Airy pattern,
:math:`\lambda` is the wavelength (in nm), and
:math:`\mathit{NA}` is the Numerical Aperture.

The Airy profile can be approximated by a Gaussian profile. The equivalent Gaussian profile is created using a standard deviation of 1.323 times the Airy width. The PSF Calculator will show a plot of the Airy profile (blue) and the corresponding Gaussian profile (red; :numref:`Figure %s <fig_psf_calculator_profile_plot>`). This is interactively updated when the parameters for the calculator are changed:

.. _fig_psf_calculator_profile_plot:
.. figure:: images/psf_calculator_profile_plot.png
    :align: center
    :figwidth: 80%

    PSF Calculator PSF profile plot

Note that the Gaussian is a good approximation until the tails of the Airy pattern. The Airy pattern contains waves of decreasing power out to infinity which are not modelled by the Gaussian.

The calculator allows for additional adjustments to be made to the calculated Gaussian standard deviation. To account for optical aberrations in the microscope the width is allowed to be wider by a set proportionality factor. The Gaussian standard deviation (:math:`s`) can then be adjusted (:math:`s_a`) for an accurate representation on square pixels (of width :math:`a`) using the following formula:

.. math::

    s_{a}=\sqrt{s^{2}+a^{2}/12}

The Gaussian Half-Width at Half-Maxima (HWHM) is calculated from the standard deviation by multiplying by :math:`\sqrt{2\ast \ln (2)}=1.177`.

The following table describes the parameters of the plugin. The calculated properties are updated dynamically.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Pixel pitch (|micro|\ m)
     - The camera pixel size.

   * - Magnification
     - The objective magnification.

   * - Beam Expander
     - Any addition magnification.

   * - Pixel pitch (nm)
     - The calculated image pixel size.

   * - Wavelength (nm)
     - The wavelength of light used for the image (lambda).

   * - Numerical Aperture
     - The objective numerical aperture (NA).

   * - Proportionality Factor
     - The proportionality factor (set to 1 to match the Gaussian to the Airy profile).

   * - Adjust for square pixels
     - Perform square pixel adjustment (set to **false** to match the Gaussian to the Airy profile).

   * - Airy Width (nm)
     - The calculated PSF Airy width in nanometres.

   * - Airy Width (pixels)
     - The calculated PSF Airy width in pixels.

   * - StdDev (nm)
     - The calculated PSF Gaussian standard deviation in nanometres.

   * - StdDev (pixels)
     - The calculated PSF Gaussian standard deviation in pixels.

   * - HWHM (pixels)
     - The calculated PSF Gaussian HWHM in pixels.


Note that the first three fields are only used to calculate the image pixel pitch. If this is already known then it can be entered into the Pixel pitch (|micro|\ m) field (note the use of micrometres and not nanometres) and the ``Magnification`` and ``Beam expander`` can be set to 1. The ``Pixel pitch`` in nanometres is then used to convert the calculated widths in nanometers to pixel dimensions.

Clicking ``OK`` will save the PSF standard deviation in pixels to the global properties. This will be used in the ``Peak Fit`` plugin.

Please contact us if you have feedback on the calculated width from the plugin verses your measured PSF using quantum dots (or other single-point light sources) on calibration images.

.. index:: psf estimator

PSF Estimator
-------------

A plugin that estimates the PSF using a test image. The fit configuration is the same as in the ``Fit Configuration`` plugin with extra parameters provided to control the estimation. The ``PSF Estimator`` dialog is show in :numref:`Figure %s <fig_psf_estimator_dialog>`. Note that a second dialog will collect parameters specific for the selected ``Fit solver``.

.. _fig_psf_estimator_dialog:
.. figure:: images/psf_estimator_dialog.png
    :align: center
    :figwidth: 80%

    PSF Estimator dialog

The estimator uses the starting configuration to fit N peaks taken from randomly selected frames in the image stack. The averages of the fitted parameters are then used as the start parameters to perform fitting again. This iterates until the Gaussian parameters do not significantly change. The parameters controlling the estimation are described below.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description.

   * - Number of Peaks
     - The number of fitted peaks to use to estimate the Gaussian parameters. The parameters are estimated by averaging across all the fitted peaks.

   * - p-Value
     - The p-value to use for significance testing, i.e. are the parameters the same using a Student's T-test at the given significance.

   * - Updates preferences
     - If selected the plugin will update the global configuration with the calculated PSF values.

   * - Log progress
     - Log progress of the estimator to the ``ImageJ`` log window.

   * - Iterate
     - When the PSF parameters have converged and a ``Free`` fitting option was chosen a test is done to determine if the angle or Y-width are significant. If not significant the estimator will ignore the insignificant parameter and restart using a simpler PSF. The order of iterations is:

       ``Free`` > ``Free Circular`` > ``Circular``

       Note these statistics often don't work so unless you expect astigmatism you can choose to start with a ``Circular`` Gaussian and just find the estimated widths.

   * - Show histograms
     - Show a histogram of the estimated parameters from the final fitting run. A histogram is shown for each parameter. These can be used to verify the mean of the parameter distribution is a suitable estimate for the parameter.

   * - Histogram bins
     - The number of bins to use on the histograms. Set to zero to auto-scale the bin width.


Note that the estimator may not find any peaks if the fitting parameters are badly configured. The estimator can be reset to defaults by holding down the ``Control`` key when running the plugin. The default values are shown below:

.. list-table::
   :widths: 40 40
   :header-rows: 1

   * - Parameter
     - Value

   * - Initial StdDev0,1
     - 1

   * - Initial Angle
     - 0

   * - Spot filter type
     - Single

   * - Spot filter
     - Mean

   * - Smoothing
     - 1.3

   * - Search width
     - 1

   * - Border
     - 1

   * - Fitting width
     - 3

   * - Fit Solver
     - Least Squares Estimator

   * - *Fit criteria*
     - *Least-squared error*

   * - *Significant digits*
     - *5*

   * - *Coord delta*
     - *0.0001*

   * - *Lambda*
     - *10*

   * - *Max iterations*
     - *20*

   * - Fit function
     - Circular

   * - Fail limit
     - 3

   * - Include neighbours
     - True

   * - Neighbour height
     - 0.3

   * - Residuals threshold
     - 1

   * - Shift factor
     - 1

   * - Signal strength
     - 0

   * - Min photons
     - 30

   * - Width factor
     - 2


.. index:: mean-variance test

Mean-Variance Test
------------------

The ``Mean-Variance Test`` plugin can be used to calculate the gain and read noise of the microscope Charged Coupled Device (CCD) camera. The plugin requires a set of calibration images. A single-image mode is available but will provide less information on the camera.

.. index:: multiple input images

Multiple Input Images
~~~~~~~~~~~~~~~~~~~~~

When run the plugin will present a folder selection dialog. The folder should contain a set of calibration images. All the images should be taken of the same view with the camera in the same gain mode.

At least one image should be taken with no exposure time. This is the image the camera records when no light has been registered on the sensor and is called the bias image.

The remaining images should be a representative series of different exposures. The purpose is to analyse how the image noise varies with exposure time. In order for the analysis to be valid no images should saturate the camera bit-depth. E.g. for a 12-bit camera all images should have pixel values below :math:`2^{12}-1 = 4095`.

All the images in the folder are opened and processed by the plugin. Each image must contain at least 2 frames. If the filename contains a valid integer delimited by a space or a period character (``.``) then this will be taken as the exposure time. Otherwise an arbitrary exposure time is used, either zero for the first image (alphabetically sorted) or 9999 for the rest.

.. index:: analysis

Analysis
~~~~~~~~

If all the images are valid (contain at least 2 frames) then the plugin will perform the mean-variance test. The average value of the bias images is used as the bias. Each image is then analysed in turn. The mean of each frame is computed. Then a pairwise difference image (i.e. one frame subtracted from the other) is computed for all-vs-all frames. The variance of the difference image is recorded and used to approximate the camera gain:

.. math::

    \mathit{gain}=\frac{\mathit{variance}}{\mathit{mean}-\mathit{bias}}

This is recorded in a summary table. A graph is then produced of the mean verses the variance. This data is fitted with a straight line. The gradient of the line is the camera gain. The read noise of the camera is computed as:

.. math::

    \mathit{read\:noise}=\frac{\sqrt{\mathit{bias\:variance}}}{\mathit{gain}}

If the bias has multiple difference images then the average bias variance is used to calculate the read noise.

.. index:: output

Output
~~~~~~

The plugin produces a summary table of the analysis for each pair of frames. The table shows the following data:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Column
     - Description

   * - Image
     - The source image.

   * - Exposure
     - The image exposure. This is the first integer number delimited by a space or period in the image title or, if no number can be found in the image title, zero for the first image and 9999 for the others.

   * - Slice1
     - The first frame (slice) used from the image.

   * - Slice2
     - The second frame (slice) used from the image.

   * - Mean1
     - The mean of slice 1.

   * - Mean2
     - The mean of slice 2.

   * - Mean
     - The mean of both slices.

   * - Variance
     - The variance of the difference image.

   * - Gain
     - The gain estimate:

       :math:`\mathit{gain}=\frac{\mathit{variance}}{\mathit{mean}-\mathit{bias}}`.


The plugin will produce a plot of the mean-variance data as show in :numref:`Figure %s <fig_mean_variance_plot>`. The plot will show the best fit line in red. If the data points with the highest mean lie well under the line it is possible that these images had saturated pixel values and should be removed from the input data set.

.. _fig_mean_variance_plot:
.. figure:: images/mean_variance_plot.png
    :align: center
    :figwidth: 80%

    Mean-variance plot produced by the Mean-Variance Test plugin.

    The best fit line is shown in red.

The plugin reports the final calculated gain and read noise to the ``ImageJ`` log, e.g.

.. code-block:: text

    Mean Variance Test
    Directory = /images/CameraCalibration/CameraGain-2-EmGain-0/
    Bias = 515.4 +/- 7.4 (ADU)
    Variance = -21.78 + 0.1557 * mean
    Read Noise = 47.53 (e-)
    Gain = 1 / 6.422 (ADU/e-)

The parameters for the best fit line are shown as ``Variance = a + b * mean``. The parameter *b* is the gain. The read noise is shown in electrons. The units for the gain are Analogue-to-Digital Unit (ADU) per electron.

Note that the gain can be expressed as electrons per ADU and so the output shows the gain using 1 over the reciprocal of the fit parameter to allow comparison with manufacturer gain values. E.g. In the example above 1 / 6.422 = 1 / (1 / 0.1557) and the gain would be 6.422 e-/ADU.

.. index:: single image mode

Single Image Mode
~~~~~~~~~~~~~~~~~

The plugin can be run using a single image. Single image mode cannot compute the camera bias or read noise and the gain values are not as accurate as the full test using multiple images.

Hold the ``Shift`` key down when running the plugin and the analysis will be performed on the currently active image. The image must have more than one slice to allow difference images to be computed and should be a white light image with a constant uniform exposure across the image field, i.e. no significant image features.

In single-image mode the plugin will compute the pairwise comparison of consecutive frames in the image and for each pair compute the approximate camera gain:

.. math::

    \mathit{gain}=\frac{\mathit{variance}}{\mathit{mean}-\mathit{bias}}

The bias must be provided since there is no input bias image; the plugin will ask the user to input the camera bias. The results will be displayed in a table as described above.

The plugin provides a plot of gain verses slice and a histogram of the values. These can be used to determine if the gain is constant throughout the image and so is a good estimate.

.. index:: mean-variance test (em-ccd)

Mean-Variance Test (EM-CCD)
---------------------------

This plugin is similar to the ``Mean-Variance Test`` plugin but is used on images taken using an Electron Multiplying (EM) CCD camera. An EM-CCD camera uses a multiplication device to increase the number of electrons that are extracted from the imaging sensor before the electrons are counted. The average number of electrons output from the multiplying device for each input electron is a constant known as the EM-gain. The plugin will compute the EM-gain of the camera using a set of calibration images. A single-image mode is available but will provide less information on the camera.

The analysis can only be performed if the gain for the camera in non-EM mode is already known. If the ``Mean-Variance Test`` plugin has been used to calculate the gain in the same ``ImageJ`` session then the value will be stored in memory. If the camera gain is not known then using a value of 1 will allow the plugin to run and the output EM-gain will be the total gain of the system.

.. index:: multiple input images

Multiple Input Images
~~~~~~~~~~~~~~~~~~~~~

Input images requirements are the same as the ``Mean-Variance Test`` plugin: images should be taken of the same view using different exposure times. Each image must have at least two frames. All images must be taken with the camera in the same gain mode and EM-gain mode. A bias image (zero exposure) must be provided.

If all the images are valid the plugin will show a dialog asking for the camera gain (:numref:`Figure %s <fig_mean_var_test_em_gain_dialog>`). This will remember the last entered value or the value computed by the ``Mean-Variance Test`` plugin.

.. _fig_mean_var_test_em_gain_dialog:
.. figure:: images/mean_var_test_em_gain_dialog.png
    :align: center
    :figwidth: 80%

    EM-gain dialog of the Mean-Variance Test (EM-CCD) plugin

.. index:: analysis

Analysis
~~~~~~~~

The images are analysed as per the
``Mean-Variance Test``
plugin. However the analysis of the difference image is used to approximate the camera EM-gain:

.. math::

    \mathit{EM\:gain}=\frac{\mathit{variance}}{(\mathit{mean}-\mathit{bias})\:(2\times\mathit{gain})}

This is recorded in a summary table. A graph is then produced of the mean verses the variance. This data is fitted with a straight line. The gradient of the line is the EM-gain multiplied by twice the camera gain therefore the EM-gain can be computed as:

.. math::

    \mathit{EM\:gain}=\frac{\mathit{gradient}}{2\times\mathit{gain}}

.. index:: output

Output
~~~~~~

The plugin summary table and mean-variance plot are the same as the ``Mean-Variance Test`` plugin. The final calculated EM-gain and total gain is reported to the ``ImageJ`` log, e.g.

.. code-block:: text

    Mean Variance Test
    Directory = /images/CameraCalibration/CameraGain-2-EmGain-250/
    Bias = 512.3 +/- 13.15 (ADU)
    Variance = -36550.0 + 79.66 * mean
    Read Noise = 0.3301 (e-)
    Gain = 1 / 6.422 (ADU/e-)
    EM-Gain = 255.8
    Total Gain = 39.83 (ADU/e-)

The total gain is the EM-gain multiplied by the camera gain. As can be seen from comparison of the analysis results with and without the EM-mode the use of EM amplification dramatically reduces the camera read noise and greatly enhances the pixel values (ADUs) produced per electron. This allows images of weak photon signals to be made, for example in single-molecule light microscopy.

The total gain can be used to convert the ADUs into photons if the camera quantum efficiency (QE) is known. The QE states how many photons are converted into an electron charge when they hit the camera sensor; the QE units are electrons per photon (e-/photon). This can be provided by the camera manufacturer and is dependent on the wavelength of light. The photon signal is therefore:

.. math::

    \mathit{Photons}=\frac{\mathit{ADUs}}{\mathit{total\:gain}\times\mathit{QE}}

The total gain multiplied by the QE is known as the system gain. The system gain is used as an input parameter in the ``Peak Fit`` plugin to convert the pixel values into photons.

.. index:: single image mode

Single Image Mode
~~~~~~~~~~~~~~~~~

The plugin can be run using a single image. Single image mode cannot compute the camera bias or read noise and the gain values are not as accurate as the full test using multiple images.

Hold the ``Shift`` key down when running the plugin and the analysis will be performed on the currently active image. The image must have more than one slice to allow difference images to be computed and should be a white light image with a constant uniform exposure across the image field, i.e. no significant image features.

In single-image mode the plugin will compute the pairwise comparison of consecutive frames in the image and for each pair compute the approximate camera gain:

.. math::

    \mathit{EM\:gain}=\frac{\mathit{variance}}{(\mathit{mean}-\mathit{bias})\:(2\times\mathit{gain})}

The bias must be provided since there is no input bias image; the plugin will ask the user to input the camera bias and camera gain. Using a camera gain of 1 will calculate the total gain of the system. The results will be displayed in a table as described above.

The plugin provides a plot of gain verses slice and a histogram of the values. These can be used to determine if the gain is constant throughout the image and so is a good estimate.

.. index:: em-gain analysis

EM-Gain Analysis
----------------

Analyses a white light image from an EM-CCD camera, construct a histogram of pixel intensity and fit the histogram to obtain the bias, EM-gain, read noise and photons per pixel (see Ulbrich & Isacoff (2007) Supplementary Information).

.. index:: em-ccd probability model

EM-CCD Probability Model
~~~~~~~~~~~~~~~~~~~~~~~~

The ``EM-Gain Analysis`` plugin uses an analysis that assumes that the EM-CCD camera has three main sources of noise:

#.  Photon shot noise occurs when light is emitted from an object. Although the average rate of light from an object is constant for a given time, e.g. 30 photons/second, each photon will arrive at a different time and the gaps between them will vary. This results in a different number of photons counted each second. This noise follows a Poisson distribution with a mean of the average photon emission rate.

#.  The photons are converted to electrons on the camera sensor. These electrons are then multiplied in the Electron Multiplication (EM) gain register. This multiplication increases the number of electrons to be read and reduces the relative size of any error introduced when reading the value. However the EM-gain process is random and introduces noise that is modelled using a Gamma distribution with a shape parameter equal to the number of input electrons and the scale parameter equal to the gain.

#.  Read noise occurs when the values stored on the camera chip for each pixel are read and converted to numbers. This noise follows a Gaussian distribution with mean zero and variable standard deviation.

The probability of observing a pixel value given an input number of photons is therefore a convolution of a Poisson, Gamma and Gaussian distribution. The convolution of the Poisson and Gamma distribution can be expressed as:

.. math::

    G_{p,m}(c)=\operatorname{e}^{-p}\delta(c)+\sqrt{\frac{p}{\mathit{cm}}}\operatorname{e}^{-{\frac{c}{m}}-p}\mathit{BesselI}_{1}(2\sqrt{\frac{\mathit{cp}}{m}})

where
:math:`p` is the average number of photons,
:math:`m` is the EM-gain multiplication factor,
:math:`c` is the observed pixel count,
:math:`\delta(c)` is the Dirac delta function (1 when c=0, 0 otherwise),
:math:`\mathit{BesselI}_1` is the modified Bessel function of the 1\ :sup:`st` kind, and
:math:`G_{p,m}(c)` is the probability of observing the pixel count c.

The output of this is subsequently convolved numerically (no algebraic solution exists) with a Gaussian function with standard deviation equal to the camera read noise and mean equal to the camera bias.

.. index:: camera bias

Camera Bias
^^^^^^^^^^^

Note that in order to observe the read noise of the camera a bias (offset) is added to the camera pixel values. This allows a pixel to record negative read noise on very low counts which would not be possible using unsigned integer values as no value below zero is allowed. The bias for the camera is set by the manufacturer and is set at a value far greater than the expected read noise of the system, e.g. 100, 400, 500 or 1000 for a read noise of 3-30 ADUs (Analogue to Digital Units, or pixel values).

.. index:: input image

Input image
~~~~~~~~~~~

The plugin requires a white light image where each pixel has been exposed to the same number of photons. This can be produced by imaging without a sample and instead using white paper in front of the objective so that images are evenly illuminated. The light can be adjusted by varying the exposure time and different calibration performed by using different camera gain settings.

The input image is used to construct a histogram of the pixel values that are observed for the given camera settings and background number of photons. This is then fit using the Poisson-Gamma-Gaussian probability mass function.

Ideally the input image should provide a minimum of 1,000,000 pixels, for example 16 frames of a 256x256 pixel image. This level of pixels is required to construct an even histogram that adequately samples the probability mass function. The pixels should have the same mean, i.e. a constant mean across the field of view. If it is not possible to achieve a constant mean across the field, for example in the instance of a gradient in the illumination, then the plugin will support rectangular ROI crops of the image. However the number of pixels should reach the minimum limit to construct a good histogram.

If the minimum pixel limit is not reached the plugin will log a warning but will continue to analyse the image.

.. index:: parameters

Parameters
~~~~~~~~~~

The following parameters can be configured:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Bias estimate
     - The initial estimate for the camera bias. The bias may be obtained from the camera manufacturer's specifications. A guess can be made by selecting the darkest part of the image, taking the mean and rounding (usually down) to the nearest hundred.

   * - Gain estimate
     - The initial estimate for the total gain of the camera. The total gain may be obtained from the camera manufacturer's specifications. A good guess would be 25-50.

   * - Noise estimate
     - The initial estimate for the camera read noise. The read noise in electrons may be obtained from the camera manufacturer's specifications. This will have to be converted to ADUs by applying the camera gain (not the total gain). A good guess would be 3-10.

   * - Show approximation
     - Show on the final output plot a function that approximates the convolution of the Poisson-Gamma distribution with a Gaussian distribution.

       This approximate PMF is used to model the EM-Gain when performing Maximum Likelihood Estimation fitting within the ``Peak Fit`` plugin.


Note that the plugin will remember the last values that were fitted for the bias, gain and noise estimates. Thus an initial guess can be used, the image analysed and then the plugin repeated with updates to the estimates if appropriate to refine the fit.

.. index:: simulation mode

Simulation Mode
~~~~~~~~~~~~~~~

Instead of using an input image to create a histogram of pixel values, it is possible to simulate pixel values by generating a Poisson-Gamma-Gaussian random variable. To run the plugin in simulation mode hold down the ``Shift`` key when running the plugin. The following additional parameters will be available:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Simulate
     - Check this box to simulate the histogram of pixel values.

   * - Bias
     - The camera bias for the simulation.

   * - Gain
     - The total gain for the simulation.

   * - Noise
     - The read noise for the simulation.

   * - Photons
     - The average number of photons per pixel for the simulation.

   * - Samples
     - The number of samples for the simulation.

   * - Sample PDF
     - Check this to generate the Probability Mass Function (PMF) using the provided parameters. Then sample randomly from within this PMF.

       The default is to generate a random Poisson sample using the average photon number, then use this to generate a Gamma sample from the photon count and then generate a Gaussian sample from the amplified photon count.


Simulation mode can be used to see if the fitting process is working given the expected parameters for bias, gain, noise and photons.

.. index:: results

Results
~~~~~~~

The plugin will create a histogram of the pixel values and attempt to fit it using the Poisson-Gamma-Gaussian PMF. The estimated and fitted parameters are written to the ``ImageJ`` log.

The histogram of pixel values, fitted PMF and the fit parameters are shown on a plot (:numref:`Figure %s <fig_em_gain_analysis_histogram_fit>`).

.. _fig_em_gain_analysis_histogram_fit:
.. figure:: images/em_gain_analysis_histogram_fit.png
    :align: center
    :figwidth: 80%

    EM-Gain Analysis histogram of pixel values and the computed fit

The values for the gain, bias and noise should be constant for different background photon levels. This can be evaluated using different input calibration images. The parameters can be used within the ``Peak Fit`` plugin to perform Maximum Likelihood Estimation modelling the camera noise of the EM-CCD camera.

.. index:: em-gain pmf

EM-Gain PMF
-----------

Displays a plot of the probability mass function (PMF) of the expected value of a pixel on an EM-CCD camera given an average number of photons. The form of the PMF is a convolution of a Poisson, Gamma and Gaussian distribution. See section :numref:`{number}: {name} <calibration_plugins:EM-CCD Probability Model>` for more details.

A fast approximation for the PMF is computed for comparison with the real PMF. This is created by analytically calculating the PMF of a Poisson-Gamma distribution and then approximating a convolution with a Gaussian distribution. The method for this approximation is taken from the supplementary Python software provided by Mortensen, *et al* (2010). They used this approximation when fitting the images of single fluorophores in TIRF (Total Internal Reflection Fluorescence) images taken with an EM-CCD camera. A second plot showing the difference between the real PMF and the approximation is displayed. This allows investigation of any situation where the approximation is not appropriate for modelling the PMF. It is rare for the approximation to differ by more than 1%.

The plugin has the following parameters:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Gain
     - The total gain for EM-CCD camera.

   * - Noise
     - The camera read noise.

   * - Photons
     - The average number of photons per pixel for the simulation.

   * - Show approximation
     - Show on the PMF plot the approximation function.

       Note: This approximate PMF is used to model the EM-Gain when performing Maximum Likelihood Estimation fitting within the ``Peak Fit`` plugin.

   * - Remove head
     - Set a limit on the initial cumulative probability to remove from the plot. This allows removing the start of the curve where the convolution of the Poisson-Gamma distribution with the Gaussian is incomplete.

   * - Remove tail
     - Set a limit on the final cumulative probability to remove from the plot. This allows removing the tail of the curve where the convolution of the Poisson-Gamma distribution with the Gaussian is incomplete. It also allows removing the long tail which can take up a large amount of the plot width.

   * - Relative delta
     - Check this to show the difference between the actual PMF and the approximate PMF as a relative score. The default is absolute.


Examples of the PMF are shown in :numref:`Table %s <table_em_gain_pmf>`. The PMF is skewed for low photons with a spike at c=0 blurred by the Gaussian read noise. Increasing photon counts return a shape more characteristic of a Poisson distribution. For this reason it is possible to use a simple Poisson model for the camera noise when performing Maximum Likelihood Estimation, i.e. ignoring the effect of EM-gain and read noise, if the number of photons within the localisation is large. This is an option available within the
``Peak Fit``
plugin and allows much faster fitting since the Poisson PMF (a) can be evaluated much faster than the Poisson-Gamma-Gaussian PMF; and (b) has an analytical derivative allowing gradient based fitting methods.

.. _table_em_gain_pmf:
.. list-table:: Example EM gain probability mass function (PMF) plots
    :align: center
    :width: 80

    * - |em_gain_pmf_1_png|
    * - |em_gain_pmf_2_png|
    * - |em_gain_pmf_3_png|
    * - The magenta line on the plot shows the position of the average number of photons after the gain has been applied.

.. |em_gain_pmf_1_png| image:: images/em_gain_pmf_1.png
.. |em_gain_pmf_2_png| image:: images/em_gain_pmf_2.png
.. |em_gain_pmf_3_png| image:: images/em_gain_pmf_3.png


.. index:: diffusion rate test

Diffusion Rate Test
-------------------

The ``Diffusion Rate Test`` plugin will simulate molecule diffusion and fit a graph of mean-squared displacement to determine the diffusion coefficient. This is a test plugin to show that the simulated diffusion in the ``Create Data`` plugin generates correct moving particles.

When a molecule is diffusing it can move in any direction. The total distance it moves and the track it took may not be visible due to the speed of movement. However the diffusion of particles in a single dimension can be modelled as a population. If the squared distances from the origin after a set time are plotted as a histogram they can be modelled using a Gaussian curve. The average distance the particles will move is zero and the variance of the Gaussian curve will be the mean-squared displacement (MSD). This can be expressed by unit time. The MSD is proportional to the diffusion coefficient (:math:`D`). The relationship for a single-dimension is MSD = :math:`2D`. This increases to :math:`4D` and :math:`6D` for two and three dimensional distances (since the diffusion in each dimension is independent).

.. index:: grid walk simulation

Grid Walk simulation
~~~~~~~~~~~~~~~~~~~~

Since the MSD in a single dimension is equal to :math:`2D`, the mean-distance a particle moves will be :math:`\sqrt{2D}`. This step size can be used to simulate diffusion using a grid walk. At each step a particle can move forward or backwards by the step size :math:`s`. If the direction is random then the population of particles will have an average displacement of zero, a mean displacement of the step size:math:`s`, and a mean squared displacement (MSD) of :math:`s^2 = 2D`. Multi-dimension diffusion is done by simulating the movement in each dimension separately.

.. index:: random move simulation

Random Move simulation
~~~~~~~~~~~~~~~~~~~~~~

Diffusion can also be simulated by moving particles on a random vector. The distance moved should be sampled from a Gaussian distribution with a mean of zero and a standard deviation of :math:`\sqrt{\mathit{MSD}}`. This is :math:`2D`, :math:`4D` or :math:`6D` for 1, 2 or 3 dimensions respectively.

However the unit vector must be directed in a random orientation. For one dimension this is either forward or backward. For higher dimensions a random vector can be produced by sampling the movement in each dimension from a Gaussian distribution with mean zero and standard deviation 1. This vector is normalised to unit length.

The generation of the unit vector and the movement distance can be combined into a single stage. The random displacement is produced by sampling each dimension from a Gaussian distribution with mean zero and standard deviation of :math:`\sqrt{2D}`. This is the equivalent of 1-dimension diffusion in 3 independent dimensions.

.. index:: confined diffusion

Confined Diffusion
~~~~~~~~~~~~~~~~~~

Particles may not be able to freely move in any direction, for example when they collide with a barrier. The ``Diffusion Rate Test`` plugin allows particles to be confined in a sphere. In this case the diffusion step is calculated and if the step would move the particle outside the sphere the move is rejected. Attempts are made to move the particle a set number of times until successful otherwise the particle coordinates are not updated. This simulation produces good results when the average step size is at least an order of magnitude less than the sphere radius (so allowing many steps inside the sphere to be valid) and the ``Randon Move`` simulation is used.

.. index:: analysis

Analysis
~~~~~~~~

The
``Diffusion Rate Test``
plugin simulates the random diffusion of many particles over a period of time. Each diffusion path is then analysed. The plugin has the following parameters:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Pixel pitch (nm)
     - The pixel size for the simulation.

   * - Seconds
     - The duration of the simulation.

   * - Steps per second
     - The number of diffusion steps the particle makes per second.

   * - Particles
     - The number of particles to simulate.

   * - Diffusion rate
     - The diffusion coefficient (D).

   * - Use grid walk
     - If **true** then simulate diffusion using a grid walk, otherwise use a random move. The grid walk simulation is approximately 3 times faster.

   * - Use confinement
     - If **true** then use a sphere to confine the particle movement.

   * - Confinement attempts
     - The number of times to attempt a confined move.

   * - Confinement radius
     - The radius of the confinement sphere.

   * - Fit N
     - When using confined diffusion only fit the first N points of the MSD plot.

   * - Show example
     - Show an example image of a diffusion path.

   * - Magnification
     - The magnification of the example image. The pixels will represent (pixel pitch) / magnification nanometres.


.. index:: output

Output
~~~~~~

.. index:: msd plot

MSD plot
^^^^^^^^

For each particle the plugin will compute the squared displacement from the origin over the time course of the simulation. Distances are computed in 2D and 3D. The plugin will plot the mean-squared distance against time for the population as shown in :numref:`Figure %s <fig_diffusion_rate_msd_plot>`. The 2D and 3D MSD data are then fit using a linear regression. The gradient of the fit can be used to calculate the diffusion coefficient by dividing by 4 or 6 respectively.

.. _fig_diffusion_rate_msd_plot:
.. figure:: images/diffusion_rate_msd_plot.png
    :align: center
    :figwidth: 80%

    Mean-squared distance (MSD) plot

    The plot shows the 2D (black) and 3D (magenta) MSD with fitted line. The upper and lower bounds for 2D MSD are shown (blue).

If confined diffusion is performed the MSD will reach a natural upper limit. This will result in a plateau of the MSD plot as shown in :numref:`Figure %s <fig_diffusion_rate_msd_plot_confined>`. In this case only the initial diffusion of the particles will be unconstrained. The analysis should therefore fit the initial linear section of the MSD plot. If the confinement radius is too small there may be no linear section to the MSD curve.

.. _fig_diffusion_rate_msd_plot_confined:
.. figure:: images/diffusion_rate_msd_plot_confined.png
    :align: center
    :figwidth: 80%

    Mean-squared distance (MSD) plot for confined diffusion

    The plot shows the 2D (black) and 3D (magenta) MSD with fitted line using the initial linear section of data. The upper and lower bounds for 2D MSD are shown (blue).

**Note:** The asymptote of the curve for confined diffusion should be defined by the average distance to the centre of a random distribution of particles within a sphere. This can be computed using the distance from the centre of all the points in a sphere divided by the number of points in a sphere. The surface area (*SA*) of a sphere is equal to the number of points at distance r from the centre. So :math:`\mathit{SA} \times r` is the sum of the distances of points at distance r from the centre. If this is integrated from zero to *R* it produces the sum of all distances from any point within a sphere of radius *R*. The number of points is the volume (*V*) of the sphere. The sum of the distances divided by the number of points is the average distance to the centre, therefore:

.. math::

    \frac{\int ^{R}\mathit{SA} \times r}{V} dr &= \frac{\int ^{R}4\pi r^{2} \times r}{4\pi R^{3}/3} dr \\
    &= \frac{\int ^{R} r^{3}}{R^{3}/3} dr \\
    &= \frac{R^{4}}{R^{3}/3} \\
    &= \frac{3R}{4}

Thus the mean-distance to the centre for particles in a sphere is 0.75 *R*. This can be used to check that the confined simulation is performing as a true random diffusion within a sphere.

.. index:: diffusion example

Diffusion example
^^^^^^^^^^^^^^^^^

If the ``Show example`` option was selected the plugin will show an image of the track of a single particle. The track is shown on a black background. The track is initialised at a value of 32 and ends with a value of 255. The movement can thus be followed using a colour lookup table (LUT), e.g. ``Image > Lookup Tables > Fire``.

The plugin will also show a plot of the displacement of the particle over time. The red line shows the X displacement and the blue shows the Y displacement.

.. index:: analysis results

Analysis results
^^^^^^^^^^^^^^^^

The fitting analysis results are output to the ``ImageJ`` log window, e.g.

.. code-block:: text

    Diffusion Rate Test : D = 1.0 um^2/sec, Precision = 0.0 nm
    Mean-displacement per dimension = 1414.0 nm/sec
    Simulation step-size = 44.72 nm
    Raw data D=1.0 um^2/s, Precision = 0.0 nm, N=22000, step=0.001 s, mean=0.004034 um^2, MSD = 4.034 um^2/s
    2D Diffusion rate = 1.013 um^2 / sec (50.22 ms)
    3D Diffusion rate = 1.175 um^2 / sec (50.22 ms)

The input diffusion coefficient is shown for reference, the units are |micro|\ m\ :sup:`2`/sec. This is converted to the expected mean-displacement in nm per second and the simulation step size (in nm). This will allow the user to experiment with the radius of the confinement sphere and the number of simulation steps. Remember that the step size should be less than the sphere radius when using confinement. The fitted diffusion coefficients from the 2D and 3D fitting are then shown. These should be close to the input diffusion rate.

If the simulation was performed using confinement then the final distance to the origin for each particle will be saved. The average distance will be shown along with the expected asymptote distance, i.e. the mean-distance to the centre of a sphere, which is calculated as 3/4 of the confinement radius, e.g.

.. code-block:: text

    3D asymptote distance = 702.7 nm (expected 750.00)

.. index:: memory results

Memory Results
^^^^^^^^^^^^^^

The coordinates of each diffusing particle, starting at the origin (0,0), are saved to a results dataset in memory. Each consecutive step of the same particle is given a new frame and particles are allocated a unique ID. The current frame is incremented between particles so that each particle track is separated in time. This allows the results set to be used within the ``Trace Diffusion`` and ``Draw Clusters`` plugins to verify their functionality.

.. index:: extra options

Extra options
~~~~~~~~~~~~~

Hold the ``Shift`` key down when running the plugin to activate extra options. The following options are available and are described in the following sections:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Aggregate steps
     - Create a second dataset by averaging N consecutive positions into a single location.

   * - MSD analysis
     - Specify the maximum number of steps used to perform MSD analysis. This is only relevant when the ``Aggregate steps`` parameter is above 1.

   * - Precision
     - Specify the localisation precision of positions.


.. index:: aggregate steps parameter

Aggregate steps parameter
~~~~~~~~~~~~~~~~~~~~~~~~~

The standard plugin simulates diffusion in small steps. These can be aggregated together to simulate the position of the particle in a frame taken on a camera. In this case the average position of a set of consecutive steps is calculated to aggregated the position into a frame. The mean-squared distance between frames is then reported to the ``ImageJ`` log:

.. code-block:: text

    Raw data D=1.0 um^2/s, Precision = 0.0 nm, N=22000, step=0.001 s, mean=0.004016 um^2, MSD = 4.016 um^2/s
    Aggregated data D=1.0 um^2/s, Precision=0.0 nm, N=200, step=0.1 s, mean=0.2713 um^2, MSD = 2.713 um^2/s

Note that the aggregation has the effect of reducing the mean-squared displacement for the dataset.

The aggregated data is saved into a dataset in memory.

.. index:: msd analysis parameter

MSD Analysis parameter
~~~~~~~~~~~~~~~~~~~~~~

The ``MSD Analysis`` option is available when data aggregation has been performed using the ``Aggregate steps`` parameter. When multiple small steps are aggregated into single coordinates this causes the observed MSD to be lower than the expected MSD given the diffusion coefficient. Effectively the averaging of the position of a particle within a frame has caused loss of information about the diffusion distance covered within that frame. MSD analysis allows the effect of aggregation to be analysed.

For each simulated track the position of the particle is computed as a rolling average of the coordinates using the configured number of ``Aggregated steps`` (N).

For example the first position is the average of the first N steps. The next position is computed by adding the N+1 coordinate and subtracting the 1\ :sup:`st` coordinate from the sum to create a new average.

Using the rolling average position the squared distance between each position and the *j*\ :sup:`th` position along the sequence is computed for all *j* up to the limit set by the ``MSD Analysis`` parameter.
The mean squared distance is then reported for each separation *j* to a summary table:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - **Column**
     - Description

   * - D
     - The input diffusion coefficient.

   * - Precision
     - The precision error of each average position.

       If positive the positions will be adjusted with a random error before the distances are computed. See section :numref:`{number} <calibration_plugins:Precision parameter>`.

   * - Dsim
     - The simulated diffusion coefficient (calculated using the mean-squared displacement between non-aggregated coordinates).

   * - Step
     - The size of a single simulation step in seconds.

   * - Resolution
     - The number of steps that are aggregated into a frame.

   * - Frame
     - The frame length in seconds.

   * - t
     - The separation *j* between the two rolling average positions in seconds.

   * - n
     - The separation *j* between the two rolling average positions in frames.

   * - N
     - The number of samples to compute the MSD.

   * - MSD
     - The mean-squared displacement (MSD).

   * - D
     - The observed diffusion coefficient (calculated as MSD / 4t).


.. index:: precision parameter

Precision parameter
~~~~~~~~~~~~~~~~~~~

By default the exact coordinates of a particle are used in the analysis and to create the output datasets. To simulate the results generated by a super-resolution image reconstruction method the coordinates can be reported with a random error added to each position. The error simulates the fitting precision of the super-resolution localisation method. Error is added independently to the X and Y coordinates using a Gaussian random variable with the given standard deviation.

The precision error has the effect of increasing the mean-squared displacement for the dataset, e.g.

.. code-block:: text

    Raw data D=1.0 um^2/s, Precision = 0.0nm, N=22000, step=0.001 s, mean=0.004016 um^2, MSD = 4.016 um^2/s
    …
    Raw data D=1.0 um^2/s, Precision = 30.0 nm, N=22000, step=0.001 s, mean=0.007614 um^2, MSD = 7.614 um^2/s

.. index:: trace diffusion

Trace Diffusion
---------------

The ``Trace Diffusion`` plugin will trace molecules through consecutive frames and then perform mean-squared displacement analysis to calculate a diffusion coefficient.

The plugin is similar to the ``Diffusion Rate Test`` plugin however instead of simulating particle diffusion the plugin will use an existing results set. This allows the analysis to be applied to results from fitting single-molecule images using the ``Peak Fit`` plugin.

.. index:: analysis

Analysis
~~~~~~~~

The plugin runs a tracing algorithm on the results to find localisations that occur in consecutive frames. Details of the tracing algorithm can be found in section :numref:`{number}: {name} <analysis_plugins:Trace Molecules>`. The distance threshold for the tracing algorithm can be specified but the time threshold is set to 1 frame, i.e. only continuous tracks will be extracted. Thus a pair of localisations within adjacent frames will be connected if they are within the distance threshold. In addition the plugin allows the track to be excluded if a second localisation occurs within an exclusion threshold of the first localisation. This effectively removes traces of particles that could overlap with another moving particle.

Once the tracks have been identified the tracks are filtered using a length criteria and shorter tracks discarded. Optionally the tracks can be truncated to the minimum length which ensures even sampling of particles with different track lengths. The plugin computes the mean-squared distance of each point from the origin. Optionally the plugin computes the mean-squared distance of each point from every other point in the track. These internal distances increase the number of points in the analysis. Therefore if the track is not truncated the number of internal distances at a given time separation is proportional to the track length. To prevent bias in the data towards the longer tracks the average distance for each time separation is computed per track and these are used in the population statistics. Thus each track contributes only once to the mean-displacement for a set time separation.

The mean-squared distance (MSD) per molecule is calculated using two methods. The ``all-vs-all`` method uses the sum of squared distances divided by the sum of time separation between points. The value includes the all-vs-all internal distances (if selected). The ``adjacent`` method uses the average of the squared distances between adjacent frames divided by the time delta (:math:`\Delta t`) between frames. The MSD values are expressed in |micro|\ m\ :sup:`2`/second and can be saved to file or shown in a histogram.

The average mean-squared distances for all the traces are plotted against the time separation and a best fit line is calculated. The mean-squared distances are proportional to the diffusion coefficient (*D*):

.. math::

    \mathit{MSD}(n\Delta t)=4\mathit{Dn}\Delta t+4\sigma ^{2}

where
:math:`n` is the number of separating frames,
:math:`\Delta t` is the time lag between frames, and
:math:`\sigma` is the localisation precision.
Thus the gradient of the best fit line can be used to obtain the diffusion coefficient. Note that the plugin will compute a fit with and without an explicit intercept and pick the solution with the best fit to the data (see :numref:`{number}: {name} <calibration_plugins:Selecting the best fit>`).
Note that an additional best fit line can be computed using a MSD correction factor
(see :numref:`{number}: {name} <calibration_plugins:MSD Correction>`).

.. index:: apparent diffusion coefficient

Apparent Diffusion Coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given that the localisations within each trace are subject to a fitting error, or precision (:math:`\sigma`), the apparent diffusion coefficient (:math:`D^{\star}`) can be calculated accounting for precision [Uphoff *et al*, 2013]:

.. math::

    D^{\star}=\mathit{max}(0,\frac{\mathit{MSD}}{4n\Delta t}-\frac{\sigma_{\mathit{loc}}^{2}}{n\Delta t})

The plugin thus computes the average precision for the localisations included in the analysis and can optionally report the apparent diffusion coefficient (:math:`D^{\star}`). If the average precision is above 100nm then the plugin prompts the user to confirm the precision value.

.. index:: jump distance analysis

Jump Distance analysis
~~~~~~~~~~~~~~~~~~~~~~

The jump distance is how far a particle moves is given time period. Analysis of a population of jump distances can be used to determine if the population contains molecules diffusing with one or more diffusion coefficients [Weimann *et al*, 2013]. For two dimensional Brownian motion the probability that a particle starting at the origin will be encountered within a shell of radius *r* and a width *dr* at time :math:`\Delta t` is given by:

.. math::

    p(r^{2},\Delta t)\mathit{dr}^{2}=\frac{1}{4D\Delta t}e^{-{\frac{r^{2}}{4D\Delta t}}}\mathit{dr}^{2}

This can be expanded to a mixed population of *m* species where each fraction (:math:`f_i`) has a diffusion coefficient :math:`D_i`\ :

.. math::

    p(r^{2},\Delta t)\mathit{dr}^{2}=\sum_{j=1}^{m}{\frac{f_{j}}{4D_{j}\Delta t}e^{-{\frac{r^{2}}{4D_{j}\Delta t}}}\mathit{dr}^{2}}

For the purposes of fitting the integrated distribution can be used. For a single population this is given by:

.. math::

    P(r^{2},\Delta t)=\int_{0}^{r^{2}}p(r^{2})\mathit{dr}^{2}=1-e^{-{\frac{r^{2}}{4D\Delta t}}}

The advantage of the integrated distribution is that specific histogram bin sizes are not
required to construct the cumulative histogram from the raw data. Note that the
integration holds for a mixed population of *m* species where each fraction (:math:`f_i`) has a diffusion coefficient :math:`D_i`\ :

.. math::

    P(r^{2},\Delta t)=1-\sum _{j=1}^{m}f_{j}e^{-{\frac{r^{2}}{4D_{j}\Delta t}}}

Weimmann *et al* (2013) show that fitting of the cumulative histogram of jump distances can
accurately reproduce the diffusion coefficient in single molecule simulations. The performance of the method uses an indicator :math:`\beta` expressed as the average distance a particle travels in the chosen time (*d*) divided by the average localisation precision (:math:`\sigma`):

.. math::

    \beta = d / \sigma

When :math:`\beta` is above 6 then jump distance analysis reproduces the diffusion coefficient as accurately as MSD analysis for single populations. For mixed populations of moving and stationary particles the MSD analysis fails (it cannot determine multiple diffusion coefficients) and the jump distance analysis yields accurate values when :math:`\beta` is above 6.

The ``Trace Diffusion`` plugin performs jump distance analysis using the jumps between frames that are *n* frames apart. The distances may be from the origin to the *n* th frame or may use all the available internal distances *n* frames apart. A cumulative histogram is produced of the jump distance. This is then fitted using a single population and then for mixed populations of *j* species by minimising the sum-of-squared residuals (SS) between the observed and expected curves. Alternatively the plugin can fit the jump distances directly without using a cumulative histogram. In this case the probability of each jump distance is computed using the formula for :math:`P(r^{2},\Delta t)` and the combined probability (likelihood) of the data given the model is computed. The best model fit is achieved by maximising the likelihood (maximum likelihood estimation, MLE).

When fitting multiple species the fit is rejected if:
(a) the relative difference between coefficients is smaller than a given factor; or
(b) the minimum fraction, :math:`f_i`, is less than a configured level.
If accepted the result must then be compared to the previous result to determine if increasing the number of parameters has improved the fit (see :numref:`{number}: {name} <calibration_plugins:Selecting the best fit>`).


Optimisation is performed using a fast search to maximise the score by varying each parameter in turn (Powell optimiser). In most cases this achieves convergence. However in the case that the default algorithm fails then a second algorithm is used that uses a directed random walk (CMAES optimiser). This algorithm depends on randomness and so can benefit from restarts. The plugin allows the number of restarts to be varied. For the optimisation of the sum-of-squares against the cumulative histogram a least-squares fitting algorithm (Levenberg-Marquardt or LVM optimiser) is used to improve the initial score where possible. The plugin will log messages on the success of the optimisers to the ``ImageJ`` log window. Extra information will be logged if using the ``Debug fitting`` option.

.. index:: msd correction

MSD Correction
~~~~~~~~~~~~~~

This corrects for the diffusion distance lost in the first and last frames of the track due to the representation of diffusion over the entire frame as an average coordinate.
A full explanation of the correction is provided in section :numref:`{number}: {name} <msd_correction:MSD Correction>`.

The observed MSD can be converted to the true MSD by dividing by a correction factor (*F*):

.. math::

    F=\frac{n-1/3}{n}

Where *n* is the number of frames over which the jump distance is measured (i.e. end - start).

When performing jump distance analysis it is not necessary to the correct each observed squared distance before fitting. Since the correction is a single scaling factor instead the computed diffusion coefficient can be adjusted by applying the correction factor after fitting. This allows the plugin to save the raw data to file and use for display on results plots.

If the ``MSD correction`` option is selected the plugin will compute the corrected diffusion coefficient as:

.. math::

    D_{\mathit{corr}}=D\cdot {\frac{n}{n-1/3}}

.. index:: fitting the plot of msd verses n frames

Fitting the plot of MSD verses N frames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When fitting the linear plot of MSD verses the number of frames the correction factor can be included. The observed MSD is composed of the actual MSD multiplied by the correction factor before being adjusted for the precision error:


.. math::

    \mathit{oMSD}(n\Delta t)=4D(n\Delta t)-\frac{4D(\Delta t)}{3}+4\sigma^{2}

This is still a linear fit with a new representation for the intercept that allows the intercept to be negative. To ensure the intercept is correctly bounded it is represented using the fit parameters and not just fit using a single constant C.

When performing the linear fit of the MSD verses jump distance plot, 3 equations are fitted and the results with the best information criterion is selected. The results of each fit are written to the ``ImageJ`` log. The following equations are fit:

Linear fit:

.. math::

    \mathit{oMSD}(n\Delta t)=4D(n\Delta t)

Linear fit with intercept:

.. math::

    \mathit{oMSD}(n\Delta t)=4D(n\Delta t)+4\sigma ^{2}

Linear fit with MSD corrected intercept:

.. math::

    \mathit{oMSD}(n\Delta t)=4D(n\Delta t)-\frac{4D(\Delta t)}{3}+4\sigma^{2}

Note: In each model the linear gradient is proportional to the diffusion coefficient.

.. index:: precision correction

Precision Correction
~~~~~~~~~~~~~~~~~~~~

Given that the localisations within each trace are subject to a fitting error, or precision (σ), the apparent diffusion coefficient (:math:`D^{\star}`) can be calculated accounting for precision [Uphoff *et al* , 2013]:

.. math::

    D^\star=\mathit{max}(0,\frac{\mathit{MSD}}{4n\Delta t}-\frac{\sigma _{\mathit{loc}}^{2}}{n\Delta t})

If the ``Precision correction`` option is selected the plugin will subtract the precision and report the apparent diffusion coefficient (:math:`D^{\star}`) from the jump distance analysis.

.. index:: msd and precision correction

MSD and Precision correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both the ``MSD correction`` and ``Precision correction`` can be applied to the fitted MSD to compute the corrected diffusion coefficient:

.. math::

    D=\frac{n}{n-1/3}\cdot \mathit{max}(0,\frac{\mathit{MSD}}{4n\Delta t}-\frac{\sigma _{\mathit{loc}}^{2}}{n\Delta t})

.. index:: selecting the best fit

Selecting the best fit
~~~~~~~~~~~~~~~~~~~~~~

The Bias Corrected Akaike Information Criterion (cAIC) [Hurvich & Tsai, 1989] is calculated for the fit using the log likelihood (*L*), the number of data points (*n*) and the number of parameters (*p*):

.. math::

   \mathit{AIC}=2p-2L

.. math::

    \mathit{cAIC}=\mathit{AIC}+(2p(p+1)/(n-p-1))

The corrected AIC penalises additional parameters. The model with the lowest cAIC is preferred. If a higher cAIC is obtained then increasing the number of fitted species in the mixed population has not improved the fit and so fitting is stopped. Note that when performing maximum likelihood estimation the log likelihood (*L*) is already known and is used directly to calculate the corrected AIC. When fitting the sum-of-squared residuals (SS) the log likelihood can be computed as:

.. math::

    L=-{\frac{n}{2}}\ln (2\pi )-\frac{n}{2}\ln(\sigma ^{2})-\frac{1}{2\sigma ^{2}}\mathit{SS}}

.. index:: parameters

Parameters
~~~~~~~~~~

The plugin dialog allowing the data to be selected is shown in :numref:`Figure %s <fig_trace_diffusion_dialog>`.

.. _fig_trace_diffusion_dialog:
.. figure:: images/trace_diffusion_dialog.png
    :align: center
    :figwidth: 80%

    Trace Diffusion dialog

The plugin has the following parameters:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Input
     - Specify the input results set.

   * - Distance threshold
     - The distance threshold for tracing.

   * - Distance exclusion
     - The exclusion distance. If a particle is within the distance threshold but a second particle is within the exclusion distance then the trace is discarded (due to overlapping tracks).

   * - Min trace length
     - The minimum length for a track (in time frames).

   * - Ignore ends
     - Ignore the end jumps in the track.

       If a fluorophore activated only part way through the first frame and bleaches only part way through the last frame the end jumps represent a shorter time-span than the frame interval. These jumps can optionally be ignored.

       This option requires tracks to be 2 frames longer than the ``Min trace length`` parameter.

   * - Save traces
     - Save the traces to file in the ``Peak Fit`` results format.


When all the datasets have been traced the plugin presents a second dialog to configure the diffusion analysis. The following parameters can be configured:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Truncate traces
     - Set to to true to only use the first N points specified by the ``Min trace length`` parameter.

   * - Internal distances
     - Compute the all-vs-all distances. Otherwise only compute distance from the origin.

   * - Fit length
     - Fit the first N points with a linear regression.

   * - MSD correction
     - Perform mean square distance (MSD) correction.

       This corrects for the diffusion distance lost in the first and last frames of the track due to the representation of diffusion over the entire frame as an average coordinate.

   * - Precision correction
     - Correct the fitted diffusion coefficient using the average precision of the localisations.

       Note that uncertainty in the position of localisations (fit precision) will contribute to the displacement between localisations. This can be corrected for by subtracting :math:`4s^2` from the measured squared distances with *s* the average precision of the localisations.

   * - Maximum likelihood
     - Perform jump distance fitting using maximum likelihood estimation (MLE). The default is sum-of-squared residuals (SS) fitting of the cumulative histogram of jump distances.

   * - Fit restarts
     - The number of restarts to attempt when fitting using the CMAES optimiser. A higher number produces and more robust fit solution since the best fit of all the restarts is selected.

       Note that the CMAES optimiser is only used when the default Powell optimiser fails to converge.

   * - Jump distance
     - The distance between frames to use for jump analysis.

   * - Minimum difference
     - The minimum relative difference (ratio) between fitted diffusion coefficients to accept the model. The difference is calculated by ranking the coefficient in descending order and then expressing successive pairs as a ratio. Models with coefficients too similar are rejected.

   * - Minimum fraction
     - The minimum fraction of the population that each species must satisfy. Models with species fractions below this are rejected.

   * - Minimum N
     - The minimum number of species to fit. This can be used to force fitting with a set number of species.

       This extra option is only available if the plugin is run with the ``Shift`` key held down, otherwise the default is 1.

   * - Maximum N
     - The maximum number of species to fit. In practice this number may not be achieved if adding more species does not improve the fit.

   * - Debug fitting
     - Output extra information to the ``ImageJ`` log window about the fitting process.

   * - Save trace distances
     - Save the traces to file. The file contains the per-molecule MSD and D* and the squared distance to the origin for each trace.

   * - Save raw data
     - Select this to select a results directory where the raw data will be saved. This is the data that is used to produce all the histograms and output plots.

   * - Show histograms
     - Show histograms of the trace data. If selected a second dialog is presented allowing the histograms to be chosen and the number of histogram bins to be configured.

   * - Title
     - A title to add to the results table.


.. index:: output

Output
~~~~~~

.. index:: msd verses time

MSD verses time
^^^^^^^^^^^^^^^

The plugin will plot the mean-squared distances against the time as show in :numref:`Figure %s <fig_trace_diffusion_msd_vs_time>`. The plot shows the best fit line. If the data is not linear then the diffusion of particles may be confined, for example by cellular structures when using *in vivo* image data. In this case the diffusion coefficient will be underestimated.

.. _fig_trace_diffusion_msd_vs_time:
.. figure:: images/trace_diffusion_msd_vs_time.png
    :align: center
    :figwidth: 80%

    Plot of mean-squared distance verses time produced by the Trace Diffusion plugin.

    The mean of the raw data is plotted with bars representing standard error of the mean. The best fit line is shown in magenta.

.. index:: jump distance histogram

Jump distance histogram
^^^^^^^^^^^^^^^^^^^^^^^

The plugin produces a cumulative probability histogram of the jump distance (see :numref:`Figure %s <fig_trace_diffusion_jump_distance_cumul_histogram>`). The best fit for a single species model will be shown in magenta. Any significant deviations of the histogram line from the single species fit are indicative of a multi-species population. If a multiple species model has a better fit than the single species model then it will be plotted in yellow.

.. _fig_trace_diffusion_jump_distance_cumul_histogram:
.. figure:: images/trace_diffusion_jump_distance_cumul_histogram.png
    :align: center
    :figwidth: 80%

    Jump distance cumulative probability histogram.

    The best fit for the single species model is shown in magenta.

.. index:: histograms

Histograms
^^^^^^^^^^

If the ``Show histograms`` option is selected the plugin presents a second dialog where the histograms can be configured. The number of bins in the histogram can be specified and outliers can optionally be removed. Outliers are any point more than 1.5 times the inter-quartile range above or below the upper and lower quartile boundaries. The following histograms can be chosen:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Parameter
     - Description

   * - Total signal
     - The total signal of each trace.

   * - Signal-per-frame
     - The signal-per-frame of the localisations in a trace.

   * - t-On
     - The on-time of a trace. This excludes the traces too short to be analysed.

   * - MSD/Molecule
     - The average mean-squared distance per molecule. Plots of the all-vs-all and adjacent MSD are shown.

       If the particles contain molecules moving with different diffusion rates or a fixed fraction of molecules then the histogram may be multi-modal.

   * - D*/Molecule
     - The apparent diffusion coefficient per molecule. Plots of the all-vs-all and adjacent D* are shown.


.. index:: summary table

Summary table
^^^^^^^^^^^^^

The plugin shows a summary table of the analysis results. This allows the plugin to be run with many different settings to view the effect on the calculated diffusion coefficient. The following columns are reported:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Field
     - Description

   * - Title
     - The title (specified by the ``Title`` parameter).

   * - Dataset
     - The input dataset.

   * - Exposure time
     - The dataset exposure time per frame.

   * - D-threshold
     - The distance threshold.

   * - Ex-threshold
     - The exclusion distance.

   * - Min-length
     - The minimum track length that was analysed.

   * - Ignore ends
     - True if the end jumps of tracks were ignored.

   * - Truncate
     - True if tracks were truncated to the min length.

   * - Internal
     - True if internal distance were used.

   * - Fit length
     - The number of points fitted in the linear regression.

   * - MSD corr
     - True if MSD correction was applied.

   * - S corr
     - True if precision correction was applied.

   * - MLE
     - True if maximum likelihood fitting was used.

   * - Traces
     - The number of traces analysed.

   * - s
     - The average precision of the localisations in the traces.

   * - D
     - The diffusion coefficient from MSD linear fitting.

   * - Fit s
     - The fitted precision when fitting an intercept in the MSD linear fit.

   * - Jump distance
     - The time distance used for jump analysis.

   * - N
     - The number of jumps for jump distance analysis.

   * - Beta
     - The beta parameter which is the ratio between the mean squared distance the localisation precision: :math:`\frac{\mathit{MSD}}{s^2}`.

       A beta above 6 indicates that jump distance analysis will produce reliable results [Weimann *et al*, 2013].

   * - Jump D
     - The diffusion coefficient(s) from jump analysis.

   * - Fractions
     - The fractions of each population from jump analysis.

   * - IC
     - The information criterion (IC) for the best model fit.

       Note that the IC is not comparable between the MLE or LSQ methods for fitting. It is also not comparable when the number of jumps is different. It can only be used to compare fitting the same jump distances with a different number of mobile species. This can can be controlled using the ``Minimum`` and ``Maximum N`` parameters.

   * - Total signal
     - The average total signal of each trace.

   * - Signal/frame
     - The average signal-per-frame of the localisations in a trace.

   * - t-On
     - The average on-time of a trace. This excludes the traces too short to be analysed.


The plugin will report the number of traces that were excluded using the length criteria and the fitting results to the ``ImageJ`` log. This includes details of the jump analysis with the fitting results for each model and the information criterion used to assess the best model, e.g.

.. code-block:: text

    783 Traces filtered to 117 using minimum length 5
    Linear fit (5 points) : Gradient = 2.096, D = 0.5239 um^2/s, SS = 0.047595 (2 evaluations)
    Jump Distance analysis : N = 151, Time = 6 frames (0.6 seconds). Mean Distance = 1371.0 nm, Precision = 38.55 nm, Beta = 35.57
    Estimated D = 0.4698 um^2/s
    Fit Jump distance (N=1) : D = 0.0498 um^2/s, SS = 0.433899, IC = -453.1 (12 evaluations)
    Fit Jump distance (N=2) : D = 1.655, 0.0346 um^2/s (0.1832, 0.8168), SS = 0.014680, IC = -960.3 (342 evaluations)
    Fit Jump distance (N=3) : D = 1.655, 0.0346, 0.0346 um^2/s (0.1832, 0.1204, 0.6964), SS = 0.014680, IC = -956.1 (407 evaluations)
    Coefficients are not different: 0.0346 / 0.0346 = 1.0
    Best fit achieved using 2 populations: D = 1.655, 0.0346 um^2/s, Fractions = 0.1832, 0.8168

.. index:: trace diffusion (multi)

Trace Diffusion (Multi)
-----------------------

This plugin allows the ``Trace Diffusion`` plugin to be run with multiple input datasets. Each dataset will be traced separately. The results are then combined for analysis. This allows analysis of multiple repeat experiments as if one single dataset.

When the plugin runs a dialog is presented that allows the datasets to be selected (:numref:`Figure %s <fig_trace_diffusion_multi_selection>`).

.. _fig_trace_diffusion_multi_selection:
.. figure:: images/trace_diffusion_multi_selection.png
    :align: center
    :figwidth: 80%

    Trace Diffusion (Multi) dataset selection dialog

*  Click a single result set to select or deselect.

*  Hold the ``Shift`` key to select or deselect a range of results starting from the last clicked result set.

*  Use the ``All`` or ``None`` buttons to select or deselect all the results.

*  Click the ``Cancel`` button to end the plugin.

*  Click the ``OK`` button to run the
   ``Trace Diffusion``
   plugin with the selected results.

When the ``Trace Diffusion`` plugin is executed it will not have the ``Input`` option as the results have already been selected. If multiple datasets are chosen the dataset name in the results table will be named using the first dataset plus the number of additional datasets, e.g. ``Dataset 1 + 6 others``.

Note that the plugin supports the ``ImageJ`` recorder to allow running within an ``ImageJ`` macro.
