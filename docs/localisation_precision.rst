.. index:: ! Localisation Precision

Localisation Precision
======================

The theoretical limit (precision) for fitting the signal (number of photons) and the XY coordinates (localisation) can be computed using the formulas of Thompson *et al* (2002) for the signal and Mortensen *et al* (2010) for the localisation.

Note that these formulas are derived from modelling the point spread function (PSF) as a 2D Gaussian for both the simulation and the fitting. Given that the true data will have a PSF defined by the microscope parameters these formulas only approximate the precision that can be obtained on image data.

The photon count *N* is computed using the volume of the fitted 2D Gaussian function (the signal) divided by the camera gain:

.. math::

    N=\frac{\mathit{Signal}}{\mathit{Gain}}

The background noise :math:`b^2` is estimated by the ``Peak Fit`` plugin during the fitting process either using a global noise estimate per frame or the local background level. This is also adjusted by the camera gain to provide the noise in photons.

Note that the precision of the localisation is the square root of the variance.

**Warning: Without a correctly calibrated camera gain and bias the precision estimates in the SMLM code will be inaccurate.**


.. index:: ! Signal Precision

Signal Precision
----------------

From Thompson *et al* (2002).

.. math::

    \mathit{Var}_{N}=F\times (N+\frac{4\pi s^{2}b^{2}}{a^{2}})

where

.. list-table::
    :widths: 20 80

    * - :math:`{Var}_{N}`
      - The variance of the signal when fitting a Gaussian 2D function to a Gaussian 2D PSF.

    * - :math:`F`
      - The noise scaling factor. 2 for an EM-CCD camera; 1 otherwise.

    * - :math:`S`
      - The standard deviation of the Gaussian function.

    * - :math:`N`
      - The number of photons in the localisation.

    * - :math:`b^2`
      - The expected number of photons per pixel from a background with spatially constant expectation value across the image.

    * - :math:`a`
      - The pixel size in nm.


.. index:: ! Localisation Precision for Least Squares Fitting

Localisation Precision for Least Squares Fitting
------------------------------------------------

From Mortensen *et al* (2010).

.. math::

    \mathit{Var}_{x}=F\times {\frac{s_{a}^{2}}{N}}\times
    (\frac{16}{9}+\frac{8\pi s_{a}^{2}b^{2}}{\mathit{Na}^{2}})

where

.. list-table::
    :widths: 20 80

    * - :math:`{Var}_{x}`
      - The variance of the localisation position in the X dimension when fitting a Gaussian 2D function to a Gaussian 2D PSF.

    * - :math:`F`
      - The noise scaling factor. 2 for an EM-CCD camera; 1 otherwise.

    * - :math:`s_a`
      - The standard deviation of the Gaussian function (:math:`s`) adjusted for square pixels:

        :math:`s_{a}=\sqrt{s^{2}+a^{2}/12}`

    * - :math:`N`
      - The number of photons in the localisation.

    * - :math:`b^2`
      - The expected number of photons per pixel from a background with spatially constant expectation value across the image.

    * - :math:`a`
      - The pixel size in nm.


.. index:: ! Localisation Precision for Maximum Likelihood Fitting

Localisation Precision for Maximum Likelihood Fitting
-----------------------------------------------------

From Mortensen *et al* (2010).

.. math::

    \mathit{Var}_{x}=F\times {\frac{s_{a}^{2}}{N}}\times {\frac{1}{-\int_{0}^{1}{\frac{t\ln (t)}{t+\rho }}\mathit{dt}}}

where

.. list-table::
    :widths: 20 80

    * - :math:`{Var}_{x}`
      - The variance of the localisation position in the X dimension when fitting a Gaussian 2D function to a Gaussian 2D PSF.

    * - :math:`F`
      - The noise scaling factor. 2 for an EM-CCD camera; 1 otherwise.

    * - :math:`s_a`
      - The standard deviation of the Gaussian function (:math:`s`) adjusted for square pixels:

        :math:`s_{a}=\sqrt{s^{2}+a^{2}/12}`

    * - :math:`N`
      - The number of photons in the localisation.

    * - :math:`\rho`
      - The integration factor:

        :math:`\frac{2\pi s_{a}^{2}b^{2}}{\mathit{Na}^{2}}`

    * - :math:`b^2`
      - The expected number of photons per pixel from a background with spatially constant expectation value across the image.

    * - :math:`a`
      - The pixel size in nm.

Note that since the formula for maximum likelihood fitting involves an integral with no analytic solution the formula is evaluated using numerical integration. This is slow to compute.

The formula can be used to demonstrate that for any given set of parameters the precision of maximum likelihood fitting is lower (i.e. better) than least squares fitting.
