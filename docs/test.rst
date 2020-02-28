:orphan:


Heading 1
=========


Heading 2
---------


Heading 3
~~~~~~~~~


Heading 4
^^^^^^^^^

Hello world!


Formatting
==========

**Bold**

*Italic*

:sup:`superscript` like 1\ :sup:`st`.

``inline block``

| line
| break in the same sentence.

b\ :sup:`2`

s\ :sub:`a`


Equations
=========

inline: :math:`y = ax + c`

.. math::

    y &= ax + c \\
      &= ax + c_0

.. math::

    \mathit{value\:1} = \Delta x

:math:`\mathit{SNR}=\frac{N}{\sqrt{b^2}}`


Tables
======


Style 1
-------

.. list-table::

    * - :math:`B`
      - The background level.

    * - :math:`\mathit{Signal}`
      - The total volume of the Gaussian.

    * - :math:`x_0`
      - The X centre of the Gaussian.

    * - :math:`y_0`
      - The Y centre of the Gaussian.

    * - :math:`\sigma_x`
      - The X standard deviation.

    * - :math:`\sigma_y`
      - The Y standard deviation.

    * - :math:`\theta`
      - The angle of rotation of the ellipse.


Style 2
-------

The indentation of the ``:widths:`` etc must match that of the table.

.. list-table:: Table title
    :widths: 20 80
    :header-rows: 1

    * - Parameter
      - Description

    * -  Use current calibration
      -  If selected use the current SMLM configuration.
         Otherwise run the configuration wizard.

         This option is only shown if a SMLM configuration file
         can be found. If no file is found then the configuration
         wizard is run by default.

    * -  Show table
      -  Show a table containing the localisations.

    * -  Show image
      -  Show a super-resolution image of the localisations.
         The image will be 1024 pixels on the long edge.


Style 3
-------

=====  ===========
Param  Description
=====  ===========
SD     The standard deviation of the Gaussian approximation to the Airy pattern
p      The proportionality factor.
       Using a value of 1 gives the theoretical lower bounds on the peak width.
=====  ===========


References
==========

See :ref:`calibration_plugins:PSF Estimator`.

See :numref:`{number}: {name} <fitting_plugins:Fitting Parameters>`

See :ref:`comparison_metrics:Comparison Metrics` for more details.


Code
====

.. code-block:: xml

    <xml>
        <start/>
    </xml>

.. code-block:: text

    any old thing
