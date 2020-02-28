.. index:: ! Comparison Metrics

Comparison Metrics
==================

Several plugins within the SMLM package compute matches between points. One set of points can be labelled as the actual result, the second can be labelled as the predicted results. When comparing actual and predicted points the following combinations are possible:

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 0

   * - Actual Point
     - Predicted Point
     - Classification
     - Label

   * - Present
     - Present
     - True Positive
     - *tp*

   * - Present
     - Absent
     - False Negative
     - *fn*

   * - Absent
     - Present
     - False Positive
     - *fp*

   * - Absent
     - Absent
     - True Negative
     - *tn*

The classification counts can be used to compute binary scoring statistics as described below.

**Notes**

#.  The categorisation of the actual and predicted points can be arbitrary. If the categorisation is reversed
    then some metrics will change and some are invariant.
#.  When comparing point coordinates that can only ever exist the definition of the true negative is invalid. However there are cases where a set of results have two classifications (either absent or present), and a set of predictions aim to predict those values. In this case the TN count can be obtained and used to compute scores.


.. index:: ! Recall

Recall
------

Recall measures the number of actual points that are correctly predicted. It is also known as the True Positive Rate
(TPR)
or sensitivity.

.. math::

    \mathit{recall}=\frac{\mathit{tp}}{\mathit{tp}+\mathit{fn}}

A score of 1 indicates that all the points were predicted, lower scores indicate that some points were missed.

Recall can be interpreted probabilistically as the chance that a randomly selected actual point will be predicted.


.. index:: ! Precision

Precision
---------

Precision measures the confidence of the predicted points. It is also known as the Positive Predicted Value
(PPV)
.

.. math::

    \mathit{precision}=\frac{\mathit{tp}}{\mathit{tp}+\mathit{fp}}

A score of 1 indicates that all the predicted points were correct, lower scores indicate that
some points are not correct.

Precision can be interpreted probabilistically as the chance that a randomly selected prediction is correct.


.. index:: ! Jaccard

Jaccard
-------

The Jaccard measures the similarity between two sets and is defined as the size of the intersection divided by the size of the union:

.. math::

    J(A,B)=\frac{\left|{A\cap B}\right|}{\left|{A\cup B}\right|}=\frac{\mathit{tp}}{\mathit{tp}+\mathit{fp}+\mathit{fn}}

A score of 1 indicates that the overlap is perfect. Zero indicates no overlap.


.. index:: ! F-score

F-score
-------

The precision and recall can be combined in a weighted score, the :math:`F_\beta`-measure or F-score.

.. math::

    F_{\beta }=(1+\beta ^{2})\cdot {\frac{\mathit{precision}\cdot
    \mathit{recall}}{\beta ^{2}\cdot {\mathit{precision}+\mathit{recall}}}}

The :math:`\beta` is a weighting.  The score was derived so that it measures the effectiveness of retrieval with respect to a user who associates :math:`\beta` times as much importance to recall as precision. For example :math:`F_{0.5}` puts more emphasis on precision and :math:`F_{2}` puts more emphasis on recall.

A weight of 1 produces the balanced F-score where precision and recall are weighted equally.

.. math::

    F_1={\frac{\mathit{precision}\cdot
    \mathit{recall}}{\mathit{precision}+\mathit{recall}}}


.. index:: ! FNR

FNR
---

False-negative rate:

.. math::

    \mathit{FNR}=\frac{\mathit{fn}}{\mathit{fn}+\mathit{tp}}


.. index:: ! FDR

FDR
---

False discovery rate:

.. math::

    \mathit{FDR}=1-\mathit{Precision}=\frac{\mathit{fp}}{\mathit{tp}+\mathit{fp}}


.. index:: ! TNR

TNR
---

True-negative rate:

.. math::

    \mathit{TNR}=\frac{\mathit{tn}}{\mathit{fp}+\mathit{tn}}


.. index:: ! NPV

NPV
---

Negative predictive value:

.. math::

    \mathit{NPV}=\frac{\mathit{tn}}{\mathit{tn}+\mathit{fn}}


.. index:: ! FPR

FPR
---

False-positive rate:

.. math::

    \mathit{FPR}=\frac{\mathit{fp}}{\mathit{fp}+\mathit{tn}}


.. index:: ! ACC

ACC
---

Accuracy:

.. math::

    \mathit{Accuracy}=\frac{\mathit{tp}+\mathit{tn}}{\mathit{tp}+\mathit{fp}+\mathit{tn}+\mathit{fn}}


.. index:: ! MCC

MCC
---

Matthews Correlation Coefficient:

.. math::

    \mathit{MCC}=\frac{\mathit{tp}\ast \mathit{tn}-\mathit{fp}\ast
    \mathit{fn}}{\sqrt{\mathit{tp}+\mathit{fp}\ast
    {\mathit{tp}+\mathit{fn}}\ast {\mathit{tn}+\mathit{fp}}\ast
    {\mathit{tn}\ast \mathit{fn}}}}

The Matthews Correlation Coefficient is used in machine learning as a measure of the quality of binary (two-class) classifications, introduced by biochemist Brian W. Matthews in 1975. It takes into account true and false positives and negatives and is generally regarded as a balanced measure which can be used even if the classes are of very different sizes. The MCC is in essence a correlation coefficient between the observed and predicted binary classifications; it returns a value between −1 and +1. A coefficient of +1 represents a perfect prediction, 0 no better than random prediction and −1 indicates total disagreement between prediction and observation. The statistic is also known as the phi coefficient.


.. index:: ! Informedness

Informedness
------------

.. math::

    \mathit{Informedness}=\mathit{TPR}+\mathit{TNR}-1


.. index:: ! Markedness

Markedness
----------

.. math::

    \mathit{Markedness}=\mathit{PPV}+\mathit{NPV}-1
