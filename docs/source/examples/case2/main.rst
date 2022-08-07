
.. _example_segre_silberberg:

.. include:: /references.txt

#######################
Segré-Silberberg effect
#######################

Comparison with the reference data by |PAN2002| is described, who simulated the Segré-Silberberg effect in a two-dimensional wall-bounded domain by means of the fictitious domain method.

*************
Configuration
*************

We consider a square domain (:math:`l_x = l_y = 1`), in which a circular object (whose radius is :math:`0.125`) is positioned at :math:`x = 0.4, y = 0`.
Everything is stationary at :math:`t = 0`, and the liquid (viscosity is adjusted so that the bulk Reynolds number leads :math:`2334`) and particle start to move because of a constant forcing in :math:`y` (stream-wise) direction, whose magnitude is :math:`-dp / dy = 2.337 \times 10^{-4}`.
Spatial resolutions are varied from :math:`32` to :math:`96` grids per domain size (:math:`8` to :math:`24` grids per diameter) to check the spatial convergence.

In practice, the configuration is specified as follows (:math:`32` grids per domain):

.. literalinclude:: config/exec32.sh
   :language: sh

*******
Results
*******

.. note::
   Conducted by GitHub Actions.

.. literalinclude:: data/32/ci.txt
   :language: text

.. literalinclude:: data/64/ci.txt
   :language: text

.. literalinclude:: data/96/ci.txt
   :language: text

The final flow field is shown below.
Note that the picture is transposed, i.e., in the picture, :math:`x` and :math:`y` directions in the simulation are shown as the vertical and horizontal directions, respectively.

.. image:: data/96/snapshot.png
   :width: 400

The colour denotes the size of the stream-wise velocity, while the center gray object is the immersed particle.

Lateral (wall-normal, :math:`x` direction) migration of a particle as a function of the traveling distance in the stream-wise (:math:`y`) direction is shown below.

.. image:: data/result.png
   :width: 800

