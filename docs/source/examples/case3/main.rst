
.. _example_rotating_ellipse:

.. include:: /references.txt

######################################
Rotation of an ellipse in a shear flow
######################################

Comparison with the reference data by |AIDUN1998| is described, who simulated an rotating ellipse in a shear flow.

*************
Configuration
*************

We consider a rectangular domain (:math:`l_x = 1, l_y = 2`), in which an elliptic object (whose major and minor axes are :math:`0.10` and :math:`0.05`, respectively) is positioned at :math:`x = 0.5, y = 1.0`.
An analytical Couette flow profile is given initially, while the particle is at rest with its major axis being parallel to the wall-normal direction (:math:`\theta_z \left( t = 0 \right) = 0`).
The liquid viscosity is adjusted so that the bulk Reynolds number :math:`Re = l_x \Delta U_{walls} / \nu` leads :math:`25`, which corresponds to :math:`Re_p = 1` (see |AIDUN1998| for the definition of the particle-based Reynolds number).
Three spatial resolutions, from :math:`48` to :math:`96` grids int the wall-normal direction, are considered.

In practice, the configuration is specified as follows (:math:`48` grids per domain):

.. literalinclude:: config/exec48.sh
   :language: sh

*******
Results
*******

.. note::
   Conducted by GitHub Actions.

.. literalinclude:: data/48/ci.txt
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

Angular velocity of the particle :math:`\omega_z` as a function of time is shown below.
Note that, although the reference data is obtained with :math:`160` grid points (:math:`320` lattice nodes) in the wall-normal direction, we limit the resolutions to lower values to reduce the computational cost.

.. image:: data/result.png
   :width: 800

