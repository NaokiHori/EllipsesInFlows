#################
Ellipses in Flows
#################

For circular (2D) or spherical (3D) objects, please refer to `the other library <https://github.com/NaokiHori/SimpleIBMSolver>`_.

|License|_ |LastCommit|_ |CI|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/EllipsesInFlows
.. _License: https://opensource.org/licenses/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/EllipsesInFlows/main
.. _LastCommit: https://github.com/NaokiHori/EllipsesInFlows/commits/main

.. |CI| image:: https://github.com/NaokiHori/EllipsesInFlows/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/EllipsesInFlows/actions/workflows/ci.yml

.. image:: https://github.com/NaokiHori/EllipsesInFlows/blob/main/docs/source/snapshot.png
   :width: 800
   :target: https://youtu.be/iuO5CxvAlio

********
Overview
********

This library numerically solves the motion of elliptic bodies governed by the Newton-Euler equations suspended in viscous liquid by means of the finite-difference and the immersed-boundary methods.

Please refer to `the documentation <https://naokihori.github.io/EllipsesInFlows/index.html>`_ for details.

********
Features
********

* MPI-parallelised
* Eulerian-based (no Lagrangian points) IBM
* Collision model between ellipses

**********
Dependency
**********

`Docker <https://www.docker.com>`_

***********
Quick start
***********

#. Create working directory

   .. code-block:: console

      mkdir /path/to/your/working/directory
      cd    /path/to/your/working/directory

#. Fetch source

   .. code-block:: console

      git clone https://github.com/NaokiHori/EllipsesInFlows
      cd EllipsesInFlows

#. Build

   .. code-block:: console

      docker build -t simplenavierstokessolver:latest .
      docker run -it --rm --cpuset-cpus="0-1" -u runner -v ${PWD}:/home/runner simplenavierstokessolver:latest
      make output
      make all

#. Run

   .. code-block:: console

      mpirun -n 2 ./a.out

********
Examples
********

Several examples can be found in the documentation.

#. `Migration of a circular object in a shear flow <https://naokihori.github.io/EllipsesInFlows/examples/case1/main.html>`_

#. `Segré-Silberberg effect <https://naokihori.github.io/EllipsesInFlows/examples/case2/main.html>`_

#. `Rotation of an ellipse in a shear flow <https://naokihori.github.io/EllipsesInFlows/examples/case3/main.html>`_

#. `Suspension in a plane Poiseuille flow <https://naokihori.github.io/EllipsesInFlows/examples/case4/main.html>`_

