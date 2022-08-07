############
Introduction
############

This library is built on `Simple Navier-Stokes Solver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_, but I made the following changes.

   * Rigid body motions

      In addition to the working liquid, finite-size elliptic objects are considered and their translational and rotational motions are integrated in time by means of the immersed boundary method and soft-sphere collision model.

   * No temperature coupling

      Temperature field is out of my focus in this project.

   * Uniform grid spacings in all directions

      In addition to the walls, boundary layers are formed on particle surfaces.
      Although their thicknesses are totally different (especially for high Reynolds number flows), using stretched grid in the wall-normal direction is less advantageous compared to the original single-phase Rayleigh-Bénard flows.
      Also stretched grid can break the angular momentum balances of the particle motions.
      Thus we adopt uniform grid spacings not only in :math:`y` but in :math:`x` directions.

   * Efficiency

      The most computationally-intensive procedure in this project is to the matrix transposes, and we need to do so for four times when solving the pressure Poisson equation.
      Since grid points are uniformly distributed in :math:`x` direction, we can consider signals in the wave space, which halves the number of parallel matrix transform and can boost the performance.

They are explained in this document, while other parts are omitted since they are basically identical to the original library.

.. note::

   This library is developed just for fun in my summer holiday.
   Although it is based on `a publicaton <https://www.sciencedirect.com/science/article/pii/S0045793021003716>`_, several major changes are made and thus this library is unofficial.

