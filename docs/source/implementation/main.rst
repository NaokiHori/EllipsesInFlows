
.. include:: /references.txt

##############
Implementation
##############

***************
Parallelisation
***************

Since this library is built on `a MPI-parallelised DNS solver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_, Eulerian domain is decomposed in :math:`y` direction.
Regarding particles, update procedures such as momentum exchange are conducted on the same Eulerian domain and thus the parallelisation is straightforward.
The velocity and position of the gravity centers are stored in a Lagrangian way and shared among all processes.
Since they are :math:`\mathcal{O} \left( N_p \right)` operations, the computational cost is negligibly small compared to the fluid solver :math:`\mathcal{O} \left( N_x \times N_y \right)`: proportional to the number of grid points.

Collisions (judgements and computations of repellent forces), on the other hand, require relatively heavy computational loads since we need :math:`\mathcal{O} \left( N_p^2 \right)` operations.
Although it is still much smaller than :math:`N_x \times N_y`, it can be fairly heavy to be completed by a main process when the volume fraction or the domain size is large.
Thus this library parallelises the procedure simply as the same way as the decomposing the Eulerian domain.
Namely, there are following collision pairs to be considered:

.. csv-table::
   :file: collision-table.csv

Here the left-most column :math:`n_0` and the top-most row :math:`n_1` show indices of the particles.
The other numbers show the indices of pairs, which are separated to each processor and considered independently.

In practice, we need to translate the index of collision pair :c-lang:`n` to the corresponding particle indices :math:`n_0` and :math:`n_1`, which are handled by :c-lang:`get_particle_indices`.

