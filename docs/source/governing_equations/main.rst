
.. include:: /references.txt

###################
Governing equations
###################

*******************************************
Translational and angular momentum balances
*******************************************

In addition to the Navier-Stokes equations

.. math::
   \der{u_i}{x_i} &= 0, \\
   \der{u_i}{t}
   +
   u_j \der{u_i}{x_j}
   &=
   -
   \der{p}{x_i}
   +
   \frac{1}{Re} \der{}{x_j} \der{u_i}{x_j}
   +
   a_i,

which govern the behaviour of the liquid, we consider the Newton-Euler equations to describe the motions of rigid particles (density: :math:`\rho_p`, volume: :math:`V_p`, surface: :math:`S_p` or :math:`\partial V_p`):

.. math::

   \rho_p V_p \oder{U_i}{t}
   =
   -
   \int_{\partial V_p} a_i d S_p
   +
   \oder{}{t} \int_{V_p} u_i d V_p
   +
   \left( 1 - \rho_p \right) \frac{1}{Fr^2}
   +
   F_i,

.. math::

   \oder{}{t} \left( I_m \Omega_i \right)
   =
   \int_{\partial V_p} \epsilon_{ijk} r_j a_k d S_p
   +
   \oder{}{t} \int_{V_p} \epsilon_{ijk} r_j u_k d V_p
   +
   T_i.

Note that they are already normalised by the reference scales, and :math:`Fr` is the Froude number defined as

.. math::

   Fr
   \equiv
   \frac{U}{\sqrt{gL}}.

The first equation describes the translational velocity :math:`U_i` of the gravity center :math:`X_i`, while the second equation tells how the angular velocity :math:`\Omega_i` is evolved in time.
The first terms in the right-hand-sides (including :math:`a_i`) denote the force and torque caused by the fluid-structure interactions, while the second terms are responsible for the motion of fictitious fluid inside objects (see |BREUGEM2012|).
The last terms :math:`F_i` and :math:`T_i` are the external force and torque, which take into account the gravitational acceleration and collisions between the other objects in this project.

Since I am interested in ellipses, it is worthwhile to introduce the following basic parameters and their relations:

.. math::

   \begin{cases}
      \text{major axis} & a, \\
      \text{minor axis} & b, \\
      \text{volume} \, V_p & \pi a b, \\
      \text{mass} \, m_p & \rho_p \pi a b, \\
      \text{moment of inertia} \, I_m & \frac{1}{4} \pi \left( a^2 + b^2 \right)
   \end{cases}

.. note::

   Since I limit my focus to two-dimensional objects for now, the volume and surface integrals reduce to the surface and line integrals, respectively.
   Also the moment of inertia inside the temporal derivative is a constant for each object and thus can be taken out, i.e., the left-hand-side of the equation of angular velocity leads

   .. math::

      \oder{}{t} \left( I_m \Omega_i \right)
      =
      I_m \oder{\Omega_i}{t}.

*********
Collision
*********

The lubrication force coming from the hydrodynamic interaction between particles (or a particle and a wall) tends to be underestimated especially when the resolution between the objects is insufficient, and as a result particles can penetrate (too much) to each other.
In order to avoid this, this library employs a collision model between objects, which is based on a simple spring model:

.. math::

   m \oder{U_i}{t} + k X_i = 0,

in which a spring whose spring constant is :math:`k` is assumed between objects and a corresponding force whose direction is parallel to the normal and magnitude is proportional to the penetration depth pushes away the objects to each other.
This force is only activated when a penetration is observed and try to resolve the penetration gently (in :math:`T`), which is the so-called soft-sphere collision model (|TSUJI1993|).

A general solution of the above equation with boundary conditions

.. math::

   X_i \left( t = 0 \right) = X_i \left( t = T \right) = 0

leads

.. math::

   X_i = A \sin \sqrt{\frac{k}{m}} t,

where

.. math::

   k = m \frac{\pi^2}{T^2}

and thus the repellent force as a function of the penetration depth in the normal direction :math:`x_i` is given by

.. math::

   f_i = k x_i.

Since the particles are rigid, :math:`T = 0` in theory, which is impractical within the current numerical scope (so-called event-driven collision model deals with :math:`T = 0`, which is not easy to integrate in time simultaneously as the fluid behaviour).

As a remedy, in this project, I consider a timescale

.. math::

   T = \sqrt{\frac{\rho r^3}{\sigma}} \propto We^{\frac{1}{2}} r^{\frac{3}{2}},

which is analogous to the capillary timescale for deformable objects whose :math:`\sigma` is the surface tension coefficient (and :math:`We` is a pseudo Weber number).
Note that :math:`r` is a particle-based length scale (e.g., harmonic average of the colliding particle radii).

In theory, :math:`\sigma \rightarrow \infty` (or equivalently :math:`T \rightarrow 0`) since the objects are rigid.
Giving larger :math:`T` stabilises the integration process because the repellent force becomes smaller, which is especially useful when the suspension volume fraction is relatively large (e.g., :math:`\gt 40 \%`).
However, the outcome might be unphysical because of the large penetration depths, and thus a proper timescale should be determined by the user eventually.

The above spring model does not include a dumper, indicating that the collision is perfectly elastic (restitution coefficient is :math:`1`) and no dissipation is involved.
Obviously this is not true in most cases and also elastic, dumping, and sliding motions in the tangential direction should be included (see |COSTA2015| for extensive analyses).

In this project, however, I neglect all these aspects just for simplicity, and include only springs in the normal directions to obtain *plausible* results.

An extension of the above concept to elliptic objects is not straightforward, which will be discussed later.

