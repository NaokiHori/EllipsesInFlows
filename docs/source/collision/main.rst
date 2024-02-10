#####################
Collision of Ellipses
#####################

In this part I consider the collisions between two ellipses; namely, in this picture:

.. image:: data/problem-setup.png
   :width: 400

I aim at investigating

* whether the two ellipses collide,
* if they do, how much the penetration depth :math:`\delta` is.

I assume the ellipses are completely smooth and rigid (i.e. no deformation), and :math:`\delta` is sufficiently small compared to the size of the ellipses.

.. note::

   This is a rough solution to obtain something plausible just to avoid the over-penetration between objects.
   Neither the collision detection nor the quantification of :math:`\delta` are rigorous.

.. seealso::

   `Documentation for Japanese readers <https://qiita.com/NaokiHori/items/daf3fd191d51a7e682f8>`_.

*************
Problem setup
*************

An ellipse whose center is at the origin and the major / minor axes are :math:`a_0` and :math:`b_0` is given by

.. math::

   \frac{x_e^2}{a_0^2}
   +
   \frac{y_e^2}{b_0^2}
   =
   1.

By introducing a parameter :math:`t_0`, this equation reads

.. math::

   \begin{pmatrix}
      x_e \\
      y_e
   \end{pmatrix}
   =
   \begin{pmatrix}
      a_0 \cos t_0 \\
      b_0 \sin t_0
   \end{pmatrix}.

To take into account the rotation of :math:`\theta_0` and off-centred objects whose centre is at :math:`\left( x_0, y_0 \right)`, I modify the equation as

.. math::

   \begin{pmatrix}
      x_e \\
      y_e
   \end{pmatrix}
   =
   \begin{pmatrix}
      x_0 \\
      y_0
   \end{pmatrix}
   +
   \begin{pmatrix}
      \cos \theta_0 & -\sin \theta_0 \\
      \sin \theta_0 &  \cos \theta_0
   \end{pmatrix}
   \begin{pmatrix}
      a_0 \cos t_0 \\
      b_0 \sin t_0
   \end{pmatrix}.

Hereafter subscripts :math:`0,1` are used to distinguish the two ellipses.

The collision / intersection of the two ellipses results in a quartic equation with respect to the intersecting point(s) :math:`x` (or :math:`y`).
The solution can be categorised as follows:

#. Four real solutions

   Four intersecting points

#. Two real and two imaginary solutions

   Two intersecting points

#. Four imaginary solutions

   No Collision

#. And more

We are interested in distinguishing the second and the third cases.
Note that the other cases are out of focus, since we assume that the collision depth is small enough.

Although we know a quartic formula, it is non-trivial to treat it numerically.
Also, let's say we obtain a theoretical solution; it is, however, still unclear how we define the penetration depth :math:`\delta`.

For these reasons, instead of solving quartic equations, I consider a simpler and approximated solution in this project.

**************************************
Collision check using circular fitting
**************************************

Although ellipses can rotate and be off-centered, by transforming the coordinate system:

.. math::
   \begin{pmatrix}
     x_p \\
     y_p
   \end{pmatrix}
   \leftarrow
   \begin{pmatrix}
     \cos ( -\theta) & - \sin ( -\theta) \\
     \sin ( -\theta) &   \cos ( -\theta)
   \end{pmatrix}
   \begin{pmatrix}
     x_p-x_0 \\
     y_p-y_0
   \end{pmatrix},

we are able to say one object is at the origin and does not rotate (the major axis coincides with the Cartesian :math:`x` axis).

======================
Circular approximation
======================

.. note::

   This part is largely inspired by `Simple Method for Distance to Ellipse <https://blog.chatfield.io/simple-method-for-distance-to-ellipse/>`_.

Recall that I consider situations where :math:`\delta` is much smaller than the sizes of the colliding objects.
This allows me to think that collisions between two circles fitting the ellipses in the vicinity of the colliding points well approximates the collisions between ellipses.

Collisions between two circles are much simpler; the penetration depth of two circles (not a rigorous definition but one sound candidate) leads to

.. math::

   \delta \equiv r_0 + r_1 - d,

where :math:`r_0` and :math:`r_1` are radii of the two circles, and :math:`d` is center-to-center distance.

Now the main focus is how to approximate an ellipse using a circle in a systematic way.
Actually there is an analytical circle which approximates an ellipse locally:

.. math::

   \left(
      x_c,
      y_c
   \right)
   =
   \left(
      a \left( 1 - \frac{b^2}{a^2} \right) \cos^3 t,
      b \left( 1 - \frac{a^2}{b^2} \right) \sin^3 t
   \right),

which is known as `evolute <https://en.wikipedia.org/wiki/Evolute#Evolute_of_an_ellipse>`_.
The local curvature is given by

.. math::

   \kappa
   =
   \frac{
     ab
   }{
      \sqrt{\left( a^2 \sin^2 t + b^2 \cos^2 t \right)^3}
   },

whose reciprocal is the radius of the circle fitting the original ellipse.

The only unknown is :math:`t` which gives a circle nicely approximating the ellipse locally.

Finding a nice :math:`t` is equal to find a normal vector from the point :math:`( x_p, y_p )` to the ellipse, and also is equal to find the minimum distance to the ellipse.

Although these lemmas are not proved here, the following picture gives a good indication:

.. image:: data/fit-circle.png
   :width: 400

Here,

* the black arrow connecting :math:`( x_p, y_p )` and :math:`( x_c, y_c )` gives a normal vector to the ellipse,

* the fitted circle gives a good approximation of the ellipse locally.

=========================
Collision of two ellipses
=========================

I use the above method to quantify the penetration depth :math:`\delta`, which is very simple: for the ellipse :math:`i`, the center of the fitted circle of the ellipse :math:`j` :math:`( x_{c_j}, y_{c_j} )` is used as the target point :math:`( x_{p_i}, y_{p_i} )` to fit a circle.
This process is iterated until the locations of :math:`( x_{c_i}, y_{c_i} )` converge.

When the two ellipses are colliding, the fitting circles lead

.. image:: data/fit-circles-0.png
   :width: 400

When the two ellipses are not colliding, the final state leads

.. image:: data/fit-circles-1.png
   :width: 400

Since we now know nice circles, the penetration depth is simply

.. math::

   \delta \equiv r_0 + r_1 - d,

where, again, :math:`r_i` are radii of the fitted circles, while :math:`d` is the distance between the two centers of the fitting circles.
Obviously two circles collide when :math:`\delta > 0` and do not otherwise.

**************************
Stand-alone implementation
**************************

.. literalinclude:: data/example.py
   :language: python

