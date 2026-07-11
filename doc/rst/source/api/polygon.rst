Polygon API
===========

Header
------

.. code-block:: c

   #include <blend_polygon.h>

Data Types
----------

``vertex``
   A 2-D point with ``x`` and ``y`` coordinates.

``polygon``
   A dynamic vertex array and vertex count.

Core Functions
--------------

``blend_polygon_alloc`` allocates a polygon.
``blend_polygon_free`` releases polygon storage.
``blend_polygon_read`` and ``blend_polygon_write`` handle text files.
``blend_polygon_validate`` checks that a polygon is valid and simple.
``blend_polygon_is_xy_monotone`` checks monotonicity.
``blend_polygon_xy_monotone_envelope`` and
``blend_polygon_xy_monotone_best_piecewise_envelope`` construct monotone
approximations.
