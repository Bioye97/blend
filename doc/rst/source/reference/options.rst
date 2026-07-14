Common Options
==============

``-R``
   Region. In 1-D this is ``xmin/xmax``. In 2-D this is
   ``xmin/xmax/ymin/ymax``. In 3-D this is
   ``xmin/xmax/ymin/ymax/zmin/zmax``.

``-I``
   Grid increment. Optional dimensions inherit from the first increment. 
   In 1-D, this is ``dx``. In 2-D, this is ``dx/dy``. In 3-D, this is ``dx/dy/dz``.

``-F``
   Window taper functions for each dimension. Run ``blend --show-windows`` to list available
   names. In 1-D, this is ``xfunction``. In 2-D, this is ``xfunction/yfunction``. In 3-D, this is
   ``xfunction/yfunction/zfunction``.

``-T``
   Taper ratios for beginning and ending sides of each dimension. In 1-D, this is ``rx1/rx2``. 
   In 2-D, this is ``rx1/rx2/ry1/ry2``. In 3-D, this is ``rx1/rx2/ry1/ry2/rz1/rz2``.

``-B``
   Blendfile input for multiple polygons/supports.

``-C``
   Clobber mode for polygon overlap handling: first ``f``, lowest ``l``, last ``o``, highest ``u``,
   arithmetic average ``a``, geometric average ``g``, or product ``p``.

``-M``
   Monotone conversion method for non-xy-monotone polygons. Options include ``e``, 
   ``b``, ``E``, and ``B``. Options ``e`` and ``b`` are only available for the ``monotone`` 
   module and use the regular mathematical definition of xy-monotonicity, for 
   ``nondecreasing vertices``. Options ``E`` and ``B`` are the only available options for
   the ``window2d`` and ``window3d`` modules and use a stricter definition of xy-monotonicity, 
   for ``strictly increasing vertices``.

``-N``
   Write polygons modified by ``-M`` by inserting ``_monotone`` before the
   file extension, for example ``south_america_monotone.txt`` for ``south_america.txt``.

``-V``
   Verbosity. Levels are quiet ``q``, errors only ``e``, warnings and errors ``w``,
   timings, warnings, and errors ``t``, informational messages, timings, warnings, and errors ``i``, 
   compatibility messages and all lower verbosity messages ``c``, and debug messages and all 
   lower verbosity messages ``d``.
