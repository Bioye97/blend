Common Options
==============

``-R``
   Region. In 1-D this is ``xmin/xmax``. In 2-D this is
   ``xmin/xmax/ymin/ymax``. In 3-D this is
   ``xmin/xmax/ymin/ymax/zmin/zmax``.

``-I``
   Grid increment. Optional dimensions inherit from the first increment.

``-F``
   Window taper function. Run ``blend --show-windows`` to list available
   names.

``-T``
   Taper ratios for beginning and ending sides of each dimension.

``-B``
   Blendfile input for multiple supports.

``-C``
   Clobber mode for overlap handling: first, lowest, last, highest,
   arithmetic average, geometric average, or product.

``-M``
   Monotone conversion method for non-xy-monotone polygons.

``-N``
   Write polygons modified by ``-M`` to ``<polygonfile>_monotone``.

``-V``
   Verbosity. Levels are ``q``, ``e``, ``w``, ``t``, ``i``, ``c``, and ``d``.
