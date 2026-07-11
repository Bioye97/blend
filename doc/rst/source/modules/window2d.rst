window2d
========

Synopsis
--------

.. code-block:: text

   blend window2d -R<xmin>/<xmax>/<ymin>/<ymax> -I<dx>[/<dy>]
                  [-F<xfunction>[/<yfunction>]]
                  [-T<rx1>/<rx2>/<ry1>/<ry2>]
                  [-B<blendfile>] [-C<f|l|o|u|a|g|p>] [-M<method>]

Description
-----------

``window2d`` returns weights on a 2-D grid or at query coordinates read from
standard input. Query coordinates use bilinear interpolation.

Blendfile rows contain:

.. code-block:: text

   <polygonfile> <xfunction>/<yfunction> <rx1>/<rx2>/<ry1>/<ry2>

Polygons must be valid, simple, and xy-monotone. If a polygon is not
xy-monotone, use ``-Me`` or ``-Mb`` to refine it. BLEND reports the original
and final vertex counts when refinement changes the polygon.

Options
-------

``-R``
   Set ``xmin/xmax/ymin/ymax``.

``-I``
   Set ``dx[/dy]``. If ``dy`` is omitted, ``dy = dx``.

``-F``
   Set x and y window functions.

``-T``
   Set x and y taper ratios.

``-B``
   Read polygon supports from a blendfile.

``-C``
   Set overlap handling.

``-M``
   Convert non-xy-monotone polygons using ``e`` or ``b``.
