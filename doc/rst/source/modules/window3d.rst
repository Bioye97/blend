window3d
========

Synopsis
--------

.. code-block:: text

   blend window3d -R<xmin>/<xmax>/<ymin>/<ymax>/<zmin>/<zmax>
                  -I<dx>[/<dy>[/<dz>]]
                  [-F<xfunction>[/<yfunction>[/<zfunction>]]]
                  [-T<rx1>/<rx2>/<ry1>/<ry2>/<rz1>/<rz2>]
                  [-B<blendfile>] [-C<f|l|o|u|a|g|p>] [-M<method>] [-N]

Description
-----------

``window3d`` extends ``window2d`` by applying the same xy polygon along a
vertical interval. Query coordinates use trilinear interpolation.

Blendfile rows contain:

.. code-block:: text

   <polygonfile> <zlo> <zhi> <xfunction>/<yfunction>/<zfunction> <rx1>/<rx2>/<ry1>/<ry2>/<rz1>/<rz2>

The xy polygon follows the same validity and xy-monotone rules as
``window2d``. Use ``-N`` with ``-M`` to write each modified polygon to
``<polygonfile>_monotone``.

Example
-------

.. code-block:: sh

   blend window3d -R0/10/0/10/0/5 -I0.5/0.5/0.25 -Bsupports.txt -Mb -N
