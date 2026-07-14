window3d
========

Description
-----------

``window3d`` extends ``window2d`` by applying the same xy polygon along a
vertical interval. Query coordinates use trilinear interpolation.

Blendfile rows contain:

.. code-block:: text

   <polygonfile> <zlo> <zhi> <xfunction>/<yfunction>/<zfunction> <rx1>/<rx2>/<ry1>/<ry2>/<rz1>/<rz2>

The xy polygon follows the same validity and strict boundary-assembly rules as
``window2d`` after it is snapped to the local support grid. Use ``-N`` with
``-M`` to write each modified polygon to names such as
``south_america_monotone.txt``.

Usage
-----

.. raw:: html

   <span id="r"></span>
   <span id="i"></span>
   <span id="f"></span>
   <span id="t"></span>
   <span id="b"></span>
   <span id="c"></span>
   <span id="m"></span>
   <span id="n"></span>
   <span id="v"></span>

.. blend-usage:: window3d
