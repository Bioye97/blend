window2d
========

Description
-----------

``window2d`` returns weights on a 2-D grid or at query coordinates read from
standard input. Query coordinates use bilinear interpolation.

Blendfile rows contain:

.. code-block:: text

   <polygonfile> <xfunction>/<yfunction> <rx1>/<rx2>/<ry1>/<ry2>

Polygons must be valid, simple, and strictly xy-monotone, meaning that 
the vertices should be strictly increasing as opposed to nondecreasing 
in the usual definition of xy-monotonicity. use ``-ME`` or ``-MB`` to
refine any input polygon to be strictly xy-monotone. BLEND reports the 
original and final vertex counts when refinement changes the polygon. 
Use ``-N`` with ``-M`` to write each modified polygon by inserting
``_monotone`` before the file extension, for example 
``south_america_monotone.txt`` for the monotone South America polygon 
``south_america.txt``.

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

.. blend-usage:: window2d
