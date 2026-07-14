window1d
========

Description
-----------

``window1d`` returns weights on a 1-D grid or at query coordinates read from
standard input. Query coordinates do not need to be exactly on grid points;
BLEND uses linear interpolation from neighboring grid weights.

If ``-B`` is given, window supports are read from a blendfile and ``-F`` and
``-T`` are ignored with a warning. Points outside all blendfile supports are
assigned zero weight.

Usage
-----

.. raw:: html

   <span id="r"></span>
   <span id="i"></span>
   <span id="f"></span>
   <span id="t"></span>
   <span id="b"></span>
   <span id="c"></span>
   <span id="v"></span>

.. blend-usage:: window1d
