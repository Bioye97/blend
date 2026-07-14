monotone
========

Description
-----------

``monotone`` reads a two-column polygon vertex file and checks whether the
polygon is valid, simple, and xy-monotone. If no polygon file is supplied,
vertices are read from standard input.

Without ``-M``, the module reports whether the polygon is xy-monotone. With
``-M``, non-xy-monotone polygons are modified and written to standard output.
Messages, including the original and final vertex counts, are written to
standard error.

Usage
-----

.. raw:: html

   <span id="m"></span>
   <span id="g"></span>
   <span id="v"></span>

.. blend-usage:: monotone
