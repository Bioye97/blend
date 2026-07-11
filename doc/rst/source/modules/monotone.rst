monotone
========

Synopsis
--------

.. code-block:: text

   blend monotone [<polygonfile>] [-M<method>] [-G<nx>[/<ny>]] [-V[q|e|w|t|i|c|d]]

Description
-----------

``monotone`` reads a two-column polygon vertex file and checks whether the
polygon is valid, simple, and xy-monotone. If no polygon file is supplied,
vertices are read from standard input.

Without ``-M``, the module reports whether the polygon is xy-monotone. With
``-M``, non-xy-monotone polygons are modified and written to standard output.
Messages, including the original and final vertex counts, are written to
standard error.

Options
-------

``-M<method>``, ``--monotone=<method>``
   Modify non-xy-monotone polygons. ``-Me`` uses the envelope method. ``-Mb``
   uses the highest-IoU piecewise envelope from both traversal directions.

``-G<nx>[/<ny>]``, ``--grid=<nx>[/<ny>]``
   Set the sampling grid used by ``-Mb``. If ``ny`` is omitted, ``ny = nx``.
   The default is ``256/256``.

``-V[level]``
   Set verbosity.

Examples
--------

.. code-block:: sh

   blend monotone polygon.txt
   blend monotone polygon.txt -Mb -G512/512 > polygon_monotone.txt
