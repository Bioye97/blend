Polygons
========

BLEND polygon windows require valid simple polygons. A polygon file is a text
file with two columns:

.. code-block:: text

   x y

Blank lines and comments are ignored in module input where supported.

xy-monotone requirement
-----------------------

The 2-D and 3-D window modules use row and column stepping vectors. For that
assembly to be well-defined, polygon supports must be xy-monotone under BLEND's
current boundary model.

Use:

.. code-block:: sh

   blend monotone polygon.txt

to check a polygon.

Conversion methods
------------------

``-Me``
   Full xy-monotone envelope.

``-Mb``
   Highest-IoU piecewise envelope, evaluated over a sampling grid.
