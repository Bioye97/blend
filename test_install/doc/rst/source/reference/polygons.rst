Polygons
========

BLEND polygon windows require valid simple polygons. A polygon file is a text
file with two columns:

.. code-block:: text

   x y

Blank lines and comments are ignored in module input where supported.

xy-monotone requirement
-----------------------

The 2-D and 3-D window modules use row and column stepping vectors. Before
assembly, BLEND maps each support polygon to its local support grid and snaps
vertices to grid nodes. For assembly to be well-defined, that grid-snapped
support must be xy-monotone under BLEND's current boundary model. The definition
required for the ``window2d`` and ``window3d`` modules is a ``strictly increasing``
condition as opposed to the usual ``nondecreasing`` condition in the standard 
definition of xy-monotonicity.

Use:

.. code-block:: sh

   blend monotone polygon.txt

to check a polygon for xy-monotonicity.

Conversion methods
------------------

``-Me``
   Full xy-monotone envelope, i.e., convex hull.

``-Mb``
   Highest-IoU(Intersection Over Union) piecewise envelope, evaluated over a sampling grid.

``-ME``
   Strict envelope for polygons that must satisfy BLEND's current window
   boundary assembly.

``-MB``
   Strict highest-IoU piecewise envelope for polygons that must satisfy
   BLEND's current window boundary assembly.
