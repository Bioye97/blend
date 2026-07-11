Window API
==========

Header
------

.. code-block:: c

   #include <blend_window.h>

The window API defines ``blend_window_function``, ``window``, window-name
lookup, and point-evaluator functions such as ``cosine``, ``trapezoid``, and
``window_function``.

The ``window`` structure stores support geometry, grid dimensions, taper ratios,
selected window functions, assembled stepping vectors, and the current
contribution value.
