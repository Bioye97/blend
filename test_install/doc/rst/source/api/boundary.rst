Boundary API
============

Header
------

.. code-block:: c

   #include <blend_boundary.h>

The boundary API contains the lower-level polygon boundary assembly routines
used by 2-D and 3-D windows. It defines ``permuted_vertex`` and functions for
assembling bottom, left, top, right, and corner-sector vertex groups.

Most users should work through the module commands or higher-level polygon and
window APIs.
