Reporting and Errors
====================

Headers
-------

.. code-block:: c

   #include <blend_error.h>
   #include <blend_report.h>

Errors
------

``blend_error_code`` defines ``SUCCESS``, ``FAIL``, and ``WFUNC_ERROR``.

Reporting
---------

``BLEND_Report`` writes messages according to the active verbosity level.
Verbosity levels follow the GMT-style sequence ``q``, ``e``, ``w``, ``t``,
``i``, ``c``, and ``d``.

``blend_elapsed_seconds`` returns elapsed wall-clock seconds for module timing
messages. The ``t`` verbosity level reports timings for grid output, streamed
queries, blendfile support preparation, and best-IoU monotone conversion.
