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
