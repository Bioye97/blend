docs
====

Synopsis
--------

.. code-block:: text

   blend docs [-Q] [-V[q|e|w|t|i|c|d]] <module-name> [<-option>]

Description
-----------

``docs`` reports documentation for a BLEND module. The current implementation
prints module documentation targets and usage text; the documentation tree added
in this manual is the backend that can later be connected to option-specific
lookup.

Required Arguments
------------------

``<module-name>``
   One of ``docs``, ``monotone``, ``window1d``, ``window2d``, or ``window3d``.

Optional Arguments
------------------

``-Q``
   Print the documentation target without opening a viewer.

``<-option>``
   Reserved for option-specific documentation.

``-V[level]``
   Set verbosity.
