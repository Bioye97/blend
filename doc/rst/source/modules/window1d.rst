window1d
========

Synopsis
--------

.. code-block:: text

   blend window1d -R<xmin>/<xmax> -I<dx> [-F<function>] [-T<r1>[/<r2>]]
                  [-B<blendfile>] [-C<f|l|o|u|a|g|p>] [-V[q|e|w|t|i|c|d]]

Description
-----------

``window1d`` returns weights on a 1-D grid or at query coordinates read from
standard input. Query coordinates do not need to land exactly on grid points;
BLEND uses linear interpolation from neighboring grid weights.

Options
-------

``-R<xmin>/<xmax>``, ``--region=<xmin>/<xmax>``
   Set the output domain.

``-I<dx>``, ``--increment=<dx>``
   Set the grid increment.

``-F<function>``, ``--function=<function>``
   Set the window taper function. Run ``blend --show-windows`` to list names.

``-T<r1>[/<r2>]``, ``--taper_ratio=<r1>[/<r2>]``
   Set beginning and ending taper ratios.

``-B<blendfile>``, ``--blendfile=<blendfile>``
   Read interval supports from a blendfile. Each row contains:
   ``left right function r1/r2``.

``-C<mode>``, ``--clobber=<mode>``
   Set overlap handling for blendfile supports. Modes are ``f``, ``l``, ``o``,
   ``u``, ``a``, ``g``, and ``p``.

Example
-------

.. code-block:: sh

   blend window1d -R40/60 -I0.5 -Ftrapezoid -T0.2/0.2 > weights.txt
