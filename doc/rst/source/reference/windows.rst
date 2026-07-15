Window Functions
================

BLEND window functions are point evaluators. A function is evaluated at the
requested grid point or interpolation neighbor; BLEND does not need to build
an entire window array before responding to a query.

Available Names
---------------

Run:

.. code-block:: sh

   blend --show-windows

Each page below shows the ``window1d`` command used to generate the
example for that window function.

Available Window Functions
--------------------------

.. toctree::
   :maxdepth: 1

   windows/boxcar
   windows/cosine
   windows/hann
   windows/tukey
   windows/trapezoid
   windows/linear
   windows/hamming
   windows/blackman
   windows/blackmanharris
   windows/welch
   windows/parzen
   windows/gaussian
   windows/smoothstep
   windows/smootherstep
   windows/exponential
   windows/sine
   windows/bohman
   windows/nuttall
   windows/kaiser
   windows/cauchy
   windows/quadratic
   windows/cubic
   windows/poisson
   windows/bartlett
   windows/barthann
   windows/bartletthann
   windows/exactblackman
   windows/blackmannuttall
   windows/flattop
   windows/lanczos
   windows/riesz
   windows/riemann
   windows/fejer
   windows/connes
   windows/hanningpoisson
   windows/kaiserbessel
   windows/plancktaper
   windows/quartic
   windows/quintic
   windows/septic
   windows/nonic
   windows/logistic
   windows/tanh
   windows/erf
   windows/arctan
   windows/gompertz
   windows/softsign
   windows/agnesi
   windows/inversequadratic
   windows/inversemultiquadric
   windows/powerlaw
   windows/root
   windows/circular
   windows/sech
   windows/sech2
   windows/student
   windows/laplace

Aliases
-------

``hann`` and ``tukey`` map to the cosine evaluator. ``linear`` maps to
``trapezoid``. ``normal`` maps to ``gaussian``.
