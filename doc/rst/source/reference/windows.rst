Window Functions
================

BLEND window functions are point evaluators. A function is evaluated at the
requested grid point or interpolation neighbor; BLEND does not need to build an
entire window array before answering a query.

Available names
---------------

Run:

.. code-block:: sh

   blend --show-windows

Current names include ``boxcar``, ``cosine``, ``trapezoid``, ``hamming``,
``blackman``, ``blackmanharris``, ``welch``, ``parzen``, ``gaussian``,
``smoothstep``, ``smootherstep``, ``exponential``, ``sine``, ``bohman``,
``nuttall``, ``kaiser``, ``cauchy``, ``quadratic``, ``cubic``, ``poisson``,
``bartlett``, ``barthann``, ``exactblackman``, ``blackmannuttall``,
``flattop``, ``lanczos``, ``riesz``, ``riemann``, ``fejer``, ``connes``,
``hanningpoisson``, ``kaiserbessel``, ``plancktaper``, ``quartic``,
``quintic``, ``septic``, ``nonic``, ``logistic``, ``tanh``, ``erf``,
``arctan``, ``gompertz``, ``softsign``, ``agnesi``, ``inversequadratic``,
``inversemultiquadric``, ``powerlaw``, ``root``, ``circular``, ``sech``,
``sech2``, ``student``, and ``laplace``.

Aliases
-------

``hann`` and ``tukey`` map to the cosine evaluator. ``linear`` maps to
``trapezoid``. ``normal`` maps to ``gaussian``.
