Installation
============

Dependencies
------------

To build BLEND, install:

* CMake 3.15 or newer
* A C compiler
* A POSIX shell

To build the documentation, install:

* Sphinx
* The Read the Docs Sphinx theme, provided by ``sphinx-rtd-theme``
* LaTeX and ``latexmk`` for PDF output

Install the theme in the same Python environment that provides
``sphinx-build``:

.. code-block:: sh

   python -m pip install sphinx-rtd-theme

or, with conda:

.. code-block:: sh

   conda install -c conda-forge sphinx_rtd_theme

Configuration
-------------

BLEND follows a GMT-like configuration style. Defaults live in
``cmake/ConfigDefault.cmake``. Local user choices live in
``cmake/ConfigUser.cmake``.

Create the user configuration file:

.. code-block:: sh

   cp cmake/ConfigUserTemplate.cmake cmake/ConfigUser.cmake

Common settings include:

.. code-block:: cmake

   set (CMAKE_INSTALL_PREFIX "/path/to/blend/install")
   set (BLEND_BUILD_TESTS ON)
   set (BLEND_BUILD_EXAMPLES ON)
   set (BLEND_BUILD_DOCS ON)

Build
-----

Configure from a separate build directory:

.. code-block:: sh

   mkdir build
   cd build
   cmake ..
   cmake --build .

Run Tests
---------

If ``BLEND_BUILD_TESTS`` is ``ON``:

.. code-block:: sh

   ctest

Build Documentation
-------------------

If ``BLEND_BUILD_DOCS`` is ``ON``:

.. code-block:: sh

   cmake --build . --target docs_html
   cmake --build . --target docs_pdf
   cmake --build . --target docs

The HTML manual is written to ``doc/html`` inside the build directory. The PDF
manual is copied into the HTML tree as ``_downloads/blend.pdf`` so the HTML
front page can link to it.

Install
-------

.. code-block:: sh

   cmake --build . --target install

This installs the command-line program, libraries, public headers, CMake
package files, examples, and documentation files.
