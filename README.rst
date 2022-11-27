.. |TRAVIS| image:: https://github.com/tmcclintock/FrisPy/workflows/Build%20Status/badge.svg?branch=master
	    :target: https://github.com/tmcclintock/FrisPy/actions
.. |COVERALLS| image:: https://coveralls.io/repos/github/tmcclintock/FrisPy/badge.svg?branch=master
	       :target: https://coveralls.io/github/tmcclintock/FrisPy?branch=master
.. |LICENSE| image:: https://img.shields.io/badge/License-MIT-yellow.svg
	     :target: https://opensource.org/licenses/MIT

|TRAVIS| |COVERALLS| |LICENSE|

FrisPy
======

Documentation for ``FrisPy`` package can be `found here on RTD
<https://frispy.readthedocs.io/en/latest/>`_.

This repository contains a physical model for a flying disc. Using this code,
one can simulate trajectories of discs with varying initial conditions, while
also changing the underlying physical model. This is useful for analyzing
the mechanics of a disc in terms of its design, as well as creating simulated
throws for things like disc launchers or other helpful tools.

This is a pure Python rebuild of the old FrisPy code, which included a version
of the integrator written in C for speed. Find the fast C simulation in the
`Frisbee_Simulator <https://github.com/tmcclintock/Frisbee_Simulator>`_
repository.

The earliest implementation of this model that I could find was by Sara Ann Hummel
for their 2003 Masters thesis for UC Davis.  You can find the document in full
`on this page <https://morleyfielddgc.files.wordpress.com/2009/04/hummelthesis.pdf>`_.

Installation
------------

The easiest way to install this package is with ``pip``. The PyPI package can
be viewed `here <https://pypi.org/project/frispy/>`_.

.. code-block:: bash

   pip install frispy


For developers
--------------

Development should be performed using `poetry <https://python-poetry.org/>`_ to handle
the development environment. Once poetry is installed, you can install the environment,
which will include ``frispy``:

.. code-block:: bash

   poetry install

All proceeding instructions assume you entered your virtual environment using ``poetry shell``,
otherwise prepend ``poetry run`` to all instructions.

If you intend to open a pull request, please make sure ``pre-commit`` is installed
before committing to your branch:

.. code-block:: bash

   pre-commit install

This will ensure that the code you submit is PEP8 compliant. Otherwise, CI checks will
fail before merging can be completed.

Verify your installation by running:

.. code-block:: bash

   pytest

Please report any problems you encounter on the `issues page
<https://github.com/tmcclintock/FrisPy/issues>`_. Thank you!
