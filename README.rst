.. |TRAVIS| image:: https://travis-ci.com/tmcclintock/FrisPy.svg?branch=v2
	    :target: https://travis-ci.com/tmcclintock/FrisPy
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
also changing the underlying physical modlel. This is useful for analyzing
the mechanics of a disc in terms of its design, as well as creating simulated
throws for things like disc launchers or other helpful tools.

This is a pure Python rebuild of the old FrisPy code, which included a version
of the integrator written in C for speed. To obtain a fast version of the
modeling code, either roll back to an old version or check out the
`Frisbee_Simulator <https://github.com/tmcclintock/Frisbee_Simulator>`_
repository.

Installation
------------

The easiest way to install this package is with ``pip``. The PyPI package can
be viewed `here <https://pypi.org/project/frispy/>`_.

.. code-block:: bash

   pip install frispy

To install from source, there are other steps involved.
First, you must obtain the code from Github. If you have
`git <https://git-scm.com/>`_ installed you can clone the repository from
the command line:

.. code-block:: bash

   git clone https://github.com/tmcclintock/FrisPy.git

or with the GitHub Desktop application. Once you have the code, change
into the directory and proceed.

Note, the only hard requirements for this package are ``python>=3.6``,
``numpy``, ``scipy``, and ``matplotlib`` (plotting only). Note that this package
uses the relatively recent
`scipy.integrate.solve_ivp
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp>`_
method, which may not exist in older versions of ``scipy``. If you have these
three packages, you can install this package with the ``setup.py`` file without
worrying about creating an environment.

From an Anaconda environment
----------------------------

The preferred method of installation is with
`anaconda
<https://docs.conda.io/projects/conda/en/latest/index.html>`_
You can install all the requirements into a compatible environment called
``frispy`` by running the following command:

.. code-block:: bash

   conda env create -f environment.yml

You can then install the package the usual way

.. code-block:: bash

   python setup.py install

You can also use ``pip`` to install the requirements from the
``requirements.txt`` file by running:

.. code-block:: bash

   pip install -r requirements.txt

Then follow this by using the ``setup.py`` file to install.

Testing
-------

Verify your installation by running:

.. code-block:: bash

   pytest

Please report any problems you encounter on the `issues page
<https://github.com/tmcclintock/FrisPy/issues>`_. Thank you!

Running
-------

Check out ``example.py`` to see how to run and view results.
In words, you create a disc and compute its trajectory.

.. code-block:: python

   from frispy import Disc

   disc = Disc()
   result = disc.compute_trajectory()
   times = result.times
   x, y, z = result.x, result.y, result.z

Once you have a trajectory, you can use that to create visualizations. For
instance, to plot the height of the disc against one of its horizontal
coordintes (``x``), you can run:

.. code-block:: python

   import matplotlib.pyplot as plt

   plt.plot(x, z)
   plt.show()

Soon
----

There are some big upgrades on the horizon! Stay tuned for:

- animated trajectories
- documentation
- example Jupyter notebooks
- plotting routines
