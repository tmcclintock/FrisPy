# FrisPy [![Build Status](https://travis-ci.com/tmcclintock/FrisPy.svg?branch=v2)](https://travis-ci.com/tmcclintock/FrisPy)

This repository contains a physical model for a flying disc. Using this code,
one can simulate trajectories of discs with varying initial conditions, while
also changing the underlying physical modlel. This is useful for analyzing
the mechanics of a disc in terms of its design, as well as creating simulated
throws for things like disc launchers or other helpful tools.

This is a pure Python rebuild of the old FrisPy code, which included a version
of the integrator written in C for speed. To obtain a fast version of the
modeling code, either roll back to an old version or check out the
[Frisbee_Simulator](https://github.com/tmcclintock/Frisbee_Simulator)
repository.

## Installation

The easiest way to install this package is with ``pip``:

```bash
pip install frispy
```

To install from source, there are other steps involved.
First, you must obtain the code from Github. If you have
[`git`](https://git-scm.com/) installed you can clone the repository from
the command line:
```bash
git https://github.com/tmcclintock/FrisPy.git
```
or with the GitHub Desktop application. Once you have the code, change
into the directory and proceed.

Note, the only hard requirements for this package are `python>=3.6`,
`numpy`, `scipy`, and `matplotlib` (plotting only). Note that this package
uses the relatively recent
[`scipy.integrate.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp)
method, which may not exist in older versions of `scipy`. If you have these
three packages, you can install this package with the `setup.py` file without
worrying about creating an environment.

### From an Anaconda environment
The preferred method of installation is with
[anaconda](https://docs.conda.io/projects/conda/en/latest/index.html).
You can install all the requirements into a compatible environment called
`frisby` by running the following command:
```bash
conda create env -f environment.yml
```
You can then install the package the usual way
```bash
python setup.py install
```

### With pip from the requirements file

If you prefer not to use conda, you can use
[`pip`](https://pypi.org/project/pip/) to install the requirements by
running:
```bash
pip install -r requirements.txt
```
And then you can install the package normally
```bash
python setup.py install
```

Verify your installation with
```bash
pytest
```

## Running

Check out `example.py` for an example of how to run and view results.
In words, you create a disc and compute its trajectory.
```python
from frispy import Disc

disc = Disc()
result = disc.compute_trajectory()
times = result.times
x, y, z = result.x, result.y, result.z
```
Once you have a trajectory, you can use that to create visualizations. For
instance, to plot the height of the disc against one of its horizontal
coordintes (`x`), you can run:
```python
import matplotlib.pyplot as plt

plt.plot(x, z)
plt.show()
```

## Soon

This rebuild is a work in progress. Check back soon for:

- animated trajectories
- documentation
- travis builds
- example Jupyter notebooks
- plotting routines
