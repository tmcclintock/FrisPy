# FrisPy [![Build Status](https://travis-ci.com/tmcclintock/FrisPy.svg?branch=v2)](https://travis-ci.com/tmcclintock/FrisPy)

This repository contains a physical model for a flying disc. Using this code, one can simulate trajectories of discs with varying initial conditions, while also allowing the physical model to be changed. This is useful for analyzing the mechanics of a disc in terms of its design, as well as creating simulated throws for things like disc launchers or other helpful tools.

This is a pure Python rebuild of the old FrisPy code, which included a version of the integrator written in C for speed. To obtain a fast version of the modeling code, either roll back to an old version or check out the [Frisbee_Simulator](https://github.com/tmcclintock/Frisbee_Simulator) repository.

## Installation

To install, at the command line just type:

```bash
python setup.py install
```

## Running

Check out `example.py` for an example of how to run and view results. It boils down to doing:
```python
disc = FrisPy.create_disc(filename = "initial_conditions_filename.txt")
times, trajectory = FrisPy.get_trajectory(disc)
```
and that's it.