import os

PATH_ROOT = os.path.dirname(__file__)
from setuptools import setup

import frispy  # noqa: E402

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="frispy",
    author="Tom McClintock",
    author_email="thmsmcclintock@gmail.com",
    version=frispy.__version__,
    description=frispy.__docs__,
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/tmcclintock/FrisPy",
    packages=["frispy"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
