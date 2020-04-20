from setuptools import setup

dist = setup(name="FrisPy",
             author="Tom McClintock",
             author_email="thmsmcclintock@gmail.com",
             version="2.0.0",
             description="Simulated flying discs",
             url="https://github.com/tmcclintock/FrisPy",
             packages=['FrisPy'],
             install_requires=['numpy','scipy']
)
