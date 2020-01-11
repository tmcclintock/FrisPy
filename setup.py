from setuptools import setup

dist = setup(name="FrisPy",
             author="Tom McClintock",
             author_email="thmsmcclintock@gmail.com",
             version="1.0.0",
             description="Model for a flying disc (frisbee) that produces simulated trajectories.",
             url="https://github.com/tmcclintock/FrisPy",
             packages=['FrisPy'],
             install_requires=['numpy','scipy']
)
