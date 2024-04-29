from setuptools import setup, find_packages

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

setup(name="ChemInterface", version='1.0', packages=find_packages(), install_requires=requirements,)
