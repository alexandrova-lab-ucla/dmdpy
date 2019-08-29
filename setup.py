
from setuptools import setup, find_packages

setup(
    name="dmdpy",
    version="1.0.0",
    author="Matthew Hennefarth",
    packages=find_packages(),
    package_data={'': ['.json', '*.j2', '*.config']}
)