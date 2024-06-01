#!/usr/bin/env python

"""
ibaqpy is a python package to compute IBAQ values from peptide identifications and protein database
"""

from setuptools import setup, find_packages

version = "0.0.3"


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="ibaqpy",
    version=version,
    description="Python package to compute intensity base absolute expression values",
    author="Yasset Perez-Riverol",
    author_email="ypriverol@gmail.com",
    long_description=readme(),
    long_description_content_type="text/markdown",
    keywords="Proteomics, Label-free, absolute quantification",
    url="https://github.com/bigbio/ibaqpy/",
    download_url="https://github.com/bigbio/ibaqpy/",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "pyopenms",
        "scikit-learn",
        "numpy",
        "click",
        "pandas",
        "matplotlib",
        "qnorm",
        "seaborn",
        "typing_extensions",
        "combat",
    ],
    entry_points={"console_scripts": ["ibaqpy=ibaqpy.ibaqpyc:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
