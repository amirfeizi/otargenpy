#!/usr/bin/env python
"""The setup script."""

from setuptools import setup, find_packages

with open("README.md") as f:
    readme = f.read()

setup(
    name="otargenpy",
    version="2.0.1",
    author="Amir Feizi",
    author_email="afeizi@gmail.com",
    python_requires=">=3.8",
    description="Tidy Python interface to the Open Targets Platform GraphQL API",
    long_description=readme,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/amirfeizi/otargenpy/",
    packages=find_packages(include=["otargenpy", "otargenpy.*"]),
    install_requires=[
        "requests",
        "pandas",
        "matplotlib",
        "numpy",
    ],
    extras_require={
        "dev": ["pytest>=7", "flake8"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="bioinformatics open-targets genetics gwas drug-targets graphql",
    zip_safe=False,
)
