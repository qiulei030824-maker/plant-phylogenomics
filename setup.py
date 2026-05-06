#!/usr/bin/env python3
"""Setup script for plant-phylogenomics."""

from setuptools import setup, find_packages

setup(
    name="plant-phylogenomics",
    version="0.1.0",
    description="Modular pipeline for plant phylogenomics analysis",
    author="Qiulei",
    url="https://github.com/qiulei030824-maker/plant-phylogenomics",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "biopython>=1.79",
        "pyyaml>=5.4",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
