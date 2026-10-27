#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="spatial-sanctuary",
    version="2.2.1",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={"spatial_sanctuary": ["data/*.csv", "data/*.tsv"]},
    include_package_data=True,
    install_requires=[
        "numpy>=1.22.0",
        "scipy>=1.9.0",
        "pandas>=1.5.0",
        "lifelines>=0.27.0",
        "matplotlib>=3.6.0"
    ],
    entry_points={
        "console_scripts": [
            "spatial-sanctuary=spatial_sanctuary.cli:main",
            "osimertinib-sanctuary=spatial_sanctuary.cli:main",
        ]
    },
)
