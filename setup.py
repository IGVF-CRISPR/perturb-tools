from setuptools import setup
import re
import os
import sys


setup(
    name="perturb-tools",
    version="0.0.12",
    python_requires=">3.7.0",
    author="Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard",
    author_email="mvinyard@broadinstitute.org",
    url="https://github.com/pinellolab/perturb-tools",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="perturb-tools - Analysis Framework for Pooled CRISPR Genome Editing Screens.",
    packages=[
        "perturb_tools",
        "perturb_tools._arithmetic",
        "perturb_tools._experimental_design",
        "perturb_tools._external_tools",
        "perturb_tools._framework",
        "perturb_tools._normalization",
	"perturb_tools._plotting",
	"perturb_tools._readwrite",
	"perturb_tools._utilities",
    ],
    
    install_requires=[
        "matplotlib>=3.4",
        "anndata>=0.7.1",
        "numpy>=1.19.2",
        "pandas>=1.1.2",
        "biopython>=1.79",
        "pool-sharq>=0.0.12",
	"plotly",
    "regex"
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.8",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
