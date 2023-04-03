#!/usr/bin/env python

from setuptools import setup, find_packages


setup(
    name="sccrispr_tools",
    version="0.0.1",
    description="Utility package for single-cell CRISPR screens with paired single-cell omics and guide counts",
    author=["IGVF CRISPR FG", "Jayoung Ryu"],
    author_email=["", "jayoung_ryu@g.harvard.edu"],
    url="https://github.com/IGVF-CRISPR/sccrispr-tools",
    packages=find_packages(),
)
