# API tutorial

[AnnData](https://anndata.readthedocs.io/en/latest/) can be used as the data structure that can naturally encode the CRISPR screen data where we have sgRNAs and experimental samples both with their own annotations.

`Screen` object is the shallow wrapper around AnnData with the CRISPR-specific methods including
* Log fold change calculation between conditions
* Log fold change aggregation across replicates
* Data export into analysis methods' input format

## Data Structures for CRISPR screen data
* [AnnData](anndata_demo.ipynb)
* [Screen](screen_demo.ipynb)
