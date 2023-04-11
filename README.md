# ![perturb-tools_logo](docs/images/perturb_tools_logo.svg)

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/perturb-tools.svg)](https://pypi.python.org/pypi/perturb-tools/)
[![PyPI version](https://badge.fury.io/py/perturb-tools.svg)](https://badge.fury.io/py/perturb-tools)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**perturb-tools** is an **analysis framework for pooled CRISPR genome-editing screens**. Thus far, development has focused on local (i.e., not genome-wide) tiling screens with specific phenotypic readouts though expansion of this scope is of interest. 




## Data Structure and Analysis Framework

```python
import perturb_tools as pt

```

Data format and organization of metadata surrounding a multidimensional experiment is inspired by [AnnData](https://anndata.readthedocs.io/en/stable/) and [MuData](https://mudata.readthedocs.io/en/latest/), a useful solution for the analogous organization of single-cell data.
<br></br>
<img src="docs/images/screendata.svg" width="600"/>

### The three main components of this data strcuture:

* **`screen.X`** (Numpy array)

* **`screen.obs`** (pandas DataFrame) of shape: `[conditions x condition_annotation]`

* **`screen.var`** (pandas DataFrame) of shape: `[guides x guide_annotation]`

See the [**tutorial**](notebooks/bulk/basic_api_demo.ipynb) for more information.


## Installation

**Install the development package**:
```BASH
# (1) clone this repository
git clone https://github.com/IGVF-CRISPR/perturb-tools.git

# (2) install the local project in editable mode
cd ./perturb-tools; pip install -e .
```


### Planning

We have a **[PR](https://github.com/IGVF-CRISPR/perturb-tools/pull/2)** for planning laying out modules and sub-modules for key functionalities required for the first-draft pipeline.

In general the structure could look something like the following:

```
IGVF-CRISPR/perturb-tools/
│
├── LICENSE
├── notebooks
│   ├── examples
│   
├── README.md
├── requirements.txt
├── setup.py
│
├── perturb-tools/
│   ├── __init__.py  
│   ├── _plotting/
│   │    ├── ...
│   │    
│   ├── _qc/
│   │    ├── ...
│   │     
│   ├── _external_tools/
│   │    ├── .../
│   │    
```

Adhering to the above structure (or some structure - as long as it follows somewhat closely to the [PEP8 style guide](https://www.python.org/dev/peps/pep-0008/) will enable seamless collaboration regardless of the module you're working on - this is especially important for code review. 

### Contributing

In order to stay organzed, let's all contribute through PRs. Think of a PR as the main topic you are planning to contribute (`qc` or `guide_counting`). If PRs and issues are foreign to you, just ask! The best way to learn git workflows is through doing. Every time we open a PR, we should organize sub-tasks as issues and link them to that PR. Conversations, feedback, and requested changes can all be mediated through the PRs. 

## General Analysis Steps
* See [**tutorial**](perturb_tools/notebooks) which includes:
  * Bulk screen data
    * [API tutorial with AnnData](perturb_tools/notebooks/anndata_demo.ipynb)
    * [API tutorial with Perturb-tools](perturb_tools/notebooks/bulk/basic_api_demo.ipynb)
      * Arithmetic: Calculating the mean, standard deviation, and log-fold change between/across replicates

