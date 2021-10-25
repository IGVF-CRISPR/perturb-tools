# ![perturb-tools_logo](docs/images/perturb_tools_logo.svg)

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/perturb-tools.svg)](https://pypi.python.org/pypi/perturb-tools/)
[![PyPI version](https://badge.fury.io/py/perturb-tools.svg)](https://badge.fury.io/py/perturb-tools)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Analysis Framework for Pooled CRISPR Genome Editing Screens

```python
import perturb_tools as pt

screen = pt.Screen(X)
```
```
Genome Editing Screen composed of n_guides x n_conditions = 100 x 3
  guides: 'experiment', 'sequence', 'target'
  conditions: 'drug', 'control', 'initial'
  layers: 'raw_counts', 'Log2Norm_counts'
```


## Data Structure
This format and organization of metadata surrounding a multidimensional experiment is inspired by [AnnData](https://anndata.readthedocs.io/en/stable/), a useful solution for the analogous organization of single-cell data.
<br></br>
<img src="docs/images/screen_data.png" width="700"/>
