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
Genome Editing Screen composed of: n_guides x n_conditions = 946 x 12

   guides:    'barcode', 'barcode_id', 'experiment', 'sequence', 'target_id', 'pred_ABE_edit', 'pred_CBE_edit'
   condit:    'conditions'
   condit_m:  'barcode_counts', 'unexpected_sequences'
   condit_p:  'correlation'
   layers:    'X_lognorm'
   uns:       'run_info', 'poolq3', 'metadata', 'SampleBarcodeReadCounts', 'CommonSampleBarcodeReadCounts'
```


## Data Structure
This format and organization of metadata surrounding a multidimensional experiment is inspired by [AnnData](https://anndata.readthedocs.io/en/stable/), a useful solution for the analogous organization of single-cell data.
<br></br>
<img src="docs/images/screendata.svg" width="700"/>
