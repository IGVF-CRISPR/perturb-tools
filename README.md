# ![perturb-tools_logo](images/perturb-tools_logo.png)

Analysis Framework for Pooled CRISPR Genome Editing Screens

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
