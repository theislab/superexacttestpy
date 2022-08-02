# superexacttestpy

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/workflow/status/ilibarra/superexacttestpy/Test/main
[link-tests]: https://github.com/theislab/superexacttestpy/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/superexacttestpy

Python implementation of the SuperExactTest algorithm

## Getting started

Please refer to the [documentation][link-docs]. In particular, the

-   [API documentation][link-api].

## What is superexacttestpy ?

Superextractestpy is a python reimplementation of the R package [SuperExactTest][r-package] allowing to perform tests on the statistical distribution as well as to visualize multiset intersection.

This algorithm calculates the intersection probability of a large number of genes in a genetic set with linear complexity.

### How to use it?

Import the package

```python
import superexacttestpy as stest
```

For example, we want to make the test on this fictive set:

```python
Set1 = [
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
]
Set2 = ["L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
Set3 = ["H", "I", "J", "K", "L", "M", "N", "O", "P", "Q"]

data = [Set1, Set2, Set3]
names = ["Set1", "Set2", "Set3"]

background_size = 1000
```

If you just want the df with the results, you can use the function `stest.tl.supertest()`

```python
stest.tl.supertest(data=data, n=background_size, names=names)
```

The function supertest has some optional arguments:

-   degree: the degree of the intersection you want to compute.
-   lower_tail: Let m be the number of elements shared in the sets : if True, p = P[overlap < m] and if False, p = P[overlap >= m].

If you want to get the df and plot the results, you can use the function `stest.pl.plot()`

```python
stest.pl.plot(data=data, n=background_size, names=names)
```

The function plot has some optional arguments:

-   degree: the degree of the intersection you want to compute.
-   sort_by: on what you want to sort the bars "degree" or "p_val"
-   show_count: if True, the number of genes in the intersection is shown.
-   color_p_val: if True, the bars are colored by their p-value.
-   size: tuple of the figsize
-   background_color: the color of the background of the plot.

### Side function

Some sides functions are available and will be described in the [readthedocs][link-api]

## Installation

You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, we recommend installing `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`\_.

There are several alternative options to install superexacttestpy:

<!--
1) Install the latest release of `superexacttestpy` from `PyPI <https://pypi.org/project/superexacttestpy/>`_:

```bash
pip install superexacttestpy
```
-->

1. Install the latest development version:

```bash
pip install git+https://github.com/ilibarra/superexacttestpy.git@main
```

## Release notes

See the [changelog][changelog].

## Contact

For questions and help requests, you can reach out in the [scverse discourse][scverse-discourse].
If you found a bug, please use the [issue tracker][issue-tracker].

## Citation

for cite superexactestpy, please use the following:

```bibtex
@software{superexacttest,
  author = {Ibarra, Mauger-Birocheau}},
  doi = {},
  month = {},
  title = {{superexacttest}},
  url = {https://github.com/theislab/superexacttestpy},
  year = {2022}
}
```

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/ilibarra/superexacttestpy/issues
[changelog]: https://superexacttestpy.readthedocs.io/latest/changelog.html
[link-docs]: https://superexacttestpy.readthedocs.io
[link-api]: https://superexacttestpy.readthedocs.io/latest/api.html
[r-package]: https://github.com/mw201608/SuperExactTest
