# Example

## Code snippet

You have a data set such as :

```
import superexacttestpy as stest
Set1 = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q"]
Set2 = ["L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
Set3 = ["H","I","J","K","L","M","N","O","P","Q"]
data = [Set1,Set2,Set3]
names = ["Set1","Set2","Set3"]
background_size = 1000
```

If you just want the df with the results, you can use the function `s.tl.supertest()`

```python
s.tl.supertest(data=data, n=background_size, names=names)
```

If you want to get the df and plot the results, you can use the function `s.pl.plot()`

```python
s.pl.plot(data=data, n=background_size, names=names)
```

```{image} _static/output_example/supertest.png

```

## Tutorials

```{toctree}
:hidden: true
:maxdepth: 1

notebooks/01_dpsets_function.ipynb
notebooks/02_intersection.ipynb
notebooks/03_cp_sets.ipynb
notebooks/04_msets.ipynb
notebooks/05_supertest.ipynb
notebooks/06_plot.ipynb
```
