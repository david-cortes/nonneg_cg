# Conjugate-gradient with non-negativity constraints

This is a package for minimization of differentiable functions subject to non-negativity constraints on all the variables, i.e.
```
min f(x)
s.t. x[i] >= 0  i = 1..n
```

It uses a modified Polak-Ribiere-Polyak method, which is described in the paper "A conjugate gradient type method for the nonnegative constraints optimization problems" (see reference at the end). Implementation is in C with wrappers for Python and R.

This algorithm tends to require more function evaluations than others, but fewer gradient evaluations. It always produces feasible descent directions.

## Installation

* Python:

Windows:
```
git clone https://www.github.com/david-cortes/nonneg_cg.git
cd nonneg_cg
python setup.py install
```
(Could also be installed from `pip` depending in NumPy version)

**Note: do NOT compile with MINGW in Windows**

Linux and Mac
```
pip install nonnegcg
```

* R:
```
devtools::install_github("david-cortes/nonneg_cg")
```
(Can also be installed with `install.packages("nonneg.cg")`, but it won't have multithreading support)

* C/C++:
```
git clone https://www.github.com/david-cortes/nonneg_cg.git
cd nonneg_cg
mkdir build
cd build
cmake ..
make

### for a system-wide install in linux
sudo make install
sudo ldconfig
```

## Usage

* Python: see example notebook [here](https://github.com/david-cortes/nonneg_cg/blob/master/example/nncg.ipynb)
* R: see example in the documentation (`help(minimize.nonneg.cg)`)
* C/C++: see example file [here](https://github.com/david-cortes/nonneg_cg/blob/master/example/c_rosenbrock.c) (- link with `-lnonnegcg`)

Simple example in Python
```python
import numpy as np
from scipy.optimize import rosen, rosen_der
from nonnegcg import minimize_nncg

x0  = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
res = minimize_nncg(x0, rosen, rosen_der)
res.x
```

The C function can also be used directly from Rcpp (through `R_GetCCallable`) and Cython (see `.pxd` file for the `cimport` - Python package in addition installs a `.h` header for usage directly in C by linking to the shared object in `nonnegcg.minimize_nncg.__file__`)

## Documentation

Package has only one function, documentation is available as docstring in Python (`help(nonnegcg.minimize_nncg)`), as regular documentation in R (`help(minimize.nonneg.cg)`), and as comment in the header file in C (link [here](https://github.com/david-cortes/nonneg_cg/blob/master/include/nonnegcg.h)).

## References
* Li, C. (2013). A conjugate gradient type method for the nonnegative constraints optimization problems. Journal of Applied Mathematics, 2013.
