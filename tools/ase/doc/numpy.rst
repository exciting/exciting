.. _numpy:

Numeric arrays in Python
========================

Links to NumPy's webpage:

* `Numpy and Scipy Documentation`_
* `Numpy guide book <http://www.tramy.us/numpybook.pdf>`_
* `Numpy example list`_
* `Numpy functions by category`_


.. _Numpy and Scipy Documentation: http://docs.scipy.org/doc
.. _Numpy example list: http://www.scipy.org/Numpy_Example_List_With_Doc
.. _Numpy functions by category:
                        http://www.scipy.org/Numpy_Functions_by_Category


ASE makes heavy use of an extension to Python called NumPy.  The
NumPy module defines an :term:`ndarray` type that can hold large arrays of
uniform multidimensional numeric data.  An array is similar to a
``list`` or a ``tuple``, but it is a lot more powerful and efficient.

XXX More examples from everyday ASE-life here ...

>>> import numpy as np
>>> a = np.zeros((3, 2))
>>> a[:, 1] = 1.0
>>> a[1] = 2.0
>>> a
array([[ 0.,  1.],
       [ 2.,  2.],
       [ 0.,  1.]])
>>> a.shape
(3, 2)
>>> a.ndim
2


The conventions of numpy's linear algebra package:

>>> import numpy as np
>>> 
>>> # Make a random hermitian matrix, H
>>> H = np.random.rand(6, 6) + 1.j * np.random.rand(6, 6)
>>> H = H + H.T.conj()
>>> 
>>> # Determine eigenvalues and rotation matrix
>>> eps, U = np.linalg.eigh(H)
>>> 
>>> # Sort eigenvalues
>>> sorted_indices = eps.real.argsort()
>>> eps = eps[sorted_indices]
>>> U = U[:, sorted_indices]
>>> 
>>> # Make print of numpy arrays less messy:
>>> np.set_printoptions(precision=3, suppress=True)
>>> 
>>> # Check that U diagonalizes H:
>>> print np.dot(np.dot(U.T.conj(), H), U) - np.diag(eps)
>>> print np.allclose(np.dot(np.dot(U.T.conj(), H), U), np.diag(eps))
>>> 
>>> # The eigenvectors of H are the *coloumns* of U:
>>> np.allclose(np.dot(H, U[:, 3]), eps[3] * U[:, 3])
>>> np.allclose(np.dot(H, U), eps * U)

The rules for multiplying 1D arrays with 2D arrays:

* 1D arrays and treated like shape (1, N) arrays (row vectors).
* left and right multiplications are treated identically.
* A length `m` *row* vector can be multiplied with an `n x m`
  matrix, producing the same result as if replaced by a matrix with
  `n` copies of the vector as rows.
* A length `n` *column* vector can be multiplied with an `n x m`
  matrix, producing the same result as if replaced by a matrix with
  `m` copies of the vector as columns.

Thus, for the arrays below:

>>> M = np.arange(5 * 6).reshape(5, 6) # A matrix af shape (5, 6)
>>> v5 = np.arange(5) + 10             # A vector of length 5
>>> v51 = v5[:, None]                  # A length 5 coulomn vector
>>> v6 = np.arange(6) - 12             # A vector of length 6
>>> v16 = v6[None, :]                  # A length 6 row vector

The following identities hold::

  v6 * M == v16 * M == M * v6 == M * v16 == M * v16.repeat(5, 0)
  v51 * M == M * v51 == M * v51.repeat(6, 1)

The exact same rules apply for adding and subtracting 1D arrays to /
from 2D arrays.
