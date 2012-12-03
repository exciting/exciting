import os
import numpy as np

from UserDict import DictMixin

# -------------------------------------------------------------------

class MemoryBase(object, DictMixin):
    """Virtual memory (VM) statistics of the current process
    obtained from the relevant entries in /proc/<pid>/status:
    VmPeak       Peak virtual memory size in bytes.
    VmLck        ???
    VmHWM        Peak resident set size ("high water mark") in bytes.
    VmRSS        Resident memory usage in bytes.
    VmSize       VM usage of the entire process in bytes.
    VmData       VM usage of heap in bytes.
    VmStk        VM usage of stack in bytes.
    VmExe        VM usage of exe's and statically linked libraries in bytes.
    VmLib        VM usage of dynamically linked libraries in bytes.
    VmPTE        ???

    Note that VmSize > VmData + VmStk + VmExe + VmLib due to overhead.
    """

    _scale = {'KB':1024.0, 'MB':1024.0**2}
    _keys = ('VmPeak', 'VmLck', 'VmHWM', 'VmRSS', 'VmSize', 'VmData', \
            'VmStk', 'VmExe', 'VmLib', 'VmPTE')

    def __init__(self, verbose=0):
        self.verbose = verbose
        if self.verbose>=2: print 'MemoryBase.__init__'
        object.__init__(self)
        self._values = np.empty(len(self._keys), dtype=np.float)

    def __repr__(self):
        """Return a representation of recorded VM statistics.
        x.__repr__() <==> repr(x)"""
        if self.verbose>=2: print 'MemoryBase.__repr__'
        s = object.__repr__(self)
        w = max(map(len, self._keys))
        unit = 'MB'
        for k,v in self.iteritems():
            res = '<N/A>'
            if not np.isnan(v):
                res = '%8.3f %s' % (v/self._scale[unit], unit)
            s += '\n\t' + k.ljust(w) + ': ' + res.rjust(8)
        return s

    def __len__(self):
        """Number of VM keys which have not been outdated.
        x.__len__() <==> len(x)"""
        if self.verbose>=3: print 'MemoryBase.__len__'
        return np.sum(~np.isnan(self._values))

    def __getitem__(self, key):
        """Return floating point number associated with a VM key.
        x.__getitem__(y) <==> x[y]"""
        if self.verbose>=2: print 'MemoryBase.__getitem__'
        if key not in self:
            raise KeyError(key)
        i = self.keys().index(key)
        return self._values[i]

    def __setitem__(self, key, value):
        """x.__setitem__(i, y) <==> x[i]=y"""
        if self.verbose>=2: print 'MemoryBase.__setitem__'
        raise Exception('Virtual member function.')

    def __delitem__(self, key):
        """x.__delitem__(y) <==> del x[y]"""
        if self.verbose>=2: print 'MemoryBase.__delitem__'
        raise Exception('Virtual member function.')

    def clear(self):
        """D.clear() -> None.  Remove all items from D."""
        if self.verbose>=1: print 'MemoryBase.clear'
        raise Exception('Virtual member function.')

    def update(self, other=None):
        """D.update(E) -> None.  Update D from E: for k in E.keys(): D[k] = E[k]"""
        if self.verbose>=1: print 'MemoryBase.update'
        DictMixin.update(self, other)

    def copy(self):
        """Return a shallow copy of a VM statistics instance.
        D.copy() -> a shallow copy of D"""
        if self.verbose>=1: print 'MemoryBase.copy'
        res = object.__new__(self.__class__)
        MemoryBase.__init__(res, self.verbose)
        DictMixin.update(res, self)
        return res

    def has_key(self, key): #necessary to avoid infinite recursion
        """Return boolean to indicate whether key is a supported VM key.
        D.has_key(k) -> True if D has a key k, else False"""
        if self.verbose>=3: print 'MemoryBase.has_key'
        return key in self._keys

    def keys(self):
        """Return list of supported VM keys.
        D.keys() -> list of D's keys"""
        if self.verbose>=3: print 'MemoryBase.keys'
        return list(self._keys)

    def values(self):
        """Return list of recorded VM statistics.
        D.values() -> list of D's values"""
        if self.verbose>=3: print 'MemoryBase.values'
        return list(self._values)

    def get(self, key, default=None):
        """Return floating point number associated with a VM key.
        D.get(k[,d]) -> D[k] if k in D, else d.  d defaults to None."""
        if self.verbose>=1: print 'MemoryBase.get'
        v = self[key]

        if type(default) in [int,float]:
            default = np.float_(default)
        if default is not None and not isinstance(default, np.floating):
            raise ValueError('Default value must be a floating point number.')

        if default is not None and np.isnan(v):
            return default
        else:
            return v

    def setdefault(self, key, default=None):
        """Return floating point number associated with a VM key.
        D.setdefault(k[,d]) -> D.get(k,d), also set D[k]=d if k not in D"""
        if self.verbose>=1: print 'MemoryBase.setdefault'
        v = self[key]

        if type(default) in [int,float]:
            default = np.float_(default)
        if default is not None and not isinstance(default, np.floating):
            raise ValueError('Default value must be a floating point number.')

        if default is not None and np.isnan(v):
            self[key] = default
            return default
        else:
            return v

    def pop(self, key, default=None):
        """Return floating point number for a VM key and mark it as outdated.
        D.pop(k[,d]) -> v, remove specified key and return the corresponding value
        If key is not found, d is returned if given, otherwise KeyError is raised"""

        if self.verbose>=1: print 'MemoryBase.pop'
        v = self[key]

        if type(default) in [int,float]:
            default = np.float_(default)
        if default is not None and not isinstance(default, np.floating):
            raise ValueError('Default value must be a floating point number.')

        if default is not None and np.isnan(v):
            return default
        else:
            del self[key]
            return v

    def popitem(self):
        """Return floating point number for some not-yet outdated VM key.
        D.popitem() -> (k, v), remove and return some (key, value) pair as a
        2-tuple; but raise KeyError if D is empty"""
        if self.verbose>=1: print 'MemoryBase.popitem'

        for k,v in self.iteritems():
            if not np.isnan(v):
                del self[k]
                return (k,v)
        raise KeyError

    def __add__(self, other):
        """x.__add__(y) <==> x+y"""
        if self.verbose>=1: print 'MemoryBase.__add__(%s,%s)' \
            % (object.__repr__(self), object.__repr__(other))
        res = self.copy()
        if isinstance(other, MemoryBase):
            res._values.__iadd__(other._values)
        elif type(other) in [int,float]:
            res._values.__iadd__(other)
        else:
            raise TypeError('Unsupported operand type')
        return res

    def __sub__(self, other):
        """x.__sub__(y) <==> x-y"""
        if self.verbose>=1: print 'MemoryBase.__sub__(%s,%s)' \
            % (object.__repr__(self), object.__repr__(other))
        res = self.copy()
        if isinstance(other, MemoryBase):
            res._values.__isub__(other._values)
        elif type(other) in [int,float]:
            res._values.__isub__(other)
        else:
            raise TypeError('Unsupported operand type')
        return res

    def __radd__(self, other):
        """x.__radd__(y) <==> y+x"""
        if self.verbose>=1: print 'MemoryBase.__radd__(%s,%s)' \
            % (object.__repr__(self), object.__repr__(other))
        res = self.copy()
        if isinstance(other, MemoryBase):
            res._values.__iadd__(other._values)
        elif type(other) in [int,float]:
            res._values.__iadd__(other)
        else:
            raise TypeError('Unsupported operand type')
        return res

    def __rsub__(self, other):
        """x.__rsub__(y) <==> y-x"""
        if self.verbose>=1: print 'MemoryBase.__rsub__(%s,%s)' \
            % (object.__repr__(self), object.__repr__(other))
        res = self.copy()
        res._values.__imul__(-1.0)
        if isinstance(other, MemoryBase):
            res._values.__iadd__(other._values)
        elif type(other) in [int,float]:
            res._values.__iadd__(other)
        else:
            raise TypeError('Unsupported operand type')
        return res

# -------------------------------------------------------------------

class MemoryStatistics(MemoryBase):

    def __init__(self, verbose=0):
        MemoryBase.__init__(self, verbose)
        self.update()

    def __setitem__(self, key, value):
        """Set VM key to a floating point number.
        x.__setitem__(i, y) <==> x[i]=y"""
        if self.verbose>=2: print 'MemoryStatistics.__setitem__'
        if key not in self:
            raise KeyError(key)
        if type(value) in [int,float]:
            value = np.float_(value)
        if not isinstance(value, np.floating):
            raise ValueError('Value must be a floating point number.')
        i = self.keys().index(key)
        self._values[i] = value

    def __delitem__(self, key):
        """Mark a VK key as outdated.
        x.__delitem__(y) <==> del x[y]"""
        if self.verbose>=2: print 'MemoryStatistics.__delitem__'
        if key not in self:
            raise KeyError(key)
        self[key] = np.nan

    def clear(self):
        """Mark all supported VM keys as outdated.
        D.clear() -> None.  Remove all items from D."""
        if self.verbose>=1: print 'MemoryStatistics.clear'
        self._values[:] = np.nan

    def refresh(self):
        """Refresh all outdated VM keys by reading /proc/<pid>/status."""
        if self.verbose>=1: print 'MemoryBase.refresh'

        # NB: Linux /proc is for humans; Solaris /proc is for programs!
        # TODO: Use pipe from 'prstat -p <pid>' or 'pmap -x <pid> 1 1'

        # Skip refresh if none are outdated (i.e. nan)
        if not np.isnan(self._values).any():
            if self.verbose>=2: print 'refresh: skipping...'
            return

        try:
            f = open('/proc/%d/status' % os.getpid(), 'r')
            for line in f:
                k, v = line.decode('ascii').split(':')

                # Only refresh supported keys that are outdated (i.e. nan)
                if k in self and np.isnan(self[k]):
                    t, s = v.strip().split(None, 1)
                    if self.verbose >= 2:
                        print 'refresh: k=%s, t=%s, s=%s' % (k, t, s)
                    self[k] = float(t) * self._scale[s.upper()]

            f.close()
        except (IOError, UnicodeError, ValueError):
            # Reset on error
            self.clear()

    def update(self, other=None):
        """Update VM statistics from a supplied dict, else clear and refresh.
        D.update(E) -> None.  Update D from E: for k in E.keys(): D[k] = E[k]"""
        if self.verbose>=1: print 'MemoryStatistics.update'

        # Call to update without arguments has special meaning
        if other is None:
            self.clear()
            self.refresh()
        else:
            MemoryBase.update(self, other)

    def __iadd__(self, other):
        """x.__iadd__(y) <==> x+=y"""
        if self.verbose>=1: print 'MemoryStatistics.__iadd__(%s,%s)' \
            % (object.__repr__(self), object.__repr__(other))
        if isinstance(other, MemoryBase):
            self._values.__iadd__(other._values)
        elif type(other) in [int,float]:
            self._values.__iadd__(other)
        else:
            raise TypeError('Unsupported operand type')
        return self

    def __isub__(self, other):
        """x.__isub__(y) <==> x-=y"""
        if self.verbose>=1: print 'MemoryStatistics.__isub__(%s,%s)' \
            % (object.__repr__(self), object.__repr__(other))
        if isinstance(other, MemoryBase):
            self._values.__isub__(other._values)
        elif type(other) in [int,float]:
            self._values.__isub__(other)
        else:
            raise TypeError('Unsupported operand type')
        return self

# -------------------------------------------------------------------

#http://www.eecho.info/Echo/python/singleton/
#http://mail.python.org/pipermail/python-list/2007-July/622333.html

class Singleton(object):
    """A Pythonic Singleton object."""
    def __new__(cls, *args, **kwargs):
        if '_inst' not in vars(cls):
            cls._inst = object.__new__(cls, *args, **kwargs)
            #cls._inst = super(type, cls).__new__(cls, *args, **kwargs)
        return cls._inst

class MemorySingleton(MemoryBase, Singleton):
    __doc__ = MemoryBase.__doc__ + """
    The singleton variant is immutable once it has been instantiated, which
    makes it suitable for recording the initial overhead of starting Python."""

    def __init__(self, verbose=0):
        if verbose>=1: print 'MemorySingleton.__init__'
        if '_values' not in vars(self):
            if verbose>=1: print 'MemorySingleton.__init__ FIRST!'
            # Hack to circumvent singleton immutability
            self.__class__ = MemoryStatistics
            self.__init__(verbose)
            self.__class__ = MemorySingleton

    def __setitem__(self, key, value):
        """Disabled for the singleton.
        x.__setitem__(i, y) <==> x[i]=y"""
        if self.verbose>=2: print 'MemorySingleton.__setitem__'
        raise ReferenceError('Singleton is immutable.')

    def __delitem__(self, key):
        """Disabled for the singleton.
        x.__delitem__(y) <==> del x[y]"""
        if self.verbose>=2: print 'MemorySingleton.__delitem__'
        raise ReferenceError('Singleton is immutable.')

    def clear(self):
        """Disabled for the singleton.
        D.clear() -> None.  Remove all items from D."""
        if self.verbose>=1: print 'MemorySingleton.clear'
        raise ReferenceError('Singleton is immutable.')

    def update(self):
        """Disabled for the singleton.
        D.update(E) -> None.  Update D from E: for k in E.keys(): D[k] = E[k]"""
        if self.verbose>=1: print 'MemorySingleton.update'
        raise ReferenceError('Singleton is immutable.')

    def copy(self):
        """Return a shallow non-singleton copy of a VM statistics instance.
        D.copy() -> a shallow copy of D"""
        if self.verbose>=1: print 'MemorySingleton.copy'
        # Hack to circumvent singleton self-copy
        self.__class__ = MemoryStatistics
        res = self.copy()
        self.__class__ = MemorySingleton
        return res

# Make sure singleton is instantiated
MemorySingleton()

# -------------------------------------------------------------------

# Helper functions for leak testing with NumPy arrays

def shapegen(size, ndims, ecc=0.5):
    """Return a generator of an N-dimensional array shape
    which approximately contains a given number of elements.

        size:       int or long in [1,inf[
                    The total number of elements
        ndims=3:    int in [1,inf[
                    The number of dimensions
        ecc=0.5:    float in ]0,1[
                    The eccentricity of the distribution
    """
    assert type(size) in [int,float] and size>=1
    assert type(ndims) is int and ndims>=1
    assert type(ecc) in [int,float] and ecc>0 and ecc<1

    for i in range(ndims-1):
        scale = size**(1.0/(ndims-i))
        c = round(np.random.uniform((1-ecc)*scale, 1.0/(1-ecc)*scale))
        size/=c
        yield c
    yield round(size)

def shapeopt(maxseed, size, ndims, ecc=0.5):
    """Return optimal estimate of an N-dimensional array shape
    which is closest to containing a given number of elements.

        maxseed:    int in [1,inf[
                    The maximal number of seeds to try
        size:       int or long in [1,inf[
                    The total number of elements
        ndims=3:    int in [1,inf[
                    The number of dimensions
        ecc=0.5:    float in ]0,1[
                    The eccentricity of the distribution
    """
    assert type(maxseed) is int and maxseed>=1
    assert type(size) in [int,float] and size>=1
    assert type(ndims) is int and ndims>=1
    assert type(ecc) in [int,float] and ecc>0 and ecc<1

    digits_best = np.inf
    shape_best = None
    for seed in range(maxseed):
        np.random.seed(seed)
        shape = tuple(shapegen(size, ndims, ecc))
        if np.prod(shape) == size:
            return -np.inf, shape
        digits = np.log10(abs(np.prod(shape)-size))
        if digits < digits_best:
            (digits_best, shape_best) = (digits, shape)
    return digits_best, shape_best

