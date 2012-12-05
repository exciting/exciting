import warnings

class Deprecate:
    def __init__(self, obj, name, newmodule, oldmodule='ase'):
        self.obj = obj
        self.name = name
        self.newmodule = newmodule
        self.oldmodule = oldmodule

    def __call__(self, *args, **kwargs):
        message = ('%s.%s is deprecated, use %s.%s instead' %
                   (self.oldmodule, self.name, self.newmodule, self.name))
        warnings.warn(message, DeprecationWarning, stacklevel=2)
        return self.obj(*args, **kwargs)

def _dep(method):
    def _method(self, *args):
        message = ('ase.%s is deprecated, use %s.%s instead' %
                   (self.name, self.newmodule, self.name))
        warnings.warn(message, DeprecationWarning, stacklevel=2)
        return method(self, *args)
    return _method

class DeprecatedFloat(float):
    def __new__(cls, value, name, newmodule):
        return float.__new__(cls, value)

    def __init__(self, value, name, newmodule):
        self.name = name
        self.newmodule = newmodule

    __mul__ = _dep(float.__mul__)
    __rmul__ = _dep(float.__rmul__)
    __div__ = _dep(float.__div__)
    __rdiv__ = _dep(float.__rdiv__)

class DeprecatedNumpyImport:
    def __init__(self):
        import numpy
        self.numpy = numpy

    def __getattr__(self, key):
        warnings.warn('ase.np is deprecated; use import numpy as np instead')
        return getattr(self.numpy, key)
