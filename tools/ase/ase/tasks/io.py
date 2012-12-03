import numpy as np

from ase.parallel import world

try:
    import json
except ImportError:
    json = None


if json is None:
    def dumps(obj):
        if isinstance(obj, str):
            return '"' + obj + '"'
        if isinstance(obj, (int, float)):
            return repr(obj)
        if isinstance(obj, dict):
            return '{' + ', '.join(dumps(key) + ': ' + dumps(value)
                                   for key, value in obj.items()) + '}'
        return '[' + ','.join(dumps(value) for value in obj) + ']'

    loads = eval
else:
    class NDArrayEncoder(json.JSONEncoder):
        def __init__(self):
            json.JSONEncoder.__init__(self, sort_keys=True, indent=4)

        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)
    

    dumps = NDArrayEncoder().encode
    loads = json.loads


def numpyfy(obj):
    if isinstance(obj, dict):
        return dict((key, numpyfy(value)) for key, value in obj.items())
    if isinstance(obj, list):
        try:
            obj = np.array(obj)
        except ValueError:
            obj = [numpyfy(value) for value in obj]
    return obj


def write_json(name, results):
    if world.rank == 0:
        fd = open(name, 'w')
        fd.write(dumps(results))
        fd.close()


def read_json(name):
    fd = open(name, 'r')
    results = loads(fd.read())
    fd.close()
    world.barrier()
    return numpyfy(results)
