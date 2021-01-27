import os

class Path():
    '''
    Class that helps to deal with data paths.
    In:
        path        str     initial path
    '''
    def __init__(self, path):
        self.path = path

    def __eq__(self, other):
        return self.path == other.path

    def __add__(other, self):
        return Path(os.path.join(other.path, self.path))

    def __str__(self):
        return self.path
        
    def split(self):
        pathsstr = self.path.split('/')
        paths = []
        for p in pathsstr:
            if p:
                paths.append(Path(p))
        return paths

    def lastElement(self):
        return os.path.basename(self.path)
        #return lastEl

    def update(self, newpath):
        return Path(os.path.join(self.path,newpath))
    
    def removeLastElement(self):
        paths = self.split()
        path_out = Path('')
        for i in range(0, len(paths)-1):
            path_out += paths[i]
        return path_out