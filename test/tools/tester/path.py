import os

class Path():
    # TODO add descriptions
    """
    Class that helps to deal with data paths.
    """
    def __init__(self, path:str):
        """
        :param path: 
        """
        self.path = path

    def __eq__(self, other:str):
        """
        Checks if two pathes are equal.
        :param other: 
        """
        return self.path == other.path

    def __add__(other:str, self):
        """
        Attaches path to other path.
        :param other: 
        """
        return Path(os.path.join(other.path, self.path))

    def __str__(self):
        return self.path
        
    def split(self):
        """
        Splits path in its individual components. 
        :return paths:  List of individual path components.
        """
        pathsstr = self.path.split('/')
        paths = []
        for p in pathsstr:
            if p:
                paths.append(Path(p))
        return paths

    def lastElement(self):
        """
        Returns last element of path.
        """
        return os.path.basename(self.path)
        #return lastEl

    def update(self, newpath:str):
        """
        Attaches new path to path.
        :param newpath: path to be attached
        """
        return Path(os.path.join(self.path,newpath))
    
    def removeLastElement(self):
        """
        Removes last element from path.
        """
        paths = self.split()
        path_out = Path('')
        for i in range(0, len(paths)-1):
            path_out += paths[i]
        return path_out