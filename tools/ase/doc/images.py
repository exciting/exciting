# creates: ase/ag.png ase/dft/water_divide_surf.png
from urllib import urlretrieve
def setup(app):
    pass
urlretrieve('http://wiki.fysik.dtu.dk/ase-files/ag.png', 'ase/ag.png')
urlretrieve('http://wiki.fysik.dtu.dk/ase-files/water_divide_surf.png',
            'ase/dft/water_divide_surf.png')
