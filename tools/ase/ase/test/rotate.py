from ase.utils import rotate, irotate

def test(xyz):
    a = rotate(xyz)
    ixyz = '%sx,%sy,%sz' % irotate(a)
    a2 = rotate(ixyz)
    print xyz
    print ixyz
    #print np.around(a-a2, 5)
    assert abs(a-a2).max() < 1e-10

test('10z')
test('155x,43y,190z')
test('55x,90y,190z')
test('180x,-90y,45z')
test('-180y')
test('40z,50x')
