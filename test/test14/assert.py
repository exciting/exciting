from lxml import etree
import math

eps_tol_omega = 1.e-5
eps_tol_dos = 1.e-3
eps_tol_pol = 1.e-5

root = etree.Element("report")

def mae_omega(ftest_path,fref_path):
    """ returns maximum absolute error of total Wannier spread per group
    """

    fref = open(fref_path,"r")
    oref = []; start = False
    for ln in fref.readlines():
        ln = ln.strip()
        if ln.startswith( ' Wannier functions', 1): start = True
        if start and ln.startswith( 'total:'):
            oref.append( float( ln[6:].split()[0]))
    #print( 'oref:', oref)
    fref.close()

    ftest = open(ftest_path,"r")
    otest = []; start = False
    for ln in ftest.readlines():
        ln = ln.strip()
        if ln.startswith( ' Wannier functions', 1): start = True
        if start and ln.startswith( 'total:'):
            otest.append( float( ln[6:].split()[0]))
    #print( 'otest:', otest)
    ftest.close()

    err = 0.0
    for o1, o2 in zip( oref, otest):
        err = max( err, abs( o1-o2))
    return err

def mse_dos(ftest_path,fref_path):
    """ returns mean square error of Wannier interpolated DOS obtained from tetrahedron integration
    """

    fref = open(fref_path,"r")
    ftest = open(ftest_path,"r")
    err = 0.0; n = 0
    for lr, lt in zip(fref.readlines(),ftest.readlines()):
        if lr.strip() == '': continue
        if lt.strip() == '': continue
        err += (float(lr.split()[1]) - float(lt.split()[1]))**2
        n += 1
    err = math.sqrt(err/n)
    fref.close()
    ftest.close()

    return err

def mse_pol(ftest_path,fref_path):
    """ returns mean square error of macroscopic polarization
    """

    fref = open(fref_path,"r")
    ftest = open(ftest_path,"r")
    err = 0.0; n = 0
    for lr, lt in zip(fref.readlines(),ftest.readlines()):
        if lr.strip().startswith( '#'): continue
        if lt.strip().startswith( '#'): continue
        for pr, pt in zip( lr.split(), lt.split()):
            err += (float(pr) - float(pt))**2
            n += 1
    err = math.sqrt(err/n)
    fref.close()
    ftest.close()

    return err

def add_test_xml(root, tpassed, name_txt, description_txt, directory_txt):
    test = etree.SubElement(root, "test")

    status = etree.SubElement(test, "status")
    if tpassed:
        status.text = "passed"
    else:
        status.text = "failed"

    name = etree.SubElement(test, "name")
    name.text = name_txt

    description = etree.SubElement(test, "description")
    description.text = description_txt

    directory = etree.SubElement(test, "directory")
    directory.text = directory_txt



########################################
# assert Wannier functions
########################################
tdir = "runwannier/"
ftest_path = tdir+"WANNIER_INFO.OUT"
fref_path = "reference/WANNIER_INFO.REF"

mae = mae_omega(ftest_path,fref_path)
#print(mae)
tpassed = mae < eps_tol_omega
name_txt = "Wannier functions in SiC"
description_txt = "Passes if total Wannier spread per group differs less than %f from the reference. Difference: %f"%(eps_tol_omega,mae)
directory_txt = "test14/"+tdir
add_test_xml(root, tpassed, name_txt, description_txt, directory_txt)

########################################
# assert Wannier interpolation
########################################
tdir = "runwannier/"
ftest_path = tdir+"TDOS_WANNIER.OUT"
fref_path = "reference/TDOS_WANNIER.REF"

mse = mse_dos(ftest_path,fref_path)
#print(mse)
tpassed = mse < eps_tol_dos
name_txt = "Wannier interpolated DOS in SiC"
description_txt = "Passes if mean square error of Wannier interpolated DOS obtained via tetrahedron integration is smaller than %f. MSE: %f"%(eps_tol_dos,mse)
directory_txt = "test14/"+tdir
add_test_xml(root, tpassed, name_txt, description_txt, directory_txt)

frep = open("report.xml","wb+")
frep.write(etree.tostring(root, pretty_print=True))
frep.close()

########################################
# assert macroscopic polarization
########################################
tdir = "runwannier/"
ftest_path = tdir+"POLARIZATION.OUT"
fref_path = "reference/POLARIZATION.REF"

mse = mse_pol(ftest_path,fref_path)
#print(mse)
tpassed = mse < eps_tol_pol
name_txt = "Macroscopic polarization in SiC"
description_txt = "Passes if mean square error of polarization is smaller than %f. MSE: %f"%(eps_tol_pol,mse)
directory_txt = "test14/"+tdir
add_test_xml(root, tpassed, name_txt, description_txt, directory_txt)

frep = open("report.xml","wb+")
frep.write(etree.tostring(root, pretty_print=True))
frep.close()
