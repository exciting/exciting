from lxml import etree
import math

eps_tol = 1.e-8

root = etree.Element("report")

def mse_eps(ftest_path,fref_path):
    """ Return MSE of real and imaginary parts of dielectric function files EPSILON*.OUT
    """
    ftest = open(ftest_path,"r")
    fref = open(fref_path,"r")

    err = 0
    n = 0
    for lt, lr in zip (ftest.readlines(),fref.readlines()):
        tRE = float(lt.split()[1]) 
        rRE = float(lr.split()[1]) 

        tIM = float(lt.split()[2]) 
        rIM = float(lr.split()[2]) 
    
        err = err + ( tRE - rRE )*( tRE - rRE ) + ( tIM - rIM )*( tIM - rIM )
        n = n + 1

    ftest.close()
    fref.close()

    mse = math.sqrt(err)/(1.0*n)
    return mse

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
# assert BSE
########################################
tdir = "runBSE"
ftest_path = tdir + "/EPSILON_NAR_BSEsinglet_SCRfull_OC11.OUT"
fref_path = "reference/"+tdir+"/EPSILON_NAR_BSEsinglet_SCRfull_OC11.OUT"

mse = mse_eps(ftest_path,fref_path)

tpassed = mse < eps_tol
name_txt = "Dielectric function of LiF with BSE" 
description_txt = "Passes if mean square error of real and imaginary parts of dielectric function lower than %g"%(eps_tol)
directory_txt = "test07/"+tdir
add_test_xml(root, tpassed, name_txt, description_txt, directory_txt)

########################################
# assert TDDFT with BSE derived kernel
########################################
tdir = "runtddftBSE"
ftest_path = tdir + "/EPSILON_NAR_FXCMB1_OC11_QMT001.OUT"
fref_path = "reference/"+tdir+"/EPSILON_NAR_FXCMB1_OC11_QMT001.OUT"

mse = mse_eps(ftest_path,fref_path)

tpassed = mse < eps_tol
name_txt = "Dielectric function of LiF with TDDFT using BSE derived kernel" 
description_txt = "Passes if mean square error of real and imaginary parts of dielectric function lower than %g"%(eps_tol)
directory_txt = "test07/"+tdir
add_test_xml(root, tpassed, name_txt, description_txt, directory_txt)


frep = open("report.xml","w+")
frep.write(etree.tostring(root, pretty_print=True))
frep.close()


