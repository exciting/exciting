#!/usr/bin/python
import os
import sys
import time
import glob
import trace
import tempfile

tmpdir = tempfile.mkdtemp(prefix='ase-')
os.chdir(tmpdir)

def build(email):
    if os.system('svn checkout ' +
                 'https://svn.fysik.dtu.dk/projects/ase/trunk ase') != 0:
        raise RuntimeError('Checkout of ASE failed!')
    os.chdir('ase')
    if os.system('python setup.py install --home=.') != 0:
        raise RuntimeError('Installation failed!')
    sys.path.insert(0, 'lib/python')
    from ase.test import test
    from ase.version import version

    # Run test-suite:
    stream = open('test-results.txt', 'w')
    results = test(verbosity=2, dir='ase/test', display=False, stream=stream)
    stream.close()
    if len(results.failures) > 0 or len(results.errors) > 0:
        address = email
        subject = 'ASE test-suite failed!'
        os.system('mail -s "%s" %s < %s' %
                  (subject, address, 'test-results.txt'))
        raise RuntimeError('Testsuite failed!')

    # Generate tar-file:
    assert os.system('python setup.py sdist') == 0

    if os.system('epydoc --docformat restructuredtext --parse-only ' +
                 '--name ASE ' +
                 '--url http://wiki.fysik.dtu.dk/ase ' +
                 '--show-imports --no-frames -v ase &> epydoc.out') != 0:
        raise RuntimeError('Epydoc failed!')

    epydoc_errors = open('epydoc.out').read()
    if ' Warning:' in epydoc_errors:
        sys.stderr.write(epydoc_errors)

    os.chdir('doc')
    os.mkdir('_build')
    if os.system('PYTHONPATH=%s/ase sphinx-build . _build' % tmpdir) != 0:
        raise RuntimeError('Sphinx failed!')
    os.system('cd _build; cp _static/searchtools.js .; ' +
              'sed -i s/snapshot.tar/%s.tar/ download.html' % version)

    if 1:
        if os.system('PYTHONPATH=%s/ase ' % tmpdir +
                     'sphinx-build -b latex . _build 2> error') != 0:
            raise RuntimeError('Sphinx failed!')
        os.system(
            'grep -v "WARNING: unusable reference target found" error 1>&2')
        
        os.chdir('_build')
        #os.system('cd ../..; ln -s doc/_static')
        if os.system('make ase-manual.pdf 2>&1') != 0:
            raise RuntimeError('pdflatex failed!')
    else:
        os.chdir('_build')

    assert os.system('mv ../../html epydoc;' +
                     'mv ../../dist/python-ase-%s.tar.gz .' % version) == 0
    
tarfiledir = None
if len(sys.argv) == 3:
    tarfiledir = sys.argv[2]
    try:
        os.remove(tarfiledir + '/ase-webpages.tar.gz')
    except OSError:
        pass

build(sys.argv[1])
    
if tarfiledir is not None:
    os.system('cd ..; tar czf %s/ase-webpages.tar.gz _build' % tarfiledir)
    os.system('cd; rm -r ' + tmpdir)
