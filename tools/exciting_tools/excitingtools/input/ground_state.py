"""Module for class of exciting ground state.

Ideally the input keywords (class attributes) should be parsed from the schema BUT
because excitingtools will also be available as a standalone package, one would need
to have a copy of the schema XML in excitingtools, which is kept synchronised with
the <EXCITINGROOT>/xml/.
"""
from excitingtools.input.base_class import ExcitingXMLInput


class ExcitingGroundStateInput(ExcitingXMLInput):

    # Reference: http://exciting.wikidot.com/ref:groundstate
    _valid_attributes = {'CoreRelativity', 'ExplicitKineticEnergy', 'PrelimLinSteps', 'ValenceRelativity', 'autokpt',
                         'beta0', 'betadec', 'betainc', 'cfdamp', 'chgexs', 'deband', 'dipolecorrection',
                         'dipoleposition', 'dlinengyfermi', 'do', 'energyref', 'epsband', 'epschg', 'epsengy',
                         'epsforcescf', 'epsocc', 'epspot', 'fermilinengy', 'findlinentype', 'fracinr', 'frozencore',
                         'gmaxvr', 'isgkmax', 'ldapu', 'lmaxapw', 'lmaxinr', 'lmaxmat', 'lmaxvr', 'lorecommendation',
                         'lradstep', 'maxscl', 'mixer', 'mixerswitch', 'modifiedsv', 'msecStoredSteps', 'nempty',
                         'ngridk', 'niterconvcheck', 'nktot', 'nosource', 'nosym', 'nprad', 'npsden', 'nwrite',
                         'outputlevel', 'ptnucl', 'radialgridtype', 'radkpt', 'reducek', 'rgkmax', 'scfconv', 'stype',
                         'swidth', 'symmorph', 'tevecsv', 'tfibs', 'tforce', 'tpartcharges', 'useDensityMatrix',
                         'vdWcorrection', 'vkloff', 'xctype'}

    def __init__(self, **kwargs):
        """Generate an object of ExcitingXMLInput for the groundstate attributes."""
        super().__init__('groundstate', self._valid_attributes, **kwargs)
