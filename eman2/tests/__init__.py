
from test_protocols_eman import (TestEmanBase, TestEmanAutopick,
                                 TestEmanCtfAuto, TestEmanInitialModelGroel,
                                 TestEmanInitialModelMda, TestEmanReconstruct,
                                 TestEmanRefine2D, TestEmanRefine2DBispec,
                                 TestEmanRefineEasy, TestEmanTiltValidate)



from pyworkflow.tests import DataSet

DataSet(name='tomo-em', folder='tomo-em',
        files={
               'tomo1': 'overview_wbp.em',
               'tomo2': 'overview_wbp2.em',
               'eman_coordinates': 'coordinates3Deman2'})