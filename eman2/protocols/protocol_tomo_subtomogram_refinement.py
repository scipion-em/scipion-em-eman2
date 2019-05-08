# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *              Ignacio del Cano  (idelcano@eyeseetea.com) [1]
# *
# * [1] EyeSeeTea Ltd, London, UK
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import basename, abspath

import pyworkflow.protocol.params as params
import pyworkflow.em as pwem
from pyworkflow.em import ImageHandler
from pyworkflow.em.data import Transform
from pyworkflow.em.convert import Ccp4Header
from pyworkflow.utils.path import createAbsLink
from pyworkflow.utils.properties import Message

from tomo.protocols import ProtTomoBase
from pyworkflow.protocol import STEPS_PARALLEL


class EmanProtTomoRefinement(pwem.EMProtocol, ProtTomoBase):
    """Protocol to performs a conventional iterative subtomogram averaging using the full set of particles."""
    _outputClassName = 'SubTomogramRefinement'
    _label = 'subtomogram refinement'

    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Params')
        group = form.addGroup('Input')
        group.addParam('inputSetOfTomograms', pwem.PointerParam,
                       pointerClass='SetOfTomograms',
                       label="Input set of tomograms", important=True,
                       help='Select the input set of tomograms.')
        group.addParam('iter', pwem.IntParam, default=1,
                       label='Number of iterations',
                       help='default(1)'
                            'The number of iterations to perform.')
        group.addParam('shrink', pwem.IntParam, default=1,
                       label='Shrink:',
                       help='Default=1 (no shrinking).'
                            'Optionally shrink the input'
                            'volumes by an integer amount for'
                            ' coarse alignment.')
        group.addParam('goldStandardOff', pwem.BooleanParam, default=False,
                       label="Gold Standard off",
                       help='This will PREVENT splitting the dataset'
                            'provided through --input into two groups, '
                            'and the entire dataset will be refined together.'
                            'If this parameter is NOT supplied (and thus '
                            'the refinement is "gold standard") and --ref is supplied, '
                            'two copies of the reference will be generated and'
                            ' randomphase-lowpass filtered to the resolution '
                            'specified through --refrandphase.')

        form.addParallelSection(threads=2, mpi=4)

    #--------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        pass

    #--------------- STEPS functions -----------------------

    def convertInputStep(self):
        pass

    def runMLStep(self, params):
        pass

    def createOutputStep(self):
        pass

    #--------------- INFO functions -------------------------

    def _validate(self):
        return []

    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []
