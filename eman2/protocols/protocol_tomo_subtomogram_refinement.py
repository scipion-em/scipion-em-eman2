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

import pyworkflow.protocol.params as params
import pyworkflow.em as pwem

from tomo.protocols import ProtTomoBase
from pyworkflow.protocol import STEPS_PARALLEL, join

import eman2
from eman2.constants import *

from tomo.objects import SetOfSubTomograms, SubTomogram
# Micrograph type constants for particle extraction
from pyworkflow.utils.dataset import DataSet

SAME_AS_PICKING = 0
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
        form.addParam('inputSetOfSubTomogram', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True, label='Input SetOfSubTomograms',
                      help='Select the SetOfSubTomograms.')
        form.addParam('inputRef', params.PointerParam,
                      pointerClass='Volume', allowsNull=True,
                      default=None, label='Input Ref Tomogram',
                      help='3D reference for initial model generation.'
                           'No reference is used by default.')

        group = form.addGroup('Input')
        group.addParam('niter', pwem.IntParam, default=5,
                       label='Number of iterations',
                       help='The number of iterations to perform.')
        group.addParam('mass', pwem.FloatParam, default=500.0,
                       label='Mass:',
                       help='Default=500.0.'
                            'mass')
        group.addParam('threads', pwem.IntParam, default=2,
                       label='Threads:',
                       help='Number of threads')
        group.addParam('pkeep', pwem.FloatParam, default=0.8,
                       label='Particle keep:',
                       help='Fraction of particles to keep')
        group.addParam('goldstandard', pwem.IntParam, default=-1,
                       label='GoldStandard:',
                       help='initial resolution for gold standard refinement')
        form.addParallelSection(threads=2, mpi=4)

    #--------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        self._insertFunctionStep('refinementSubtomogram')


    def runImportParticlesSqlite(cls, pattern, samplingRate):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(pwem.ProtImportParticles,
                                         importFrom=4,
                                         sqliteFile=pattern, samplingRate=samplingRate)
        cls.launchProtocol(cls.protImport)
        # check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. outputParticles is None.' % pattern)
        return cls.protImport

    def getFile(self, key):
        if key in self.filesDict:
            return join(self.path, self.filesDict[key])
        return join(self.path, key)

    #--------------- STEPS functions -----------------------
    def refinementSubtomogram(self):
        print ('---------------------------------------------------------------------------------')
        print ('refinementSubtomogram')
        print ('---------------------------------------------------------------------------------')


        print("import particles")

        self.partsFn = self.getFile(self.inputSetOfSubTomogram.get().getFileName())
        self.protImport = self.runImportParticlesSqlite(self.partsFn, 3.5)
        args = ' %s' % (
            pwem.os.getcwd() + "/" + self.inputSetOfSubTomogram.get().getFileName()).replace("subtomograms.sqlite", "extra/sptboxer_01/basename.hdf")
        print(self.inputSetOfSubTomogram.get().get())
        print ('---------------------------------------------------------------------------------')
        print(self.inputRef.get())
        if self.inputRef is not None:
            args += (' --reference=%s ' % self.inputRef.get().getFileName())

        if self.mass:
            args += (' --mass=%f' % self.mass)

        if self.niter > 1:
            args += ' --niter=%d' % self.niter

        args += ' --threads=%d' % self.threads
        args += ' --goldstandard=%d ' % self.goldstandard
        args += ' --pkeep=%f ' % self.pkeep

        args += ' --sym=c1 --maxtilt=90.0 '
        print ('---------------------------------------------------------------------------------')
        print("command: " + args)
        program = eman2.Plugin.getProgram('e2spt_refine.py')
        self.runJob(program, args,
                    cwd=self._getExtraPath())

    def convertInputStep(self):
        pass

    def runMLStep(self, params):
        pass

    def createOutputStep(self):
        pass

    #--------------- INFO functions -------------------------

    def _validate(self):
        errors = []

        if not eman2.Plugin.isNewVersion():
            errors.append('Your EMAN2 version does not support the subtomogram refinement. '
                          'Please update your installation to EMAN 2.23 or newer.')

        return errors

    def _citations(self):
        return []

    def _summary(self):
        summary = []
        summary.append("SetOfSubTomograms source: %s" % (self.inputSetOfSubTomogram.get().getFileName()))

        if self.getOutputsSize() >= 1:
            summary.append("Subtomogram refinamented")
        else:
            summary.append("Output subtomograms not ready yet.")

        return summary

    def _methods(self):
        return []
