# coding=utf-8
# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
# *              Ignacio del Cano  (idelcano@eyeseetea.com) [1]
# *              Arnau Sanchez  (arnau@eyeseetea.com) [1]
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
from types import NoneType
from glob import glob
import os
import re

from pyworkflow import utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.em as pwem
from pyworkflow.protocol import STEPS_PARALLEL

from eman2.convert import writeSetOfSubTomograms, getLastParticlesParams, updateSetOfSubTomograms
import eman2

from tomo.protocols import ProtTomoBase
from tomo.objects import AverageSubTomogram, SetOfSubTomograms


SAME_AS_PICKING = 0


class EmanProtMultiReferenceRefinement(pwem.EMProtocol, ProtTomoBase):
    """
    This protocol wraps *e2spt_classify.py* EMAN2 program.

    Protocol to performs a conventional iterative subtomogram averaging
    using the full set of particles.
    It will take a set of subtomograms (particles) and a set of Averages(reference)
    and 3D reconstruct a subtomogram.
    It also builds a set of subtomograms that contains the original particles
    plus the score, coverage and align matrix per subtomogram .
    """

    _outputClassName = 'MultiReferenceRefinement'
    _label = 'multi-reference refinement'
    OUTPUT_PREFIX = 'outputSetOf3DClassesSubTomograms'
    OUTPUT_DIR = "spt_00"

    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSetOfSubTomogram', params.PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True, label='Input SubTomograms',
                      help='Select the set of subtomograms to perform the reconstruction.')
        form.addParam('inputRef', params.PointerParam,
                      # SetOfSubtomograms for now, actually SetOfAverages
                      pointerClass='SetOfSubTomograms', allowsNull=True,
                      default=None, label='Input Ref Set of Averages',
                      help='3D reference for initial model generation.'
                           'No reference is used by default.')

        form.addSection(label='Optimization')
        form.addParam('mask', params.PointerParam,
                      label='Mask',
                      allowsNull=True,
                      pointerClass='VolumeMask',
                      help='Select a 3D Mask to be applied to the initial model')

        form.addParam('niter', pwem.IntParam, default=5,
                       label='Number of iterations',
                       help='The number of iterations to perform.')

        form.addParam('sym', params.StringParam, default='c1',
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Symmetry',
                      help='Symmetry (Default: c1')

        form.addParam('tarres', pwem.FloatParam, default=20.0,
                      label='Tarres:',
                      help='Target Resolution. Default = 20.0')

        form.addParam('mass', pwem.FloatParam, default=500.0,
                       label='Mass:',
                       help='Default=500.0.'
                            'mass')

        form.addParam('localfilter', params.BooleanParam, default=False,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Local filter',
                      help='use tophat local')

        form.addParam('threads', pwem.IntParam, default=2,
                       label='Threads:',
                       help='Number of threads')

        form.addParam('maskFile', params.PointerParam, allowsNull=True,
                       expertLevel=params.LEVEL_ADVANCED,
                       pointerClass='VolumeMask', label='Mask file',
                       help='Select the mask object')


        form.addParallelSection(threads=2, mpi=4)

    #--------------- INSERT steps functions ----------------


    def _insertAllSteps(self):
        #TODO: Get the basename.hdf from the inputSetOfSubTomogram
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('multiLevelRefinement')
        #TODO: Set and show the output
        #self._insertFunctionStep('createOutputStep')

    #--------------- STEPS functions -----------------------
    def convertInputStep(self):
        storePath = self._getExtraPath("subtomograms")
        pwutils.makePath(storePath)

        self.newFn = pwutils.removeBaseExt(list(self.inputSetOfSubTomogram.get().getFiles())[0]).split('__ctf')[0] + '.hdf'
        self.newFn = pwutils.join(storePath, self.newFn)
        writeSetOfSubTomograms(self.inputSetOfSubTomogram.get(), storePath)

    def multiLevelRefinement(self):
        """ Run the Multi-Level refinement. """
        args = ' %s' % self.newFn
        if not isinstance(self.inputRef.get(), NoneType):
            args += (' --refs=%s ' % self.newFn)
        args += ' --mass=%f' % self.mass
        args += ' --tarres=%f' % self.tarres
        args += ' --mask=%s' % self.mask
        args += ' --threads=%d' % self.threads
        args += ' --sym=%s ' % self.sym
        args += ' --path=%s ' % self.getOutputPath()
        if self.niter > 1:
            args += ' --niter=%d' % self.niter
        if self.localfilter:
            args += ' --localfilter '

        program = eman2.Plugin.getProgram('e2spt_classify.py')
        self._log.info('Launching: ' + program + ' ' + args)
        self.runJob(program, args)

    def getLastFromOutputPath(self, pattern):
        threedPaths = glob(self.getOutputPath("*"))
        imagePaths = sorted(path for path in threedPaths if re.match(pattern, os.path.basename(path)))
        if not imagePaths:
            raise Exception("No file in output directory matches pattern: %s" % pattern)
        else:
            return imagePaths[-1]

    def createOutputStep(self):
        '''
        lastImage = self.getLastFromOutputPath("threed_\d+.hdf")
        samplingRate = self.inputRef.get().getSamplingRate()
        inputSetOfSubTomograms = self.inputSetOfSubTomogram.get()

        # Output 1: Subtomogram
        averageSubTomogram = AverageSubTomogram()
        averageSubTomogram.setFileName(lastImage)
        averageSubTomogram.copyInfo(inputSetOfSubTomograms)
        averageSubTomogram.setSamplingRate(samplingRate)

        # Output 2: setOfSubTomograms
        particleParams = getLastParticlesParams(self.getOutputPath())
        outputSetOfSubTomograms = self._createSet(SetOfSubTomograms, 'subtomograms%s.sqlite', "")
        outputSetOfSubTomograms.copyInfo(inputSetOfSubTomograms)
        outputSetOfSubTomograms.setCoordinates3D(inputSetOfSubTomograms.getCoordinates3D())
        outputSetOfSubTomograms.setSamplingRate(samplingRate) # diff
        updateSetOfSubTomograms(inputSetOfSubTomograms, outputSetOfSubTomograms, particleParams)

        self._defineOutputs(outputSubTomogram=averageSubTomogram, outputSetOfSubTomograms=outputSetOfSubTomograms)
        self._defineSourceRelation(self.inputSetOfSubTomogram, averageSubTomogram)
        self._defineSourceRelation(self.inputSetOfSubTomogram, outputSetOfSubTomograms)
        '''
    def getOutputPath(self, *args):
        return self._getExtraPath(self.OUTPUT_DIR, *args)

    #--------------- INFO functions -------------------------

    @classmethod
    def isDisabled(cls):
        return not eman2.Plugin.isTomoAvailableVersion()

    def _summary(self):
        summary = []
        summary.append("Set Of SubTomograms source: %s" % (self.inputSetOfSubTomogram.get().getFileName()))

        if not isinstance(self.inputRef.get(), NoneType):
            summary.append("Referenced Tomograms source: %s" % (self.inputRef.get().getFileName()))

        if self.getOutputsSize() >= 1:
            summary.append("Subtomogram refinement Completed")
        else:
            summary.append("Subtomogram refinement not ready yet.")

        return summary

    def _methods(self):
        inputSetOfSubTomgrams = self.inputSetOfSubTomogram.get()
        return [
            "Applied refinement using e2spt_classify (stochastic gradient descent)",
            "A total of %d particles of dimensions %s were used"
            % (inputSetOfSubTomgrams.getSize(), inputSetOfSubTomgrams.getDimensions()),
        ]