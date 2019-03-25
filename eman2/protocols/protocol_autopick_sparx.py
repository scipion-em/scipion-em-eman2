# **************************************************************************
# *
# * Authors:     Airen Zaldivar (azaldivar@cnb.csic.es) [1]
# *              J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] Science for Life Laboratory, Stockholm University
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

import os

from pyworkflow.protocol.params import (IntParam, FloatParam,
                                        BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtParticlePickingAuto

import eman2
from eman2.convert import readSetOfCoordinates


class SparxGaussianProtPicking(ProtParticlePickingAuto):
    """
    Automated particle picker for SPA. Uses Sparx gaussian picker.
    For more information see http://sparx-em.org/sparxwiki/e2boxer
    """
    _label = 'sparx gaussian picker'

    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)
        self.extraParams = 'pixel_input=1:pixel_output=1'

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)
        form.addParam('boxSize', IntParam, default=100, allowsPointers=True,
                      label='Box Size', help='Box size in pixels')
        line = form.addLine('Picker range',
                            help='CCF threshold range for automatic picking')
        line.addParam('lowerThreshold', FloatParam, default='1',
                      label='Lower')
        line.addParam('higherThreshold', FloatParam, default='30',
                      label='Higher')
        form.addParam('gaussWidth', FloatParam, default='1',
                      label='Gauss Width',
                      help='Width of the Gaussian kernel used')
        form.addParam('useVarImg', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Use variance image?')
        form.addParam('doInvert', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Invert contrast?',
                      help='Picker expects particles to be white.')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertInitialSteps(self):
        initId = self._insertFunctionStep('initSparxDb',
                                          self.lowerThreshold.get(),
                                          self.higherThreshold.get(),
                                          self.boxSize.get(), self.gaussWidth.get(),
                                          self.useVarImg, self.doInvert)
        return [initId]

    # --------------------------- STEPS functions -----------------------------
    def initSparxDb(self, lowerThreshold, higherThreshold, boxSize,
                    gaussWidth, useVarImg, doInvert):
        args = {"lowerThreshold": lowerThreshold,
                "higherThreshold": higherThreshold,
                "boxSize": boxSize,
                "gaussWidth": gaussWidth,
                "useVarImg": "true" if useVarImg else "false",
                "doInvert": "true" if doInvert else "false",
                "extraParams": self.extraParams}
        params = 'demoparms --makedb=thr_low=%(lowerThreshold)s:'
        params += 'thr_hi=%(higherThreshold)s:boxsize=%(boxSize)s:'
        params += 'invert_contrast=%(doInvert)s:use_variance=%(useVarImg)s:'
        params += 'gauss_width=%(gaussWidth)s:%(extraParams)s'

        self.runJob(eman2.Plugin.getProgram('sxprocess.py'),
                    params % args, cwd=self.getCoordsDir(),
                    numberOfThreads=1)

    def _pickMicrograph(self, mic, *args):
        micFile = os.path.relpath(mic.getFileName(), self.getCoordsDir())
        params = ('--gauss_autoboxer=demoparms --write_dbbox --boxsize=%d %s'
                  % (self.boxSize, micFile))
        program = eman2.Plugin.getBoxerCommand(boxerVersion='old')

        self.runJob(program, params, cwd=self.getCoordsDir())

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------
    def getCoordsDir(self):
        return self._getExtraPath()

    def getFiles(self):
        return (self.inputMicrographs.get().getFiles() |
                ProtParticlePickingAuto.getFiles(self))

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        coordSet.setBoxSize(self.boxSize.get())
        readSetOfCoordinates(workingDir, micList, coordSet, newBoxer=False)
