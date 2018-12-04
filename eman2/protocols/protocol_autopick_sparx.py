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

from pyworkflow.protocol.params import IntParam, FloatParam, BooleanParam, PointerParam
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
        self.extraParams = ('pixel_input=1:pixel_output=1:invert_contrast'
                            '=True:use_variance=True')

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)

        form.addParam('bxSzFromCoor', BooleanParam, dafault=False,
                      label='Use an input coordinates for box size')
        form.addParam('boxSize', IntParam, default=100,
                      conditon='bxSzFromCoor==False',
                      label='Box Size', help='Box size in pixels')
        form.addParam('coordsToBxSz', PointerParam, pointerClass='SetOfCoordinates',
                      condition='bxSzFromCoor==True',
                      label='Coordinates to extract the box size.',
                      help='Coordinates to extract the box size. '
                           'It can be an empty set.')
        line = form.addLine('Picker range',
                            help='CCF threshold range for automatic picking')
        line.addParam('lowerThreshold', FloatParam, default='1',
                      label='Lower')
        line.addParam('higherThreshold', FloatParam, default='30',
                      label='Higher')

        form.addParam('gaussWidth', FloatParam, default='1',
                      label='Gauss Width',
                      help='Width of the Gaussian kernel used')

    # --------------------------- INSERT steps functions ----------------------
    def _insertInitialSteps(self):
        initId = self._insertFunctionStep('initSparxDb',
                                          self.lowerThreshold.get(),
                                          self.higherThreshold.get(),
                                          self.getBoxSize(),
                                          self.gaussWidth.get())
        return [initId]

    # --------------------------- STEPS functions -----------------------------
    def initSparxDb(self, lowerThreshold, higherThreshold, boxSize, gaussWidth):
        args = {"lowerThreshold": lowerThreshold,
                "higherThreshold": higherThreshold,
                "boxSize": boxSize,
                "gaussWidth": gaussWidth,
                "extraParams": self.extraParams}
        params = 'demoparms --makedb=thr_low=%(lowerThreshold)s:'
        params += 'thr_hi=%(higherThreshold)s:boxsize=%(boxSize)s:'
        params += 'gauss_width=%(gaussWidth)s:%(extraParams)s'

        self.runJob(eman2.Plugin.getProgram('sxprocess.py'),
                    params % args, cwd=self.getCoordsDir())

    def _pickMicrograph(self, mic, *args):
        micFile = os.path.relpath(mic.getFileName(), self.getCoordsDir())
        params = ('--gauss_autoboxer=demoparms --write_dbbox --boxsize=%d %s'
                  % (self.getBoxSize(), micFile))
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
        coordSet.setBoxSize(self.getBoxSize())
        readSetOfCoordinates(workingDir, micList, coordSet, newBoxer=False)

    def getBoxSize(self):
        if self.bxSzFromCoor:
            return self.coordsToBxSz.get().getBoxSize()
        else:
            return self.boxSize.get()
