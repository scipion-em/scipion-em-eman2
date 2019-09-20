# **************************************************************************
# *
# * Authors:     Adrian Quintana (adrian@eyeseetea.com) [1]
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


from pyworkflow import utils as pwutils
import pyworkflow.em as pwem
import pyworkflow.protocol.params as params

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram

from xmipp3.convert import writeSetOfCoordinates

import eman2
from eman2.constants import *

import xmippLib

# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1

class EmanProtTomoExtraction(pwem.EMProtocol, ProtTomoBase):
    """ Manual picker for Tomo. Uses EMAN2 e2spt_boxer.py.
    """
    _label = 'tomo extraction'
    OUTPUT_PREFIX = 'outputSetOfSubtomogram'

    @classmethod
    def isDisabled(cls):
        return not eman2.Plugin.isTomoAvailableVersion()

    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Input Coordinates", important=True,
                      pointerClass='SetOfCoordinates3D',
                      help='Select the SetOfCoordinates3D.')

        # The name for the followig param is because historical reasons
        # now it should be named better 'micsSource' rather than
        # 'micsSource', but this could make inconsistent previous executions
        # of this protocols, we will keep the name
        form.addParam('downsampleType', params.EnumParam,
                      choices=['same as picking', 'other'],
                      default=0, important=True,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Tomogram source',
                      help='By default the subtomograms will be extracted '
                           'from the tomogram used in the picking '
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide '
                           'a different tomogram to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided '
                           'tomogram and coordinates are related ')

        form.addParam('inputTomogram', params.PointerParam,
                      pointerClass='Tomogram',
                      condition='downsampleType != %s' % SAME_AS_PICKING,
                      important=True, label='Input tomogram',
                      help='Select the tomogram from which to extract.')

        form.addParam('boxSize', params.FloatParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size')

        form.addParam('downFactor', params.FloatParam, default=1.0,
                      label='Downsampling factor',
                      help='Select a value greater than 1.0 to reduce the size '
                           'of tomogram before extracting the subtomograms. '
                           'If 1.0 is used, no downsample is applied. '
                           'Non-integer downsample factors are possible. ')

        form.addSection(label='Preprocess')
        form.addParam('doInvert', params.BooleanParam, default=False,
                      label='Invert contrast?',
                      help='Invert the contrast if your tomogram is black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')

        form.addParam('doNormalize', params.BooleanParam, default=False,
                      label='Normalize subtomogram?',
                      help='Normalization processor applied to subtomograms before extraction.')
        form.addParam('normproc', params.EnumParam,
                      choices=['normalize', 'normalize.edgemean'],
                      label='Normalize method',
                      condition='doNormalize',
                      default=PROC_NORMALIZE,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Use normalize.edgemean if the particles have a'
                        ' clear solvent background (i.e., they are not part of a '
                        'larger complex or embeded in a membrane)')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _tomosOther(self):
        """ Return True if other micrographs are used for extract. """
        return self.downsampleType == OTHER

    def getInputTomogram(self):
        """ Return the tomogram associated to the SetOfCoordinates3D or
        Other micrographs. """
        if not self._tomosOther():
            return self.inputCoordinates.get().getVolumes()
        else:
            return self.inputTomogram.get()

    def _insertAllSteps(self):
        self._insertFunctionStep('writeSetOfCoordinates3D')
        self._insertFunctionStep('extractParticles')
        self._insertFunctionStep('createOutputStep')

    def readSetOfTomograms(self, workDir, tomogramsSet):

        subtomogram = SubTomogram()

        imgh = pwem.ImageHandler()
        x, y, z, n = imgh.getDimensions(workDir)
        for index in range(1, n + 1):
            subtomogram.cleanObjId()
            subtomogram.setLocation(index, workDir)
            subtomogram.setCoordinate3D(self.coordDict[index-1])
            subtomogram.setAcquisition(self.getInputTomogram().getAcquisition())
            tomogramsSet.append(subtomogram)

    def createOutputStep(self):
        suffix = self._getOutputSuffix(SetOfSubTomograms)
        self.outputSubTomogramsSet = self._createSetOfSubTomograms(suffix)
        self.outputSubTomogramsSet.setSamplingRate(self.getInputTomogram.get().getSamplingRate()*self.cshrink)
        self.outputSubTomogramsSet.setCoordinates3D(self.inputCoordinates)

        self.readSetOfTomograms(self._getExtraPath(pwutils.join('sptboxer_01', 'basename.hdf')),
                                self.outputSubTomogramsSet)

        self._defineOutputs(outputSetOfSubtomogram=self.outputSubTomogramsSet)
        self._defineSourceRelation(self.inputCoordinates, self.outputSubTomogramsSet)

    def writeSetOfCoordinates3D(self):
        self.coordsFileName = self._getExtraPath(
            pwutils.replaceBaseExt(self.getInputTomogram().getFileName(), 'coords'))
        self.coordDict = []
        out = file(self.coordsFileName, "w")
        for coord3DSet in self.inputCoordinates.get().iterCoordinates():
            out.write("%d\t%d\t%d\n" % (coord3DSet.getX(), coord3DSet.getY(), coord3DSet.getZ()))
            self.coordDict.append(coord3DSet.clone())
        out.close()

    # --------------------------- STEPS functions -----------------------------
    def extractParticles(self):
        TsTomo = self.inputTomogram.get().getSamplingRate()
        if self.inputCoordinates.get().getSamplingRate() is not None:
            TsCoord = self.inputCoordinates.get().getSamplingRate()
        else:
            TsCoord = TsTomo
        self.cshrink = float(TsCoord/(TsTomo*self.downFactor.get()))
        self.boxSize = self.boxSize.get()
        fnTomo = self.getInputTomogram().getFileName()
        if self.downFactor.get() != 1:
            I = xmippLib.Image(fnTomo)
            x,y,z,_ = I.getDimensions()
            I = I.scale(int(x*self.downFactor.get()), int(y*self.downFactor.get()), int(z*self.downFactor.get()))
            fnTomo = self._getTmpPath("downsampled_tomo.em")
            I.write(fnTomo)
        args = '%s --coords %s --boxsize %d' % (fnTomo, pwutils.replaceBaseExt(self.getInputTomogram().getFileName(), 'coords'),
            self.boxSize)
        if self.doInvert:
            args += ' --invert'

        if self.doNormalize:
            args += ' --normproc %s' % self.getEnumText('normproc')

        if self.cshrink > 1:
            args += ' --cshrink %d' % self.cshrink
        program = eman2.Plugin.getProgram('e2spt_boxer_old.py')
        self.runJob(program, args,
                    cwd=self._getExtraPath())

    def _validate(self):
        errors = []

        if not eman2.Plugin.isTomoAvailableVersion():
            errors.append('Your EMAN2 version does not support the tomo boxer. '
                          'Please update your installation to EMAN 2.3 or newer.')

        return errors

    def _methods(self):
        methodsMsgs = []

        if self.getOutputsSize() >= 1:
            msg = ("A total of %s subtomograms of size %s were extracted"
                   % (str(self.inputCoordinates.get().getSize()), self.boxSize))

            if self._tomosOther():
                msg += (" from another set of micrographs: %s"
                        % self.getObjectTag('inputMicrographs'))

            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Set of Subtomograms not ready yet")

        if self.doInvert:
            methodsMsgs.append("Inverted contrast on images.")
        if self.cshrink > 1:
            methodsMsgs.append("Coordinates multiple by factor %d."
                               % self.cshrink)
        if self.doNormalize:
            methodsMsgs.append("Particles were normalised. Using normalization method %s") % self.getEnumText('normproc')

        return methodsMsgs

    def _summary(self):
        summary = []
        summary.append("Micrographs source: %s"
                       % self.getEnumText("downsampleType"))

        if self.getOutputsSize() >= 1:
            summary.append("Particle box size: %s" % self.boxSize)
            summary.append("Subtomogram extracted: %s" %
                           self.inputCoordinates.get().getSize())
        else:
            summary.append("Output subtomograms not ready yet.")

        return summary
