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

import os

from pwem.emlib.image import ImageHandler
from pyworkflow import utils as pwutils
from pwem.protocols import EMProtocol
import pyworkflow.protocol.params as params
from pyworkflow.utils.path import moveFile, cleanPath
from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram, TomoAcquisition
import eman2
from eman2.constants import *

# Tomogram type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class EmanProtTomoExtraction(EMProtocol, ProtTomoBase):
    """ Extraction for Tomo. Uses EMAN2 e2spt_boxer_old.py."""
    _label = 'tomo extraction'
    OUTPUT_PREFIX = 'outputSetOfSubtomogram'

    @classmethod
    def isDisabled(cls):
        return not eman2.Plugin.isTomoAvailableVersion()

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Input Coordinates", important=True,
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('tomoSource', params.EnumParam,
                      choices=['same as picking', 'other'],
                      default=0,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Tomogram source',
                      help='By default the subtomograms will be extracted '
                           'from the tomogram used in the picking '
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide '
                           'a different tomogram to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided '
                           'tomogram and coordinates are related ')

        form.addParam('inputTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      condition='tomoSource != %s' % SAME_AS_PICKING,
                      label='Input tomogram',
                      help='Select the tomogram from which to extract.')

        form.addParam('boxSize', params.FloatParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size. '
                           'The wizard selects same box size as picking')

        form.addParam('downFactor', params.FloatParam, default=1.0,
                      label='Downsampling factor',
                      help='If 1.0 is used, no downsample is applied. '
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
                      help='Use normalize.edgemean if the particles have a clear solvent background '
                           '(i.e., they are not part of a larger complex or embeded in a membrane)')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('writeSetOfCoordinates3D')
        self._insertFunctionStep('extractParticles')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------

    def writeSetOfCoordinates3D(self):
        self.lines = []
        self.tomoFiles = []
        inputSet = self.getInputTomograms()
        for item in inputSet:
            coordDict = []
            tomo = item.clone()
            self.coordsFileName = self._getExtraPath(
                pwutils.replaceBaseExt(tomo.getFileName(), 'coords'))

            out = open(self.coordsFileName, "w")
            for coord3DSet in self.inputCoordinates.get().iterCoordinates():
                if os.path.basename(tomo.getFileName()) == os.path.basename(coord3DSet.getVolName()):
                    out.write("%d\t%d\t%d\n" % (coord3DSet.getX(), coord3DSet.getY(), coord3DSet.getZ()))
                    coordDict.append(coord3DSet.clone())
            if coordDict:
                self.lines.append(coordDict)
                self.tomoFiles.append(tomo.getFileName())
                self.samplingRateTomo = tomo.getSamplingRate()
            out.close()

    def extractParticles(self):
        samplingRateCoord = self.inputCoordinates.get().getSamplingRate()
        samplingRateTomo = self.getInputTomograms().getFirstItem().getSamplingRate()
        for tomo in self.tomoFiles:
            args = '%s ' % tomo
            args += "--coords % s --boxsize % d" % (pwutils.replaceBaseExt(tomo, 'coords'), self.boxSize.get())
            if self.doInvert:
                args += ' --invert'
            if self.doNormalize:
                args += ' --normproc %s' % self.getEnumText('normproc')
            self.cshrink = float(samplingRateCoord/samplingRateTomo)
            if self.cshrink > 1:
                args += ' --cshrink %d' % self.cshrink
            program = eman2.Plugin.getProgram('e2spt_boxer_old.py')
            self.runJob(program, args, cwd=self._getExtraPath())
            moveFile(self._getExtraPath(os.path.join('sptboxer_01', 'basename.hdf')),
                     self._getExtraPath(pwutils.replaceBaseExt(tomo, 'hdf')))
            cleanPath(self._getExtraPath("sptboxer_01"))

    def createOutputStep(self):
        self.outputSubTomogramsSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        self.outputSubTomogramsSet.setSamplingRate(self.getInputTomograms().getSamplingRate() / self.downFactor.get())
        self.outputSubTomogramsSet.setCoordinates3D(self.inputCoordinates)
        acquisition = TomoAcquisition()
        acquisition.setAngleMin(self.getInputTomograms().getFirstItem().getAcquisition().getAngleMin())
        acquisition.setAngleMax(self.getInputTomograms().getFirstItem().getAcquisition().getAngleMax())
        acquisition.setStep(self.getInputTomograms().getFirstItem().getAcquisition().getStep())
        self.outputSubTomogramsSet.setAcquisition(acquisition)
        for item in self.getInputTomograms().iterItems():
            for ind, tomoFile in enumerate(self.tomoFiles):
                if os.path.basename(tomoFile) == os.path.basename(item.getFileName()):
                    coordSet = self.lines[ind]
                    outputSet = self.readSetOfSubTomograms(self._getExtraPath(pwutils.replaceBaseExt(tomoFile, "hdf")),
                                            self.outputSubTomogramsSet, coordSet)

        self._defineOutputs(outputSetOfSubtomogram=outputSet)
        self._defineSourceRelation(self.inputCoordinates, outputSet)

    # --------------------------- INFO functions --------------------------------

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
                   % (str(self.inputCoordinates.get().getSize()), self.boxSize.get()))
            if self._tomosOther():
                msg += (" from another set of tomograms: %s"
                        % self.getObjectTag('inputTomogram'))
            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Set of Subtomograms not ready yet")
        if self.doInvert:
            methodsMsgs.append("Inverted contrast on images.")
        if self.downFactor.get() != 1:
            methodsMsgs.append("Subtomograms downsample by factor %d."
                               % self.downFactor.get())
        if self.doNormalize:
            methodsMsgs.append("Particles were normalised. Using normalization method %s") % \
            self.getEnumText('normproc')
        return methodsMsgs

    def _summary(self):
        summary = []
        summary.append("Tomogram source: %s"
                       % self.getEnumText("tomoSource"))
        if self.getOutputsSize() >= 1:
            summary.append("Particle box size: %s" % self.boxSize.get())
            summary.append("Subtomogram extracted: %s" %
                           self.inputCoordinates.get().getSize())
        else:
            summary.append("Output subtomograms not ready yet.")
        return summary

    # --------------------------- UTILS functions ----------------------------------

    def _tomosOther(self):
        """ Return True if other tomograms are used for extract. """
        return self.tomoSource == OTHER

    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        if self.tomoSource.get() == SAME_AS_PICKING:
            return self.inputCoordinates.get().getPrecedents()
        else:
            return self.inputTomograms.get()

    def readSetOfSubTomograms(self, workDir, outputSubTomogramsSet, coordSet):
        imgh = ImageHandler()
        x, y, z, n = imgh.getDimensions(workDir)
        for index in range(1, n + 1):
            subtomogram = SubTomogram()
            subtomogram.cleanObjId()
            subtomogram.setLocation(index, workDir)
            dfactor = self.downFactor.get()
            if dfactor != 1:
                fnSubtomo = self._getExtraPath("downsampled_subtomo%d.mrc" % index)
                ImageHandler.scaleSplines(subtomogram.getLocation(), fnSubtomo, dfactor)
                subtomogram.setLocation(fnSubtomo)
            subtomogram.setCoordinate3D(coordSet[index-1])
            outputSubTomogramsSet.append(subtomogram)
        return outputSubTomogramsSet
