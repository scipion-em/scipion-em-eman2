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
from pyworkflow.protocol.params import PointerParam
from pyworkflow.protocol.constants import STATUS_FINISHED
import pyworkflow.em as pwem
import pyworkflow.protocol.params as params

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram

import eman2
from eman2.constants import *

# Micrograph type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1

class EmanProtTomoExtraction(pwem.EMProtocol, ProtTomoBase):
    """ Manual picker for Tomo. Uses EMAN2 e2spt_boxer.py.
    """
    _label = 'tomo extraction'
    OUTPUT_PREFIX = 'outputSetOfSubtomogram'

    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)

    @classmethod
    def isDisabled(cls):
        return not eman2.Plugin.isNewVersion()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates", important=True,
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
                      label='Micrographs source',
                      help='By default the particles will be extracted '
                           'from the micrographs used in the picking '
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide '
                           'a different set of micrographs to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided '
                           'micrographs and coordinates are related '
                           'by micName or by micId. Difference in pixel size '
                           'will be handled automatically.')

        form.addParam('inputTomogram', params.PointerParam,
                      pointerClass='SetOfMicrographs',
                      condition='downsampleType != %s' % SAME_AS_PICKING,
                      important=True, label='Input tomogram',
                      help='Select the tomogram from which to extract.')

        form.addSection(label='Preprocess')
        form.addParam('doInvert', params.BooleanParam, default=None,
                      label='Invert contrast?',
                      help='Invert the contrast if your particles are black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background. Frealign (up to v9.07) requires black '
                           'particles over a white background')

        form.addParam('doNormalize', params.BooleanParam, default=True,
                      label='Normalize particles?',
                      help='Normalization processor applied to particles before extraction.')
        form.addParam('normproc', params.EnumParam,
                      choices=['normalize', 'normalize.edgemean'],
                      label='',
                      condition='doNormalize',
                      default=PROC_NORMALIZE,
                      display=params.EnumParam.DISPLAY_COMBO,
                      help='Use normalize.edgemean if the particles have a'
                        ' clear solvent background (i.e., they are not part of a '
                        'larger complex or embeded in a membrane)')

        form.addParam('doRescale', params.BooleanParam, default=False,
                      label='Coordinates Factor?',
                      help='Specifies the factor by which to multiply the'
                        'coordinates, so that they can be at the'
                        'same scale as the tomogram.')

        form.addParam('cshrink', params.IntParam,
                      validators=[params.Positive],
                      condition='doRescale',
                      label='Factor',
                      help='For example, provide --cshrink=2 if the coordinates'
                        'were determined in a binned-by-2 (shrunk-by-2)'
                        'tomogram, but you want to extract the subvolumes from'
                        'a tomogram without binning/shrinking (which should be'
                        '2x larger).')



        # cshrink, normproc, invertedY

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _micsOther(self):
        """ Return True if other micrographs are used for extract. """
        return self.downsampleType == OTHER

    def getInputTomogram(self):
        """ Return the micrographs associated to the SetOfCoordinates or
        Other micrographs. """
        if not self._micsOther():
            return self.inputCoordinates.get().getVolumes()
        else:
            return self.inputTomogram.get()

    def _insertAllSteps(self):

        self._insertFunctionStep('extractParticles')
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):

        suffix = self._getOutputSuffix(SetOfSubTomograms)
        self.outputSubTomogramsSet = self._createSetOfSubTomograms(suffix)
        subtomogram=SubTomogram()
        subtomogram.setFileName(self._getExtraPath(pwutils.join('sptboxer_01', 'basename.hdf')))
        subtomogram.setSamplingRate(self.getInputTomogram().getSamplingRate())
        self.outputSubTomogramsSet.append(subtomogram)

        self._defineOutputs(outputSetOfSubtomogram=self.outputSubTomogramsSet)
        self._defineSourceRelation(self.inputCoordinates, self.outputSubTomogramsSet)


    def writeSetOfCoordinates3D(self):
        self.coordsFileName = self._getExtraPath(
            pwutils.replaceBaseExt(self.getInputTomogram().getFileName(), 'coords'))

        out = file(self.coordsFileName, "w")
        for coord3DSet in self.inputCoordinates.get().iterCoordinates():
            out.write("%d\t%d\t%d\n" % (coord3DSet.getX(), coord3DSet.getY(), coord3DSet.getZ()))
        out.close()

    # --------------------------- STEPS functions -----------------------------
    def extractParticles(self):

        self.writeSetOfCoordinates3D()
        self.boxSize = self.inputCoordinates.get().getBoxSize()

        args = '%s --coords %s --boxsize %d' % (
            self.getInputTomogram().getFileName(), pwutils.replaceBaseExt(self.getInputTomogram().getFileName(), 'coords'),
            self.boxSize)

        if self.doInvert:
            args += ' --invert'

        if self.doNormalize:
            args += ' --normproc '

        program = eman2.Plugin.getProgram('e2spt_boxer.py')
        self.runJob(program, args,
                    cwd=self._getExtraPath())

    def _validate(self):
        errors = []

        if not eman2.Plugin.isNewVersion():
            errors.append('Your EMAN2 version does not support the tomo boxer. '
                          'Please update your installation to EMAN 2.21 or newer.')

        return errors


    def _methods(self):
        methodsMsgs = []


        if self.getStatus() == STATUS_FINISHED:
            msg = ("A total of %d particles of size %d were extracted"
                   % (self.outputSubTomogramsSet.getDim(), self.boxSize))

            if self._micsOther():
                msg += (" from another set of micrographs: %s"
                        % self.getObjectTag('inputMicrographs'))

            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)

        # if self.doInvert:
        #     methodsMsgs.append("Inverted contrast on images.")
        # if self._doDownsample():
        #     methodsMsgs.append("Particles downsampled by a factor of %0.2f."
        #                        % self.downFactor)
        # if self.doNormalize:
        #     methodsMsgs.append("Particles were normalised.")

        return methodsMsgs


    def _summary(self):
        summary = []
        summary.append("Micrographs source: %s"
                       % self.getEnumText("downsampleType"))
        summary.append("Particle box size: %d" % self.boxSize)

        if not hasattr(self, 'outputSubTomogramsSet'):
            summary.append("Output subtomograms not ready yet.")
        else:
            summary.append("Subtomograms extracted: %d" %
                           self.outputSubTomogramsSet.getDim())

        return summary

