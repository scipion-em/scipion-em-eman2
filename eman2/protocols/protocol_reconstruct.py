# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from enum import Enum

from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        EnumParam, StringParam, BooleanParam,
                                        LEVEL_ADVANCED)
from pyworkflow.constants import PROD
from pyworkflow.utils.path import cleanPattern, makePath
from pwem.objects.data import Volume
from pwem.protocols import ProtReconstruct3D

from .. import Plugin
from ..convert import writeSetOfParticles
from ..constants import (RECON_FOURIER, FOURIER_GAUSS2, KEEP_PERCENTAGE,
                         KEEP_STDDEV, KEEP_ABSQUAL)


class outputs(Enum):
    outputVolume = Volume


class EmanProtReconstruct(ProtReconstruct3D):
    """
    This protocol wraps *e2make3d.py* EMAN2 program.

    Reconstructs 3D volumes using a set of 2D images.
    Euler angles are extracted from the 2D image headers
    and symmetry is imposed. Several reconstruction methods
    are available. The fourier method is the default and
    recommended reconstructor.
    """

    _label = 'reconstruct'
    _devStatus = PROD
    _possibleOutputs = outputs

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {
            'partSet': 'sets/inputSet.lst',
            'partFlipSet': 'sets/inputSet__ctf_flip.lst',
            'volume': self._getExtraPath('volume.hdf'),
        }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles",
                      pointerCondition='hasAlignmentProj',
                      help='Select the input images from the project.')
        form.addParam('skipctf', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Skip ctf estimation?',
                      help='Use this if you want to skip running e2ctf.py. '
                           'It is not recommended to skip this step unless CTF '
                           'estimation was already done with EMAN2.')
        form.addParam('useE2make3d', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Use old e2make3d?',
                      help='Use the traditional e2make3d program instead of '
                           'the new e2make3dpar program.')
        form.addParam('numberOfIterations', IntParam, default=2,
                      condition='useE2make3d',
                      label='Number of iterations:',
                      help='Set the number of iterations. Iterative '
                           'reconstruction improves the overall normalization '
                           'of the 2D images as they are inserted into the '
                           'reconstructed volume, and allows for the '
                           'exclusion of the poorer quality images.')
        form.addParam('symmetry', StringParam, default='c1',
                      label='Symmetry group',
                      help='Set the symmetry; if no value is given then the '
                           'model is assumed to have no symmetry. \n'
                           'Choices are: *i, c, d, tet, icos, or oct* \n'
                           'See https://blake.bcm.edu/emanwiki/EMAN2/Symmetry \n'
                           'for a detailed description of symmetry in Eman.')
        line = form.addLine('Padding to Reconstruct: ',
                            expertLevel=LEVEL_ADVANCED,
                            help='Will zero-pad images to the specifed size '
                                 '(x,y) or (x,x) prior to reconstruction. '
                                 'If not specified no padding occurs.')
        line.addParam('padX', IntParam, default=0, label='X ')
        line.addParam('padY', IntParam, default=0, label='Y ')

        line = form.addLine('Dimensions Volume: ',
                            expertLevel=LEVEL_ADVANCED,
                            help='Defines the dimensions (x,y,z) or (x,x,x) '
                                 'of the reconstructed volume. If omitted, '
                                 'implied value based on padded 2D images '
                                 'is used. ')
        line.addParam('dimVolX', IntParam, default=0, label='X')
        line.addParam('dimVolY', IntParam, default=0, label='Y')
        line.addParam('dimVolZ', IntParam, default=0, label='Z')

        line = form.addLine('Dimensions to Write Volume: ',
                            expertLevel=LEVEL_ADVANCED,
                            help='Defines the dimensions (x,y,z) or (x,x,x) '
                                 'of the final volume written to disk, if '
                                 'omitted, size will be based on unpadded '
                                 'input size. ')
        line.addParam('dimWriteVolX', IntParam, default=0, label='X')
        line.addParam('dimWriteVolY', IntParam, default=0, label='Y')
        line.addParam('dimWriteVolZ', IntParam, default=0, label='Z')
        form.addParam('reconstructionMethod', EnumParam,
                      choices=['back_projection', 'fourier', 'fourier_iter',
                               'fouriersimple2D', 'nn4', 'nn4_ctf',
                               'nn4_ctf_rect', 'nn4_ctfw', 'nn4_ctfws', 'nn4_rect',
                               'nnSSNR', 'nnSSNR_ctf', 'real_median', 'wiener_fourier'],
                      label="Reconstruction Method:", default=RECON_FOURIER,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Reconstructor to use. See e2help.py reconstructors '
                           '-v 9. Default is fourier:mode=gauss_2.')
        form.addParam('fourierMode', EnumParam,
                      condition="reconstructionMethod==1 or reconstructionMethod==2 or reconstructionMethod==12",
                      choices=['nearest_neighbor', 'gauss_2', 'gauss_3',
                               'gauss_5', 'gauss_5_slow', 'gypergeom_5',
                               'experimental'],
                      label="Mode to Fourier method:", default=FOURIER_GAUSS2,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Fourier pixel insertion mode. See e2help.py reconstructors '
                           'fourier -v 9. Default mode is gauss_2.')
        form.addParam('keepSense', EnumParam, expertLevel=LEVEL_ADVANCED,
                      choices=['percentage', 'standard deviation',
                               'absolute quality'],
                      label="Sense of keep:", default=KEEP_PERCENTAGE,
                      display=EnumParam.DISPLAY_COMBO,
                      help="If *percentage* is selected, *keep* parameter "
                           "will be interpreted as a percentage. It is the "
                           "default option.\nIf *standard deviation* is "
                           "selected, *keep* parameter will be interpreted "
                           "as a standard deviation coefficient instead of "
                           "as a percentage.\n"
                           "If *absolute quality* is selected, *keep* "
                           "parameter will refer to the absolute quality "
                           "of the class-average, not a local quality "
                           "relative to other similar sized classes.")
        form.addParam('keep', FloatParam, default=1.0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Fraction of slices to keep",
                      help='The fraction of slices to keep, in fraction,'
                           ' based on quality scores (1.0 = use all slices).')
        form.addParam('doNotAutoWt', BooleanParam, default=False,
                      label='Disable automatic weighting?',
                      help='This argument turns automatic weighting off '
                           'causing all images to be weighted by 1. If this '
                           'argument is False images inserted into the '
                           'reconstructed volume are weighted by the number '
                           'of particles that contributed to them (i.e. as '
                           'in class averages), which is extracted from the '
                           'image header.')
        form.addParam('extraParams', StringParam,
                      expertLevel=LEVEL_ADVANCED,
                      default='',
                      label='Additional parameters',
                      help="In this box command-line arguments may be "
                           "provided that are not generated by the GUI. "
                           "See e2make3dpar.py -h.")

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertImagesStep')
        self._insertFunctionStep('reconstructVolumeStep', self._prepareParams())
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertImagesStep(self):
        partSet = self.inputParticles.get()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)
        if not self.skipctf:
            program = Plugin.getProgram('e2ctf.py')
            acq = partSet.getAcquisition()

            args = " --voltage %3d" % acq.getVoltage()
            args += " --cs %f" % acq.getSphericalAberration()
            args += " --ac %f" % (100 * acq.getAmplitudeContrast())
            if not partSet.isPhaseFlipped():
                args += " --phaseflip"
            args += " --computesf --apix %f " % partSet.getSamplingRate()
            args += " --allparticles --autofit --curdefocusfix --storeparm -v 8"
            args += " --threads=%d" % self.numberOfThreads.get()
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfThreads=1)

        program = Plugin.getProgram('e2buildsets.py')
        args = " --setname=inputSet --allparticles"
        self.runJob(program, args, cwd=self._getExtraPath(), numberOfThreads=1)

    def reconstructVolumeStep(self, args):
        """ Run the EMAN program to reconstruct a volume. """
        cleanPattern(self._getFileName("volume"))
        if self.useE2make3d:
            program = Plugin.getProgram('e2make3d.py')
        else:
            program = Plugin.getProgram('e2make3dpar.py')
        self.runJob(program, args, cwd=self._getExtraPath(), numberOfThreads=1)

    def createOutputStep(self):
        partSet = self.inputParticles.get()
        vol = Volume()
        vol.setFileName(self._getFileName("volume"))
        vol.copyInfo(partSet)
        self._defineOutputs(**{outputs.outputVolume.name: vol})
        self._defineSourceRelation(self.inputParticles, vol)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        if not self.useE2make3d and self.reconstructionMethod != RECON_FOURIER:
            errors.append('e2make3dpar.py program can use only Fourier method '
                          'for reconstruction!')

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Input images: %s" % self.getObjectTag('inputParticles'))
            summary.append("Output volume: %s" % self.getObjectTag('outputVolume'))
        return summary

    # --------------------------- UTILS functions -----------------------------

    def _prepareParams(self):
        args = "--input %(imgsFn)s --output %(outputVol)s --sym %(sym)s"

        if self.useE2make3d:
            args += " --iter %(numberOfIterations)d "
            args += " --recon %(reconsMethod)s"
        else:
            args += " --mode %s" % self.getEnumText('fourierMode')
            args += " --threads=%d" % self.numberOfThreads.get()

        if self.extraParams.hasValue():
            args += ' ' + self.extraParams.get()

        reconsMethod = self.getEnumText('reconstructionMethod')
        if (reconsMethod == 'fourier' or reconsMethod == 'fourier_plane' or
                reconsMethod == 'fouriersimple2D' or
                reconsMethod == 'wiener_fourier'):
            reconsMethod = reconsMethod + ':mode=' + self.getEnumText('fourierMode')

        params = {'imgsFn': self._getParticlesStack(),
                  'outputVol': self._getBaseName("volume"),
                  'numberOfIterations': self.numberOfIterations.get(),
                  'sym': self.symmetry.get(),
                  'reconsMethod': reconsMethod
                  }

        args %= params

        if self.padX.get() > 0:
            if self.padY.get() <= 0 or self.padX.get() == self.padY.get():
                args += " --pad %d" % self.padX.get()
            else:
                args += " --pad %d,%d" % (self.padX.get(), self.padY.get())

        if self.dimVolX.get() > 0:
            if ((self.dimVolY.get() <= 0 and self.dimVolZ.get() <= 0) or
                    (self.dimVolY.get() == self.dimVolX.get() and
                     self.dimVolZ.get() == self.dimVolX.get())):
                args += " --padvol %d" % self.dimVolX.get()
            else:
                args += " --padvol %d,%d,%d" % (self.dimVolX.get(),
                                                self.dimVolY.get(),
                                                self.dimVolZ.get())

        if self.dimWriteVolX.get() > 0:
            if ((self.dimWriteVolY.get() <= 0 and self.dimWriteVolZ.get() <= 0) or
                    (self.dimWriteVolY.get() == self.dimWriteVolX.get() and
                     self.dimWriteVolZ.get() == self.dimWriteVolX.get())):
                args += " --outsize %d" % self.dimWriteVolX.get()
            else:
                args += " --outsize %d,%d,%d" % (self.dimWriteVolX.get(),
                                                 self.dimWriteVolY.get(),
                                                 self.dimWriteVolZ.get())

        if self.keepSense == KEEP_STDDEV:
            args += "  --keep %f --keepsig" % self.keep.get()
        elif self.keepSense == KEEP_ABSQUAL:
            args += "  --keep %f --keepabs" % self.keep.get()

        if self.keep.get() != 1.0 and self.keepSense == KEEP_PERCENTAGE:
            args += " --keep %f" % self.keep.get()

        if self.doNotAutoWt:
            args += " --no_wt"

        return args

    def _getBaseName(self, key):
        """ Remove the folders and return the file from the filename. """
        return os.path.basename(self._getFileName(key))

    def _getParticlesStack(self):
        if not self.inputParticles.get().isPhaseFlipped() and not self.skipctf:
            return self._getFileName("partFlipSet")
        else:
            return self._getFileName("partSet")
