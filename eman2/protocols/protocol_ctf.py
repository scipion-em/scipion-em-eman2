# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pyworkflow.protocol.params import (FloatParam, EnumParam,
                                        BooleanParam)
from pyworkflow.constants import PROD
import pyworkflow.utils as pwutils
from pwem.objects.data import CTFModel, SetOfParticles
from pwem.protocols import ProtProcessParticles

from .. import Plugin
from ..constants import HIRES, INVAR_AUTO
from ..convert import writeSetOfParticles, iterLstFile, jsonToCtfModel


class EmanProtCTFAuto(ProtProcessParticles):
    """
    This protocol wraps *e2ctf_auto.py* EMAN2 program.
    It automates the CTF fitting and structure factor
    generation process.
    """

    _label = 'ctf auto'
    _devStatus = PROD
    _possibleOutputs = {
        'outputParticles_flip_fullRes': SetOfParticles,
        'outputParticles_flip_invar': SetOfParticles,
        'outputParticles_flip_lp5': SetOfParticles,
        'outputParticles_flip_lp7': SetOfParticles,
        'outputParticles_flip_lp12': SetOfParticles,
        'outputParticles_flip_lp20': SetOfParticles
    }

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """
        myDict = {
            'partSet': self._getExtraPath('sets/all.lst'),
            'partSetFlipInvar': self._getExtraPath('sets/all__ctf_flip_invar.lst'),
            'partSetFlipFullRes': self._getExtraPath('sets/all__ctf_flip_fullres.lst'),
            'partSetFlipLp5': self._getExtraPath('sets/all__ctf_flip_lp5.lst'),
            'partSetFlipLp7': self._getExtraPath('sets/all__ctf_flip_lp7.lst'),
            'partSetFlipLp12': self._getExtraPath('sets/all__ctf_flip_lp12.lst'),
            'partSetFlipLp20': self._getExtraPath('sets/all__ctf_flip_lp20.lst'),
        }
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineProcessParams(self, form):
        form.addParam('type', EnumParam,
                      choices=['hires', 'midres', 'lores'],
                      label='type', default=HIRES,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Performs CTF processing targeting '
                           'different resolution:\n'
                           '*hires* - 2-6 Angstrom\n*midres* - 7-15 Angstrom\n'
                           '*lores* - 15-30 Angstrom')
        line = form.addLine('Defocus search range (microns)',
                            help='Select _minimum_ and _maximum_ values for '
                                 'defocus search range (in microns). Underfocus'
                                 ' is represented by a positive number.')
        line.addParam('minDefocus', FloatParam, default=0.6,
                      label='Min')
        line.addParam('maxDefocus', FloatParam, default=4.,
                      label='Max')
        form.addParam('fromScratch', BooleanParam, default=False,
                      label='Fit from scratch',
                      help='Force refitting of CTF from scratch, ignoring any '
                           'previous fits.')
        form.addParam('astig', BooleanParam, default=False,
                      label='Estimate astigmatism',
                      help='Includes astigmatism in automatic fitting.')
        form.addParam('phaseEst', BooleanParam, default=False,
                      label='Estimate phase shift',
                      help='Include phase/amplitude contrast in CTF '
                           'estimation. For use with hole-less phase plates.')

        form.addSection(label='Advanced')
        form.addParam('extrapad', BooleanParam, default=False,
                      label='Extra padding',
                      help='If particles were boxed more tightly than EMAN '
                           'requires, this will add some extra padding.')
        form.addParam('invarType', EnumParam,
                      choices=['auto', 'bispec', 'harmonic'],
                      label='Invariant type', default=INVAR_AUTO,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Which type of invariants to generate')
        form.addParam('highDensity', BooleanParam, default=False,
                      label='High density ',
                      help='If particles are very close together, this will '
                           'interfere with SSNR estimation. '
                           'If set uses an alternative strategy, '
                           'but may over-estimate SSNR.')
        form.addParam('invert', BooleanParam, default=False,
                      label='Invert contrast',
                      help='Invert the contrast of the particles in output '
                           'files (default false)')
        form.addParam('constBfact', FloatParam, default=-1.0,
                      label='Constant B-factor',
                      help='Set B-factor to a fixed value, negative value '
                           'enables autofitting.')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertImagesStep')
        args = self._prepareParams()
        self._insertFunctionStep('runCTFStep', args)
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertImagesStep(self):
        partSet = self._getInputParticles()
        partAlign = partSet.getAlignment()
        storePath = self._getExtraPath("particles")
        pwutils.makePath(storePath)
        writeSetOfParticles(partSet, storePath, alignType=partAlign)

    def runCTFStep(self, args):
        """ Run the EMAN e2ctf_auto.py program. """
        program = Plugin.getProgram('e2ctf_auto.py')
        self.runJob(program, args, cwd=self._getExtraPath(),
                    numberOfThreads=1)

    def createOutputStep(self):
        inputSet = self._getInputParticles(pointer=True)
        outputSets = self._getOutputSets()
        outputs = {}

        for key, fn in outputSets.items():
            outputSet = self._createSetOfParticles(suffix='_%s' % key)
            outputSet.copyInfo(inputSet.get())
            outputSet.setIsPhaseFlipped(True)
            outputSet.setHasCTF(True)
            outputSet.copyItems(inputSet.get(),
                                updateItemCallback=self._updateCTF,
                                itemDataIterator=iterLstFile(self._getFileName(fn)))
            newPix = self._getNewPixSize(outputSet.getDimensions()[0])
            outputSet.setSamplingRate(newPix)

            summary = self.getSummary(key)
            outputSet.setObjComment(summary)
            if key == 'FL':
                outputName = 'outputParticles_flip_fullRes'
            elif key == 'bispec':
                outputName = 'outputParticles_flip_invar'
            else:
                outputName = 'outputParticles_flip_lp%s' % key

            outputs[outputName] = outputSet

        self._defineOutputs(**outputs)
        for _, out in self.iterOutputAttributes(SetOfParticles):
            self._defineSourceRelation(inputSet, out)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        partSet = self._getInputParticles()
        if partSet.isPhaseFlipped():
            errors.append('Input particles are already phase-flipped. '
                          'Please provide original raw particle images.')
        if not self.fromScratch and not partSet.hasCTF():
            errors.append('Input particles have no CTF information, '
                          'please select _Fit from scratch_ option.')
        if partSet.getAcquisition() is None:
            errors.append('Acquisition information missing for input '
                          'particles, you cannot estimate CTF!')

        return errors

    def _summary(self):
        summary = []

        if self.hasAttribute('outputParticles_flip_invar'):
            summary.append('CTF estimation on particles completed, '
                           'produced filtered particles and bispectra.')
        else:
            summary.append('Output is not ready yet.')

        return summary

    def getSummary(self, key):
        if key == 'FL':
            return 'Phase flipped, full resolution'
        elif key == 'bispec':
            return 'Bispectra footprints computed from high pass filtered normalized particles'
        else:
            return "Phase flipped, low-pass filtered to %d A" % int(key)

    # --------------------------- UTILS functions -----------------------------

    def _prepareParams(self):
        partSet = self._getInputParticles()
        acq = partSet.getAcquisition()
        args = "--%s" % self.getEnumText('type')
        args += " --voltage %3d" % acq.getVoltage()
        args += " --cs %0.3f" % acq.getSphericalAberration()
        args += " --ac %0.2f" % (100 * acq.getAmplitudeContrast())
        args += " --apix %0.3f" % partSet.getSamplingRate()

        if self.fromScratch:
            args += " --fromscratch"
        else:
            args += " --curdefocusfix"
        if self.astig:
            args += " --astigmatism"
        if self.phaseEst:
            args += " --phaseplate"
        if self.extrapad:
            args += " --extrapad"
        if self.highDensity:
            args += " --highdensity"
        if self.invert:
            args += " --invert"

        args += " --invartype %s" % self.getEnumText('invarType')
        args += " --constbfactor %0.2f --defocusmin %0.2f --defocusmax %0.2f" % (
            self.constBfact.get(),
            self.minDefocus.get(),
            self.maxDefocus.get())
        args += " --threads %d" % self.numberOfThreads.get()
        args += " --minqual 0"

        return args

    def _getInputParticles(self, pointer=False):
        return self.inputParticles if pointer else self.inputParticles.get()

    def _updateCTF(self, item, row):
        fileName = self._getExtraPath(row[1])
        item.setLocation(row[0], fileName)
        if not item.hasCTF():
            item.setCTF(CTFModel())
        jsonToCtfModel(fileName, item.getCTF())

    def _getOutputSets(self):
        protType = self.getEnumText('type')
        outputs = {}
        if protType == 'hires':
            outputs.update({'FL': 'partSetFlipFullRes',
                            '12': 'partSetFlipLp12',
                            '5': 'partSetFlipLp5'})
        elif protType == 'midres':
            outputs.update({'FL': 'partSetFlipFullRes',
                            '20': 'partSetFlipLp20',
                            '7': 'partSetFlipLp7'})
        else:  # lores
            outputs.update({'20': 'partSetFlipLp20',
                            '12': 'partSetFlipLp12'})

        outputs['bispec'] = 'partSetFlipInvar'

        return outputs

    def _getNewPixSize(self, newBox):
        # calculates new pix size for binned particles
        inputParts = self.inputParticles.get()
        oldDimX = inputParts.getDimensions()[0]
        oldPixSize = inputParts.getSamplingRate()
        newPixSize = float(oldDimX) / newBox * oldPixSize
        return newPixSize
