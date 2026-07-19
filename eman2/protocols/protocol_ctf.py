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
    Performs automated contrast transfer function estimation and particle
    preprocessing for cryo-EM single-particle workflows using EMAN2 tools.
    The protocol estimates microscope defocus parameters, evaluates image
    quality, and generates multiple particle outputs optimized for different
    downstream reconstruction and classification strategies.

    AI Generated:

    CTF Auto (EmanProtCTFAuto) — User Manual
        Overview

        The CTF Auto protocol provides an automated workflow for estimating
        and correcting the contrast transfer function (CTF) of cryo-EM
        particle images. In cryo-electron microscopy, accurate CTF estimation
        is one of the most important preprocessing steps because the microscope
        optics strongly modulate image contrast and attenuate structural
        information at different spatial frequencies. Correct estimation and
        compensation of these effects are essential before classification,
        refinement, or high-resolution reconstruction.

        This protocol is designed to simplify CTF processing by combining
        parameter estimation, particle conditioning, and generation of
        processed particle sets into a single workflow. It is particularly
        useful in large-scale cryo-EM processing pipelines where rapid and
        reproducible preprocessing is required.

        Biological Purpose and Typical Applications

        In practical cryo-EM studies, biological users employ this protocol
        after particle extraction and before extensive downstream analysis.
        The protocol prepares particles for classification, ab initio model
        generation, refinement, and structural interpretation by improving
        image consistency and compensating for microscope-induced phase
        distortions.

        Different biological projects may require different preprocessing
        strategies. High-resolution studies targeting atomic interpretation
        often prioritize preservation of fine structural detail, whereas
        medium- or low-resolution analyses may focus on robust particle
        classification or identification of major conformational states.
        This protocol supports these different goals by generating several
        processed particle outputs tailored to distinct resolution ranges.

        Inputs and Experimental Requirements

        The protocol requires a particle dataset together with acquisition
        metadata describing the imaging conditions. Accurate values for
        microscope voltage, spherical aberration, amplitude contrast, and
        pixel size are essential because they directly influence CTF
        estimation quality.

        The input particles should generally represent raw, non-phase-flipped
        images. Using particles that have already undergone CTF correction
        can lead to inaccurate estimation and biologically misleading results.
        Similarly, datasets lacking proper acquisition metadata are unsuitable
        for reliable CTF analysis.

        In many cryo-EM workflows, previously estimated CTF parameters may
        already exist. The protocol can either refine existing estimations or
        perform a complete estimation from scratch. Refinement is typically
        faster and useful when prior estimates are trustworthy, while full
        re-estimation is preferable for problematic datasets or after major
        preprocessing changes.

        Resolution-Oriented Processing Modes

        The protocol supports several operating modes corresponding to
        different target resolution regimes. These modes are intended to adapt
        preprocessing behavior to the biological objective of the study.

        The high-resolution mode is optimized for projects aiming at near-
        atomic or atomic detail. It attempts to preserve high-frequency
        information and is most appropriate for well-behaved datasets with
        high signal quality and stable microscope conditions.

        The medium-resolution mode is often appropriate for routine structural
        studies, intermediate refinement stages, and conformational analysis.
        It balances robustness and detail preservation and is frequently used
        during exploratory classification workflows.

        The low-resolution mode emphasizes stability and noise suppression.
        It is particularly useful for heterogeneous datasets, challenging
        membrane proteins, flexible assemblies, or early-stage exploratory
        analyses where major structural features are more important than fine
        detail.

        Defocus Search and Optical Estimation

        One of the most biologically important parameters in CTF analysis is
        the defocus search range. The protocol allows users to define the
        expected minimum and maximum defocus values in microns. Appropriate
        ranges improve estimation reliability and reduce the risk of incorrect
        fitting.

        In many practical experiments, underfocus values depend strongly on
        imaging strategy and specimen thickness. Narrow search ranges are
        generally preferable when acquisition conditions are well controlled,
        while broader ranges may be necessary for heterogeneous collections
        or historical datasets.

        The protocol can also estimate astigmatism, which models directional
        distortions introduced by imperfections in the microscope optics.
        Astigmatism estimation is often important for high-resolution work,
        although it may increase computational complexity and occasionally
        destabilize fitting for low-quality data.

        For experiments using phase plates, the protocol can estimate phase
        shifts in addition to conventional defocus parameters. This capability
        is especially relevant for hole-less phase plate imaging strategies
        where contrast formation differs substantially from conventional
        defocus-based imaging.

        Particle Conditioning and Filtering

        Beyond parameter estimation, the protocol generates multiple processed
        particle datasets designed for different analytical purposes. These
        outputs include phase-flipped particles at full resolution as well as
        low-pass filtered variants targeting specific spatial frequency ranges.

        Full-resolution phase-flipped particles are generally appropriate for
        refinement and reconstruction workflows where preservation of high-
        frequency information is essential. Low-pass filtered particles are
        often more stable for classification and exploratory analyses because
        they suppress noisy high-resolution components that may interfere with
        alignment or clustering.

        The protocol also generates invariant or bispectral representations of
        particles. These representations can be valuable in classification and
        orientation-independent analyses because they emphasize robust global
        image characteristics rather than fine alignment-sensitive details.

        Additional preprocessing options help adapt the workflow to difficult
        experimental conditions. Extra padding may improve behavior when
        particles were boxed too tightly. High-density handling strategies are
        useful when neighboring particles interfere with signal estimation.
        Contrast inversion can also be applied when imaging conditions produce
        inverted particle appearance.

        B-Factor Handling and Signal Interpretation

        The protocol supports both automatic and fixed B-factor estimation.
        Biologically, the B-factor describes the attenuation of high-frequency
        signal and is related to specimen disorder, beam-induced motion,
        radiation damage, and experimental noise.

        Automatic estimation is generally appropriate for routine processing
        because it adapts to the characteristics of the dataset. Fixed values
        may be useful in specialized workflows where consistent filtering
        behavior across multiple datasets is required.

        Outputs and Biological Interpretation

        After completion, the protocol produces several particle datasets with
        distinct preprocessing characteristics. These outputs are intended to
        support different downstream cryo-EM tasks rather than represent a
        single universally optimal result.

        Phase-flipped full-resolution particles are usually the preferred
        input for high-resolution refinement workflows. Low-pass filtered
        particles are often more appropriate for heterogeneous classification,
        ab initio reconstruction, or difficult alignment problems. Bispectral
        representations are useful for invariant-based analyses and specialized
        classification approaches.

        Because each output emphasizes different aspects of the signal, users
        should interpret results in the context of the biological question.
        Aggressive filtering may improve stability while simultaneously hiding
        subtle conformational variability or fine structural features.

        Practical Recommendations

        For most standard cryo-EM projects, medium-resolution processing with
        automatic parameter estimation provides a reliable starting point.
        High-resolution mode is most beneficial when particle quality and
        microscope stability are excellent and when downstream refinement aims
        for near-atomic interpretation.

        When working with noisy, heterogeneous, or flexible biological
        assemblies, low-pass filtered outputs often improve classification
        stability and reveal major conformational differences more clearly.
        Users should visually inspect representative particles and estimated
        CTF parameters to confirm that fitting remains biologically plausible.

        Astigmatism and phase-shift estimation should generally be enabled
        only when justified by the acquisition strategy or microscope behavior.
        Over-parameterization can occasionally reduce robustness in low-signal
        datasets.

        Final Perspective

        Automated CTF estimation is a foundational step in cryo-EM image
        processing because it directly influences all subsequent structural
        analyses. Reliable estimation and carefully chosen preprocessing
        strategies improve particle consistency, facilitate classification,
        and enhance reconstruction quality. Successful use of this protocol
        depends not only on computational settings but also on understanding
        the biological specimen, imaging conditions, and ultimate structural
        objectives of the experiment.
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
        self._insertFunctionStep('convertImagesStep', needsGPU=False)
        args = self._prepareParams()
        self._insertFunctionStep('runCTFStep', args, needsGPU=False)
        self._insertFunctionStep('createOutputStep', needsGPU=False)

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
