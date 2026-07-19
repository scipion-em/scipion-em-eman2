# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (josue.gomez-blanco@mcgill.ca) [1]
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pyworkflow.object import String
from pyworkflow.constants import PROD
from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import getExt
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.protocol.params import BooleanParam, IntParam, StringParam
from pwem.protocols import ProtParticlePicking

from .. import Plugin
from ..convert import readSetOfCoordinates


class EmanProtBoxing(ProtParticlePicking):
    """
    Provides a semi-automated particle picking environment for single particle
    analysis workflows using the EMAN2 boxing framework. The protocol is intended
    to help users identify and extract particle coordinates from cryo-EM
    micrographs through an interactive graphical interface that combines manual
    supervision with neural-network-assisted particle detection.

    AI Generated:

    EMAN Boxer Particle Picking (EmanProtBoxing) - User Manual
        Overview

        The EMAN Boxer protocol is designed for the interactive selection of
        particle coordinates from cryo-EM micrographs prior to particle extraction
        and downstream reconstruction workflows. In single particle analysis,
        particle picking is one of the most biologically influential preprocessing
        steps because the quality of the selected particles directly determines
        the quality of all subsequent analyses, including classification,
        refinement, and three-dimensional reconstruction.

        This protocol provides a semi-automated environment in which the user can
        combine manual expertise with machine-assisted picking strategies. Rather
        than relying entirely on automated detection, the workflow encourages
        iterative visual inspection and refinement, allowing the user to maintain
        biological control over the particle selection process.

        Inputs and General Workflow

        The protocol requires a set of input micrographs representing the raw or
        preprocessed cryo-EM images from which particles will be selected. These
        micrographs should ideally possess good contrast, accurate calibration,
        and limited contamination. While automated methods can tolerate moderate
        levels of noise, extremely poor micrograph quality or strong ice
        contamination can significantly reduce picking reliability.

        During execution, the protocol launches an interactive graphical
        interface where particles can be selected manually, automatically, or
        through iterative neural-network-assisted training. The user may begin by
        selecting representative particles and contaminants, allowing the system
        to learn the appearance of true particles relative to background noise or
        artifacts.

        In practical biological workflows, users often alternate between manual
        correction and automated prediction until the selected coordinates
        adequately represent the particle population across the dataset.

        Particle and Box Size Selection

        Two biologically important parameters are the particle size and the box
        size. The particle size should correspond approximately to the largest
        visible dimension of the molecular complex in the micrograph. Accurate
        estimation improves particle centering and helps the picking algorithms
        distinguish true particles from contaminants or ice features.

        The box size defines the region extracted around each selected particle.
        A box that is too small may truncate peripheral structural information,
        while a box that is too large introduces unnecessary background noise and
        increases computational cost. In most cryo-EM workflows, the selected box
        should comfortably enclose the particle together with a margin of solvent
        around it.

        Biological users should verify that the chosen dimensions remain
        appropriate across the full dataset, particularly when particles display
        multiple orientations or variable conformations.

        Interactive and Neural-Network-Assisted Picking

        One of the strengths of the protocol is its support for iterative
        learning-based picking strategies. The user can provide examples of good
        particles, bad particles, and background regions, enabling the system to
        improve particle discrimination progressively.

        This approach is especially valuable for challenging datasets containing
        heterogeneous particle populations, low contrast, preferred orientations,
        or strong contamination. In such cases, fully automated methods may
        generate large numbers of false positives or miss biologically relevant
        particles entirely.

        The interactive nature of the workflow allows the user to continuously
        evaluate the biological validity of the detected particles. Careful human
        supervision remains essential because even sophisticated automated methods
        may preferentially select contaminants, carbon edges, ice crystals, or
        aggregated particles if not properly trained.

        Device Selection and Computational Considerations

        The protocol allows execution on either CPU or GPU hardware. GPU-based
        execution is particularly advantageous during neural-network training and
        prediction because it substantially accelerates iterative optimization.
        For large cryo-EM datasets, GPU acceleration can reduce picking time from
        hours to minutes.

        Nevertheless, CPU execution remains suitable for smaller projects,
        exploratory analyses, or environments without dedicated accelerators.
        Biological users should balance computational speed with resource
        availability within their local infrastructure.

        Coordinate System Considerations

        Different microscopy image formats occasionally exhibit differences in
        coordinate orientation, particularly along the vertical axis. The protocol
        therefore provides an option to invert Y coordinates when necessary.
        This becomes important when importing coordinates into downstream Scipion
        workflows or external cryo-EM software packages.

        Users should always inspect the resulting particle positions visually
        after picking to confirm that the coordinates correctly overlay the
        biological particles in the original micrographs.

        Outputs and Their Interpretation

        The main output of the protocol is a validated set of particle
        coordinates associated with the input micrographs. These coordinates
        define the particle centers that will later be used for particle
        extraction and refinement.

        Biologically meaningful outputs depend not only on the quantity of picked
        particles but also on their quality and diversity. Overly aggressive
        picking strategies may include contaminants and reduce reconstruction
        quality, whereas excessively conservative picking may discard rare but
        important structural states.

        Users should therefore evaluate both particle purity and structural
        diversity before proceeding to downstream classification and refinement.

        Practical Recommendations

        In routine cryo-EM workflows, it is often advisable to begin with a small
        subset of representative micrographs and manually inspect the initial
        picks carefully. Once reliable particle examples are established,
        automated or neural-network-assisted picking can be expanded to the full
        dataset.

        For heterogeneous or low-contrast datasets, repeated cycles of manual
        correction and retraining usually produce substantially better results
        than fully unsupervised picking. It is also recommended to inspect picks
        near carbon edges, contamination regions, and crowded particle areas,
        where false positives are most common.

        Users should verify that particles remain centered and fully enclosed by
        the selected box size before committing to large-scale extraction.

        Final Perspective

        Particle picking is not simply a technical preprocessing step but a
        biologically critical selection process that determines which molecular
        observations contribute to the final reconstruction. Careful supervision,
        realistic particle sizing, and iterative validation are essential for
        obtaining reliable cryo-EM datasets suitable for high-quality structural
        interpretation.
    """
    _label = 'boxer'
    _devStatus = PROD

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {'goodRefsFn': self._getExtraPath('info/boxrefs.hdf'),
                  'badRefsFn': self._getExtraPath('info/boxrefsbad.hdf'),
                  'bgRefsFn': self._getExtraPath('info/boxrefsbg.hdf'),
                  'nnetFn': self._getExtraPath('nnet_pickptcls.hdf'),
                  'nnetClFn': self._getExtraPath('nnet_classify.hdf'),
                  'trainoutFn': self._getExtraPath('trainout_pickptcl.hdf'),
                  'trainoutClFn': self._getExtraPath('trainout_classify.hdf')
                  }
        self._updateFilenamesDict(myDict)

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtParticlePicking._defineParams(self, form)
        form.addParam('boxSize', IntParam, default=-1,
                      label='Box size (px)',
                      allowsPointers=True,
                      help="Box size in pixels.")
        form.addParam('particleSize', IntParam, default=-1,
                      label='Particle size (px)',
                      help="Longest axis of particle in pixels (diameter, "
                           "not radius).")
        form.addParam('device', StringParam, default='cpu',
                      label='Device',
                      help='For Convnet training only.\n'
                           'Pick a device to use. Choose from cpu, '
                           'gpu, or gpuX (X=0,1,...) when multiple '
                           'gpus are available. Default is cpu.')

        form.addParam('invertY', BooleanParam, default=False,
                      label='Invert Y coordinates',
                      help='In some cases, using dm3 or tiff Y coordinates '
                           'must be flipped. Check output and activate this'
                           ' if needed.')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self.inputMics = self.inputMicrographs.get()
        micList = [os.path.relpath(mic.getFileName(),
                                   self.getCoordsDir()) for mic in self.inputMics]

        self._params = {'inputMics': ' '.join(micList)}
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    # --------------------------- STEPS functions -----------------------------
    def launchBoxingGUIStep(self):
        # Print the eman version, useful to report bugs
        self.runJob(Plugin.getProgram('e2version.py'), '')
        # Program to execute and it arguments
        program = Plugin.getProgram('e2boxer.py')
        arguments = " --apix=%(pixSize)f --boxsize=%(boxSize)d"
        arguments += " --ptclsize=%(ptclSize)d --gui --threads=%(thr)d --no_ctf"

        acq = self.inputMics.getAcquisition()
        arguments += " --voltage %d" % acq.getVoltage()
        arguments += " --cs %0.2f" % acq.getSphericalAberration()
        arguments += " --ac %0.2f" % (100 * acq.getAmplitudeContrast())

        self._params.update({
            'pixSize': self.inputMics.getSamplingRate(),
            'boxSize': self.boxSize.get(),
            'ptclSize': self.particleSize.get(),
            'thr': self.numberOfThreads.get()
        })
        arguments += " --device=%s" % self.device.get()

        arguments += " %(inputMics)s"

        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(program, arguments % self._params, cwd=self.getCoordsDir())

        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self._createOutput(self.getCoordsDir())

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        return errors

    def _warnings(self):
        warnings = []
        firstMic = self.inputMicrographs.get().getFirstItem()
        fnLower = firstMic.getFileName().lower()

        ext = getExt(fnLower)

        if ext in ['.tif', '.dm3'] and not self.invertY.get():
            warnings.append(
                'We have seen a flip in Y when using %s files in EMAN2' % ext)
            warnings.append(
                'The generated coordinates may or may not be valid in Scipion.')
            warnings.append(
                'TIP: Activate "Invert Y coordinates" if you find it wrong.')
        return warnings

    # --------------------------- UTILS functions -----------------------------
    def getCoordsDir(self):
        return self._getExtraPath()

    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMics, coordSet,
                             self.invertY.get())
