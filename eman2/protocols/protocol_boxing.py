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
    """ Semi-automated particle picker for SPA. Uses EMAN2 e2boxer.py.
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
