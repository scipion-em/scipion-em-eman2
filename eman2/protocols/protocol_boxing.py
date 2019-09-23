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

from pyworkflow.object import String
from pyworkflow.utils.properties import Message
from pyworkflow.utils.path import join, getExt
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.protocol.params import (BooleanParam, IntParam,
                                        StringParam)

import eman2
from eman2.convert import loadJson, readSetOfCoordinates


class EmanProtBoxing(ProtParticlePicking):
    """ Semi-automated particle picker for SPA. Uses EMAN2 e2boxer.py.
    """
    _label = 'boxer'

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

        form.addParam('useNewBoxer', BooleanParam,
                      default=True,
                      label='Use new e2boxer?',
                      help='Requires EMAN version 2.21 or newer')
        form.addParam('boxSize', IntParam, default=-1,
                      condition='useNewBoxer',
                      label='Box size (px)',
                      help="Box size in pixels.")
        form.addParam('particleSize', IntParam, default=-1,
                      condition='useNewBoxer',
                      label='Particle size (px)',
                      help="Longest axis of particle in pixels (diameter, "
                           "not radius).")
        if self._isVersion23():
            form.addParam('device', StringParam, default='cpu',
                          condition='useNewBoxer',
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
        self.runJob(eman2.Plugin.getProgram('e2version.py'), '')
        useNewBoxer = self._useNewBoxer()
        # Program to execute and it arguments
        boxerVersion = 'old' if not useNewBoxer else 'new'
        boxer = eman2.Plugin.getBoxerCommand(boxerVersion=boxerVersion)
        program = eman2.Plugin.getProgram(boxer)
        arguments = ''

        if useNewBoxer:
            arguments = " --apix=%(pixSize)f --boxsize=%(boxSize)d"
            arguments += " --ptclsize=%(ptclSize)d --gui --threads=%(thr)d --no_ctf"
            self._params.update({
                'pixSize': self.inputMics.getSamplingRate(),
                'boxSize': self.boxSize.get(),
                'ptclSize': self.particleSize.get(),
                'thr': self.numberOfThreads.get()
            })
            if self._isVersion23():
                arguments += " --device=%s" % self.device.get()

        arguments += " %(inputMics)s"

        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(program, arguments % self._params, cwd=self.getCoordsDir())

        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self.check_gauss()
            self._createOutput(self.getCoordsDir())

    def check_gauss(self):
        if self._useNewBoxer():
            # gauss picker is not implemented for new e2boxer
            pass
        else:
            # Function to check if gaussian algorithm was used to pick
            # and if so ask user if she wants to perform an automatic
            # picking for the remaining micrographs
            gaussJsonFile = join("e2boxercache", "gauss_box_DB.json")
            # Check if gauss json file exists and load it
            if os.path.exists(gaussJsonFile):
                jsonGaussDict = loadJson(gaussJsonFile)
                gaussParsDict = None
                micList = [os.path.relpath(mic.getFileName(),
                                           self.getCoordsDir()) for mic in self.inputMics]
                # Go over the list of input micrographs and see if
                # gaussian was used to pick any of them
                for mic in micList:
                    if mic in jsonGaussDict:
                        gaussParsDict = jsonGaussDict[mic]
                        break
                if gaussParsDict is not None:
                    # If found ask user if she wats to perform an automatic gaussian
                    # picking for the rest of mics
                    self._params['boxSize'] = gaussParsDict['boxsize']
                    # Run sxprocess.py to store parameters
                    program = eman2.Plugin.getProgram("sxprocess.py")
                    argsList = ["'%s'=%s:" % (key, val) for (key, val) in gaussParsDict.iteritems()]
                    args = 'demoparms --makedb ' + "".join(argsList)
                    # Remove last ":" to avoid error
                    args = args[:-1]
                    # Run the command with formatted parameters
                    self._log.info('Launching: ' + program + ' ' + args)
                    self.runJob(program, args)
                    # Now run e2boxer.py with stored parameters
                    arguments = "--gauss_autoboxer=demoparms --write_dbbox "
                    arguments += " --boxsize=%(boxSize)s " + "%(inputMics)s"
                    boxer = eman2.Plugin.getBoxerCommand(boxerVersion='old')
                    program = eman2.Plugin.getProgram(boxer)
                    self._log.info('Launching: ' + program + ' ' + arguments % self._params)
                    self.runJob(program, arguments % self._params, cwd=self.getCoordsDir())

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
                             self.invertY.get(), newBoxer=self._useNewBoxer())

    def _useNewBoxer(self):
        return True if self.useNewBoxer else False

    def _isVersion23(self):
        return eman2.Plugin.isVersion('2.3')
