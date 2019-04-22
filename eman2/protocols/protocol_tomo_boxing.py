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

from pyworkflow.object import String

from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutils

import eman2
from eman2.convert import loadJson


class EmanProtTomoBoxing(EMProtocol):
    """ Semi-automated particle picker for SPA. Uses EMAN2 e2boxer.py.
    """
    _label = 'tomo boxer'
    OUTPUT_PREFIX = 'output3DCoordinates'

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, label="Volume to compare", important=True,
                      pointerClass='Volume',
                      help='This volume will be compared to the reference volume.')

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # TODO
        # input should be volume or set of volumes?
        # should we ask for an option to hide file menu?
        # clean up eman tmp files?
        # it looks like we can have different boxsizes for different classes, do we want to support it?
        # review input files? what is missing?
        # inverseY?
        self.inputVol = self.inputVolume.get().getFileName()
        self._params = {'inputVol': self.inputVol}
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    # TODO
    # Extracted from protocol_particles_picking
    def __getOutputSuffix(self):
        """ Get the name to be used for a new output.
        For example: outputCoordiantes7.
        It should take into account previous outputs
        and number with a higher value.
        """
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(SetOf3DCoordinates):
            suffix = attrName.replace(self.OUTPUT_PREFIX, '')
            try:
                counter = int(suffix)
            except:
                counter = 1 # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter+1) if maxCounter > 0 else '' # empty if not output

    def _createOutput(self, outputDir):
        inputVol = self.inputVolume.get()
        suffix = self.__getOutputSuffix()
        outputName = self.OUTPUT_PREFIX + suffix

        coord3DSet = self._createSetOfCoordinates(inputVol, suffix)


        # TODO
        # Extracted from eman2 convert
        # Assumes single volume
        # read boxSize from info/project.json
        jsonFnbase = pwutils.join(outputDir, 'info', 'extra-%s_info.json' % pwutils.removeBaseExt(self.inputVolume.get().getFileName()))
        jsonBoxDict = loadJson(jsonFnbase)
        size = int(jsonBoxDict["class_list"]["0"]["boxsize"])

        if jsonBoxDict.has_key("boxes_3d"):
            boxes = jsonBoxDict["boxes_3d"]

            for box in boxes:
                x, y, z = box[:3]

                coord = Coordinate_3D()
                coord.setPosition(x, y, z)
                coord.setVolume(inputVol)
                coord3DSet.append(coord)

        coord3DSet.setBoxSize(size)
        coord3DSet.setObjComment(self.getSummary(coord3DSet))

        outputs = {outputName: coord3DSet}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inputVolume, coord3DSet)

    # --------------------------- STEPS functions -----------------------------
    def launchBoxingGUIStep(self):

        program = eman2.Plugin.getProgram("e2spt_boxer.py")

        arguments = "%(inputVol)s"
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        self.runJob(program, arguments % self._params)

        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self._leaveDir()  # going back to project dir
            self._createOutput(self.getWorkingDir())
    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []


        return errors

    # def _warnings(self):
    #     warnings = []
    #     firstMic = self.inputMicrographs.get().getFirstItem()
    #     fnLower = firstMic.getFileName().lower()
    #
    #     ext = getExt(fnLower)
    #
    #     if ext in ['.tif', '.dm3'] and not self.invertY.get():
    #         warnings.append(
    #             'We have seen a flip in Y when using %s files in EMAN2' % ext)
    #         warnings.append(
    #             'The generated coordinates may or may not be valid in Scipion.')
    #         warnings.append(
    #             'TIP: Activate "Invert Y coordinates" if you find it wrong.')
    #     return warnings

    # --------------------------- UTILS functions -----------------------------
    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        EMProtocol._runSteps(self, startIndex)

    def getSummary(self, coord3DSet):
        summary = []
        summary.append("Summary")
        return "\n".join(summary)

