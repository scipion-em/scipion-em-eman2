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
import json
import fileinput
import sys
from os.path import exists

from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.protocol.params import BooleanParam, PointerParam
from pyworkflow import utils as pwutils

import eman2
from eman2.convert import loadJson, writeJson

from tomo.protocols import ProtTomoPicking
from tomo.objects import Coordinate3D, SetOfCoordinates3D


class EmanProtTomoBoxing(ProtTomoPicking):
    """ Manual picker for Tomo. Uses EMAN2 e2spt_boxer.py.
    """
    _label = 'tomo boxer'

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)

    @classmethod
    def isDisabled(cls):
        return not eman2.Plugin.isTomoAvailableVersion()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)

        form.addParam('inMemory', BooleanParam, default=False,
                      label='Read in Memory',
                      help='This will read the entire tomogram into memory.'
                           'Much faster, but you must have enough ram.')
        form.addParam('inputCoordinates', PointerParam, label="Input Coordinates",
                      allowsNull=True, pointerClass='SetOfCoordinates3D',
                      help='Select the SetOfCoordinates3D.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self.inputTomo = self.inputTomogram.get()
        self._params = {'inputTomogram': self.inputTomo.getFileName()}
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', self._params, interactive=True)

    def _readCoordinates3D(self, box, tomo, coord3DSet):
        x, y, z = box[:3]
        coord = Coordinate3D()
        coord.setPosition(x, y, z)
        coord.setVolume(tomo)
        coord3DSet.append(coord)

    def _readSetOfCoordinates3D(self, jsonBoxDict, tomo, coord3DSetDict):
        if jsonBoxDict.has_key("boxes_3d"):
            boxes = jsonBoxDict["boxes_3d"]

            for box in boxes:
                classKey = box[5]
                coord3DSet = coord3DSetDict[classKey]
                coord3DSet.enableAppend()

                self._readCoordinates3D(box, tomo, coord3DSet)

    def _createOutput(self, outputDir):
        jsonFnbase = pwutils.join(outputDir, 'info',
                                  'extra-%s_info.json'
                                  % pwutils.removeBaseExt(self.inputTomo.getFileName()))
        jsonBoxDict = loadJson(jsonFnbase)

        # Create a Set of 3D Coordinates per class
        coord3DSetDict = {}
        coord3DMap = {}
        for key, classItem in jsonBoxDict["class_list"].iteritems():
            index = int(key)
            suffix = self._getOutputSuffix(SetOfCoordinates3D)
            coord3DSet = self._createSetOfCoordinates3D(self.inputTomo, suffix)
            coord3DSet.setBoxSize(int(classItem["boxsize"]))
            coord3DSet.setName(classItem["name"])
            coord3DSet.setVolumes(self.inputTomogram)

            name = self.OUTPUT_PREFIX + suffix
            args = {}
            args[name] = coord3DSet
            self._defineOutputs(**args)
            self._defineSourceRelation(self.inputTomogram, coord3DSet)

            coord3DSetDict[index] = coord3DSet
            coord3DMap[index] = name

        # Populate Set of 3D Coordinates with 3D Coordinates
        self._readSetOfCoordinates3D(jsonBoxDict, self.inputTomo, coord3DSetDict)

        # Update Outputs
        for index, coord3DSet in coord3DSetDict.iteritems():
            coord3DSet.setObjComment(self.getSummary(coord3DSet))
            self._updateOutputSet(coord3DMap[index], coord3DSet, state=coord3DSet.STREAM_CLOSED)


    # --------------------------- STEPS functions -----------------------------
    def launchBoxingGUIStep(self, tomo):
        inputCoor = self.inputCoordinates.get()
        if inputCoor is not None:
            infoDir = 'info'
            fnInputCoor = 'extra-%s_info.json' % pwutils.removeBaseExt(self.inputTomo.getFileName())
            pathInputCoor = pwutils.join(infoDir, fnInputCoor)
            if not exists(pathInputCoor):
                pwutils.makePath(infoDir)
            f = open(pathInputCoor, 'w')
            initFile = '{\n"boxes_3d": [\n\n],\n"class_list": {\n"0": {\n"boxsize": 32,\n"name": "particles_00"\n}\n}\n}'
            f.write(initFile)
            f.close()
            linCoor = '[%d, %d, %d, "manual", 0.0, 0]' % (inputCoor.getFirstItem().getX(),
                                                          inputCoor.getFirstItem().getY(),
                                                          inputCoor.getFirstItem().getZ())
            r = open(pathInputCoor, "r")
            contents = r.readlines()
            r.close()
            contents.insert(2, linCoor)
            w = open(pathInputCoor, "w")
            contents = "".join(contents)
            w.write(contents)
            w.close()
            for coor in inputCoor.iterCoordinates():
                linCoor = '[%d, %d, %d, "manual", 0.0, 0],\n' % (coor.getX(), coor.getY(), coor.getZ())
                r = open(pathInputCoor, "r")
                contents = r.readlines()[1:]
                r.close()
                contents.insert(2, linCoor) # Do not read first line (already read above) # change to 3
                w = open(pathInputCoor, "w")
                contents = "".join(contents)
                w.write(contents)
                w.close()

        program = eman2.Plugin.getProgram("e2spt_boxer.py")
        arguments = "%(inputTomogram)s"
        if self.inMemory:
            arguments += " --inmemory"
        self._log.info('Launching: ' + program + ' ' + arguments % tomo)
        self.runJob(program, arguments % tomo)
        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self._leaveDir()  # going back to project dir
            self._createOutput(self.getWorkingDir())

    def _validate(self):
        errors = []

        if not eman2.Plugin.isTomoAvailableVersion():
            errors.append('Your EMAN2 version does not support the tomo boxer. '
                          'Please update your installation to EMAN 2.3 or newer.')
        return errors

    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        ProtTomoPicking._runSteps(self, startIndex)

    def getMethods(self, output):
        msg = 'User picked %d particles ' % output.getSize()
        msg += 'with a particle size of %s.' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.inputTomogram is None:
            return ['Input tomogram not available yet.']

        methodsMsgs.append("Input tomogram %s of dims %s."
                           % (self.getObjectTag('inputTomogram'),
                              str(self.inputTomogram.get().getDim())))

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                msg = self.getMethods(output)
                methodsMsgs.append("%s: %s" % (self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

    def getSummary(self, coord3DSet):
        summary = []
        summary.append("Number of particles picked: %s" % coord3DSet.getSize())
        summary.append("Particle size: %s" % coord3DSet.getBoxSize())
        return "\n".join(summary)

    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append("Output 3D Coordinates not ready yet.")

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

