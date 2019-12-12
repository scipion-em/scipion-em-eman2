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

from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.protocol.params import BooleanParam, PointerParam, LEVEL_ADVANCED
from pyworkflow import utils as pwutils
from pyworkflow.utils import importFromPlugin

import eman2
from eman2.convert import loadJson, coordinates2json, readSetOfCoordinates3D
from eman2.constants import TOMO_NEEDED_MSG

ProtTomoPicking = importFromPlugin("tomo.protocols", "ProtTomoPicking", errorMsg=TOMO_NEEDED_MSG)
SetOfCoordinates3D = importFromPlugin("tomo.objects", "SetOfCoordinates3D", errorMsg=TOMO_NEEDED_MSG)


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
                      label='Read in Memory', expertLevel=LEVEL_ADVANCED,
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

    def _createOutput(self, outputDir):
        jsonFnbase = pwutils.join(outputDir, 'info',
                                  'extra-%s_info.json'
                                  % pwutils.removeBaseExt(self.inputTomo.getFileName()))
        jsonBoxDict = loadJson(jsonFnbase)

        coord3DSetDict = {}
        coord3DMap = {}
        setTomograms = self.getParentSet()
        for key, classItem in jsonBoxDict["class_list"].iteritems():
            index = int(key)
            suffix = self._getOutputSuffix(SetOfCoordinates3D)
            coord3DSet = self._createSetOfCoordinates3D(self.inputTomo, suffix)
            coord3DSet.setBoxSize(int(classItem["boxsize"]))
            coord3DSet.setName(classItem["name"])
            coord3DSet.setVolumes(setTomograms)
            coord3DSet.setSamplingRate(setTomograms.getSamplingRate())

            name = self.OUTPUT_PREFIX + suffix
            args = {}
            args[name] = coord3DSet
            self._defineOutputs(**args)
            self._defineSourceRelation(setTomograms, coord3DSet)

            coord3DSetDict[index] = coord3DSet
            coord3DMap[index] = name

        # Populate Set of 3D Coordinates with 3D Coordinates
        readSetOfCoordinates3D(jsonBoxDict, coord3DSetDict, self.inputTomo)

        # Update Outputs
        for index, coord3DSet in coord3DSetDict.iteritems():
            coord3DSet.setObjComment(self.getSummary(coord3DSet))
            self._updateOutputSet(coord3DMap[index], coord3DSet, state=coord3DSet.STREAM_CLOSED)

    # --------------------------- STEPS functions -----------------------------
    def _createInputCoordsFile(self):
        cwd = os.getcwd()
        infoDir = pwutils.join(cwd, 'info')
        self._leaveWorkingDir()
        fnInputCoor = 'extra-%s_info.json' % pwutils.removeBaseExt(self.inputTomo.getFileName())
        pathInputCoor = pwutils.join(infoDir, fnInputCoor)
        if not os.path.exists(pathInputCoor):
            pwutils.makePath(infoDir)
        return pathInputCoor

    def launchBoxingGUIStep(self, tomo):
        inputCoor = self.inputCoordinates.get()
        if inputCoor is not None:
            pathInputCoor = self._createInputCoordsFile()
            coordinates2json(pathInputCoor, inputCoor)
            self._enterWorkingDir()

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

    def getParentSet(self):
        parentObj = self.inputTomogram.getObjValue()
        return getattr(parentObj, 'outputTomograms')

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
