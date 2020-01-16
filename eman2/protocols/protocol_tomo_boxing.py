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

import eman2
from eman2.convert import loadJson, coordinates2json, readSetOfCoordinates3D
from eman2.viewers.views_tkinter_tree import EmanDialog

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider


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
        # Copy input coordinates to Extra Path
        if self.inputCoordinates.get():
            self._insertFunctionStep('copyInputCoords')

        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    def _createOutput(self):
        coord3DSetDict = {}
        coord3DMap = {}
        setTomograms = self.inputTomograms.get()
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coord3DSet = self._createSetOfCoordinates3D(setTomograms, suffix)
        coord3DSet.setName("tomoCoord")
        coord3DSet.setVolumes(setTomograms)
        coord3DSet.setSamplingRate(setTomograms.getSamplingRate())
        first = True
        for tomo in setTomograms.iterItems():
            jsonFnbase = pwutils.join(self._getExtraPath(),
                                      'extra-%s_info.json'
                                      % pwutils.removeBaseExt(tomo.getFileName()))
            if not os.path.isfile(jsonFnbase):
                continue
            jsonBoxDict = loadJson(jsonFnbase)

            for key, classItem in jsonBoxDict["class_list"].iteritems():
                index = int(key)
                if first:
                    coord3DSet.setBoxSize(int(classItem["boxsize"]))
                    first = False

                name = self.OUTPUT_PREFIX + suffix
                args = {}
                args[name] = coord3DSet
                coord3DSetDict[index] = coord3DSet
                coord3DMap[index] = name

                # Populate Set of 3D Coordinates with 3D Coordinates
                readSetOfCoordinates3D(jsonBoxDict, coord3DSetDict, tomo.clone())

        self._defineOutputs(**args)
        self._defineSourceRelation(setTomograms, coord3DSet)

        # Update Outputs
        for index, coord3DSet in coord3DSetDict.iteritems():
            coord3DSet.setObjComment(self.getSummary(coord3DSet))
            self._updateOutputSet(coord3DMap[index], coord3DSet, state=coord3DSet.STREAM_CLOSED)

    # --------------------------- STEPS functions -----------------------------
    def copyInputCoords(self):
        setTomograms = self.inputTomograms.get()
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        for tomo in setTomograms.iterItems():
            coord3DSet = self._createSetOfCoordinates3D(setTomograms, suffix)
            coord3DSet.setName("tomoCoord")
            coord3DSet.setVolumes(setTomograms)
            coord3DSet.setBoxSize(self.inputCoordinates.get().getBoxSize())
            coord3DSet.setSamplingRate(setTomograms.getSamplingRate())
            tomoName = os.path.basename(tomo.getFileName())
            for coord in self.inputCoordinates.get().iterCoordinates():
                if tomoName in coord.getVolName():
                    coord3DSet.append(coord)
            if not coord3DSet.isEmpty():
                fnInputCoor = 'extra-%s_info.json' % pwutils.removeBaseExt(tomo.getFileName())
                pathInputCoor = pwutils.join(self._getExtraPath(), fnInputCoor)
                coordinates2json(pathInputCoor, coord3DSet)
                pwutils.cleanPattern("tomoCoord*")

    def launchBoxingGUIStep(self):

        tomoList = [tomo.clone() for tomo in self.inputTomograms.get().iterItems()]

        tomoProvider = TomogramsTreeProvider(tomoList, self._getExtraPath(), "eman")

        self.dlg = EmanDialog(None, self._getExtraPath(), provider=tomoProvider, inMemory=self.inMemory.get(),)

        # Open dialog to request confirmation to create output
        import Tkinter as tk
        frame = tk.Frame()
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
            self._createOutput()

    def _validate(self):
        errors = []

        if not eman2.Plugin.isTomoAvailableVersion():
            errors.append('Your EMAN2 version does not support the tomo boxer. '
                          'Please update your installation to EMAN 2.3 or newer.')
        return errors

    def getMethods(self, output):
        msg = 'User picked %d particles ' % output.getSize()
        msg += 'with a particle size of %s.' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.inputTomograms is None:
            return ['Input tomogram not available yet.']

        methodsMsgs.append("Input tomograms imported of dims %s." %(
                              str(self.inputTomograms.get().getDim())))

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
