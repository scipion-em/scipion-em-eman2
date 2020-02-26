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

from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.protocol.params import BooleanParam, PointerParam, LEVEL_ADVANCED

import eman2
from eman2.convert import setCoords2Jsons, jsons2SetCoords
from eman2.viewers.views_tkinter_tree import EmanDialog

from tomo.protocols import ProtTomoPicking
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
        jsons2SetCoords(self, self.inputTomograms.get(), self._getExtraPath())

    # --------------------------- STEPS functions -----------------------------
    def copyInputCoords(self):
        setCoords2Jsons(self.inputTomograms.get(), self.inputCoordinates.get(), self._getExtraPath())

    def launchBoxingGUIStep(self):

        tomoList = [tomo.clone() for tomo in self.inputTomograms.get().iterItems()]

        tomoProvider = TomogramsTreeProvider(tomoList, self._getExtraPath(), "json")

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
