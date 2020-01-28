# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import pyworkflow.viewer as pwviewer
import pyworkflow.em.viewers.views as vi
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils

import tomo.objects

from eman2.convert import setCoords2Jsons, jsons2SetCoords


class EmanDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.SetOfCoordinates3D
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return vi.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        if issubclass(cls, tomo.objects.SetOfCoordinates3D):
            from eman2.viewers.views_tkinter_tree import EmanDialog
            from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider
            from tomo.objects import SetOfCoordinates3D

            suffix = self.protocol._getOutputSuffix(SetOfCoordinates3D)
            suffix = int(suffix)-1
            prefix = self.protocol.OUTPUT_PREFIX

            if suffix > 1:
                suffix = str(suffix)
                name = prefix + suffix
            else:
                suffix = ''
                name = prefix

            check = pwutils.removeBaseExt(obj.getFileName())
            # FIXME Now it doesn't open one viewer per output but can't choose one explicitly
            if filter(str.isdigit, check) == suffix:
                outputCoords = getattr(self.protocol, name)

                tomoList = [item.clone() for item in outputCoords.getPrecedents().iterItems()]

                path = self.protocol._getTmpPath()

                tomoProvider = TomogramsTreeProvider(tomoList, path, 'json',)

                setCoords2Jsons(outputCoords.getPrecedents(), outputCoords, path)

                setView = EmanDialog(self._tkRoot, path, provider=tomoProvider)

                import Tkinter as tk
                frame = tk.Frame()
                if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
                    jsons2SetCoords(self.protocol, outputCoords.getPrecedents(), path)

        return views
