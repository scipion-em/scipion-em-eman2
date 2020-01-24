# **************************************************************************
# *
# * Authors:     David Herreros
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

import os, threading

from pyworkflow import utils as pwutils
from pyworkflow.utils.process import runJob
from pyworkflow.gui.dialog import ToolbarListDialog
from pyworkflow.utils.path import moveFile, cleanPath, copyFile

import eman2


class EmanDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an Eman subprocess from a list of Tomograms.
    """

    def __init__(self, parent, path, **kwargs):
        self.path = path
        self.provider = kwargs.get("provider", None)
        self.inMemory = kwargs.get("inMemory", None)
        self.dir = os.getcwd()
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                    allowsEmptySelection=False,
                                    itemDoubleClick=self.doubleClickOnTomogram,
                                    **kwargs)

    def refresh_gui(self):
        if self.proc.isAlive():
            self.after(1000, self.refresh_gui)
        else:
            outFile = 'extra-%s_info.json' % pwutils.removeBaseExt(self.tomo.getFileName())
            moveFile((os.path.join(self.path, "info", outFile)), os.path.join(self.path, outFile))
            cleanPath(os.path.join(self.path, "info"))
            self.tree.update()


    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchEmanForTomogram, args=(self.inMemory, self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def lanchEmanForTomogram(self, inMemory, tomo):
        pathCoor = self._moveCoordsToInfo(tomo)

        program = eman2.Plugin.getProgram("e2spt_boxer.py")
        arguments = "%s" % tomo.getFileName()
        if inMemory:
            arguments += " --inmemory"
        runJob(None, program, arguments, env=eman2.Plugin.getEnviron(), cwd=self.path)

    def _moveCoordsToInfo(self, tomo):
        cwd = os.getcwd()
        infoDir = pwutils.join(os.path.abspath(self.path), 'info')
        fnCoor = 'extra-%s_info.json' % pwutils.removeBaseExt(tomo.getFileName())
        pathCoor = os.path.join(infoDir, fnCoor)
        if os.path.exists(os.path.join(self.path, fnCoor)):
            pwutils.makePath(infoDir)
            copyFile(os.path.join(self.path, fnCoor), pathCoor)
        return pathCoor