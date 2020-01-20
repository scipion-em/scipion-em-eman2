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

import os

from pyworkflow.utils.properties import Message
from pyworkflow import utils as pwutils
from pyworkflow.protocol.params import (PointerParam, IntParam,
                                        StringParam, FloatParam, LEVEL_ADVANCED)
from pyworkflow.utils.path import moveFile

import eman2
from eman2.convert import loadJson, readSetOfCoordinates3D

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D


class EmanProtTomoTempMatch(ProtTomoPicking):
    """
    This protocol wraps *e2spt_tempmatch.py* EMAN2 program.

    It will perform a sweep of an initial volume against a tomogram
    to find correlation peaks and extract the corresponding subtomogram
    coordinates

    """

    _label = 'template matching'

    def __init__(self, **args):
        ProtTomoPicking.__init__(self, **args)

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input tomograms", important=True,
                      help='Specify tomograms containing reference-like particles\n'
                           'to be extracted')
        form.addParam('ref', PointerParam, pointerClass="Volume",
                      label='Reference volume',
                      help='Specify a 3D volume')
        form.addParam('nptcl', IntParam, default=500,
                      label='Number of particles',
                      help='Maximum number of particles per tomogram')
        form.addParam('dthr', FloatParam, default=-1,
                      label='Distance threshold',
                      help='', expertLevel=LEVEL_ADVANCED)
        form.addParam('vthr', FloatParam, default=2.0,
                      label='Value threshold',
                      help='', expertLevel=LEVEL_ADVANCED)
        form.addParam('delta', FloatParam, default=30.0,
                      label='Delta angle',
                      help='', expertLevel=LEVEL_ADVANCED)
        form.addParam('sym', StringParam, default='c1',
                      label='Point-group symmetry',
                      help='')
        form.addParam('boxSize', FloatParam, important=True, label='Box size',
                      help="The wizard selects same box size as reference size")

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('preprocess')
        self._insertFunctionStep('tempMatchStep')
        self._insertFunctionStep("createOutputStep")


    # --------------------------- STEPS functions -----------------------------

    def preprocess(self):
        program = eman2.Plugin.getProgram("e2proc3d.py")
        setDim = self.inputSet.get().getDim()
        if max(setDim) > 1000:
            sizeThreshold = max(self.inputSet.get().getDim())
        else:
            sizeThreshold = 1000
        if (setDim[0] < sizeThreshold) or (setDim[1] < sizeThreshold):
            for tomo in self.inputSet.get():
                self.runJob(program, '%s %s --clip=%d,%d,%d' % (tomo.getFileName(), tomo.getFileName(),
                                                             sizeThreshold, sizeThreshold, tomo.getDim()[2]),
                            env=eman2.Plugin.getEnviron())

    def tempMatchStep(self):

        self.box = self.boxSize.get()

        volFile = os.path.abspath(self.ref.get().getFileName())
        params = ""

        for tomo in self.inputSet.get():
            params = params + " %s" % tomo.getFileName()

        params = params + " --reference=%s --nptcl=%d --dthr=%f --vthr=%f --delta=%f --sym=%s " \
                          "--rmedge --rmgold --boxsz=%d" % (volFile, self.nptcl.get(), self.dthr.get(),
                          self.vthr.get(), self.delta.get(), self.sym.get(), self.box)

        program = eman2.Plugin.getProgram("e2spt_tempmatch.py")

        self.runJob(program, params, cwd=os.path.abspath(self._getTmpPath()),
                    env=eman2.Plugin.getEnviron())

        # Move output files to Extra Path
        moveFile(self._getTmpPath("ccc.hdf"), self._getExtraPath("particles" + ".hdf"))

        for tomo in self.inputSet.get():
            tomoName = os.path.basename(tomo.getFileName())
            tomoName = os.path.splitext(tomoName)[0]
            tomoCoord = "extra-" + tomoName + "_info.json"
            moveFile(self._getTmpPath(os.path.join("info", tomoCoord)),
                     self._getExtraPath("extra-" + tomoName + "_info.json"))

    # --------------------------- UTILS functions ------------------------------
    def createOutputStep(self):
        # Create a Set of 3D Coordinates per class
        coord3DSetDict = {}
        coord3DMap = {}
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coord3DSet = self._createSetOfCoordinates3D(self.inputSet.get(), suffix)
        coord3DSet.setBoxSize(self.box)
        coord3DSet.setName("tomoCoord")
        coord3DSet.setVolumes(self.inputSet.get())
        coord3DSet.setSamplingRate(self.inputSet.get().getSamplingRate())

        for tomo in self.inputSet.get():
            inputTomo = tomo.clone()
            tomoName = os.path.basename(tomo.getFileName())
            tomoName = os.path.splitext(tomoName)[0]

            jsonFnbase = pwutils.join(self._getExtraPath(), 'extra-%s_info.json' % tomoName)
            jsonBoxDict = loadJson(jsonFnbase)

            for key, classItem in jsonBoxDict["class_list"].iteritems():
                index = int(key)
                coord3DSetDict[index] = coord3DSet
                name = self.OUTPUT_PREFIX + suffix
                coord3DMap[index] = name
                args = {}
                args[name] = coord3DSet
                # Populate Set of 3D Coordinates with 3D Coordinates
                readSetOfCoordinates3D(jsonBoxDict, coord3DSetDict, inputTomo)

        self._defineOutputs(**args)
        self._defineSourceRelation(self.inputSet.get(), coord3DSet)

        # Update Outputs
        for index, coord3DSet in coord3DSetDict.iteritems():
            coord3DSet.setObjComment(self.getSummary(coord3DSet))
            self._updateOutputSet(coord3DMap[index], coord3DSet, state=coord3DSet.STREAM_CLOSED)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    def getSummary(self, coord3DSet):
        summary = []
        summary.append("Number of particles picked: %s" % coord3DSet.getSize())
        summary.append("Particle size: %s" % coord3DSet.getBoxSize())
        return "\n".join(summary)

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output 3D Coordinates not ready yet.")

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:* \n %s " % (key, output.getObjComment()))
        else:
            summary.append(Message.TEXT_NO_OUTPUT_CO)
        return summary

    def _methods(self):
        tomos = self.inputSet.get()
        return [
            "Subtomogram coordinates obtained with e2spt_tempmatch.py",
            "A total of %d tomograms of dimensions %s were used"
               % (tomos.getSize(), tomos.getDimensions()),
        ]