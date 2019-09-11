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

import os
from glob import glob

from tomo.objects import SetOfTomograms, Tomogram
from pyworkflow import utils as pwutils
from pyworkflow.protocol.params import (PointerParam, IntParam,
                                        BooleanParam, LEVEL_ADVANCED,
                                        StringParam, FloatParam)
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Volume
from pyworkflow.utils.path import moveFile, cleanPath, makePath

import eman2
from eman2.convert import loadJson

from tomo.protocols import ProtTomoPicking
from tomo.objects import Coordinate3D, SetOfCoordinates3D


class EmanProtTempMatch(ProtTomoPicking):
    """
    This protocol wraps *e2spt_tempmatch.py* EMAN2 program.

    It will perform a sweep of an initial volume against a tomogram
    to find correlation peaks and extract those particle coordinates

    """

    _label = 'template matching'

    def __init__(self, **args):
        ProtTomoPicking.__init__(self, **args)

    # --------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam,
                      pointerClass='Tomogram',
                      label="Input tomograms", important=True,
                      help='Specify tomograms containing reference-like particles\n'
                           'to be extracted')
        form.addParam('ref', PointerParam, pointerClass="Volume",
                      label='Reference volume',
                      help='Specify a 3D volume')
        form.addParam('nptcl', IntParam, default=500,
                      label='Number of particles',
                      help='Maximum number of particles')
        form.addParam('dthr', FloatParam, default=16.0,
                      label='Distance threshold',
                      help='')
        form.addParam('vthr', FloatParam, default=2.0,
                      label='Value threshold',
                      help='')
        form.addParam('delta', FloatParam, default=30.0,
                      label='Delta angle',
                      help='')
        form.addParam('sym', StringParam, default='c1',
                      label='Point-group symmetry',
                      help='')
        form.addParam('boxSize', IntParam, default=126, label='Box size')

        form.addParallelSection(threads=8, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._insertFunctionStep("convertInput")

        input = self.inputSet.get()

        deps = []
        if isinstance(input,SetOfTomograms):
            for tomo in input.iterItems():
                copyTomo = tomo.clone()
                deps.append(self._insertFunctionStep('tempMatchStep', copyTomo,  prerequisites=[]))
        else:
                deps.append(self._insertFunctionStep('tempMatchStep', input, prerequisites=[]))


    # --------------------------- STEPS functions -----------------------------

    def convertInput(self):
        img = ImageHandler()
        makePath(self._getTmpPath("Vol"))
        img.convert(self.ref.get(), self._getTmpPath(os.path.join("Vol","ref_vol.mrc")))

    def tempMatchStep(self,inputTomo):

        tomoName = inputTomo.getFileName()

        makePath(self._getTmpPath("tmpDir"))

        volFile = os.path.abspath(self._getTmpPath(os.path.join("Vol","ref_vol.mrc")))

        params = " %s --reference=%s --nptcl=%d --dthr=%f --vthr=%f --delta=%f --sym=%s --rmedge --rmgold --boxsz=%d" %(
                 tomoName, volFile, self.nptcl.get(), self.dthr.get(), self.vthr.get(), self.delta.get(), self.sym.get(),
                 self.boxSize.get())

        self.runJob(eman2.Plugin.getTemplateCommand(), params, cwd=os.path.abspath(self._getTmpPath("tmpDir")),
                    env=eman2.Plugin.getEnviron())

        #Move output files to Extra Path
        baseTomo = os.path.basename(tomoName)
        baseTomo = os.path.splitext(baseTomo)[0]
        self.nameTomo = baseTomo


        moveFile(self._getTmpPath(os.path.join("tmpDir","ccc.hdf")),self._getExtraPath(baseTomo+".hdf"))
        moveFile(self._getTmpPath(os.path.join("tmpDir","info","extra-"+baseTomo+"_info.json")),
                 self._getExtraPath("extra-"+baseTomo+"_info.json"))
        cleanPath(self._getTmpPath("tmpDir"))

        self.createOutput(inputTomo)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions ------------------------------
    def getSummary(self, coord3DSet):
        summary = []
        summary.append("Number of particles picked: %s" % coord3DSet.getSize())
        summary.append("Particle size: %s" % coord3DSet.getBoxSize())
        return "\n".join(summary)

    def _readSetOfCoordinates3D(self, jsonBoxDict, coord3DSetDict, inputTomo):
        if jsonBoxDict.has_key("boxes_3d"):
            boxes = jsonBoxDict["boxes_3d"]

            for box in boxes:
                classKey = box[5]
                coord3DSet = coord3DSetDict[classKey]
                coord3DSet.enableAppend()

                self._readCoordinates3D(box, coord3DSet, inputTomo)

    def _readCoordinates3D(self, box, coord3DSet, inputTomo):
        x, y, z = box[:3]
        coord = Coordinate3D()
        coord.setPosition(x, y, z)
        coord.setVolume(inputTomo)
        coord3DSet.append(coord)


    def createOutput(self, inputTomo):

        tomoName = inputTomo.getFileName()
        baseTomo = os.path.basename(tomoName)
        baseTomo = os.path.splitext(baseTomo)[0]

        jsonFnbase = pwutils.join(self._getExtraPath(), 'extra-%s_info.json' % baseTomo)
        jsonBoxDict = loadJson(jsonFnbase)

        # Create a Set of 3D Coordinates per class
        coord3DSetDict = {}
        coord3DMap = {}

        for key, classItem in jsonBoxDict["class_list"].iteritems():
            index = int(key)
            suffix = self._getOutputSuffix(SetOfCoordinates3D)
            coord3DSet = self._createSetOfCoordinates3D(inputTomo, suffix)

            coord3DSet.setBoxSize(self.boxSize.get())
            coord3DSet.setName(baseTomo)

            name = self.OUTPUT_PREFIX + suffix
            args = {}
            args[name] = coord3DSet
            self._defineOutputs(**args)
            self._defineSourceRelation(inputTomo, coord3DSet)

            coord3DSetDict[index] = coord3DSet
            coord3DMap[index] = name

        # Populate Set of 3D Coordinates with 3D Coordinates
        self._readSetOfCoordinates3D(jsonBoxDict, coord3DSetDict, inputTomo)

        # Update Outputs
        for index, coord3DSet in coord3DSetDict.iteritems():
            coord3DSet.setObjComment(self.getSummary(coord3DSet))
            coord3DSet.setVolumes(inputTomo)
            self._updateOutputSet(coord3DMap[index], coord3DSet, state=coord3DSet.STREAM_CLOSED)