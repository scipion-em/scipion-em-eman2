# **************************************************************************
# *
# *  Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *  Authors:     J.M. de la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] SciLifeLab, Stockholm University
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

import pyworkflow.em as em
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, IntParam,
                                        StringParam, BooleanParam)
from pyworkflow.utils import createLink, copyTree, removeBaseExt


import eman2
import eman2.convert
from eman2.constants import *


class EmanProtEvalRefine(em.ProtAnalysis3D):
    """
    This protocol wraps *e2evalrefine.py* EMAN2 program.

    This program performs various assessments of
    e2refine_easy (or similar) runs, and operates in
    one of several possible modes.
"""
    _label = 'evaluate refine'

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {
            'mask': self._getExtraPath('mask.hdf'),
            'partSetGood': self._getExtraPath('sets/pf01_good_inputSet.lst'),
            'partSetBad': self._getExtraPath('sets/pf01_bad_inputSet.lst'),
            'partSetGood2': self._getExtraPath('sets/pf01_good_inputSet__ctf_flip.lst'),
            'partSetBad2': self._getExtraPath('sets/pf01_bad_inputSet__ctf_flip.lst'),
            #'partSetGoodInvar': self._getExtraPath('sets/pf01_good_inputSet__ctf_flip_invar.lst'),
            #'partSetBadInvar': self._getExtraPath('sets/pf01_bad_inputSet__ctf_flip_invar.lst'),
        }
        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputProt', PointerParam,
                      label='Choose e2refine_easy protocol',
                      pointerClass='EmanProtRefine',
                      help='Select a completed EMAN refinement protocol.')
        form.addParam('evalPtcls', BooleanParam, default=True,
                      label='Evaluate particle-map agreement',
                      help='Evaluates the particle-map agreement. '
                           'This may be used to identify bad particles.')
        form.addParam('evalCls', BooleanParam, default=False,
                      label='Evaluate class-average-projection agreement')
        form.addParam('evalAniso', IntParam, default=-1,
                      label='Class number to evaluate anisotropy',
                      help='Assesses the amount and direction of any '
                           'magnification anisotropy present in a raw '
                           'data set by considering particles in a '
                           'range of different orientations. Works '
                           'best with large particles. Specify a '
                           'class-average number for a class '
                           'containing many particles (within a '
                           'particular refine_xx folder). It is '
                           'a good idea to compare results among '
                           'classes in different highly occupied '
                           'orientations.')
        form.addParam('includeProj', BooleanParam, default=False,
                      label='Write out projections?',
                      help='If specified with --evalptclqual, projections '
                           'will be written to disk for easy comparison.')
        form.addParam('iter', IntParam, default=-1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Iteration',
                      help='If a refine_XX folder is being used, '
                           'this selects a particular refinement '
                           'iteration. _-1_ means the last complete '
                           'iteration is used.')
        form.addParam('evalClsDetail', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Generate individual FRC curves?',
                      help='If specified with --evalclassqual, will '
                           'generate individual FRC curves for '
                           'each class average in the even subset')
        form.addParam('inputMask', PointerParam,
                      expertLevel=LEVEL_ADVANCED,
                      pointerClass='VolumeMask', allowsNull=True,
                      label='Input 3D mask',
                      help='Mask to be used to focus --evalptclqual '
                           'and other options. May be useful for '
                           'separating heterogeneous data.')
        form.addParam('extraParams', StringParam, default='',
                      expertLevel=LEVEL_ADVANCED,
                      label='Additional arguments:',
                      help='In this box command-line arguments may be provided '
                           'that are not generated by the GUI. This may be '
                           'useful for testing developmental options and/or '
                           'expert use of the program. \n'
                           'The command "e2evalrefine.py -h" will print a list '
                           'of possible options.')

        form.addParallelSection(threads=4, mpi=0)

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self.convertInputStep()
        self._insertFunctionStep('createLinksStep')
        args = self._prepareParams()
        self._insertFunctionStep('evalStep', args)
        self._insertFunctionStep('createOutputStep')

    # -------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        if self.inputMask.hasValue():
            mask = self.inputMask.get()
            orig = os.path.relpath(mask.getFileName(), self._getExtraPath())
            maskBase = os.path.basename(self._getFileName("mask"))
            args = "%s %s --apix=%0.3f" % (orig, maskBase,
                                           mask.getSamplingRate())
            self.runJob(eman2.Plugin.getProgram('e2proc3d.py'), args,
                        cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)

    def createLinksStep(self):
        prot = self._inputProt()
        presResDir = prot._getExtraPath("refine_01")
        currResDir = self._getExtraPath("refine_01")
        prevPartDir = prot._getExtraPath("particles")
        currPartDir = self._getExtraPath("particles")
        prevSetsDir = prot._getExtraPath("sets")
        currSetsDir = self._getExtraPath("sets")

        createLink(presResDir, currResDir)
        createLink(prevPartDir, currPartDir)
        copyTree(prevSetsDir, currSetsDir)

    def evalStep(self, args):
        """ Run the EMAN program to evalrefine. """
        program = eman2.Plugin.getProgram('e2evalrefine.py')
        # mpi and threads are handled by EMAN itself
        try:
            self.runJob(program, args, cwd=self._getExtraPath(),
                        numberOfMpi=1, numberOfThreads=1)
        except Exception:
            print("Seems that e2evalrefine.py has failed, although files were "
                  "generated.")

    def createOutputStep(self):
        skipCtf = self._inputProt().skipctf.get()
        inputSet = self._inputProt()._getInputParticles()

        def _nextItem(outputIter):
            try:
                item = next(outputIter)
                r = None if item is None else item
            except StopIteration:
                r = None
            return r

        def _createOutputTuple(label, fileKey, fileKey2):
            outputSet = self._createSetOfParticles(suffix='_%s' % label)
            outputSet.copyInfo(inputSet)
            outputSet.setIsPhaseFlipped(True)
            outputSet.setHasCTF(True)
            outputSet.setSamplingRate(inputSet.getSamplingRate())
            fn = fileKey if skipCtf else fileKey2
            outputIter = eman2.convert.iterLstFile(self._getFileName(fn))
            return outputSet, outputIter, _nextItem(outputIter)

        goodSet, goodIter, goodFirst = _createOutputTuple('good', 'partSetGood', 'partSetGood2')
        badSet, badIter, badFirst = _createOutputTuple('bad', 'partSetBad', 'partSetBad2')

        for _, part in eman2.convert.iterParticlesByMic(inputSet):
            loc, fn = part.getIndex(), part.getCoordinate().getMicName()
            if loc == goodFirst[0] and removeBaseExt(fn) == removeBaseExt(goodFirst[1]):
                goodSet.append(part)
                goodFirst = _nextItem(goodIter)
            elif loc == badFirst[0] and removeBaseExt(fn) == removeBaseExt(badFirst[1]):
                badSet.append(part)
                badFirst = _nextItem(badIter)
            else:
                raise Exception("Particle %d@%s (id=%d) not found in any set"
                                % (part.getIndex(), part.getFileName(), part.getObjId()))

        self._defineOutputs(outputParticlesGood=goodSet,
                            outputParticlesBad=badSet)
        self._defineSourceRelation(inputSet, goodSet)
        self._defineSourceRelation(inputSet, badSet)

    # -------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []

        return summary

    def _methods(self):
        methods = []

        return methods

    # -------------------------- UTILS functions ------------------------------
    def _prepareParams(self):
        args = " --threads %d" % self.numberOfThreads.get()

        if self.evalPtcls:
            args += " --evalptclqual"

        if self.evalCls:
            args += " --evalclassqual"

        if self.evalAniso.get() != -1:
            args += " --anisitropy %d" % self.evalAniso.get()

        if self.includeProj:
            args += " --includeprojs"

        if self.iter.get() != -1:
            args += " --iter %d" % self.iter.get()

        if self.evalClsDetail:
            args += " --evalclassdetail"

        if self.inputMask.hasValue():
            args += " --mask %s" % os.path.basename(self._getFileName("mask"))

        if self.extraParams.hasValue():
            args += " " + self.extraParams.get()

        args += " refine_01"

        return args

    def _inputProt(self):
        return self.inputProt.get()
