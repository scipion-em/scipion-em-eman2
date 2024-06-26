# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.protocol.params import (IntParam, FloatParam,
                                        EnumParam, PointerParam,
                                        StringParam, USE_GPU,
                                        GPU_LIST, BooleanParam)
from pyworkflow.utils import makePath, createLink
from pyworkflow.constants import PROD
from pwem.protocols import ProtParticlePickingAuto

from .. import Plugin
from ..convert import readSetOfCoordinates, convertReferences
from ..constants import AUTO_CONVNET, AUTO_GAUSS


class EmanProtAutopick(ProtParticlePickingAuto):
    """ Automated particle picker for SPA. Uses EMAN2 (versions 2.2+) e2boxer.py
    """
    _label = 'boxer auto'
    _devStatus = PROD

    def _createFilenameTemplates(self):
        """ Centralize the names of the files. """

        myDict = {'goodRefsFn': self._getExtraPath('info/boxrefs.hdf'),
                  'badRefsFn': self._getExtraPath('info/boxrefsbad.hdf'),
                  'bgRefsFn': self._getExtraPath('info/boxrefsbg.hdf'),
                  'nnetFn': self._getExtraPath('nnet_pickptcls.hdf'),
                  'nnetClFn': self._getExtraPath('nnet_classify.hdf'),
                  'trainoutFn': self._getExtraPath('trainout_pickptcl.hdf'),
                  'trainoutClFn': self._getExtraPath('trainout_classify.hdf')
                  }
        self._updateFilenamesDict(myDict)

    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtParticlePickingAuto._defineParams(self, form)
        form.addHidden(USE_GPU, BooleanParam, default=False,
                       label="Use GPU?",
                       help="Set to Yes if you want to run Neural Net "
                            "boxer on GPU. Default is CPU.")
        form.addHidden(GPU_LIST, StringParam, default='0',
                       label="Choose GPU ID",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on.\n"
                            "Eman boxer can use only one GPU.")
        form.addParam('boxSize', IntParam, default=128,
                      allowsPointers=True,
                      label='Box size (px)',
                      help="Box size in pixels. See https://eman2.org/BoxSize")
        form.addParam('particleSize', IntParam, default=100,
                      label='Particle size (px)',
                      help="Longest axis of particle in pixels (diameter, "
                           "not radius).")
        form.addParam('boxerMode', EnumParam,
                      choices=['local search', 'by ref', 'neural net', 'gauss'],
                      label="Autopicker mode:", default=AUTO_CONVNET,
                      display=EnumParam.DISPLAY_COMBO,
                      help="Choose autopicker mode:\n\n"
                           " _local search_ - Heavily downsamples the "
                           "particles and references, and actually performs "
                           "a 2-D alignment of each putative particle to each "
                           "reference to identify the best particles in the "
                           "image. In theory this should produce fewer "
                           "false positives. Reference requirements are "
                           "similar to By Ref, though a smaller number of "
                           "references may be fine.\n"
                           " _by ref_ - This is a classic reference based "
                           "particle picker. To use it, you need to have "
                           "high quality good references in all possible "
                           "3-D orientations. The algorithm will do in-plane "
                           "rotation of the references and cross-correlate "
                           "to look for peaks. It is recommended that you "
                           "use projections of a 3-D map (low resolution) "
                           "as references.\n"
                           " _neural net_ - This is the most accurate boxer "
                           "by far, both in terms of false positives and "
                           "false negatives, when trained properly. It is "
                           "based on a pair of neural networks, one to "
                           "discriminate between putative particles and "
                           "the background, and a second to discriminate "
                           "between real particles and contamination or "
                           "other high contrast non-particles, which is "
                           "why it has two thresholds.\n"
                           " _gauss_ - This is a simple and fast "
                           "reference-free picker, which provides a simple "
                           "solution for easy particle picking cases, "
                           "where the particles have good contrast and are "
                           "monodisperse. This may work well for things "
                           "like viruses or ribosomes. It is very fast "
                           "and since it requires no references, it's "
                           "easy to try. It likely won't work well for "
                           "most projects. It was ported from the old "
                           "boxer program by a volunteer (Vadim Kotov).")
        form.addParam('threshold', FloatParam, default=5.0,
                      label='Threshold',
                      condition='boxerMode!=%d' % AUTO_GAUSS)
        form.addParam('threshold2', FloatParam, default=-5.0,
                      condition='boxerMode==%d' % AUTO_CONVNET,
                      label='Threshold2')
        form.addParam('gaussLow', FloatParam, default=1.,
                      condition='boxerMode==%d' % AUTO_GAUSS,
                      label='Threshold low')
        form.addParam('gaussHigh', FloatParam, default=2.,
                      condition='boxerMode==%d' % AUTO_GAUSS,
                      label='Threshold high')
        form.addParam('gaussWidth', FloatParam, default=1.,
                      condition='boxerMode==%d' % AUTO_GAUSS,
                      label='Gaussian width')

        form.addSection('References')
        form.addParam('boxerProt', PointerParam,
                      pointerClass='EmanProtBoxing',
                      condition='boxerMode==%d' % AUTO_CONVNET,
                      label='Previous e2boxer protocol',
                      help='Provide previously executed e2boxer protocol '
                           'that has all 3 types of references and '
                           'pre-trained neural network.')
        form.addParam('goodRefs', PointerParam,
                      pointerClass='SetOfAverages',
                      condition='boxerMode<%d' % AUTO_CONVNET,
                      label="Good references",
                      help="Good particle references.")

        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertInitialSteps(self):
        self._createFilenameTemplates()
        initId = self._insertFunctionStep('convertInputStep')
        return [initId]

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        goodRefs = self.goodRefs.get() if self.goodRefs.hasValue() else None
        boxerProt = self.boxerProt.get() if self.boxerProt.hasValue() else None
        storePath = self._getExtraPath("info")
        makePath(storePath)

        if goodRefs is not None:
            convertReferences(goodRefs, self._getFileName('goodRefsFn'))

        if boxerProt is not None:
            boxerProt._createFilenameTemplates()
            keys = ['goodRefsFn', 'badRefsFn', 'bgRefsFn',
                    'nnetFn', 'nnetClFn',
                    'trainoutFn', 'trainoutClFn']

            for fn in keys:
                if os.path.exists(boxerProt._getFileName(fn)):
                    createLink(boxerProt._getFileName(fn),
                               self._getFileName(fn))

    def _pickMicrograph(self, mic, *args):
        micFile = os.path.relpath(mic.getFileName(), self.getCoordsDir())
        params = " --apix=%f --no_ctf" % self.inputMicrographs.get().getSamplingRate()
        params += " --boxsize=%d" % self.boxSize.get()
        params += " --ptclsize=%d" % self.particleSize.get()
        params += " --threads=%d" % self.numberOfThreads.get()

        modes = ['auto_local', 'auto_ref', 'auto_convnet', 'auto_gauss']
        params += " --autopick=%s" % modes[self.boxerMode.get()]

        if self.boxerMode.get() == AUTO_GAUSS:
            params += ":gauss_width=%0.3f:thr_low=%0.3f:thr_high=%0.3f:boxsize=%d" % (
                self.gaussWidth.get(), self.gaussLow.get(),
                self.gaussHigh.get(), self.boxSize.get())
        else:
            params += ":threshold=%0.2f" % self.threshold.get()

        if self.boxerMode.get() == AUTO_CONVNET:
            params += ":threshold2=%0.2f" % self.threshold2.get()

            if self.useGpu:
                params += " --device=gpu%s" % self.gpuList.get().strip()
            else:
                params += " --device=cpu"

        params += ' %s' % micFile
        program = Plugin.getProgram('e2boxer.py')

        self.runJob(program, params, cwd=self.getCoordsDir())

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []

        if self.useGpu and (self.boxerMode.get() != AUTO_CONVNET):
            errors.append("You can use GPU only for neural net picker!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def getCoordsDir(self):
        return self._getExtraPath()

    def getFiles(self):
        return (self.inputMicrographs.get().getFiles() |
                ProtParticlePickingAuto.getFiles(self))

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        coordSet.setBoxSize(self.boxSize.get())
        readSetOfCoordinates(workingDir, micList, coordSet)
